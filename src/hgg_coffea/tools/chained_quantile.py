import json
from typing import Any, Dict, List, Optional

import awkward
import numpy
import xgboost

from hgg_coffea.tools.bdt_loader import load_bdt


class wrapped_xgb:
    def __init__(
        self,
        model: xgboost.Booster,
        scale: Optional[float],
        center: Optional[float],
        variables: Optional[List[str]],
    ) -> None:
        self._model = model
        self._scale = scale
        self._center = center
        self._variables = variables

    def __call__(self, array: numpy.ndarray) -> numpy.ndarray:
        asDMatrix = xgboost.DMatrix(array)
        out = self._model.predict(asDMatrix)
        return out * (self._scale or 1.0) + (self._center or 0.0)

    @property
    def scale(self) -> Optional[float]:
        return self._scale

    @property
    def center(self) -> Optional[float]:
        return self._center

    @property
    def variables(self) -> Optional[List[str]]:
        return self._variables


def create_evaluator(
    weights: str,
    scale: Optional[float] = None,
    center: Optional[float] = None,
    variables: Optional[List[str]] = None,
    **kwargs: Dict[Any, Any],
) -> xgboost.Booster:
    model = load_bdt(weights)
    if model is None:
        raise RuntimeError(f"Could not load {weights}, check warnings!")
    return wrapped_xgb(model=model, scale=scale, center=center, variables=variables)


class ChainedQuantileRegression:
    def __init__(
        self,
        corrections_summary: str,
        SS_variables: List[str],
        correctShowerShapes: bool = True,
        correctIsolations: bool = True,
        correctPreshower: bool = False,
    ):
        with open(corrections_summary) as f:
            cq_config = json.load(f)

        self.transforms: Dict[str, Any] = {}

        for vargroup, varlist in cq_config.items():
            if not correctIsolations and vargroup == "isolations":
                continue
            if not correctShowerShapes and vargroup == "shower_shapes":
                continue
            self.transforms[vargroup] = {}
            for varname, regions in varlist.items():
                self.transforms[vargroup][varname] = {}
                for region, steps in regions.items():
                    if vargroup == "isolations":
                        for stepname, config in steps.items():
                            keys = list(config.keys())
                            if "variables" in keys:
                                config["variables"] = [
                                    var.split(":=")[-1].strip()
                                    for var in config["variables"]
                                ]
                            if "weights_data" in keys and "weights_mc" in keys:
                                if (
                                    f"{stepname}_data"
                                    not in self.transforms[vargroup][varname]
                                ):
                                    self.transforms[vargroup][varname][
                                        f"{stepname}_data"
                                    ] = {}
                                if (
                                    f"{stepname}_mc"
                                    not in self.transforms[vargroup][varname]
                                ):
                                    self.transforms[vargroup][varname][
                                        f"{stepname}_mc"
                                    ] = {}
                                config["weights"] = config["weights_data"]
                                self.transforms[vargroup][varname][f"{stepname}_data"][
                                    region
                                ] = create_evaluator(**config)
                                config["weights"] = config["weights_mc"]
                                self.transforms[vargroup][varname][f"{stepname}_mc"][
                                    region
                                ] = create_evaluator(**config)
                            else:
                                if stepname not in self.transforms[vargroup][varname]:
                                    self.transforms[vargroup][varname][stepname] = {}
                                self.transforms[vargroup][varname][stepname][
                                    region
                                ] = create_evaluator(**config)
                    elif vargroup == "shower_shapes":
                        self.transforms[vargroup][varname][region] = create_evaluator(
                            **steps
                        )
                    else:
                        raise Exception(
                            f"{vargroup} is not understood by ChainedQuantileRegression"
                        )
        self.ssvs = [ssv.split(":=")[-1].strip() for ssv in SS_variables]

    def apply_shower_shapes(
        self,
        photons: awkward.Array,
        rho: awkward.Array,
        isEB: numpy.ndarray,
        isEE: numpy.ndarray,
    ) -> awkward.Array:
        xforms = self.transforms["shower_shapes"]
        photons["uncorr_r9"] = photons.r9
        photons["uncorr_s4"] = photons.s4
        photons["uncorr_sieie"] = photons.sieie
        photons["uncorr_sieip"] = photons.sieip
        photons["uncorr_etaWidth"] = photons.etaWidth
        photons["uncorr_phiWidth"] = photons.phiWidth

        # Now for the xgboosty bits.
        # r9
        irho = self.ssvs.index("fixedGridRhoAll")
        stack_vars = [awkward.to_numpy(photons[name]) for name in self.ssvs[:irho]]
        stack_vars.append(awkward.to_numpy(rho))
        stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in self.ssvs[irho + 1 :]]
        )
        eval_vars = numpy.column_stack(stack_vars)
        eval_vars_eb = eval_vars[isEB]
        eval_vars_ee = eval_vars[~isEB]

        for var, xform in xforms.items():
            npy = awkward.to_numpy(photons[var])
            npy[isEB] = npy[isEB] + xform["EB"](eval_vars_eb)
            npy[isEE] = npy[isEE] + xform["EE"](eval_vars_ee)
            photons[var] = npy

        return photons

    def apply_photon_isolation(
        self,
        photons: awkward.Array,
        rho: awkward.Array,
        isEB: numpy.ndarray,
        isEE: numpy.ndarray,
    ) -> awkward.Array:
        xforms = self.transforms["isolations"]["phoIso"]

        clf_mc = xforms["peak_tail_clfs_mc"]
        clf_data = xforms["peak_tail_clfs_data"]
        p2t = xforms["peak2tail"]
        morphing = xforms["morphing"]

        photons["uncorr_pfPhoIso03"] = photons.pfPhoIso03

        # clfs (input variables are the same)
        clf_vars = clf_mc["EB"].variables
        irho = clf_vars.index("fixedGridRhoAll")
        clf_stack_vars = [awkward.to_numpy(photons[name]) for name in clf_vars[:irho]]
        clf_stack_vars.append(awkward.to_numpy(rho))
        clf_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in clf_vars[irho + 1 :]]
        )
        clf_eval_vars = numpy.column_stack(clf_stack_vars)
        # conversion from T[-1,1] to probability [0,1]
        p_tail_data = numpy.ones((clf_eval_vars.shape[0],))
        p_tail_mc = numpy.ones_like(p_tail_data)

        p_tail_data[isEB] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_data["EB"](clf_eval_vars[isEB])) - 1.0)
        )
        p_tail_data[isEE] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_data["EE"](clf_eval_vars[isEE])) - 1.0)
        )
        p_tail_mc[isEB] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_mc["EB"](clf_eval_vars[isEB])) - 1.0)
        )
        p_tail_mc[isEE] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_mc["EE"](clf_eval_vars[isEE])) - 1.0)
        )

        p_peak_data = 1 - p_tail_data
        p_peak_mc = 1 - p_tail_mc

        migration = numpy.random.uniform(size=clf_eval_vars.shape[0])
        pfPhoIso = awkward.to_numpy(photons.pfPhoIso03)

        p_move_to_tail = (p_tail_data - p_tail_mc) / p_peak_mc
        p_move_to_peak = (p_peak_data - p_peak_mc) / p_tail_mc

        # peak2tail
        to_tail = (
            (pfPhoIso == 0) & (p_tail_data > p_tail_mc) & (migration < p_move_to_tail)
        )
        to_peak = (
            (pfPhoIso > 0) & (p_peak_data > p_peak_mc) & (migration <= p_move_to_peak)
        )

        p2t_vars = p2t["EB"].variables
        irho = p2t_vars.index("fixedGridRhoAll")
        irnd = p2t_vars.index("peak2tail_rnd")
        p2t_stack_vars = [awkward.to_numpy(photons[name]) for name in p2t_vars[:irho]]
        p2t_stack_vars.append(awkward.to_numpy(rho))
        p2t_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in p2t_vars[irho + 1 : irnd]]
        )
        # https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/plugins/DifferentialPhoIdInputsCorrector.cc#L301
        p2t_stack_vars.append(
            numpy.random.uniform(low=0.01, high=0.99, size=clf_eval_vars.shape[0])
        )
        p2t_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in p2t_vars[irnd + 1 :]]
        )

        p2t_eval_vars = numpy.column_stack(p2t_stack_vars)
        if numpy.any(isEB & to_tail):
            pfPhoIso[isEB & to_tail] = p2t["EB"](p2t_eval_vars[isEB & to_tail])
        if numpy.any(isEE & to_tail):
            pfPhoIso[isEE & to_tail] = p2t["EE"](p2t_eval_vars[isEE & to_tail])
        pfPhoIso[to_peak] = 0

        # update photon for morph
        photons["pfPhoIso03"] = pfPhoIso

        # morphing
        needs_morph = pfPhoIso > 0
        morph_vars = morphing["EB"].variables
        irho = morph_vars.index("fixedGridRhoAll")
        morph_stack_vars = [
            awkward.to_numpy(photons[name]) for name in morph_vars[:irho]
        ]
        morph_stack_vars.append(awkward.to_numpy(rho))
        morph_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in morph_vars[irho + 1 :]]
        )

        morph_eval_vars = numpy.column_stack(morph_stack_vars)
        pfPhoIso[isEB & needs_morph] = pfPhoIso[isEB & needs_morph] + morphing["EB"](
            morph_eval_vars[isEB & needs_morph]
        )
        pfPhoIso[isEE & needs_morph] = pfPhoIso[isEE & needs_morph] + morphing["EE"](
            morph_eval_vars[isEE & needs_morph]
        )

        photons["pfPhoIso03"] = pfPhoIso

        return photons

    def apply_charged_isolation(
        self,
        photons: awkward.Array,
        rho: awkward.Array,
        isEB: numpy.ndarray,
        isEE: numpy.ndarray,
    ) -> awkward.Array:
        xforms = self.transforms["isolations"]["chIso"]

        clf_mc = xforms["peak_tail_clfs_mc"]
        clf_data = xforms["peak_tail_clfs_data"]
        p2t = xforms["chIso_peak2tail"]
        p2t_worst = xforms["chIsoWorst_peak2tail"]
        morphing = xforms["chIso_morphing"]
        morphing_worst = xforms["chIsoWorst_morphing"]

        photons["uncorr_pfChargedIsoPFPV"] = photons.pfChargedIsoPFPV
        photons["uncorr_pfChargedIsoWorstVtx"] = photons.pfChargedIsoWorstVtx

        # clfs (input variables are the same)
        clf_vars = clf_mc["EB"].variables
        irho = clf_vars.index("fixedGridRhoAll")
        clf_stack_vars = [awkward.to_numpy(photons[name]) for name in clf_vars[:irho]]
        clf_stack_vars.append(awkward.to_numpy(rho))
        clf_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in clf_vars[irho + 1 :]]
        )
        clf_eval_vars = numpy.column_stack(clf_stack_vars)
        # ---Charge isolations
        #  ----------------+-------------------------+
        #  ChIsoWorst tail | 01         | 11         |
        #  ChIsoWorst peak | 00         | X          |
        #  ----------------+-------------------------+
        #                  | ChIso peak | ChIso tail |
        #                  +-------------------------+
        probs_data = numpy.ones((clf_eval_vars.shape[0], 3))
        probs_mc = numpy.ones_like(probs_data)

        # [[00, 01, 11], ...]
        probs_data[isEB] = clf_data["EB"](clf_eval_vars[isEB])
        probs_data[isEE] = clf_data["EE"](clf_eval_vars[isEE])
        probs_mc[isEB] = clf_mc["EB"](clf_eval_vars[isEB])
        probs_mc[isEE] = clf_mc["EE"](clf_eval_vars[isEE])

        migration = numpy.random.uniform(size=clf_eval_vars.shape[0])
        migration_subcat = numpy.random.uniform(size=clf_eval_vars.shape[0])
        pfChgIso = awkward.to_numpy(photons.pfChargedIsoPFPV)
        pfChgIsoWorst = awkward.to_numpy(photons.pfChargedIsoWorstVtx)

        can_migrate = probs_mc > probs_data
        should_migrate = migration[:, None] < (1 - probs_data / probs_mc)

        # peak2tail
        to_00 = (
            (pfChgIso == 0)
            & (pfChgIsoWorst == 0)
            & can_migrate[:, 0]
            & should_migrate[:, 0]
        )
        to_01 = (
            (pfChgIso == 0)
            & (pfChgIsoWorst > 0)
            & can_migrate[:, 1]
            & should_migrate[:, 1]
        )
        to_11 = (
            (pfChgIso > 0)
            & (pfChgIsoWorst > 0)
            & can_migrate[:, 2]
            & should_migrate[:, 2]
        )

        p2t_vars = p2t["EB"].variables
        irho = p2t_vars.index("fixedGridRhoAll")
        irnd = p2t_vars.index("peak2tail_chIso_rnd")
        p2t_stack_vars = [awkward.to_numpy(photons[name]) for name in p2t_vars[:irho]]
        p2t_stack_vars.append(awkward.to_numpy(rho))
        p2t_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in p2t_vars[irho + 1 : irnd]]
        )
        # https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/plugins/DifferentialPhoIdInputsCorrector.cc#L301
        p2t_stack_vars.append(
            numpy.random.uniform(low=0.01, high=0.99, size=clf_eval_vars.shape[0])
        )
        p2t_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in p2t_vars[irnd + 1 :]]
        )
        # worst
        p2t_worst_vars = p2t_worst["EB"].variables
        irho_worst = p2t_worst_vars.index("fixedGridRhoAll")
        irnd_worst = p2t_worst_vars.index("peak2tail_chIsoWorst_rnd")
        p2t_worst_stack_vars = [
            awkward.to_numpy(photons[name]) for name in p2t_worst_vars[:irho_worst]
        ]
        p2t_worst_stack_vars.append(awkward.to_numpy(rho))
        p2t_worst_stack_vars.extend(
            [
                awkward.to_numpy(photons[name])
                for name in p2t_worst_vars[irho_worst + 1 : irnd_worst]
            ]
        )
        # https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/plugins/DifferentialPhoIdInputsCorrector.cc#L301
        p2t_worst_stack_vars.append(
            numpy.random.uniform(low=0.01, high=0.99, size=clf_eval_vars.shape[0])
        )
        p2t_worst_stack_vars.extend(
            [
                awkward.to_numpy(photons[name])
                for name in p2t_worst_vars[irnd_worst + 1 :]
            ]
        )

        def get_z(cat1: int, cat2: int) -> numpy.ndarray:
            return (probs_data[:, cat1] - probs_mc[:, cat1]) / (
                probs_mc[:, cat2] - probs_data[:, cat2]
            )

        p2t_eval_vars = numpy.column_stack(p2t_stack_vars)
        p2t_worst_eval_vars = numpy.column_stack(p2t_worst_stack_vars)
        # 00
        #  00 -> 01
        to_00_01 = to_00 & (~can_migrate[:, 1]) & can_migrate[:, 2]
        to_00_01_EB = isEB & to_00_01
        to_00_01_EE = isEE & to_00_01
        if numpy.any(to_00_01_EB):
            pfChgIsoWorst[to_00_01_EB] = p2t_worst["EB"](
                p2t_worst_eval_vars[to_00_01_EB]
            )
        if numpy.any(to_00_01_EE):
            pfChgIsoWorst[to_00_01_EE] = p2t_worst["EE"](
                p2t_worst_eval_vars[to_00_01_EE]
            )
        #  00 -> 11
        to_00_11 = to_00 & can_migrate[:, 1] & (~can_migrate[:, 2])
        to_00_11_EB = isEB & to_00_11
        to_00_11_EE = isEE & to_00_11
        if numpy.any(to_00_11_EB):
            pfChgIso[to_00_11_EB] = p2t["EB"](p2t_eval_vars[to_00_11_EB])
            pfChgIsoWorst[to_00_11_EB] = p2t_worst["EB"](
                p2t_worst_eval_vars[to_00_11_EB]
            )
        if numpy.any(to_00_11_EE):
            pfChgIso[to_00_11_EE] = p2t["EE"](p2t_eval_vars[to_00_11_EE])
            pfChgIsoWorst[to_00_11_EE] = p2t_worst["EE"](
                p2t_worst_eval_vars[to_00_11_EE]
            )
        #  00 -> either
        to_00_either = to_00 & (~can_migrate[:, 1]) & (~can_migrate[:, 2])
        to_00_either_EB = isEB & to_00_either
        to_00_either_EE = isEE & to_00_either
        which_01_11 = migration_subcat <= get_z(1, 0)
        if numpy.any(to_00_either_EB & which_01_11):
            pfChgIsoWorst[to_00_either_EB & which_01_11] = p2t_worst["EB"](
                p2t_worst_eval_vars[to_00_either_EB & which_01_11]
            )
        if numpy.any(to_00_either_EE & which_01_11):
            pfChgIsoWorst[to_00_either_EE & which_01_11] = p2t_worst["EE"](
                p2t_worst_eval_vars[to_00_either_EE & which_01_11]
            )
        if numpy.any(to_00_either_EB & (~which_01_11)):
            pfChgIso[to_00_either_EB & (~which_01_11)] = p2t["EB"](
                p2t_eval_vars[to_00_either_EB & (~which_01_11)]
            )
            pfChgIsoWorst[to_00_either_EB & (~which_01_11)] = p2t_worst["EB"](
                p2t_worst_eval_vars[to_00_either_EB & (~which_01_11)]
            )
        if numpy.any(to_00_either_EE & (~which_01_11)):
            pfChgIso[to_00_either_EE & (~which_01_11)] = p2t["EE"](
                p2t_eval_vars[to_00_either_EE & (~which_01_11)]
            )
            pfChgIsoWorst[to_00_either_EE & (~which_01_11)] = p2t_worst["EE"](
                p2t_worst_eval_vars[to_00_either_EE & (~which_01_11)]
            )

        # 01
        #  01 -> 00
        to_01_00 = to_01 & (~can_migrate[:, 0]) & can_migrate[:, 2]
        pfChgIso[to_01_00] = 0
        pfChgIsoWorst[to_01_00] = 0
        #  01 -> 11
        to_01_11 = to_01 & can_migrate[:, 0] & (~can_migrate[:, 2])
        to_01_11_EB = isEB & to_01_11
        to_01_11_EE = isEE & to_01_11
        if numpy.any(to_01_11_EB):
            pfChgIso[to_01_11_EB] = p2t["EB"](p2t_eval_vars[to_01_11_EB])
        if numpy.any(to_01_11_EE):
            pfChgIso[to_01_11_EE] = p2t["EE"](p2t_eval_vars[to_01_11_EE])
        #  01 -> either
        to_01_either = to_01 & (~can_migrate[:, 0]) & (~can_migrate[:, 2])
        to_01_either_EB = isEB & to_01_either
        to_01_either_EE = isEE & to_01_either
        which_00_11 = migration_subcat <= get_z(0, 1)
        if numpy.any(to_01_either_EB & (~which_00_11)):
            pfChgIso[to_01_either_EB & (~which_00_11)] = p2t["EB"](
                p2t_eval_vars[to_01_either_EB & (~which_00_11)]
            )
        if numpy.any(to_01_either_EE & (~which_00_11)):
            pfChgIso[to_01_either_EE & (~which_00_11)] = p2t["EE"](
                p2t_eval_vars[to_01_either_EE & (~which_00_11)]
            )
        pfChgIsoWorst[to_01_either & which_01_11] = 0

        # 11
        #  11 -> 00
        to_11_00 = to_11 & (~can_migrate[:, 0]) & can_migrate[:, 1]
        pfChgIso[to_01_00] = 0
        pfChgIsoWorst[to_11_00] = 0
        #  11 -> 01
        to_11_01 = to_11 & can_migrate[:, 0] & (~can_migrate[:, 1])
        pfChgIsoWorst[to_11_01] = 0
        #  11 -> either
        to_11_either = to_00 & (~can_migrate[:, 0]) & (~can_migrate[:, 1])
        which_00_01 = migration_subcat <= get_z(0, 2)
        pfChgIso[to_11_either] = 0
        pfChgIsoWorst[to_11_either & which_00_01] = 0

        # update photon for morph
        photons["pfChargedIsoPFPV"] = pfChgIso
        photons["pfChargedIsoWorstVtx"] = pfChgIsoWorst

        # tail morphing PV
        needs_morph = pfChgIso > 0
        morph_vars = morphing["EB"].variables
        irho = morph_vars.index("fixedGridRhoAll")
        morph_stack_vars = [
            awkward.to_numpy(photons[name]) for name in morph_vars[:irho]
        ]
        morph_stack_vars.append(awkward.to_numpy(rho))
        morph_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in morph_vars[irho + 1 :]]
        )

        morph_eval_vars = numpy.column_stack(morph_stack_vars)
        pfChgIso[isEB & needs_morph] = pfChgIso[isEB & needs_morph] + morphing["EB"](
            morph_eval_vars[isEB & needs_morph]
        )
        pfChgIso[isEE & needs_morph] = pfChgIso[isEE & needs_morph] + morphing["EE"](
            morph_eval_vars[isEE & needs_morph]
        )

        photons["pfChargedIsoPFPV"] = pfChgIso

        # tail morphing worst vertex
        needs_morph_worst = pfChgIsoWorst > 0
        morph_worst_vars = morphing_worst["EB"].variables
        irho = morph_worst_vars.index("fixedGridRhoAll")
        morph_worst_stack_vars = [
            awkward.to_numpy(photons[name]) for name in morph_worst_vars[:irho]
        ]
        morph_worst_stack_vars.append(awkward.to_numpy(rho))
        morph_worst_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in morph_worst_vars[irho + 1 :]]
        )

        morph_worst_eval_vars = numpy.column_stack(morph_worst_stack_vars)
        pfChgIsoWorst[isEB & needs_morph_worst] = pfChgIsoWorst[
            isEB & needs_morph_worst
        ] + morphing_worst["EB"](morph_worst_eval_vars[isEB & needs_morph_worst])
        pfChgIsoWorst[isEE & needs_morph_worst] = pfChgIsoWorst[
            isEE & needs_morph_worst
        ] + morphing_worst["EE"](morph_worst_eval_vars[isEE & needs_morph_worst])

        photons["pfChargedIsoWorstVtx"] = pfChgIsoWorst

        return photons

    def apply(self, events: awkward.Array) -> awkward.Array:
        # We're going to work in flattened data within this
        # function. Less mind-bending.
        photons = events.Photon
        rho = awkward.ones_like(photons.pt) * events.fixedGridRhoAll
        counts = awkward.num(photons, axis=1)
        photons = awkward.flatten(photons)
        rho = awkward.flatten(rho)

        isEB = awkward.to_numpy(numpy.abs(photons.eta) < 1.5)
        isEE = ~isEB

        if "shower_shapes" in self.transforms:
            photons = self.apply_shower_shapes(photons, rho, isEB, isEE)

        if "isolations" in self.transforms:
            for key in self.transforms["isolations"].keys():
                if "phoIso" == key:
                    photons = self.apply_photon_isolation(photons, rho, isEB, isEE)
                if "chIso" == key:
                    photons = self.apply_charged_isolation(photons, rho, isEB, isEE)

        return awkward.unflatten(photons, counts)
