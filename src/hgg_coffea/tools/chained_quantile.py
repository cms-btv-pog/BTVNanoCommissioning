import json
from typing import Any, Dict, List, Optional

import awkward
import numpy
import xgboost


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
        return self._model.predict(asDMatrix) * (self._scale or 1.0) + (
            self._center or 0.0
        )

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
    try:
        model = xgboost.Booster()
        model.load_model(weights)
        return wrapped_xgb(model=model, scale=scale, center=center, variables=variables)
    except xgboost.core.XGBoostError:
        raise ValueError(f"could not find: {weights}")


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
            self.transforms[vargroup] = {}
            for varname, regions in varlist.items():
                for region, steps in regions.items():
                    if not correctIsolations and vargroup == "isolations":
                        continue
                    if not correctShowerShapes and vargroup == "shower_shapes":
                        continue
                    self.transforms[vargroup][varname] = {}
                    if vargroup == "isolations":
                        self.transforms[vargroup][varname][region] = {}
                        for stepname, config in steps.items():
                            keys = list(config.keys())
                            if "variables" in keys:
                                config["variables"] = [
                                    var.split(":=")[-1].strip()
                                    for var in config["variables"]
                                ]
                            if "weights_data" in keys and "weights_mc" in keys:
                                config["weights"] = config["weights_data"]
                                self.transforms[vargroup][varname][region][
                                    f"{stepname}_data"
                                ] = create_evaluator(**config)
                                config["weights"] = config["weights_mc"]
                                self.transforms[vargroup][varname][region][
                                    f"{stepname}_mc"
                                ] = create_evaluator(**config)
                            else:
                                self.transforms[vargroup][varname][region][
                                    stepname
                                ] = create_evaluator(**config)
                    elif vargroup == "shower_shapes":
                        self.transforms[vargroup][varname][region] = []
                        self.transforms[vargroup][varname][region].append(
                            create_evaluator(**steps)
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

        for var, xform in xforms:
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
        clf_vars = clf_mc.variables
        irho = clf_vars.index("fixedGridRhoAll")
        clf_stack_vars = [awkward.to_numpy(photons[name]) for name in clf_vars[:irho]]
        clf_stack_vars.append(awkward.to_numpy(rho))
        clf_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in clf_vars[irho + 1 :]]
        )
        clf_eval_vars = numpy.column_stack(clf_stack_vars)
        # conversion from T[-1,1] to probability [0,1]
        p_tail_data = numpy.ones(size=clf_eval_vars.shape[0])
        p_tail_mc = numpy.ones_like(p_tail_data)

        p_tail_data[isEB] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_data["isEB"](clf_eval_vars[isEB])) - 1.0)
        )
        p_tail_data[isEE] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_data["isEE"](clf_eval_vars[isEE])) - 1.0)
        )
        p_tail_mc[isEB] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_mc["isEB"](clf_eval_vars[isEB])) - 1.0)
        )
        p_tail_mc[isEE] = 1.0 / (
            1.0 + numpy.sqrt(2.0 / (1.0 + clf_mc["isEE"](clf_eval_vars[isEE])) - 1.0)
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

        p2t_vars = p2t.variables
        irho = p2t.variables.index("fixedGridRhoAll")
        irnd = p2t.variables.index("peak2tail_rnd")
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
        pfPhoIso[isEB & to_tail] = p2t["EB"](p2t_eval_vars[isEB & to_tail])
        pfPhoIso[isEE & to_tail] = p2t["EE"](p2t_eval_vars[isEE & to_tail])
        pfPhoIso[to_peak] = 0

        # update photon for morph
        photons["pfPhoIso03"] = pfPhoIso

        # morphing
        needs_morph = pfPhoIso > 0
        morph_vars = morphing.variables
        irho = clf_vars.index("fixedGridRhoAll")
        morph_stack_vars = [
            awkward.to_numpy(photons[name]) for name in morph_vars[:irho]
        ]
        morph_stack_vars.append(awkward.to_numpy(rho))
        clf_stack_vars.extend(
            [awkward.to_numpy(photons[name]) for name in morph_vars[irho + 1 :]]
        )

        morph_eval_vars = numpy.column_stack(morph_stack_vars)
        pfPhoIso[isEB & needs_morph] = pfPhoIso[needs_morph] + morphing["EB"](
            morph_eval_vars
        )
        pfPhoIso[isEE & needs_morph] = pfPhoIso[needs_morph] + morphing["EE"](
            morph_eval_vars
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
        return photons

    def apply(self, events: awkward.Array) -> awkward.Array:
        # We're going to work in flattened data within this
        # function. Less mind-bending.
        photons = events.Photon
        rho = awkward.ones_list(photons.pt) * events.fixedGridRhoAll
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
