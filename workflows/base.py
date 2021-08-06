import coffea
from coffea import hist, processor
import numpy as np
import xgboost as xg
import awkward as ak
import vector
import pandas as pd
import functools as ft
import operator as op
import os
import shutil as shu
import pathlib as pl
import sys
import warnings

vector.register_awkward()


class HggBaseProcessor(processor.ProcessorABC):
    def __init__(
        self,
        metaconditions,
        do_systematics,
        apply_trigger,
        output_location,
        taggers,
        trigger_group,
        analysis,
    ):
        self.meta = metaconditions
        self.do_systematics = do_systematics
        self.apply_trigger = apply_trigger
        self.output_location = output_location
        self.trigger_group = trigger_group
        self.analysis = analysis

        # diphoton preselection cuts
        self.min_pt_photon = 25.0
        self.min_pt_lead_photon = 35.0
        self.min_mvaid = -0.9
        self.max_sc_eta = 2.5
        self.gap_barrel_eta = 1.4442
        self.gap_endcap_eta = 1.566
        self.max_hovere = 0.08
        self.min_full5x5_r9 = 0.8
        self.max_chad_iso = 20.0
        self.max_chad_rel_iso = 0.3

        self.taggers = []
        if taggers is not None:
            self.taggers = taggers
            self.taggers.sort(key=lambda x: x.priority)

        self.prefixes = {"pho_lead": "lead", "pho_sublead": "sublead"}

        # initialize diphoton mva
        try:
            self.diphoton_mva = xg.Booster()
            self.diphoton_mva.load_model(self.meta["flashggDiPhotonMVA"]["weightFile"])
        except xg.core.XGBoostError:
            warnings.warn(
                f"SKIPPING diphoton_mva, could not find: {self.meta['flashggDiPhotonMVA']['weightFile']}"
            )
            self.diphoton_mva = None

    def photon_preselection(self, photons):
        photon_abs_eta = np.abs(photons.eta)
        return photons[
            (photons.pt > self.min_pt_photon)
            & (photon_abs_eta < self.max_sc_eta)
            & (
                (photon_abs_eta < self.gap_barrel_eta)
                | (photon_abs_eta > self.gap_endcap_eta)
            )
            & (photons.mvaID > self.min_mvaid)
            & (photons.hoe < self.max_hovere)
            & (
                (photons.r9 > self.min_full5x5_r9)
                | (photons.pfRelIso03_chg < self.max_chad_iso)
                | (photons.pfRelIso03_chg / photons.pt < self.max_chad_rel_iso)
            )
        ]

    def diphoton_list_to_pandas(self, diphotons):
        output = pd.DataFrame()
        for field in ak.fields(diphotons):
            prefix = self.prefixes.get(field, "")
            if len(prefix) > 0:
                for subfield in ak.fields(diphotons[field]):
                    output[f"{prefix}_{subfield}"] = ak.to_numpy(
                        diphotons[field][subfield]
                    )
            else:
                output[field] = ak.to_numpy(diphotons[field])
        return output

    def dump_pandas(self, pddf, fname, location, subdirs=[]):
        xrd_prefix = "root://"
        pfx_len = len(xrd_prefix)
        xrootd = False
        if xrd_prefix in location:
            try:
                import XRootD
                import XRootD.client

                xrootd = True
            except ImportError:
                raise ImportError(
                    "Install XRootD python bindings with: conda install -c conda-forge xroot"
                )
        local_file = (
            os.path.abspath(os.path.join(".", fname))
            if xrootd
            else os.path.join(".", fname)
        )
        subdirs = "/".join(subdirs) if xrootd else os.path.sep.join(subdirs)
        destination = (
            location + subdirs + f"/{fname}"
            if xrootd
            else os.path.join(location, os.path.join(subdirs, fname))
        )
        pddf.to_parquet(local_file)
        if xrootd:
            copyproc = XRootD.client.CopyProcess()
            copyproc.add_job(local_file, destination)
            copyproc.prepare()
            copyproc.run()
            client = XRootD.client.FileSystem(
                location[: location[pfx_len:].find("/") + pfx_len]
            )
            status = client.locate(
                destination[destination[pfx_len:].find("/") + pfx_len + 1 :],
                XRootD.client.flags.OpenFlags.READ,
            )
            assert status[0].ok
            del client
            del copyproc
        else:
            dirname = os.path.dirname(destination)
            if not os.path.exists(dirname):
                pl.Path(dirname).mkdir(parents=True, exist_ok=True)
            shu.copy(local_file, destination)
            assert os.path.isfile(destination)
        pl.Path(local_file).unlink()

    def process_extra(self, events):
        raise NotImplementedError

    def process(self, events):

        # data or monte carlo?
        data_kind = "mc" if "GenPart" in ak.fields(events) else "data"

        # met filters
        met_filters = self.meta["flashggMetFilters"][data_kind]
        filtered = ft.reduce(
            op.and_,
            (events.Flag[metfilter.split("_")[-1]] for metfilter in met_filters),
        )

        triggered = ak.ones_like(filtered)
        if self.apply_trigger:
            triggers = self.meta["TriggerPaths"][self.trigger_group][self.analysis]
            triggered = ft.reduce(
                op.or_, (events.HLT[trigger[4:-1]] for trigger in triggers)
            )

        # apply met filters and triggers to data
        events = events[filtered & triggered]

        # modifications to photons

        # photon preselection
        photons = self.photon_preselection(events.Photon)
        # sort photons in each event descending in pt
        # make descending-pt combinations of photons
        photons = photons[ak.argsort(photons.pt, ascending=False)]
        diphotons = ak.combinations(photons, 2, fields=["pho_lead", "pho_sublead"])
        # the remaining cut is to select the leading photons
        # the previous sort assures the order
        diphotons = diphotons[diphotons["pho_lead"].pt > self.min_pt_lead_photon]

        # now turn the diphotons into candidates with four momenta and such
        diphoton_4mom = diphotons["pho_lead"] + diphotons["pho_sublead"]
        diphotons["pt"] = diphoton_4mom.pt
        diphotons["eta"] = diphoton_4mom.eta
        diphotons["phi"] = diphoton_4mom.phi
        diphotons["mass"] = diphoton_4mom.mass
        diphotons["charge"] = diphoton_4mom.charge
        diphotons = ak.with_name(diphotons, "PtEtaPhiMCandidate")

        # sort diphotons by pT
        diphotons = diphotons[ak.argsort(diphotons.pt, ascending=False)]

        # baseline modifications to diphotons
        diphotons = self.add_diphoton_mva(diphotons, events)

        # set diphotons as part of the event record
        events["diphotons"] = diphotons

        # workflow specific processing
        events = self.process_extra(events)

        # run taggers on the events list with added diphotons
        # the shape here is ensured to be broadcastable
        for tagger in self.taggers:
            diphotons["_".join([tagger.name, str(tagger.priority)])] = tagger(events)

        # if there are taggers to run, arbitrate by them first
        if len(self.taggers):
            counts = ak.num(diphotons.pt, axis=1)
            flat_tags = np.stack(
                (
                    ak.flatten(diphotons["_".join([tagger.name, str(tagger.priority)])])
                    for tagger in self.taggers
                ),
                axis=1,
            )
            tags = ak.from_regular(ak.unflatten(flat_tags, counts), axis=2)
            winner = ak.min(tags[tags != 0], axis=2)
            diphotons["best_tag"] = winner

            # lowest priority is most important (ascending sort)
            # leave in order of diphoton pT in case of ties (stable sort)
            sorted = ak.argsort(diphotons.best_tag, stable=True)
            diphotons = diphotons[sorted]

        diphotons = ak.firsts(diphotons)

        # annotate diphotons with event information
        diphotons["event"] = events.event
        diphotons["lumi"] = events.luminosityBlock
        diphotons["run"] = events.run

        # drop events without a preselected diphoton candidate
        # drop events without a tag, if there are tags
        if len(self.taggers):
            diphotons = diphotons[
                ~(ak.is_none(diphotons) | ak.is_none(diphotons.best_tag))
            ]
        else:
            diphotons = diphotons[~ak.is_none(diphotons)]

        if self.output_location is not None:
            df = self.diphoton_list_to_pandas(diphotons)
            fname = (
                events.behavior["__events_factory__"]._partition_key.replace("/", "_")
                + ".parquet"
            )
            subdirs = []
            if "dataset" in events.metadata:
                subdirs.append(events.metadata["dataset"])
            self.dump_pandas(df, fname, self.output_location, subdirs)

        return {}

    def postprocess(self, accumulant):
        raise NotImplementedError

    def add_diphoton_mva(self, diphotons, events):
        if self.diphoton_mva is None:
            return diphotons

        var_order = self.meta["flashggDiPhotonMVA"]["inputs"]

        bdt_vars = {}

        bdt_vars["dipho_leadIDMVA"] = diphotons.pho_lead.mvaID
        bdt_vars["dipho_subleadIDMVA"] = diphotons.pho_sublead.mvaID
        bdt_vars["dipho_leadEta"] = diphotons.pho_lead.eta
        bdt_vars["dipho_subleadEta"] = diphotons.pho_sublead.eta
        bdt_vars["dipho_lead_ptoM"] = diphotons.pho_lead.pt / diphotons.mass
        bdt_vars["dipho_sublead_ptoM"] = diphotons.pho_sublead.pt / diphotons.mass

        def calc_displacement(photons, events):
            x = photons.x_calo - events.PV.x
            y = photons.y_calo - events.PV.y
            z = photons.z_calo - events.PV.z
            return ak.zip({"x": x, "y": y, "z": z}, with_name="Vector3D")

        v_lead = calc_displacement(diphotons.pho_lead, events)
        v_sublead = calc_displacement(diphotons.pho_sublead, events)

        p_lead = v_lead.unit() * diphotons.pho_lead.energyRaw
        p_lead["energy"] = diphotons.pho_lead.energyRaw
        p_lead = ak.with_name(p_lead, "Momentum4D")
        p_sublead = v_sublead.unit() * diphotons.pho_sublead.energyRaw
        p_sublead["energy"] = diphotons.pho_sublead.energyRaw
        p_sublead = ak.with_name(p_sublead, "Momentum4D")

        sech_lead = 1.0 / np.cosh(p_lead.eta)
        sech_sublead = 1.0 / np.cosh(p_sublead.eta)
        tanh_lead = np.cos(p_lead.theta)
        tanh_sublead = np.cos(p_sublead.theta)

        cos_dphi = np.cos(p_lead.deltaphi(p_sublead))

        numerator_lead = sech_lead * (
            sech_lead * tanh_sublead - tanh_lead * sech_sublead * cos_dphi
        )
        numerator_sublead = sech_sublead * (
            sech_sublead * tanh_lead - tanh_sublead * sech_lead * cos_dphi
        )

        denominator = (
            1.0 - tanh_lead * tanh_sublead - sech_lead * sech_sublead * cos_dphi
        )

        add_reso = (
            0.5
            * (-np.sqrt(2.0) * events.BeamSpot.sigmaZ / denominator)
            * (numerator_lead / p_lead.mag + numerator_sublead / p_sublead.mag)
        )

        dEnorm_lead = diphotons.pho_lead.energyErr / diphotons.pho_lead.energy
        dEnorm_sublead = diphotons.pho_sublead.energyErr / diphotons.pho_sublead.energy

        sigma_m = 0.5 * np.sqrt(dEnorm_lead ** 2 + dEnorm_sublead ** 2)
        sigma_wv = np.sqrt(add_reso ** 2 + sigma_m ** 2)

        vtx_prob = ak.full_like(sigma_m, 0.999)  # !!!! placeholder !!!!

        bdt_vars["CosPhi"] = cos_dphi
        bdt_vars["vtxprob"] = vtx_prob
        bdt_vars["sigmarv"] = sigma_m
        bdt_vars["sigmawv"] = sigma_wv

        counts = ak.num(diphotons, axis=-1)
        bdt_inputs = np.column_stack(
            [ak.to_numpy(ak.flatten(bdt_vars[name])) for name in var_order]
        )
        tempmatrix = xg.DMatrix(bdt_inputs, feature_names=var_order)
        scores = self.diphoton_mva.predict(tempmatrix)

        for var, arr in bdt_vars.items():
            if "dipho" not in var:
                diphotons[var] = arr

        diphotons["bdt_score"] = ak.unflatten(scores, counts)

        return diphotons
