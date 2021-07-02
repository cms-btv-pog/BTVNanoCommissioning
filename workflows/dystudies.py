import coffea
from coffea import hist, processor
import numpy as np
import awkward as ak
import pandas as pd
import functools as ft
import operator as op
import os
import shutil as shu
import pathlib as pl

class DYStudiesProcessor(processor.ProcessorABC):
    def __init__(self, metaconditions, do_systematics=False, apply_trigger=False, output_location=None):
        self.meta = metaconditions
        self.do_systematics = do_systematics
        self.apply_trigger = apply_trigger
        self.output_location = output_location
        self.trigger_group = ".*DoubleEG.*"
        self.analysis = "mainAnalysis"

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

        self.prefixes = {"0": "lead", "1": "sublead"}
        
    def photon_preselection(self, photons):
        photon_abs_eta = np.abs(photons.eta)
        return photons[  (photons.pt > self.min_pt_photon)
                       & (photon_abs_eta < self.max_sc_eta)
                       & ((photon_abs_eta < self.gap_barrel_eta) | (photon_abs_eta > self.gap_endcap_eta))
                       & (photons.mvaID > self.min_mvaid)
                       & (photons.hoe < self.max_hovere)
                       & (  (photons.r9 > self.min_full5x5_r9)
                          | (photons.pfRelIso03_chg < self.max_chad_iso)
                          | (photons.pfRelIso03_chg/photons.pt < self.max_chad_rel_iso))]

    def diphoton_list_to_pandas(self, diphotons):
        output = pd.DataFrame()
        for field in ak.fields(diphotons):
            prefix = self.prefixes.get(field, "")
            if len(prefix) > 0:
                for subfield in ak.fields(diphotons[field]):
                    output[f"{prefix}_{subfield}"] = ak.to_numpy(diphotons[field][subfield])
            else:
                output[field] = ak.to_numpy(diphotons[field])
        return output

    def dump_pandas(self, pddf, fname, location, subdirs=[]):
        xrd_prefix = 'root://'
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
        local_file = os.path.join(".", fname)
        subdirs = "/".join(subdirs) if xrootd else os.path.sep.join(subdirs)
        destination = locations + subdirs + f"/{fname}" if xrootd else os.path.join(location, os.path.join(subdirs, fname))
        pddf.to_parquet(local_file)
        if xrootd:
            pfx_len = len(xrd_prefix)
            client = XRootD.client.FileSystem(location[:location[pfx_len:].find('/') + pfx_len])
            status = client.copy(local_file, destination)
            assert status[0].ok
        else:
            dirname = os.path.dirname(destination)
            if not os.path.exists(dirname):
                pl.Path(dirname).mkdir(parents=True, exist_ok=True)
            shu.copy(local_file, destination)
            assert os.path.isfile(destination)
        
            
    
    def process(self, events):

        # data or monte carlo?
        data_kind = "mc" if "GenPart" in ak.fields(events) else "data"

        # met filters
        met_filters = self.meta["flashggMetFilters"][data_kind]
        filtered = ft.reduce(op.and_, (events.Flag[metfilter.split("_")[-1]] for metfilter in met_filters))

        triggered = ak.ones_like(filtered)
        if self.apply_trigger:
            triggers = self.meta["TriggerPaths"][self.trigger_group][self.analysis]
            triggered = ft.reduce(op.or_, (events.HLT[trigger[4:-1]] for trigger in triggers))
        
        # apply met filters and triggers to data
        events = events[filtered & triggered]

        # photon preselection
        photons = self.photon_preselection(events.Photon)
        # sort photons in each event descending in pt
        # make descending-pt combinations of photons
        photons = photons[ak.argsort(photons.pt, ascending=False)]
        diphotons = ak.combinations(photons, 2)
        # the remaining cut is to select the leading photons
        # the previous sort assures the order
        diphotons = diphotons[diphotons["0"].pt > self.min_pt_lead_photon]

        # now turn the diphotons into candidates with four momenta and such
        diphoton_4mom = diphotons["0"] + diphotons["1"]
        diphotons["pt"] = diphoton_4mom.pt
        diphotons["eta"] = diphoton_4mom.eta
        diphotons["phi"] = diphoton_4mom.phi
        diphotons["mass"] = diphoton_4mom.mass
        diphotons = ak.with_name(diphotons, "PtEtaPhiMCandidate")
        
        # arbitrate diphotons
        diphotons = diphotons[ak.argsort(diphotons.pt, ascending=False)]
        diphotons = ak.firsts(diphotons)

        # annotate diphotons with event information
        diphotons["event"] = events.event
        diphotons["lumi"] = events.luminosityBlock
        diphotons["run"] = events.run
        
        # drop events without a preselected diphoton candidate
        diphotons = diphotons[~ak.is_none(diphotons)]

        if self.output_location is not None:
            df = self.diphoton_list_to_pandas(diphotons)
            fname = events.behavior["__events_factory__"]._partition_key.replace("/", "_") + ".parquet"
            subdirs = []
            if 'dataset' in events.metadata:
                subdirs.append(events.metadata["dataset"])
            self.dump_pandas(df, fname, self.output_location, subdirs)
        
        return {
            
        }

    def postprocess(self):
        pass
