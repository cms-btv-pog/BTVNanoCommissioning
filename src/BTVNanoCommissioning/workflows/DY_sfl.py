import awkward as ak
import numpy as np
from coffea import processor

from BTVNanoCommissioning.helpers.func import PFCand_link, dump_lumi, update
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.correction import (
    common_shifts,
    load_lumi,
    load_SF,
    weight_manager,
)
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histo_writter,
    histogrammer,
)
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    MET_filters,
    btag_wp_dict,
    ele_mvatightid,
    jet_id,
    mu_idiso,
)


class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
        selectionModifier="DYM",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        self.selMod = selectionModifier
        # Load corrections
        self.SF_map = load_SF(self._year, self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        events = missing_branch(events)
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        isMu = False
        if "DYM" in self.selMod:
            triggers = ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"]
            isMu = True
        elif "DYE" in self.selMod:
            triggers = ["Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
        else:
            raise ValueError(self.selMod, "is not a valid selection modifier.")

        if self.noHist:
            output = {}
        else:
            output = histogrammer(
                jet_fields=events.Jet.fields,
                obj_list=["posl", "negl", "dilep", "jet0"],
                hist_collections=["common", "fourvec", "DY_sfl"],
                include_m=isMu,
                year=self._year,
                campaign=self._campaign,
            )
        if shift_name is None:
            output["sumw"] = len(events) if isRealData else ak.sum(events.genWeight)
        ####################
        #    Selections    #
        ####################

        # Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        # HLT
        req_trig = HLT_helper(events, triggers)
        req_metfilter = MET_filters(events, self._campaign)

        # Muon cuts
        dilep_mu = events.Muon[(events.Muon.pt > 20) & mu_idiso(events, self._campaign)]

        # Electron cuts
        dilep_ele = events.Electron[
            (events.Electron.pt > 25) & ele_mvatightid(events, self._campaign)
        ]
        if isMu:
            thisdilep = dilep_mu
            otherdilep = dilep_ele
        else:
            thisdilep = dilep_ele
            otherdilep = dilep_mu

        # dilepton
        pos_dilep = thisdilep[thisdilep.charge > 0]
        neg_dilep = thisdilep[thisdilep.charge < 0]
        req_pl = ak.count(pos_dilep.pt, axis=1) >= 1
        req_nl = ak.count(neg_dilep.pt, axis=1) >= 1
        req_dilep_chrg = ak.num(thisdilep.charge) >= 2
        req_otherdilep_chrg = ak.num(otherdilep.charge) == 0
        req_dilep = ak.fill_none(
            req_pl & req_nl & req_dilep_chrg & req_otherdilep_chrg,
            False,
            axis=-1,
        )

        pl_iso = ak.all(
            events.Jet.metric_table(pos_dilep) > 0.4, axis=2, mask_identity=True
        )
        nl_iso = ak.all(
            events.Jet.metric_table(neg_dilep) > 0.4, axis=2, mask_identity=True
        )
        jet_sel = ak.fill_none(
            jet_id(events, self._campaign) & pl_iso & nl_iso,
            False,
            axis=-1,
        )

        pos_dilep = ak.pad_none(pos_dilep, 1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep, 1, axis=1)

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81) & (dilep_mass.mass < 101) & (dilep_mass.pt > 15)
        )

        # Jet cuts
        pl_iso = ak.all(
            events.Jet.metric_table(pos_dilep[:, 0]) > 0.4, axis=2, mask_identity=True
        )
        nl_iso = ak.all(
            events.Jet.metric_table(neg_dilep[:, 0]) > 0.4, axis=2, mask_identity=True
        )
        event_jet = events.Jet[
            ak.fill_none(
                jet_id(events, self._campaign) & pl_iso & nl_iso,
                False,
                axis=-1,
            )
        ]
        req_jets = ak.count(event_jet.pt, axis=1) >= 1

        # store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            jet_sel == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        event_level = ak.fill_none(
            req_lumi & req_trig & req_dilep & req_dilepmass & req_jets & req_metfilter,
            False,
        )
        if len(events[event_level]) == 0:
            if self.isArray:
                array_writer(
                    self,
                    events[event_level],
                    events,
                    None,
                    ["nominal"],
                    dataset,
                    isRealData,
                    empty=True,
                )
            return {dataset: output}
        ####################
        # Selected objects #
        ####################
        sposmu = pos_dilep[event_level][:, 0]
        snegmu = neg_dilep[event_level][:, 0]
        sz = sposmu + snegmu

        sel_jet = event_jet[event_level][:, 0]

        sel_mu = ak.concatenate([sposmu, snegmu])
        smu = ak.zip(
            {
                b: ak.Array(np.reshape(sel_mu[b].to_numpy(), (len(sposmu[b]), 2)))
                for b in sposmu.fields
            }
        )
        # Jets passed tagger working point
        if not isRealData:
            flav = ak.where(
                (event_jet.hadronFlavour == 0) & (event_jet.partonFlavour == 0),
                event_jet.hadronFlavour + 1,
                event_jet.hadronFlavour,
            )
            flav = ak.where(event_jet.partonFlavour == 21, flav + 6, flav)

        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level]
        pruned_ev["AllSelJet"] = event_jet[event_level]  # untouched by histo_writter

        for tagger, tag_obj in btag_wp_dict[f"{self._year}_{self._campaign}"].items():
            for stringency, wp in tag_obj["b"].items():
                if stringency == "No":
                    continue
                mask_postag = event_jet[f"btag{tagger}B"] > wp
                mask_negtag = event_jet[f"btagNeg{tagger}B"] > wp
                postag_jet = event_jet[mask_postag][event_level]
                negtag_jet = event_jet[mask_negtag][event_level]
                key = f"{tagger}{stringency}"
                pruned_ev[f"{key}_postag_jet"] = postag_jet
                pruned_ev[f"{key}_negtag_jet"] = negtag_jet
                if isRealData:
                    pruned_ev[f"{key}_postag_jet", "flavor"] = ak.zeros_like(
                        postag_jet.pt,
                        dtype=int,
                    )
                    pruned_ev[f"{key}_negtag_jet", "flavor"] = ak.zeros_like(
                        negtag_jet.pt,
                        dtype=int,
                    )
                else:
                    pruned_ev[f"{key}_postag_jet", "flavor"] = flav[mask_postag][
                        event_level
                    ]
                    pruned_ev[f"{key}_negtag_jet", "flavor"] = flav[mask_negtag][
                        event_level
                    ]
        if isMu:
            pruned_ev["MuonPlus"] = sposmu
            pruned_ev["MuonMinus"] = snegmu
            pruned_ev["posl"] = sposmu
            pruned_ev["negl"] = snegmu
            pruned_ev["SelMuon"] = smu
        else:
            pruned_ev["ElectronPlus"] = sposmu
            pruned_ev["ElectronMinus"] = snegmu
            pruned_ev["posl"] = sposmu
            pruned_ev["negl"] = snegmu
            pruned_ev["SelElectron"] = smu
        pruned_ev["dilep"] = sz
        pruned_ev["dilep", "pt"] = pruned_ev.dilep.pt
        pruned_ev["dilep", "eta"] = pruned_ev.dilep.eta
        pruned_ev["dilep", "phi"] = pruned_ev.dilep.phi
        pruned_ev["dilep", "mass"] = pruned_ev.dilep.mass
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        pruned_ev["dr_mu1jet"] = sposmu.delta_r(sel_jet)
        pruned_ev["dr_mu2jet"] = snegmu.delta_r(sel_jet)
        if isRealData:
            pruned_ev["AllSelJet", "flavor"] = ak.zeros_like(
                event_jet[event_level].pt, dtype=int
            )
        else:
            pruned_ev["AllSelJet", "flavor"] = flav[event_level]

        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            pruned_ev["PFCands"] = PFCand_link(events, event_level, jetindx)

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        # Configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]
        # Configure histograms
        if not self.noHist:
            output = histo_writter(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )
        # Output arrays
        if self.isArray:
            array_writer(
                self, pruned_ev, events, weights, systematics, dataset, isRealData
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
