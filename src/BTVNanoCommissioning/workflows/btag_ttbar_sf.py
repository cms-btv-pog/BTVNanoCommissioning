# https://git.rwth-aachen.de/3pia/cms_analyses/jet_tagging_sf/-/blob/master/recipes/selection3.py?ref_type=heads

from functools import reduce
from operator import and_, or_

# import collections
import awkward as ak
import numpy as np

# import os
# import uproot
from coffea import processor

# user helper function
from BTVNanoCommissioning.helpers.func import (
    # PFCand_link,
    # uproot_writeable,
    dump_lumi,
    update,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.array_writer import array_writer

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    common_shifts,
    load_lumi,
    load_SF,
    weight_manager,
)

## load histograms & selctions for this workflow
from BTVNanoCommissioning.utils.histogrammer import histo_writter, histogrammer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    ele_mvatightid,
    jet_id,
    mu_idiso,
)


def make_p4(obj):
    return ak.zip(
        {
            "pt": obj.pt,
            "eta": obj.eta,
            "phi": obj.phi,
            "mass": obj.mass,
        },
        with_name="PtEtaPhiMCandidate",
    )


def min_dr(particles):
    di_particles = ak.combinations(
        particles,
        n=2,
        replacement=False,
        axis=1,
        fields=["p0", "p1"],
    )
    return ak.min(
        make_p4(di_particles.p0).delta_r(make_p4(di_particles.p1)),
        axis=-1,
        mask_identity=False,
    )


def reduce_and(*conditions):
    return reduce(and_, conditions)


def reduce_or(*conditions):
    return reduce(or_, conditions)


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
        selectionModifier=None,
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        self.channel = selectionModifier
        self._workflow = "btag_ttbar_sf"
        if self.channel in ("mumu", "ee", "emu"):
            self._workflow += f"_{self.channel}"
        ## Load corrections
        self.SF_map = load_SF(year=self._year, campaign=self._campaign)

        if self._year == "2024":
            ParT_name = "btagUParTAK4B"
        else:
            ParT_name = "btagRobustParTAK4B"

        self.b_tagger_config = {
            "btagDeepFlavB": {
                "loose": 0.0614,
                "medium": 0.3196,
                "tight": 0.73,
            },
            "btagPNetB": {
                "loose": 0.1,
                "medium": 0.5,
                "tight": 0.9,
            },
            ParT_name: {
                "loose": 0.0897,
                "medium": 0.451,
                "tight": 0.8604,
            },
        }

    @property
    def accumulator(self):
        return self._accumulator

    ## Apply corrections on momentum/mass on MET, Jet, Muon
    def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, shift_name):
        """Selection following
        https://git.rwth-aachen.de/3pia/cms_analyses/jet_tagging_sf/-/blob/master/recipes/selection3.py?ref_type=heads
        """
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        # print(f"processing dataset {dataset} ({isRealData=}) with {len(events)=}")
        ######################
        #  Create histogram  #
        ######################
        output = (
            {}
            if self.noHist
            else histogrammer(
                events,
                workflow=self._workflow,
                year=self._year,
                campaign=self._campaign,
            )
        )
        if shift_name is None:
            if isRealData:
                output["sumw"] = len(events)
            else:
                output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        ## HLT
        channels = ("ee", "emu", "mumu")
        triggers = {
            "ee": [
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            ],
            "emu": [
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            ],
            "mumu": [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            ],
        }

        req_trig = {
            channel: HLT_helper(events, triggers[channel]) for channel in channels
        }

        # load the objects
        muons = events.Muon
        electrons = events.Electron
        jets = events.Jet
        met = events.MET

        ##### Add some selections
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2

        muon_sel = (muons.pt > 15) & (mu_idiso(events, self._campaign))
        good_muons = muons[muon_sel]
        # Electron cut
        ele_sel = (electrons.pt > 15) & (ele_mvatightid(events, self._campaign))
        good_electrons = electrons[ele_sel]

        good_leptons = ak.with_name(
            ak.concatenate([good_muons, good_electrons], axis=1), "PtEtaPhiMCandidate"
        )

        # lepton cuts
        lep_tpls = ak.combinations(
            good_leptons[..., :2], n=2, replacement=False, axis=-1, fields=["l0", "l1"]
        )
        dilep = ak.firsts(lep_tpls.l0 + lep_tpls.l1)
        dilep_mass = ak.fill_none(dilep.mass, np.nan)
        zll_mass_min = dilep_mass >= 12.0
        zll_mass_cut = abs(dilep_mass - 91.1876) <= 10.0
        zll_mass_veto = abs(dilep_mass - 91.1876) >= 10.0
        z_pt_cut = ak.fill_none(dilep.pt, np.nan) > 10.0
        ll_opp_charge = ak.fill_none(dilep.charge, np.nan) == 0

        dr_ll = min_dr(good_leptons[..., :2])
        dr_ll_cut = dr_ll > 0.2

        # leading lepton pt cuts
        leading_muon_pt = ak.fill_none(ak.firsts(good_muons.pt), np.nan)
        leading_electron_pt = ak.fill_none(ak.firsts(good_electrons.pt), np.nan)
        leading_muon_pt_cut = leading_muon_pt > 25.0
        leading_electron_pt_cut = leading_electron_pt > 25.0

        # configure triggers
        trigger_ee = reduce_and(req_trig["ee"], leading_electron_pt_cut)
        trigger_emu = reduce_and(
            req_trig["emu"], reduce_or(leading_electron_pt_cut, leading_muon_pt_cut)
        )
        trigger_mumu = reduce_and(req_trig["mumu"], leading_muon_pt_cut)

        # configure channels
        ch_mumu = reduce_and(
            trigger_mumu,
            ak.num(good_muons) == 2,
            ak.num(good_electrons) == 0,
            zll_mass_min,
            ll_opp_charge,
            dr_ll_cut,
        )
        ch_ee = reduce_and(
            trigger_ee,
            ak.num(good_electrons) == 2,
            ak.num(good_muons) == 0,
            zll_mass_min,
            ll_opp_charge,
            dr_ll_cut,
        )
        ch_emu = reduce_and(
            trigger_emu,
            ak.num(good_electrons) == 1,
            ak.num(good_muons) == 1,
            zll_mass_min,
            ll_opp_charge,
            dr_ll_cut,
        )
        assert not ak.any(reduce_and(ch_mumu, ch_ee, ch_emu))
        ch_incl = reduce_or(ch_mumu, ch_ee, ch_emu)

        ## Jet cuts
        jet_sel = jet_id(events, self._campaign)
        good_jets = jets[jet_sel]
        # TODO: clean jets, leptons
        clean_good_jets = good_jets
        # TODO: jet veto map?

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(ak.local_index(events.Jet.pt), jet_sel)
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]
        ## Other cuts
        met_pt_cut = met.pt > 30

        ## Apply all selections
        n_jets = ak.num(clean_good_jets)
        two_jets = n_jets == 2
        # ge_two_jets = n_jets >= 2

        mht_x = ak.sum(clean_good_jets.x, axis=-1) + ak.sum(good_leptons.x, axis=-1)
        mht_y = ak.sum(clean_good_jets.y, axis=-1) + ak.sum(good_leptons.y, axis=-1)
        mht = np.sqrt(mht_x**2 + mht_y**2)
        z_diamond = (
            (dilep_mass < (65.5 + 3 * mht / 8))
            | (dilep_mass > (108.0 - mht / 4))
            | (dilep_mass < (79.0 - 3 * mht / 4))
            | (dilep_mass > (99.0 + mht / 2))
        )

        lepton_region_cut_HF = reduce_and(met_pt_cut, zll_mass_veto)
        lepton_region_cut_LF = reduce_and(
            ~met_pt_cut, zll_mass_cut, z_pt_cut, ~z_diamond
        )

        # combine all cuts
        # TODO: met filter
        common_cuts = [req_lumi, two_jets]  # met_filter

        if self.channel == "mumu":
            event_level = reduce_and(*common_cuts, ch_mumu)
        elif self.channel == "ee":
            event_level = reduce_and(*common_cuts, ch_ee)
        elif self.channel == "emu":
            event_level = reduce_and(*common_cuts, ch_emu)
        else:
            # raise ValueError(f"Invalid channel: {self.channel}")
            event_level = reduce_and(*common_cuts, ch_incl)
        # Skip empty events
        if len(events[event_level]) == 0:
            if self.isArray:
                array_writer(
                    self,
                    events[event_level],
                    events,
                    "nominal",
                    dataset,
                    isRealData,
                )
            return {dataset: output}

        ####################
        # Selected objects # : Pruned objects with reduced event_level
        ####################
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = clean_good_jets[event_level][:, :2]
        # get number of jets distribution without cut on number of jets
        # does not work shape wise
        pruned_ev["njet"] = ak.count(clean_good_jets[event_level].pt, axis=1)

        if self.channel == "mumu":
            pruned_ev["SelMuon"] = good_muons[event_level][:, :2]
        elif self.channel == "ee":
            pruned_ev["SelElectron"] = good_electrons[event_level][:, :2]
        elif self.channel == "emu":
            pruned_ev["SelMuon"] = good_muons[event_level][:, 0]
            pruned_ev["SelElectron"] = good_electrons[event_level][:, 0]
        else:
            pruned_ev["SelMuon"] = good_muons[event_level][:, :2]
            pruned_ev["SelElectron"] = good_electrons[event_level][:, :2]

        # pruned_ev["ch_mumu"] = ch_mumu[event_level]
        # pruned_ev["ch_ee"] = ch_ee[event_level]
        # pruned_ev["ch_emu"] = ch_emu[event_level]
        # if "PFCands" in events.fields:
        #     pruned_ev.PFCands = PFCand_link(events, event_level, jetindx)

        ####################
        #  Tag and Probe   #
        ####################

        for tagger_name in self.b_tagger_config.keys():
            for i_jet in range(2):
                pruned_ev[f"{tagger_name}_region_HF_jet{i_jet}"] = np.zeros(
                    len(pruned_ev.njet), dtype=bool
                )
                pruned_ev[f"{tagger_name}_region_LF_jet{i_jet}"] = np.zeros(
                    len(pruned_ev.njet), dtype=bool
                )

        for i_tag_jet, i_probe_jet in [(0, 1), (1, 0)]:
            for tagger_name, tagger_config in self.b_tagger_config.items():
                b_score_btag = ak.fill_none(
                    ak.firsts(
                        pruned_ev["SelJet"][:, i_tag_jet : i_tag_jet + 1][tagger_name]
                    ),
                    np.nan,
                )

                tag_jet_cut = {
                    "HF": b_score_btag >= tagger_config["medium"],
                    "LF": b_score_btag <= tagger_config["loose"],
                }

                is_HF = reduce_and(tag_jet_cut["HF"], lepton_region_cut_HF[event_level])
                is_LF = reduce_and(tag_jet_cut["LF"], lepton_region_cut_LF[event_level])
                assert not ak.any(is_HF & is_LF)
                pruned_ev[f"{tagger_name}_region_HF_jet{i_probe_jet}"] = is_HF
                pruned_ev[f"{tagger_name}_region_LF_jet{i_probe_jet}"] = is_LF

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)

        # configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]

        if not isRealData:
            pruned_ev["weight"] = weights.weight()
            for ind_wei in weights.weightStatistics.keys():
                pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                    include=[ind_wei]
                )

        # Configure histograms
        if not self.noHist:
            output = histo_writter(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )

        # Output arrays
        if self.isArray:
            array_writer(
                self,
                pruned_ev,
                events,
                systematics[0],
                dataset,
                isRealData,
            )

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
