import awkward as ak
import numpy as np
import correctionlib
import os
from coffea import processor
from coffea.analysis_tools import Weights, PackedSelection

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    common_shifts,
    weight_manager,
    reweighting,
)

# user helper function
from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    uproot_writeable,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch

## load histograms & selctions for this workflow
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histogrammer,
)
from BTVNanoCommissioning.utils.histogramming.histograms.qgtag import (
    qg_writer,
    trijet_writer,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    MET_filters,
)


class NanoProcessor(processor.ProcessorABC):
    """Gluon-enriched trijet selection.

    Follows Gallicchio & Schwartz, JHEP 10 (2011) 103 [arXiv:1104.1175]: the
    softest of three hard, well-separated jets is the gluon-enriched probe, and
    the single kinematic variable

        disc = |eta_j3| - |eta_j1 - eta_j2|

    reproduces almost all of the discrimination power of a full multivariate
    analysis (their figure 15). The gluon-like tail is at *negative* values, so
    the ``*_gluon`` selection modifiers cut ``disc < 0``.
    """

    def __init__(
        self,
        year="2022",
        campaign="Summer22",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
        selectionModifier="PFJet",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ## Load corrections
        self.SF_map = load_SF(self._year, self._campaign)
        self.selectionModifier = selectionModifier

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        events = missing_branch(events)
        sumws = reweighting(events, self.isSyst)
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), sumws, name)
            for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, sumws, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        ## HLT
        # The gluon-enrichment cut is applied on top of the same trigger scheme,
        # so strip the suffix before resolving the trigger menu.
        selection_mode = self.selectionModifier.replace("_gluon", "")
        gluon_enriched = self.selectionModifier.endswith("_gluon")

        if selection_mode == "ZB":
            triggers = {
                "ZeroBias": [15, 60],
            }
            ptmin, ptmax = 15, 80
        elif selection_mode == "PFJet":
            triggers = {
                "PFJet40": [60, 86],
                "PFJet60": [86, 110],
                "PFJet80": [110, 170],
                "PFJet140": [170, 236],
                "PFJet200": [236, 305],
                "PFJet260": [305, 375],
                "PFJet320": [375, 460],
                "PFJet400": [460, 575],
                "PFJet500": [575, 1e7],
            }
            ptmin, ptmax = 60, 1e7
        elif selection_mode == "DiPFJetAve":
            triggers = {
                "DiPFJetAve40": [40, 60],
                "DiPFJetAve60": [60, 80],
                "DiPFJetAve80": [80, 140],
                "DiPFJetAve140": [140, 200],
                "DiPFJetAve200": [200, 260],
                "DiPFJetAve260": [260, 320],
                "DiPFJetAve320": [320, 400],
                "DiPFJetAve400": [400, 500],
                "DiPFJetAve500": [500, 1e7],
            }
            ptmin, ptmax = 40, 1e7
        else:
            raise ValueError(
                self.selectionModifier, "is not a valid selection modifier."
            )

        # Sort the jets by pt (the corrections applied above can reorder them)
        events["Jet"] = events.Jet[ak.argsort(events.Jet.pt, axis=1, ascending=False)]

        req_metfilter = MET_filters(events, self._campaign)

        ##### Add some selections
        ## Jet cuts
        # |eta| < 2.5 for all three jets, matching the acceptance of the paper
        # (the discriminant is built out of the rapidities, so jets outside the
        # tracker acceptance are not usable as probes anyway).
        jet_sel = jet_id(events, self._campaign, max_eta=2.5, min_pt=20)

        if self._year == "2016":
            jet_puid = events.Jet.puId >= 1
        elif self._year in ["2017", "2018"]:
            jet_puid = events.Jet.puId >= 4
        else:
            jet_puid = ak.ones_like(jet_sel)

        jet_sel = jet_sel & jet_puid

        # Index with the selection rather than ak.mask: masking leaves the failing
        # jets in place as None, so slots 0/1/2 kept pointing at the leading
        # *unselected* jets and the event was then dropped by None propagation
        # instead of falling back on the leading selected jets.
        event_jet = ak.pad_none(events.Jet[jet_sel], 4)
        njet = ak.count(event_jet.pt, axis=1)

        # Validate the paths against the sample (raises if none of them exist).
        # The returned OR is unused: each path is pt-binned individually below.
        HLT_helper(events, list(triggers.keys()))

        # NB: `events.HLT[trg] = ...` writes into a temporary copy of the HLT
        # record and is silently discarded, leaving the pt binning with no effect.
        # Keep the binned decisions in a local dict instead.
        trig_pass = {}
        for trg in triggers:
            if not hasattr(events.HLT, trg):
                continue
            trig_pass[trg] = ak.fill_none(
                events.HLT[trg]
                & (event_jet[:, 0].pt >= triggers[trg][0])
                & (event_jet[:, 0].pt < triggers[trg][1]),
                False,
            )

        req_trig = np.zeros(len(events), dtype="bool")
        for trg_pass in trig_pass.values():
            req_trig = req_trig | trg_pass

        def _req(mask):
            # PackedSelection needs a plain boolean array: the padded jet slots
            # make every cut option-typed, and a missing jet must fail the cut.
            return ak.to_numpy(ak.fill_none(mask, False))

        j1, j2, j3, j4 = (
            event_jet[:, 0],
            event_jet[:, 1],
            event_jet[:, 2],
            event_jet[:, 3],
        )

        # The discriminant of arXiv:1104.1175: |eta_j3| - |eta_j1 - eta_j2|
        disc = np.abs(j3.eta) - np.abs(j1.eta - j2.eta)

        # Build selections with PackedSelection for cleaner tracking and cutflow
        selections = PackedSelection()
        selections.add("lumi", _req(req_lumi))
        selections.add("metfilter", _req(req_metfilter))
        selections.add("trigger", _req(req_trig))
        selections.add("jets", _req(njet > 2))
        # pT cut on *all three* jets, i.e. the softest jet sets the pt bin
        selections.add("softjet", _req(j3.pt >= ptmin))
        selections.add("leadjet", _req(j1.pt < ptmax))
        # dR > 1.0 between any two of the three jets (paper's generation cut):
        # keeps the jets well separated so the quark/gluon assignment is
        # meaningful and the discriminant is not smeared by overlapping jets.
        selections.add(
            "drjj",
            _req(
                (j1.delta_r(j2) > 1.0) & (j1.delta_r(j3) > 1.0) & (j2.delta_r(j3) > 1.0)
            ),
        )
        # Exclusive 3-jet topology: any further jet must be soft compared with
        # the two leading ones (same handle as the dijet workflow).
        selections.add(
            "excl3jet",
            _req(
                ak.where(
                    njet > 3,
                    j4.pt / (0.5 * (j1 + j2).pt) < 0.15,
                    ak.ones_like(req_trig, dtype=bool),
                )
            ),
        )
        if "GenVtx_z" in events.fields:
            selections.add("vtx", _req(np.abs(events.GenVtx_z - events.PV_z) < 0.2))
        else:
            selections.add("vtx", ak.ones_like(events.run, dtype=bool))
        # Gluon-enriched working point: the gluon tail of the discriminant sits
        # at negative values (figure 14 of the paper).
        if gluon_enriched:
            selections.add("gluon_disc", _req(disc < 0.0))

        cuts = [
            "lumi",
            "metfilter",
            "trigger",
            "jets",
            "softjet",
            "leadjet",
            "drjj",
            "excl3jet",
            "vtx",
        ]
        if gluon_enriched:
            cuts.append("gluon_disc")
        event_level = selections.all(*cuts)

        ##<==== finish selection

        ######################
        #  Create histogram  # : Get the histogram dict from `histogrammer`
        ######################

        output = {}

        if not self.noHist:
            output = histogrammer(
                jet_fields=events.Jet.fields,
                obj_list=[],
                hist_collections=["qgtag"],
                axes_collections=["qgtag"],
                is_trijet=True,
            )

        if shift_name is None:
            output["sumw"] = sumws["sumw"]
            if not isRealData and self.isSyst:
                if "LHEPdfWeight" in events.fields:
                    output["PDF_sumwUp"] = sumws["PDF_sumwUp"]
                    output["PDF_sumwDown"] = sumws["PDF_sumwDown"]
                    output["aS_sumwUp"] = sumws["aS_sumwUp"]
                    output["aS_sumwDown"] = sumws["aS_sumwDown"]
                    output["PDFaS_sumwUp"] = sumws["PDFaS_sumwUp"]
                    output["PDFaS_sumwDown"] = sumws["PDFaS_sumwDown"]
                if "LHEScaleWeight" in events.fields:
                    output["muR_sumwUp"] = sumws["muR_sumwUp"]
                    output["muR_sumwDown"] = sumws["muR_sumwDown"]
                    output["muF_sumwUp"] = sumws["muF_sumwUp"]
                    output["muF_sumwDown"] = sumws["muF_sumwDown"]
                if "PSWeight" in events.fields:
                    if len(events.PSWeight[0]) == 4:
                        output["ISR_sumwUp"] = sumws["ISR_sumwUp"]
                        output["ISR_sumwDown"] = sumws["ISR_sumwDown"]
                        output["FSR_sumwUp"] = sumws["FSR_sumwUp"]
                        output["FSR_sumwDown"] = sumws["FSR_sumwDown"]

        event_level = ak.fill_none(event_level, False)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        # Skip empty events -
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

        ##===>  Ntuplization  : store custom information
        ####################
        # Selected objects # : Pruned objects with reduced event_level
        ####################
        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]

        # Built from the *selected* jets, so the stored objects are the ones the
        # cuts above were evaluated on.
        pruned_sel_jet = event_jet[event_level]
        pruned_ev["SelJet"] = pruned_sel_jet[:, :3]
        # j3, the softest of the three: the gluon-enriched probe jet
        pruned_ev["SoftJet"] = pruned_sel_jet[:, 2]
        pruned_ev["LeadJet"] = pruned_sel_jet[:, 0]
        pruned_ev["SubleadJet"] = pruned_sel_jet[:, 1]

        # Cross-object quantities have to be assigned explicitly on pruned_ev
        pruned_ev["trijet_deta12"] = np.abs(
            pruned_sel_jet[:, 0].eta - pruned_sel_jet[:, 1].eta
        )
        pruned_ev["trijet_abseta3"] = np.abs(pruned_sel_jet[:, 2].eta)
        # |eta_j3| - |eta_j1 - eta_j2| : the composite discriminant of the paper
        pruned_ev["trijet_disc"] = (
            pruned_ev["trijet_abseta3"] - pruned_ev["trijet_deta12"]
        )
        pruned_ev["trijet_mass"] = (
            pruned_sel_jet[:, 0] + pruned_sel_jet[:, 1] + pruned_sel_jet[:, 2]
        ).mass
        pruned_ev["njet"] = ak.count(pruned_sel_jet.pt, axis=1)

        ## <========= end: store custom objects

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        if isRealData:
            if self._year == "2022":
                run_num = "355374_362760"
            elif self._year == "2023":
                run_num = "366727_370790"
            elif self._year == "2024":
                run_num = "378985_386951"
            elif self._year == "2025":
                run_num = "391658_398860"
            else:
                raise NotImplementedError(
                    f"Prescale weights not available for data in {self._year}."
                )

            pruned_ev["psweight"] = np.zeros(len(pruned_ev))
            for trigger in trig_pass:
                # Check if the prescale weight file exists for the given trigger and year
                psfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{trigger}_year{self._year}.json"
                if not os.path.isfile(psfile):
                    psfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{trigger}_run{run_num}.json"
                    if not os.path.isfile(psfile):
                        raise NotImplementedError(
                            f"Prescale weights not available for {trigger} in {self._year}. Please run `scripts/dump_prescale.py`."
                        )
                pseval = correctionlib.CorrectionSet.from_file(psfile)
                thispsweight = pseval["prescaleWeight"].evaluate(
                    pruned_ev.run,
                    f"HLT_{trigger}",
                    ak.values_astype(pruned_ev.luminosityBlock, np.float32),
                )
                # Use the pt-binned decision: the raw HLT bits overlap, so the
                # lowest (most prescaled) path would otherwise claim every event.
                pruned_ev["psweight"] = ak.where(
                    (trig_pass[trigger][event_level]) & (pruned_ev["psweight"] == 0),
                    thispsweight,
                    pruned_ev["psweight"],
                )
            weights.add("psweight", pruned_ev["psweight"])

        # Configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]

        # Configure histograms
        if not self.noHist:
            output = qg_writer(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )
            output = trijet_writer(pruned_ev, output, weights, systematics, self.isSyst)
        # Output arrays
        if self.isArray:
            othersData = [
                "SV_*",
                "PV_npvs",
                "PV_npvsGood",
                "Rho_*",
                "run",
                "luminosityBlock",
            ]
            for trigger in trig_pass:
                othersData.append(f"HLT_{trigger}")
            array_writer(
                self,
                pruned_ev,
                events,
                weights,
                systematics,
                dataset,
                isRealData,
                othersData=othersData,
            )

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
