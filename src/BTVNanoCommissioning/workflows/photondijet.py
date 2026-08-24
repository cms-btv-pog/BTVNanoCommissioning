import awkward as ak
import os
import numpy as np
import correctionlib
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
    photondijet_writer,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    MET_filters,
)

# Quark-enriched working point on eta_gamma*eta_j1 + dR_gamma,j2. Both terms are
# smaller for quark probes (figures 9 and 10 of arXiv:1104.1175), so quark-like
# is the *low* side. The value is a starting point read off those figures; the
# discriminant is stored per event, so tighter points can be chosen offline.
QUARK_DISC_CUT = 2.0


class NanoProcessor(processor.ProcessorABC):
    """Quark-enriched photon+2jet selection.

    Follows Gallicchio & Schwartz, JHEP 10 (2011) 103 [arXiv:1104.1175]. Where
    the `QG_photonjet` workflow is their gamma+1jet sample - whose quark purity
    saturates around 88% because a 2-body final state leaves no kinematic
    handles - this is their gamma+2jet sample: the *softer* jet is the probe,
    and the single variable

        disc = eta_gamma * eta_j1 + dR(gamma, j2)

    recovers essentially the full 9-input BDT performance (their figure 11).
    Quark-like probes sit at low values: the photon is collinear with a quark
    probe (the q -> q gamma singularity, there is no g -> g gamma vertex), and
    the photon and the harder jet are widely separated in eta.
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
        selectionModifier="",
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

    ## Apply corrections on momentum/mass on MET, Jet, Muon
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

        # "softprobe": drop the paper's requirement that both jets share the
        # photon pt scale, keeping a softer (and much more abundant) probe.
        # "quark": cut on the discriminant for a quark-enriched sample.
        modifiers = [m for m in self.selectionModifier.split("_") if m]
        for m in modifiers:
            if m not in ["softprobe", "quark"]:
                raise ValueError(
                    self.selectionModifier, "is not a valid selection modifier."
                )
        common_ptscale = "softprobe" not in modifiers
        quark_enriched = "quark" in modifiers

        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        ## HLT - same photon menu as the gamma+1jet workflow
        if self._year == "2022" or self._year == "2023":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 50],
                "Photon50EB_TightID_TightIso": [50, 75],
                "Photon75EB_TightID_TightIso": [75, 90],
                "Photon90EB_TightID_TightIso": [90, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        elif self._year == "2024":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 50],
                "Photon50EB_TightID_TightIso": [50, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        elif self._year == "2025":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 40],
                "Photon40EB_TightID_TightIso": [40, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        else:
            raise ValueError(self._year, "is not a valid year.")

        # Sort the jets by pt (the corrections applied above can reorder them)
        events["Jet"] = events.Jet[ak.argsort(events.Jet.pt, axis=1, ascending=False)]

        req_metfilter = MET_filters(events, self._campaign)

        ##### Add some selections
        ## Jet cuts
        # |eta| < 2.5, the acceptance of the paper: the discriminant is built out
        # of rapidities, and the flavour of an HF jet is not usable as a probe.
        jet_sel = jet_id(events, self._campaign, max_eta=2.5, min_pt=20)

        if self._year == "2016":
            jet_puid = events.Jet.puId >= 1
        elif self._year in ["2017", "2018"]:
            jet_puid = events.Jet.puId >= 4
        else:
            jet_puid = ak.ones_like(jet_sel)

        jet_sel = jet_sel & jet_puid

        ## Photon cuts
        photon_sel = (
            (events.Photon.cutBased == 3)
            & (events.Photon.hoe < 0.02148)
            & (events.Photon.r9 > 0.94)
            & (events.Photon.r9 < 1.0)
            & (np.abs(events.Photon.eta) < 1.3)
        )

        # Index with the selection rather than ak.mask: masking leaves the failing
        # candidates in place as None, so slot 0 remains the leading *unselected*
        # candidate and the event is then dropped by None propagation instead of
        # falling back on the leading selected candidate.
        event_ph = ak.pad_none(events.Photon[photon_sel], 1)
        event_jet = ak.pad_none(events.Jet[jet_sel], 3)
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
                & (event_ph[:, 0].pt >= triggers[trg][0])
                & (event_ph[:, 0].pt < triggers[trg][1]),
                False,
            )

        req_trig = np.zeros(len(events), dtype="bool")
        # The pt scale of the event is set by the photon trigger bin it falls in;
        # the paper applies its pt cut to *all* jets at that same scale.
        ptmin_ev = np.zeros(len(events))
        for trg, passed in trig_pass.items():
            req_trig = req_trig | passed
            ptmin_ev = np.where(
                ak.to_numpy(passed) & (ptmin_ev == 0), triggers[trg][0], ptmin_ev
            )

        def _req(mask):
            # PackedSelection needs a plain boolean array: the padded object slots
            # make every cut option-typed, and a missing object must fail the cut.
            return ak.to_numpy(ak.fill_none(mask, False))

        photon, j1, j2, j3 = (
            event_ph[:, 0],
            event_jet[:, 0],
            event_jet[:, 1],
            event_jet[:, 2],
        )

        # eta_gamma*eta_j1 + dR(gamma, j2): the composite variable of the paper.
        # Note both terms use the *harder* jet for the eta product and the
        # *softer* (probe) jet for the distance to the photon.
        etaprod = photon.eta * j1.eta
        drgj2 = photon.delta_r(j2)
        disc = etaprod + drgj2

        # Build selections with PackedSelection for cleaner tracking and cutflow
        selections = PackedSelection()
        selections.add("lumi", _req(req_lumi))
        selections.add("metfilter", _req(req_metfilter))
        selections.add("trigger", _req(req_trig))
        selections.add("photon", _req(ak.count(event_ph.pt, axis=1) > 0))
        selections.add("jets", _req(njet > 1))
        # dR > 1.0 between the two jets and dR > 0.5 between the photon and each
        # of them (the paper's generation cuts): well separated objects, and no
        # photon/jet overlap for either the probe or the eta-product jet.
        selections.add("drjj", _req(j1.delta_r(j2) > 1.0))
        selections.add(
            "drgj", _req((photon.delta_r(j1) > 0.5) & (photon.delta_r(j2) > 0.5))
        )
        if common_ptscale:
            # pt cut on *both* jets at the scale of the photon pt bin
            selections.add("jetpt", _req((j1.pt >= ptmin_ev) & (j2.pt >= ptmin_ev)))
        else:
            selections.add("jetpt", _req(j2.pt >= 20))
        # Exclusive 2-jet topology: any further jet must be soft compared with
        # the two probing ones (same handle as the dijet workflow).
        selections.add(
            "excl2jet",
            _req(
                ak.where(
                    njet > 2,
                    j3.pt / (0.5 * (j1 + j2).pt) < 0.15,
                    ak.ones_like(req_trig, dtype=bool),
                )
            ),
        )
        if "GenVtx_z" in events.fields:
            selections.add("vtx", _req(np.abs(events.GenVtx_z - events.PV_z) < 0.2))
        else:
            selections.add("vtx", ak.ones_like(events.run, dtype=bool))
        if quark_enriched:
            selections.add("quark_disc", _req(disc < QUARK_DISC_CUT))

        cuts = [
            "lumi",
            "metfilter",
            "trigger",
            "photon",
            "jets",
            "drjj",
            "drgj",
            "jetpt",
            "excl2jet",
            "vtx",
        ]
        if quark_enriched:
            cuts.append("quark_disc")
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
                is_photondijet=True,
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

        # Built from the *selected* collections, so the stored objects are the
        # ones the cuts above were evaluated on.
        pruned_sel_jet = event_jet[event_level]
        pruned_ev["Tag"] = event_ph[event_level][:, 0]
        pruned_ev["Tag", "pt"] = pruned_ev["Tag"].pt
        pruned_ev["Tag", "eta"] = pruned_ev["Tag"].eta
        pruned_ev["Tag", "phi"] = pruned_ev["Tag"].phi
        pruned_ev["SelJet"] = pruned_sel_jet[:, :2]
        pruned_ev["LeadJet"] = pruned_sel_jet[:, 0]
        # j2, the softer of the two: the quark-enriched probe
        pruned_ev["SoftJet"] = pruned_sel_jet[:, 1]

        # Cross-object quantities have to be assigned explicitly on pruned_ev
        pruned_ev["photondijet_etaprod"] = etaprod[event_level]
        pruned_ev["photondijet_drgj2"] = drgj2[event_level]
        # eta_gamma*eta_j1 + dR(gamma, j2)
        pruned_ev["photondijet_disc"] = (
            pruned_ev["photondijet_etaprod"] + pruned_ev["photondijet_drgj2"]
        )
        pruned_ev["photondijet_mass"] = (
            pruned_sel_jet[:, 0] + pruned_sel_jet[:, 1] + pruned_ev["Tag"]
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
            output = photondijet_writer(
                pruned_ev, output, weights, systematics, self.isSyst
            )
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
