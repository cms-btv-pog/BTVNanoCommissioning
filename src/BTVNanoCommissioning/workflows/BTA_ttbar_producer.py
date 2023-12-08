import collections
import awkward as ak
import numpy as np
import uproot
from coffea import processor
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.helpers.BTA_helper import (
    BTA_ttbar_HLT_chns,
    to_bitwise_trigger,
)
from BTVNanoCommissioning.helpers.func import update
from BTVNanoCommissioning.utils.correction import load_SF, JME_shifts
import os


## Based on coffea_array_producer.ipynb from Congqiao
class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=True,
        noHist=False,
        chunksize=75000,
    ):
        self._year = year
        self._campaign = campaign
        self.chunksize = chunksize

        self.name = name
        self.SF_map = load_SF(self._campaign)

        ### Custom initialzations for BTA_ttbar workflow ###

        # note: in BTA TTbarSelectionProducer it says "disable trigger selection in MC"
        # for consistency, we will disable both trigger selection in MC and data here
        self.do_trig_sel = False

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        events = missing_branch(events)
        shifts = []

        if "JME" in self.SF_map.keys():
            shifts = JME_shifts(
                shifts, self.SF_map, events, self._campaign, isRealData, False, True
            )
        else:
            shifts = [
                ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
            ]

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        events = missing_branch(events)

        # basic variables
        basic_vars = {
            "Run": events.run,
            "Evt": events.event,
            "LumiBlock": events.luminosityBlock,
            "rho": events.fixedGridRhoFastjetAll,
            "npvs": events.PV.npvs,
            "npvsGood": events.PV.npvsGood,
        }
        if not isRealData:
            basic_vars["nPU"] = events.Pileup.nPU
            basic_vars["nPUtrue"] = events.Pileup.nTrueInt
            basic_vars["pthat"] = events.Generator.binvar

        ###########
        #   HLT   #
        ###########
        triggers = [t[0] for t in BTA_ttbar_HLT_chns]

        checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
        if ak.all(checkHLT == False):
            raise ValueError("HLT paths:", triggers, " are all invalid in", dataset)
        elif ak.any(checkHLT == False):
            print(
                np.array(triggers)[~checkHLT],
                " not exist in",
                dataset,
                ", please check!",
            )

        # pass_trigger dim: (n_evt, n_trig), append False for non-exist triggers
        pass_trig = np.array(
            [
                getattr(events.HLT, _trig, np.zeros(len(events), dtype="bool"))
                for _trig in triggers
            ]
        ).T

        # get trigger pass for each channel (e, m, em)
        pass_trig_chn = {
            chn: np.zeros(len(events), dtype="bool")
            for chn in ["e", "m", "ee", "mm", "em"]
        }
        for i, (trig, chn) in enumerate(
            BTA_ttbar_HLT_chns
        ):  # loop over e, m, em chanenl
            pass_trig_chn[chn] |= pass_trig[:, i]

        # convert to bitwise trigger (total number of triggers is less than 32)
        basic_vars[f"ttbar_trigWord"] = to_bitwise_trigger(
            pass_trig, ak.ArrayBuilder()
        ).snapshot()[:, 0]

        #########################
        #  Lep/di-lep channels  #
        #########################
        # select offline electrons and muons (BTA TTbarSelectionProducer definition)
        eles = events.Electron
        eles = eles[
            (eles.pt > 20)
            & (abs(eles.deltaEtaSC + eles.eta) < 2.4)
            & eles.convVeto  # pass conversion veto
            & (eles.cutBased >= 3)  # pass cut-based medium ID
        ]

        muons = events.Muon
        muons = muons[
            (muons.pt > 20)
            & (abs(muons.eta) < 2.4)
            & muons.tightId  # pass cut-based tight ID
            & (muons.pfRelIso04_all < 0.12)  # muon isolation cut
        ]

        # assign channels: 13, 11, 13*13, 13*11, 11*11
        # according to TTbarSelectionProducer::AssignChannel

        # pairing eletrons and muons
        pair_em = ak.cartesian([eles, muons])  # dim: (event, N[ele+mu pair])
        sumPt_em = ak.fill_none(ak.max(pair_em["0"].pt + pair_em["1"].pt, axis=1), 0)
        # the index to the max pT of (ele, mu) pair
        idx_em = ak.singletons(
            ak.argmax(pair_em["0"].pt + pair_em["1"].pt, axis=1)
        )  # dim: (event, 1 or 0)

        # pairing electrons
        pair_ee = ak.combinations(eles, 2)
        sumPt_ee = ak.fill_none(ak.max(pair_ee["0"].pt + pair_ee["1"].pt, axis=1), 0)
        idx_ee = ak.singletons(ak.argmax(pair_ee["0"].pt + pair_ee["1"].pt, axis=1))

        # pairing muons
        pair_mm = ak.combinations(muons, 2)
        sumPt_mm = ak.fill_none(ak.max(pair_mm["0"].pt + pair_mm["1"].pt, axis=1), 0)
        idx_mm = ak.singletons(ak.argmax(pair_mm["0"].pt + pair_mm["1"].pt, axis=1))

        # criteria for each channel (m, e, em, mm, ee), according to TTbarSelectionProducer::AssignChannel
        n_eles, n_muons = ak.num(eles), ak.num(muons)
        criteria = {
            "m": (n_eles == 0) & (n_muons == 1),
            "e": (n_eles == 1) & (n_muons == 0),
            "em": (sumPt_em > sumPt_mm) & (sumPt_em > sumPt_ee),
            "mm": (sumPt_mm > sumPt_em) & (sumPt_mm > sumPt_ee),
            "ee": (sumPt_ee > sumPt_em) & (sumPt_ee > sumPt_mm),
        }
        # if required, should also pass corresponding trigger for each channel (doTrigSel=True in TTbarSelectionProducer_cfi, but it was again set to False in the code)
        # note that here we set default do_trig_sel to False
        if self.do_trig_sel:
            for chn in criteria.keys():
                criteria[chn] &= pass_trig_chn[chn]

        # finally, assign the channel of each event, in the order of m, e, em, mm, ee (according to BTA TTbarSelectionProducer)
        # - first, assign channel m or e, if only one lepton is presented
        # - then, assign di-lep channels, based on the larger sumPt
        zeros = ak.zeros_like(events.run, dtype=int)
        chsel = zeros
        chsel = ak.where(
            criteria["m"] & (chsel == 0),
            ak.fill_none((-13) * ak.firsts(muons).charge, 0),
            chsel,
        )
        chsel = ak.where(
            criteria["e"] & (chsel == 0),
            ak.fill_none((-11) * ak.firsts(eles).charge, 0),
            chsel,
        )
        chsel = ak.where(
            criteria["em"] & (chsel == 0),
            ak.fill_none(
                (-13)
                * (-11)
                * ak.firsts(pair_em[idx_em])["0"].charge
                * ak.firsts(pair_em[idx_em])["1"].charge,
                0,
            ),
            chsel,
        )
        chsel = ak.where(
            criteria["mm"] & (chsel == 0),
            ak.fill_none(
                (-13)
                * (-13)
                * ak.firsts(pair_mm[idx_mm])["0"].charge
                * ak.firsts(pair_mm[idx_mm])["1"].charge,
                0,
            ),
            chsel,
        )
        chsel = ak.where(
            criteria["ee"] & (chsel == 0),
            ak.fill_none(
                (-11)
                * (-11)
                * ak.firsts(pair_ee[idx_ee])["0"].charge
                * ak.firsts(pair_ee[idx_ee])["1"].charge,
                0,
            ),
            chsel,
        )
        basic_vars["ttbar_chan"] = chsel

        #######################
        #  Electrons & muons  #
        #######################
        # the stored electrons and muons, depending on the assigned channel
        eles_pass = ak.where(
            criteria["em"],
            pair_em[idx_em]["0"],
            ak.where(
                criteria["mm"],
                eles[:, :0],
                ak.where(
                    criteria["ee"],
                    ak.concatenate(
                        [pair_ee[idx_ee]["0"], pair_ee[idx_ee]["1"]], axis=1
                    ),
                    eles,
                ),
            ),
        )

        muons_pass = ak.where(
            criteria["em"],
            pair_em[idx_em]["1"],
            ak.where(
                criteria["mm"],
                ak.concatenate([pair_mm[idx_mm]["0"], pair_mm[idx_mm]["1"]], axis=1),
                ak.where(criteria["ee"], muons[:, :0], muons),
            ),
        )

        lep_arrays = ak.zip(
            {
                "pt": ak.concatenate([eles_pass.pt, muons_pass.pt], axis=1),
                "eta": ak.concatenate([eles_pass.eta, muons_pass.eta], axis=1),
                "phi": ak.concatenate([eles_pass.phi, muons_pass.phi], axis=1),
                "m": ak.concatenate([eles_pass.mass, muons_pass.mass], axis=1),
                "ch": ak.concatenate([eles_pass.charge, muons_pass.charge], axis=1),
                "id": ak.concatenate(
                    [abs(eles_pass.pdgId), abs(muons_pass.pdgId)], axis=1
                ),
            }
        )
        if not isRealData:
            lep_arrays["gid"] = ak.fill_none(
                ak.concatenate(
                    [eles_pass.matched_gen.pdgId, muons_pass.matched_gen.pdgId],
                    axis=1,
                ),
                0,
            )

        if not isRealData:
            ############
            #  Genlep  #
            ############
            # Genlep: same with nominal BTA workflow
            _fix = lambda x: ak.fill_none(x, 0)
            is_lep = (
                lambda p: (abs(_fix(p.pdgId)) == 11)
                | (abs(_fix(p.pdgId)) == 13)
                | (abs(_fix(p.pdgId)) == 15)
            )
            is_WZ = lambda p: (abs(_fix(p.pdgId)) == 23) | (abs(_fix(p.pdgId)) == 24)
            is_heavy_hadron = lambda p, pid: (abs(_fix(p.pdgId)) // 100 == pid) | (
                abs(_fix(p.pdgId)) // 1000 == pid
            )

            sel = (
                is_lep(events.GenPart)
                & (events.GenPart.hasFlags("isLastCopy"))
                & (events.GenPart.pt > 3.0)
            )  # requires pT > 3 GeV
            genlep = events.GenPart[sel]

            # trace parents up to 4 generations (from BTA code)
            genlep_pa1G = genlep.parent
            genlep_pa2G = genlep.parent.parent
            genlep_pa3G = genlep.parent.parent.parent
            genlep_pa4G = genlep.parent.parent.parent.parent
            istau = abs(genlep_pa1G.pdgId) == 15
            isWZ = is_WZ(genlep_pa1G) | is_WZ(genlep_pa2G)
            isD = is_heavy_hadron(genlep_pa1G, 4) | is_heavy_hadron(genlep_pa2G, 4)
            isB = (
                is_heavy_hadron(genlep_pa1G, 5)
                | is_heavy_hadron(genlep_pa2G, 5)
                | is_heavy_hadron(genlep_pa3G, 5)
                | is_heavy_hadron(genlep_pa4G, 5)
            )

            Genlep = ak.zip(
                {
                    "pT": genlep.pt,
                    "eta": genlep.eta,
                    "phi": genlep.phi,
                    "pdgID": genlep.pdgId,
                    "status": genlep.status,
                    "mother": ak.fill_none(
                        ak.where(isB | isD, 5 * isB + 4 * isD, 10 * istau + 100 * isWZ),
                        0,
                    ),
                }
            )

            ###################
            #  Gen particles  #
            ###################
            # Gen particles of interest in the BTA ttbar workflow
            genpart = events.GenPart[events.GenPart.hasFlags("isHardProcess")]
            gen_channel = ak.prod(
                genpart[(abs(genpart.pdgId) == 11) | (abs(genpart.pdgId) == 13)].pdgId,
                axis=1,
            )

            gen_arrays = ak.zip(
                {
                    "pt": genpart.pt,
                    "eta": genpart.eta,
                    "phi": genpart.phi,
                    "m": genpart.mass,
                    "id": genpart.pdgId,
                }
            )

            #################
            #  Gen weights  #
            #################
            # top pT reweighting, according to BTA code
            # the numbers are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_dat
            genttbar = events.GenPart[
                (events.GenPart.hasFlags("isLastCopy"))
                & (abs(events.GenPart.pdgId) == 6)
            ]  # is last copy and is t/tbar
            gentop = ak.firsts(genttbar[genttbar.pdgId == 6])
            genantitop = ak.firsts(genttbar[genttbar.pdgId == -6])
            basic_vars["ttbar_ptweight"] = ak.fill_none(
                np.sqrt(
                    np.exp(0.0615 - 0.0005 * np.minimum(gentop.pt, 500))
                    * np.exp(0.0615 - 0.0005 * np.minimum(genantitop.pt, 500))
                ),
                1.0,
            )

            # store weights according to BTA code
            # note: in original BTA code, all weights are appended in a same variable ttbar_w. Here we separate them according to the logic in NanoAOD
            basic_vars["ttbar_w"] = events.genWeight
            lhe_pdf_w_arrays = getattr(events, "LHEPdfWeight", None)
            lhe_scale_w_arrays = getattr(events, "LHEScaleWeight", None)
            ps_w_arrays = getattr(events, "PSWeight", None)

        ###############
        #     Jet     #
        ###############
        # Jet: same with nominal BTA workflow
        jet = events.Jet[
            (events.Jet.pt > 20.0) & (abs(events.Jet.eta) < 2.5)
        ]  # basic selection
        zeros = ak.zeros_like(jet.pt, dtype=int)
        Jet = ak.zip(
            {
                # basic kinematics
                "pt": jet.pt,
                "eta": jet.eta,
                "phi": jet.phi,
                "mass": jet.mass,
                "uncorrpt": jet.pt_raw,
                # jet ID/pileup ID // !!!
                "looseID": jet.jetId >= 2,
                "tightID": jet.jetId >= 4,
                "tightlepvetoID": jet.jetId >= 6,
                # pileup ID (essentially userInt('puId106XUL18Id') for UL18)
                # PU ID for Run 3 is not ready
                # 'pileup_tightID':ak.values_astype((jet.puId & (1 << 0) > 0) | (jet.pt > 50.), int)  ,
                # 'pileup_mediumID': ak.values_astype((jet.puId & (1 << 1) > 0) | (jet.pt > 50.), int) ,
                # 'pileup_looseID': ak.values_astype((jet.puId & (1 << 2) > 0) | (jet.pt > 50.), int) ,
                # taggers/vars // JP/JBP to be calibrated.. !!!
                "area": jet.area,
                # Deep Jet
                "DeepFlavourBDisc": jet.btagDeepFlavB,
                "DeepFlavourCvsLDisc": jet.btagDeepFlavCvL,
                "DeepFlavourCvsBDisc": jet.btagDeepFlavCvB,
            }
        )
        if isRealData:
            Jet["vetomap"] = jet.veto

        if not isRealData:
            Jet["partonFlavour"] = jet.partonFlavour
            Jet["hadronFlavour"] = jet.hadronFlavour
            # corrected flavours defined in BTA
            Jet["flavour"] = ak.where(
                jet.hadronFlavour != 0,
                jet.hadronFlavour,
                ak.where(
                    (abs(jet.partonFlavour) == 4) | (abs(jet.partonFlavour) == 5),
                    zeros,
                    jet.partonFlavour,
                ),
            )

            # genJet pT
            genJetIdx = ak.where(
                (jet.genJetIdx < ak.num(events.GenJet)) & (jet.genJetIdx != -1),
                jet.genJetIdx,
                zeros,
            )  # in case the genJet index out of range

            Jet["genpt"] = ak.where(
                jet.genJetIdx != -1, events.GenJet[genJetIdx].pt, -99
            )

            # gen-level jet cleaning aginst prompt leptons
            genlep_prompt = genlep[(Genlep.mother != 0) & (Genlep.mother % 10 == 0)]
            overlapped = ak.zeros_like(jet.pt, dtype=bool)

            for lep_pid, dr_thres in [(11, 0.2), (13, 0.2), (15, 0.3)]:
                lep = genlep_prompt[abs(genlep_prompt.pdgId) == lep_pid]
                pair = ak.cartesian(
                    [jet, lep], axis=1, nested=True
                )  # dim: (event, jet, lep)
                dr_min = ak.min(
                    pair["0"].delta_r(pair["1"]), axis=-1
                )  # dim: (event, jet)
                dr_min = ak.fill_none(dr_min, 99.0)  # for 0 lepton case, assign dr=99
                overlapped = overlapped | (dr_min < dr_thres)
            Jet["flavourCleaned"] = ak.where(
                overlapped, zeros - 999, Jet["flavour"]
            )  # fill -999 for jets with lepton overlapping

            # jet cleaning aginst pileup
            Jet["flavourCleaned"] = ak.where(
                (Jet["genpt"] < 8.0) & (Jet["genpt"] > 0), zeros, Jet["flavourCleaned"]
            )

        # define Jet_clean: special for BTA ttbar workflow which cleans jets against selected leptons
        overlapped = ak.zeros_like(jet.pt, dtype=bool)
        for lep in [eles, muons]:
            pair = ak.cartesian(
                [jet, lep], axis=1, nested=True
            )  # dim: (event, jet, lep)
            dr_min = ak.min(pair["0"].delta_r(pair["1"]), axis=-1)  # dim: (event, jet)
            dr_min = ak.fill_none(dr_min, 99.0)  # for 0 lepton case, assign dr=99
            overlapped = overlapped | (dr_min < 0.4)

        Jet_clean = Jet[overlapped == False]  # remove jets overlapping with leptons

        ###############
        #     MET     #
        ###############
        if "Run3" in self._campaign:
            basic_vars["ttbar_met_pt"] = events.PuppiMET.pt
            basic_vars["ttbar_met_phi"] = events.PuppiMET.phi
        else:
            basic_vars["ttbar_met_pt"] = events.MET.pt
            basic_vars["ttbar_met_phi"] = events.MET.phi

        ###############
        #  Selection  #
        ###############
        # apply selections on events based on BTA TTbarSelectionProducer
        passLepSel = chsel != 0
        passJetSel = (
            ((abs(chsel) == 11) | (abs(chsel) == 13)) & (ak.num(Jet_clean) >= 4)
        ) | ((abs(chsel) > 13) & (ak.num(Jet_clean) >= 1))
        passMetSel = events.PuppiMET.pt > 0

        # and the channel selection, configured in TTbarSelectionFilter
        # https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements/blob/10_6_X/python/TTbarSelectionFilter_cfi.py#L4
        passChannelSel = (chsel == -11 * 11) | (chsel == -11 * 13) | (chsel == -13 * 13)

        # final event selection
        passEvent = passJetSel & passLepSel & passMetSel & passChannelSel

        ###############
        #  Write root #
        ###############
        if ak.any(passEvent) == False:
            return {dataset: 0}

        output = {
            **{k: v[passEvent] for k, v in basic_vars.items()},
            "Jet": Jet_clean[passEvent],
            "ttbar_lep": lep_arrays[passEvent],
        }
        if not isRealData:
            output["ttbar_gen"] = gen_arrays[passEvent]
            if lhe_pdf_w_arrays is not None:
                output["ttbar_lhe_pdf_w"] = lhe_pdf_w_arrays[passEvent]
            if lhe_scale_w_arrays is not None:
                output["ttbar_lhe_scale_w"] = lhe_scale_w_arrays[passEvent]
            if ps_w_arrays is not None:
                output["ttbar_ps_w"] = ps_w_arrays[passEvent]

        # customize output file name: <dataset>_<nanoAOD file name>_<chunk index>.root
        fname = f"{dataset}/{events.metadata['filename'].split('/')[-1].replace('.root','')}_{int(events.metadata['entrystop']/self.chunksize)}.root"
        os.system(f"mkdir -p {dataset}")
        with uproot.recreate(fname) as fout:
            output_root = {}
            for bname in output.keys():
                if not output[bname].fields:
                    output_root[bname] = ak.packed(ak.without_parameters(output[bname]))
                else:
                    b_nest = {}
                    for n in output[bname].fields:
                        b_nest[n] = ak.packed(ak.without_parameters(output[bname][n]))
                    output_root[bname] = ak.zip(b_nest)
            fout["btagana/ttree"] = output_root
            if isRealData:
                fout["sumw"] = {"total_events": ak.Array([len(events)])}
            else:
                fout["sumw"] = {
                    "total_events": ak.Array([len(events)]),
                    "total_pos_events": ak.Array([ak.sum(events.genWeight > 0)]),
                    "total_neg_events": ak.Array([-1.0 * ak.sum(events.genWeight < 0)]),
                    "total_wei_events": ak.Array([ak.sum(events.genWeight)]),
                    "total_poswei_events": ak.Array(
                        [ak.sum(events.genWeight[events.genWeight > 0.0])]
                    ),
                    "total_negwei_events": ak.Array(
                        [ak.sum(events.genWeight[events.genWeight < 0.0])]
                    ),
                }

        return {dataset: len(events)}

    def postprocess(self, accumulator):
        return accumulator
