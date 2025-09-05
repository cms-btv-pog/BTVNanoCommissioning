import collections
import numpy as np, awkward as ak

from coffea import processor
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)

from BTVNanoCommissioning.helpers.func import update, dump_lumi
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper, jet_id, MET_filters,
    mu_idiso, ele_mvatightid,
    btag_wp,
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
        selectionModifier="mu",  # "mu" or "el"
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        if selectionModifier not in ["el", "mu"]:
            raise ValueError(f"Invalid selectionModifier: {selectionModifier}")
        self.channel = selectionModifier
        self.SF_map = load_SF(self._year, self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)
        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        output = {} if self.noHist else histogrammer(events, "sf_ttsemilep_tnp")

        if shift_name is None:
            output["sumw"] = len(events) if isRealData else ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################
        # Lumi mask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        # HLT
        triggers = ["IsoMu24", "Mu50"] if self.channel == "mu" else ["Ele32_WPTight_Gsf", "Ele50_CaloIdVT_GsfTrkIdT_PFJet165"]
        req_trig = HLT_helper(events, triggers)

        # MET filters (normalized campaign key)
        req_metf = MET_filters(events, self._campaign)

        # Veto leptons
        veto_mu = events.Muon[(events.Muon.pt > 15) & (abs(events.Muon.eta) < 2.4) & ak.fill_none(events.Muon.looseId, False)]
        veto_el = events.Electron[
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.4)
            & ~((abs(events.Electron.eta) > 1.4442) & (abs(events.Electron.eta) < 1.566))
            & (ak.fill_none(events.Electron.cutBased >= 2, False)
               | ak.fill_none(getattr(events.Electron, "mvaIso_WP90", False), False))
        ]

        # Lepton selection 
        if self.channel == "mu":
            mu_mask = (events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & mu_idiso(events, self._campaign)
            iso_lep = events.Muon[mu_mask]
            req_lepveto = (ak.num(veto_mu, axis=1) == 1) & (ak.num(veto_el, axis=1) == 0)
        else:
            el_mask = (events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4) & ele_mvatightid(events, self._campaign)
            iso_lep = events.Electron[el_mask]
            req_lepveto = (ak.num(veto_mu, axis=1) == 0) & (ak.num(veto_el, axis=1) == 1)
        req_lep = ak.num(iso_lep, axis=1) == 1

        # delta R cleaning vs veto leptons (keep jets near selected iso lepton)
        dr_mu = events.Jet.metric_table(veto_mu)
        dr_el = events.Jet.metric_table(veto_el)
        clean_mu = ak.all(dr_mu > 0.4, axis=2, mask_identity=True)
        clean_el = ak.all(dr_el > 0.4, axis=2, mask_identity=True)

        base_jet_mask = jet_id(events, self._campaign)
        jet_mask = ak.fill_none(base_jet_mask & clean_mu & clean_el, False, axis=-1)
        event_jet = events.Jet[jet_mask]
        req_jets = ak.num(event_jet) >= 4

        # b-tag medium requirement
        # TODO: check if DeepFlav is adequate
        bmask_M = btag_wp(event_jet, self._year, self._campaign, tagger="DeepFlav", borc="b", wp="M")
        nbtagM  = ak.sum(ak.fill_none(bmask_M, False), axis=1)
        req_btags = nbtagM >= 2

        # Mask for event selection
        evmask = req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets & req_btags
        if not ak.any(evmask):
            pruned_ev = events[evmask]
            weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
            if not self.noHist:
                output = histo_writter(pruned_ev, output, weights, ["nominal"], self.isSyst, self.SF_map)
            return {dataset: output}

        # Object selection
        ev = events[evmask]
        jets = event_jet[evmask]
        bmask_M = bmask_M[evmask]
        lep = ak.firsts(iso_lep[evmask])

        # MET four-vector
        met = ak.zip(
            {
                "x": ev.MET.pt * np.cos(ev.MET.phi),
                "y": ev.MET.pt * np.sin(ev.MET.phi),
                "z": ak.zeros_like(ev.MET.pt),
                "t": ev.MET.pt,
            },
            with_name="FourVector",
        )

        #############################
        # ttbar reconstruction (experimental)
        #############################
        # TODO: check if equivalent to otto's code

        # Cap at N leading jets to keep combos in check
        NMAX = 6
        jets = jets[:, :NMAX]
        bmask_M = bmask_M[:, :NMAX]

        # Indices per event
        idx = ak.local_index(jets, axis=1)

        # All ordered (BLep, BHad) pairs with bl != bh
        b_pairs = ak.cartesian({"bl": idx, "bh": idx}, axis=1, nested=True)
        b_pairs = b_pairs[b_pairs.bl != b_pairs.bh]

        # All unordered (Ja, Jb) pairs for W (ja < jb)
        w_pairs = ak.combinations(idx, 2, replacement=False, axis=1)

        # Cross product → all candidate assignments
        comb = ak.cartesian({"b": b_pairs, "w": w_pairs}, axis=1, nested=True)

        # Keep only candidates with 4 distinct indices
        distinct = (
            (comb.b.bl != comb.w["0"])
            & (comb.b.bl != comb.w["1"])
            & (comb.b.bh != comb.w["0"])
            & (comb.b.bh != comb.w["1"])
        )
        comb = comb[distinct]

        has_cand = ak.num(comb.b.bl, axis=1) > 0

        # Gather jets for each role
        iBL, iBH, iJa, iJb = comb.b.bl, comb.b.bh, comb.w["0"], comb.w["1"]
        BL = ak.take(jets, iBL, axis=1)
        BH = ak.take(jets, iBH, axis=1)
        JA = ak.take(jets, iJa, axis=1)
        JB = ak.take(jets, iJb, axis=1)

        # Broadcast lepton & MET to candidate shape
        lep_cand = ak.broadcast_arrays(lep, BL)[0]
        met_cand = ak.broadcast_arrays(met, BL)[0]

        # Lepton subtraction for BLep if ΔR(BL, lep) < 0.4
        dphi = ak.where(
            (BL.phi - lep_cand.phi) > np.pi, BL.phi - lep_cand.phi - 2*np.pi,
            ak.where((BL.phi - lep_cand.phi) < -np.pi, BL.phi - lep_cand.phi + 2*np.pi, BL.phi - lep_cand.phi)
        )
        dR = ak.sqrt((BL.eta - lep_cand.eta) ** 2 + dphi ** 2)
        overlap = dR < 0.4

        BL_px = ak.where(overlap, BL.px - lep_cand.px, BL.px)
        BL_py = ak.where(overlap, BL.py - lep_cand.py, BL.py)
        BL_pz = ak.where(overlap, BL.pz - lep_cand.pz, BL.pz)
        BL_e  = ak.where(overlap, BL.energy - lep_cand.energy, BL.energy)

        # Vectorized neutrino solution (W-mass constraint), choose root closer to m_top^lep
        px_l, py_l, pz_l, e_l = lep_cand.px, lep_cand.py, lep_cand.pz, lep_cand.energy
        px_n, py_n = met_cand.x, met_cand.y

        mW, mT = 80.4, 172.5
        muW = (mW**2) / 2.0 + px_l * px_n + py_l * py_n
        a = (e_l**2 - pz_l**2)
        b = -2.0 * muW * pz_l
        c = e_l**2 * (px_n**2 + py_n**2) - muW**2
        disc = b**2 - 4.0 * a * c

        sqrt_disc = ak.where(disc > 0, ak.sqrt(disc), 0.0)
        pz_n_1 = (-b + sqrt_disc) / (2.0 * a)
        pz_n_2 = (-b - sqrt_disc) / (2.0 * a)

        e_n_1 = ak.sqrt(px_n**2 + py_n**2 + pz_n_1**2)
        e_n_2 = ak.sqrt(px_n**2 + py_n**2 + pz_n_2**2)

        # Build both TLep masses, choose root closer to mT
        tlep1_e  = e_l + e_n_1 + BL_e
        tlep1_px = px_l + px_n + BL_px
        tlep1_py = py_l + py_n + BL_py
        tlep1_pz = pz_l + pz_n_1 + BL_pz
        m_tlep1  = ak.sqrt(ak.maximum(0.0, tlep1_e**2 - (tlep1_px**2 + tlep1_py**2 + tlep1_pz**2)))

        tlep2_e  = e_l + e_n_2 + BL_e
        tlep2_px = px_l + px_n + BL_px
        tlep2_py = py_l + py_n + BL_py
        tlep2_pz = pz_l + pz_n_2 + BL_pz
        m_tlep2  = ak.sqrt(ak.maximum(0.0, tlep2_e**2 - (tlep2_px**2 + tlep2_py**2 + tlep2_pz**2)))

        choose_1 = (abs(m_tlep1 - mT) < abs(m_tlep2 - mT))
        pz_n = ak.where(choose_1, pz_n_1, pz_n_2)
        e_n  = ak.where(choose_1, e_n_1, e_n_2)

        # Final TLep 4-vector + guard window (Otto’s Dnu proxy)
        tlep_e  = ak.where(choose_1, tlep1_e, tlep2_e)
        tlep_px = ak.where(choose_1, tlep1_px, tlep2_px)
        tlep_py = ak.where(choose_1, tlep1_py, tlep2_py)
        tlep_pz = ak.where(choose_1, tlep1_pz, tlep2_pz)
        m_tlep  = ak.where(choose_1, m_tlep1, m_tlep2)
        valid_dnu = (m_tlep > 100.0) & (m_tlep < 240.0)

        # Hadronic side masses
        W_px = JA.px + JB.px
        W_py = JA.py + JB.py
        W_pz = JA.pz + JB.pz
        W_e  = JA.energy + JB.energy
        m_whad = ak.sqrt(ak.maximum(0.0, W_e**2 - (W_px**2 + W_py**2 + W_pz**2)))

        th_px = BH.px + W_px
        th_py = BH.py + W_py
        th_pz = BH.pz + W_pz
        th_e  = BH.energy + W_e
        m_thad = ak.sqrt(ak.maximum(0.0, th_e**2 - (th_px**2 + th_py**2 + th_pz**2)))

        # χ² per candidate
        sig_w, sig_t = 30.0, 40.0
        chi2 = ((m_whad - mW) / sig_w) ** 2 + ((m_thad - mT) / sig_t) ** 2 + ((m_tlep - mT) / sig_t) ** 2
        chi2 = ak.where(valid_dnu, chi2, np.inf)

        # Pick best candidate per event (take first if degenerate minima)
        best_mask = chi2 == ak.min(chi2, axis=1, initial=np.inf)
        pick = lambda arr: ak.firsts(arr[best_mask])
        best_BL, best_BH = pick(BL), pick(BH)
        best_JA, best_JB = pick(JA), pick(JB)

        best_nu = ak.zip(
            {"x": ak.firsts(px_n[best_mask]),
             "y": ak.firsts(py_n[best_mask]),
             "z": ak.firsts(pz_n[best_mask]),
             "t": ak.firsts(e_n[best_mask])},
            with_name="FourVector",
        )
        best_tlep = ak.zip(
            {"x": ak.firsts((tlep_px)[best_mask]),
             "y": ak.firsts((tlep_py)[best_mask]),
             "z": ak.firsts((tlep_pz)[best_mask]),
             "t": ak.firsts((tlep_e )[best_mask])},
            with_name="FourVector",
        )
        best_thad = ak.zip(
            {"x": ak.firsts((th_px)[best_mask]),
             "y": ak.firsts((th_py)[best_mask]),
             "z": ak.firsts((th_pz)[best_mask]),
             "t": ak.firsts((th_e )[best_mask])},
            with_name="FourVector",
        )
        best_chi2 = ak.firsts(chi2[best_mask])

        # Category based on DeepFlav-M of the chosen b jets
        bl_isM = ak.firsts(ak.take(bmask_M, iBL, axis=1)[best_mask])
        bh_isM = ak.firsts(ak.take(bmask_M, iBH, axis=1)[best_mask])
        twom   = ak.values_astype(bl_isM, np.int32) + ak.values_astype(bh_isM, np.int32)
        best_cat = ak.where(twom == 2, "RES", ak.where(twom == 1, "ONEB", "NONE"))

        # Assemble pruned_ev
        pruned_ev = ev[has_cand]
        pruned_ev = ak.with_field(pruned_ev, lep[has_cand], "sel_lep")
        pruned_ev = ak.with_field(pruned_ev, jets[has_cand], "sel_jets")
        pruned_ev = ak.with_field(pruned_ev, ak.sum(ak.fill_none(bmask_M[has_cand], False), axis=1), "nbtagM")
        pruned_ev = ak.with_field(pruned_ev, best_BL[has_cand], "blep")
        pruned_ev = ak.with_field(pruned_ev, best_BH[has_cand], "bhad")
        pruned_ev = ak.with_field(pruned_ev, best_JA[has_cand], "ja")
        pruned_ev = ak.with_field(pruned_ev, best_JB[has_cand], "jb")
        pruned_ev = ak.with_field(pruned_ev, best_nu[has_cand], "nu")
        pruned_ev = ak.with_field(pruned_ev, best_tlep[has_cand], "tlep")
        pruned_ev = ak.with_field(pruned_ev, best_thad[has_cand], "thad")
        pruned_ev = ak.with_field(pruned_ev, best_chi2[has_cand], "chi2")
        pruned_ev = ak.with_field(pruned_ev, best_cat[has_cand], "ttcat")
        pruned_ev = ak.with_field(
            pruned_ev,
            ak.zip(
                {"x": best_tlep.x + best_thad.x,
                 "y": best_tlep.y + best_thad.y,
                 "z": best_tlep.z + best_thad.z,
                 "t": best_tlep.t + best_thad.t},
                with_name="FourVector",
            ),
            "tt",
        )

        ####################
        #     Output       #
        ####################
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        systematics = [shift_name] if shift_name is not None else ["nominal"] + list(weights.variations)

        if not self.noHist:
            output = histo_writter(pruned_ev, output, weights, systematics, self.isSyst, self.SF_map)
        if self.isArray:
            array_writer(self, pruned_ev, events, weights, systematics, dataset, isRealData)

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
