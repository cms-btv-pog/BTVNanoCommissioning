import collections
import numpy as np, awkward as ak
import hist as Hist

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
    btag_wp, wp_dict
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
        tag_tagger = "UParTAK4",
        tag_wp = "M",
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
        self._wp_table = wp_dict(year, campaign)
        self._all_taggers = sorted(self._wp_table.keys())
        self.tag_tagger = tag_tagger
        self.tag_wp = tag_wp

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
        triggers_mu = ["IsoMu24", "Mu50"]
        triggers_el = ["Ele32_WPTight_Gsf", "Ele50_CaloIdVT_GsfTrkIdT_PFJet165"]
        triggers = triggers_mu if self.channel == "mu" else triggers_el
        req_trig = HLT_helper(events, triggers)

        # MET filters
        req_metf = MET_filters(events, self._campaign)

        # Loose (veto) masks
        mu_loose_mask = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < 2.4)
            & ak.fill_none(events.Muon.looseId, False)
        )
        el_loose_mask = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.4)
            & ~((abs(events.Electron.eta) > 1.4442) & (abs(events.Electron.eta) < 1.566))
            & (ak.fill_none(events.Electron.cutBased >= 2, False)
            | ak.fill_none(getattr(events.Electron, "mvaIso_WP90", False), False))
        )

        # Tight lepton selection
        if self.channel == "mu":
            mu_tight_mask = (events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & mu_idiso(events, self._campaign)
            iso_lep = events.Muon[mu_tight_mask]
            # exactly one tight muon, zero loose electrons
            req_lepveto = (ak.num(events.Muon[mu_loose_mask], axis=1) == 1) & (ak.num(events.Electron[el_loose_mask], axis=1) == 0)
            # other loose muons (exclude the tight one) for cleaning
            other_mu = events.Muon[mu_loose_mask & ~mu_tight_mask]
            other_el = events.Electron[el_loose_mask]
        else:
            el_tight_mask = (events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4) & ele_mvatightid(events, self._campaign)
            iso_lep = events.Electron[el_tight_mask]
            req_lepveto = (ak.num(events.Muon[mu_loose_mask], axis=1) == 0) & (ak.num(events.Electron[el_loose_mask], axis=1) == 1)
            other_mu = events.Muon[mu_loose_mask]
            other_el = events.Electron[el_loose_mask & ~el_tight_mask]

        req_lep = ak.num(iso_lep, axis=1) == 1

        # ΔR cleaning vs *other* leptons (not the selected tight one)
        dr_mu = events.Jet.metric_table(other_mu)
        dr_el = events.Jet.metric_table(other_el)
        all_true_jets = ak.ones_like(events.Jet.pt, dtype=bool)
        has_o_mu = ak.num(other_mu, axis=1) > 0
        has_o_el = ak.num(other_el, axis=1) > 0
        clean_mu = ak.where(
            has_o_mu,
            ak.all(dr_mu > 0.4, axis=-1, mask_identity=True),
            all_true_jets
        )
        clean_el = ak.where(
            has_o_el,
            ak.all(dr_el > 0.4, axis=-1, mask_identity=True),
            all_true_jets
        )

        base_jet_mask = jet_id(events, self._campaign)
        jet_mask = ak.fill_none(base_jet_mask & clean_mu & clean_el, False, axis=-1)
        event_jet = events.Jet[jet_mask]
        req_jets = ak.num(event_jet, axis=1) == 4

        # require at least one tag jet
        bmask_tag = btag_wp(event_jet, self._year, self._campaign, tagger=self.tag_tagger,   borc="b", wp=self.tag_wp)
        nbtag_tag = ak.sum(ak.fill_none(bmask_tag, False), axis=1)
        req_btags = nbtag_tag >= 1

        # --- cutflow helper ---
        if not self.noHist and "cutflow" not in output:
            cf_axis = Hist.axis.StrCategory([], name="step", growth=True)
            output["cutflow"] = Hist.Hist(cf_axis, Hist.storage.Weight())
        def _cf(step, mask):
            if self.noHist:
                return
            if isinstance(mask, ak.Array):
                mask = ak.to_numpy(mask)
            mask = np.asarray(mask, dtype=bool)
            output["cutflow"].fill(step, weight=float(mask.sum()))
        def _expand_full(full_len, base_sel, per_sel):
            """per_sel has length sum(base_sel); expand back to length full_len"""
            arr = np.zeros(full_len, dtype=bool)
            arr[base_sel] = ak.to_numpy(per_sel)
            return arr

        _cf("all", np.ones(len(events), dtype=bool))
        _cf("lumi", req_lumi)
        _cf("trig", req_lumi & req_trig)
        _cf("metf", req_lumi & req_trig & req_metf)
        _cf("lepveto", req_lumi & req_trig & req_metf & req_lepveto)
        _cf("tightlep", req_lumi & req_trig & req_metf & req_lepveto & req_lep)
        _cf("jets==4", req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets)
        _cf("bjets>=1", req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets & req_btags)

        # Event mask
        evmask = req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets & req_btags
        if not ak.any(evmask):
            return {dataset: output}

        # Objects after event selection mask
        ev = events[evmask]
        jets = event_jet[evmask]
        bmask_tag = ak.fill_none(bmask_tag[evmask], False)
        lep = ak.firsts(iso_lep[evmask])

        # MET 4-vector
        met = ak.zip(
            {"x": ev.MET.pt * np.cos(ev.MET.phi),
             "y": ev.MET.pt * np.sin(ev.MET.phi),
             "z": ak.zeros_like(ev.MET.pt),
             "t": ev.MET.pt},
            with_name="FourVector",
        )

        # -------- ttbar reconstruction --------

        # Limit to N leading jets for speed 
        # (keep this for flexibility, in case the jet multiplicity cut is loosened)
        NMAX = 4
        jets = jets[:, :NMAX]
        bmask_tag = bmask_tag[:, :NMAX]

        # Indices
        idx = ak.local_index(jets, axis=1)

        # (bl,bh) ordered pairs, bl != bh ; (ja,jb) combinations
        b_pairs = ak.cartesian({"bl": idx, "bh": idx}, axis=1, nested=False)
        b_pairs = b_pairs[b_pairs.bl != b_pairs.bh]
        w_pairs = ak.combinations(idx, 2, axis=1, replacement=False)

        # all assignments per event
        comb = ak.cartesian({"b": b_pairs, "w": w_pairs}, axis=1, nested=False)

        # keep only distinct indices
        distinct = (
            (comb.b.bl != comb.w["0"])
            & (comb.b.bl != comb.w["1"])
            & (comb.b.bh != comb.w["0"])
            & (comb.b.bh != comb.w["1"])
        )
        comb = comb[distinct]

        has_cand = ak.num(comb.b.bl, axis=1) > 0
        if not ak.any(has_cand):
            return {dataset: output}

        # Cutflow: tt candidate found (expand back to full length)
        _cf("tt cand found", _expand_full(len(events), ak.to_numpy(evmask), has_cand))

        # Role assignment:
        # - BL: b-quark from leptonic top decay
        # - BH: b-quark from hadronic top decay
        # - JA: first jet from hadronic W decay (label: "A")
        # - JB: second jet from hadronic W decay (label: "B")
        iBL, iBH, iJa, iJb = comb.b.bl, comb.b.bh, comb.w["0"], comb.w["1"]
        BL = jets[iBL]
        BH = jets[iBH]
        JA = jets[iJa]
        JB = jets[iJb]

        # Broadcast lepton & MET
        lep_eta = ak.broadcast_arrays(lep.eta, BL.eta)[0]
        lep_phi = ak.broadcast_arrays(lep.phi, BL.phi)[0]
        lep_px  = ak.broadcast_arrays(lep.px , BL.px )[0]
        lep_py  = ak.broadcast_arrays(lep.py , BL.py )[0]
        lep_pz  = ak.broadcast_arrays(lep.pz , BL.pz )[0]
        lep_e   = ak.broadcast_arrays(lep.energy, BL.energy)[0]

        met_x = ak.broadcast_arrays(met.x, BL.px)[0]
        met_y = ak.broadcast_arrays(met.y, BL.px)[0]

        # BLep lepton-subtraction if ΔR < 0.4
        dphi = ak.where(
            (BL.phi - lep_phi) > np.pi, BL.phi - lep_phi - 2*np.pi,
            ak.where((BL.phi - lep_phi) < -np.pi, BL.phi - lep_phi + 2*np.pi, BL.phi - lep_phi)
        )
        dR = np.sqrt((BL.eta - lep_eta) ** 2 + dphi ** 2)
        overlap = dR < 0.4
        BL_px = ak.where(overlap, BL.px - lep_px, BL.px)
        BL_py = ak.where(overlap, BL.py - lep_py, BL.py)
        BL_pz = ak.where(overlap, BL.pz - lep_pz, BL.pz)
        BL_e  = ak.where(overlap, BL.energy - lep_e, BL.energy)

        # Neutrino pz from W mass (choose root closer to mT)
        px_l, py_l, pz_l, e_l = lep_px, lep_py, lep_pz, lep_e
        px_n, py_n = met_x, met_y
        mW, mT = 80.4, 172.5
        muW = (mW**2)/2.0 + px_l*px_n + py_l*py_n
        a = (e_l**2 - pz_l**2)
        b = -2.0*muW*pz_l
        c = e_l**2*(px_n**2 + py_n**2) - muW**2
        disc = b**2 - 4.0*a*c
        sqrt_disc = ak.where(disc > 0, np.sqrt(disc), 0.0)
        pz_n_1 = (-b + sqrt_disc)/(2.0*a)
        pz_n_2 = (-b - sqrt_disc)/(2.0*a)
        e_n_1 = np.sqrt(px_n**2 + py_n**2 + pz_n_1**2)
        e_n_2 = np.sqrt(px_n**2 + py_n**2 + pz_n_2**2)

        # tlep masses for both roots
        tlep1_e  = e_l + e_n_1 + BL_e
        tlep1_px = px_l + px_n + BL_px
        tlep1_py = py_l + py_n + BL_py
        tlep1_pz = pz_l + pz_n_1 + BL_pz
        m_tlep1 = np.sqrt(np.maximum(0.0, tlep1_e**2 - (tlep1_px**2 + tlep1_py**2 + tlep1_pz**2)))
        tlep2_e  = e_l + e_n_2 + BL_e
        tlep2_px = px_l + px_n + BL_px
        tlep2_py = py_l + py_n + BL_py
        tlep2_pz = pz_l + pz_n_2 + BL_pz
        m_tlep2 = np.sqrt(np.maximum(0.0, tlep2_e**2 - (tlep2_px**2 + tlep2_py**2 + tlep2_pz**2)))

        choose_1 = (abs(m_tlep1 - mT) < abs(m_tlep2 - mT))
        pz_n = ak.where(choose_1, pz_n_1, pz_n_2)
        e_n  = ak.where(choose_1, e_n_1, e_n_2)
        tlep_e  = ak.where(choose_1, tlep1_e,  tlep2_e)
        tlep_px = ak.where(choose_1, tlep1_px, tlep2_px)
        tlep_py = ak.where(choose_1, tlep1_py, tlep2_py)
        tlep_pz = ak.where(choose_1, tlep1_pz, tlep2_pz)
        m_tlep  = ak.where(choose_1, m_tlep1,  m_tlep2)
        valid_dnu = (m_tlep > 100.0) & (m_tlep < 240.0)

        # Hadronic side
        W_px = JA.px + JB.px
        W_py = JA.py + JB.py
        W_pz = JA.pz + JB.pz
        W_e  = JA.energy + JB.energy
        m_whad = np.sqrt(np.maximum(0.0, W_e**2  - (W_px**2  + W_py**2  + W_pz**2)))
        th_px = BH.px + W_px
        th_py = BH.py + W_py
        th_pz = BH.pz + W_pz
        th_e  = BH.energy + W_e
        m_thad = np.sqrt(np.maximum(0.0, th_e**2 - (th_px**2 + th_py**2 + th_pz**2)))

        # χ² & best candidates
        sig_w, sig_t = 30.0, 40.0
        chi2 = ((m_whad - mW)/sig_w)**2 + ((m_thad - mT)/sig_t)**2 + ((m_tlep - mT)/sig_t)**2
        chi2 = ak.where(valid_dnu, chi2, np.inf)
        best_idx = ak.where(has_cand, ak.argmin(chi2, axis=1), -1)
        li       = ak.local_index(chi2, axis=1)
        best_mask = (li == best_idx) & (best_idx >= 0)

        best_BL = ak.firsts(BL[best_mask])
        best_BH = ak.firsts(BH[best_mask])
        best_JA = ak.firsts(JA[best_mask])
        best_JB = ak.firsts(JB[best_mask])
        W_best = ak.zip(
            {
                "x": ak.firsts((JA.x + JB.x)[best_mask]),
                "y": ak.firsts((JA.y + JB.y)[best_mask]),
                "z": ak.firsts((JA.z + JB.z)[best_mask]),
                "t": ak.firsts((JA.t + JB.t)[best_mask]),
            },
            with_name="FourVector",
        )
        best_nu = ak.zip(
            {"x": ak.firsts(px_n[best_mask]),
             "y": ak.firsts(py_n[best_mask]),
             "z": ak.firsts(pz_n[best_mask]),
             "t": ak.firsts(e_n[best_mask])},
            with_name="FourVector",
        )
        best_tlep = ak.zip(
            {"x": ak.firsts(tlep_px[best_mask]),
             "y": ak.firsts(tlep_py[best_mask]),
             "z": ak.firsts(tlep_pz[best_mask]),
             "t": ak.firsts(tlep_e [best_mask])},
            with_name="FourVector",
        )
        best_thad = ak.zip(
            {"x": ak.firsts(th_px[best_mask]),
             "y": ak.firsts(th_py[best_mask]),
             "z": ak.firsts(th_pz[best_mask]),
             "t": ak.firsts(th_e [best_mask])},
            with_name="FourVector",
        )
        best_chi2 = ak.firsts(chi2[best_mask])

        # helper for FourVector mass / pt (NumPy ops are fine on Awkward arrays)
        def _four_mass(v):
            return np.sqrt(np.maximum(0.0, v.t * v.t - (v.x * v.x + v.y * v.y + v.z * v.z)))

        def _four_pt(v):
            return np.sqrt(v.x * v.x + v.y * v.y)

        # masses & pTs (use the “best” candidates)
        tlep_mass = _four_mass(best_tlep)
        thad_mass = _four_mass(best_thad)
        whad_mass = _four_mass(W_best)

        tlep_pt = _four_pt(best_tlep)
        thad_pt = _four_pt(best_thad)

        # simple ΔR helper (objects here have eta/phi)
        def _dR(a, b):
            dphi = ak.where((a.phi - b.phi) > np.pi, a.phi - b.phi - 2*np.pi,
                ak.where((a.phi - b.phi) < -np.pi, a.phi - b.phi + 2*np.pi, a.phi - b.phi))
            return np.sqrt((a.eta - b.eta) ** 2 + dphi ** 2)

        dr_lep_blep = _dR(lep,     best_BL)
        dr_lep_bhad = _dR(lep,     best_BH)
        dr_ja_jb    = _dR(best_JA, best_JB)

        drs = ak.stack([
            _dr(best_BH, best_JA), _dr(best_BH, best_JB), _dr(best_BH, best_BL),
            _dr(best_BL, best_JA), _dr(best_BL, best_JB)
        ], axis=1)
        mindr = ak.min(drs, axis=1)

        # tag and probe code flags (the actual logic is in the histogrammer)
        had_tag = ak.fill_none(ak.firsts(bmask_tag[iBH][best_mask]), False)
        lep_tag = ak.fill_none(ak.firsts(bmask_tag[iBL][best_mask]), False)
        had_pt_ok = ak.fill_none(ak.firsts(BH.pt[best_mask]) >= 30.0, False)
        lep_pt_ok = ak.fill_none(ak.firsts(BL.pt[best_mask]) >= 30.0, False)
        tnp_had_fill = had_pt_ok & lep_tag & (mindr >= 0.8)
        tnp_lep_fill = lep_pt_ok & had_tag & (mindr >= 0.8)

        had_pt_for_tnp = ak.firsts(BH.pt[best_mask])
        lep_pt_for_tnp = ak.firsts(BL.pt[best_mask])

        # Assemble pruned_ev (aliases for histogrammer)
        pruned_ev = ev[has_cand]
        pruned_ev = ak.with_field(pruned_ev, lep[has_cand], "sel_lep")
        pruned_ev = ak.with_field(pruned_ev, jets[has_cand], "sel_jets")
        pruned_ev = ak.with_field(pruned_ev, best_BL[has_cand], "blep")
        pruned_ev = ak.with_field(pruned_ev, best_BH[has_cand], "bhad")
        pruned_ev = ak.with_field(pruned_ev, best_JA[has_cand], "ja")
        pruned_ev = ak.with_field(pruned_ev, best_JB[has_cand], "jb")
        pruned_ev = ak.with_field(pruned_ev, best_nu[has_cand], "nu")
        pruned_ev = ak.with_field(pruned_ev, best_tlep[has_cand], "tlep")
        pruned_ev = ak.with_field(pruned_ev, best_thad[has_cand], "thad")
        pruned_ev = ak.with_field(pruned_ev, best_chi2[has_cand], "chi2")

        # map aliases for histogrammer
        pruned_ev = ak.with_field(pruned_ev, pruned_ev.sel_jets[:, :4], "SelJet")  # trim to 4
        if self.channel == "mu":
            pruned_ev = ak.with_field(pruned_ev, pruned_ev.sel_lep, "SelMuon")
        else:
            pruned_ev = ak.with_field(pruned_ev, pruned_ev.sel_lep, "SelElectron")
        pruned_ev = ak.with_field(pruned_ev, ev.MET[has_cand], "MET")
        pruned_ev = ak.with_field(pruned_ev, ak.num(pruned_ev.SelJet, axis=1), "njet")
        pruned_ev = ak.with_field(pruned_ev, tlep_mass[has_cand], "tlep_mass")
        pruned_ev = ak.with_field(pruned_ev, thad_mass[has_cand], "thad_mass")
        pruned_ev = ak.with_field(pruned_ev, whad_mass[has_cand], "whad_mass")

        pruned_ev = ak.with_field(pruned_ev, tlep_pt[has_cand],  "tlep_pt")
        pruned_ev = ak.with_field(pruned_ev, thad_pt[has_cand],  "thad_pt")

        pruned_ev = ak.with_field(pruned_ev, dr_lep_blep[has_cand], "dr_lep_blep")
        pruned_ev = ak.with_field(pruned_ev, dr_lep_bhad[has_cand], "dr_lep_bhad")
        pruned_ev = ak.with_field(pruned_ev, dr_ja_jb   [has_cand], "dr_ja_jb")

        # attach tag-and-probe payload
        pruned_ev = ak.with_field(pruned_ev, self._year,     "run_year")
        pruned_ev = ak.with_field(pruned_ev, self._campaign, "run_campaign")
        pruned_ev = ak.with_field(pruned_ev, tnp_had_fill[has_cand], "tnp_had_fill")
        pruned_ev = ak.with_field(pruned_ev, tnp_lep_fill[has_cand], "tnp_lep_fill")
        pruned_ev = ak.with_field(pruned_ev, had_pt_for_tnp[has_cand], "tnp_had_pt")
        pruned_ev = ak.with_field(pruned_ev, lep_pt_for_tnp[has_cand],    "tnp_lep_pt")


        ####################
        #     Output       #
        ####################
        # If still empty (extremely rare), skip weights/histos
        if len(pruned_ev) == 0:
            return {dataset: output}

        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        systematics = [shift_name] if shift_name is not None else ["nominal"] + list(weights.variations)

        if not self.noHist:
            output = histo_writter(pruned_ev, output, weights, systematics, self.isSyst, self.SF_map)
        if self.isArray:
            array_writer(self, pruned_ev, events, weights, systematics, dataset, isRealData)

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
