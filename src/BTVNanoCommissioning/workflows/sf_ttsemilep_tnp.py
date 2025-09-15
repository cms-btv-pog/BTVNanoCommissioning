import collections
import numpy as np, awkward as ak
import hist as Hist
from coffea import processor

from BTVNanoCommissioning.utils.correction import (
    load_lumi, load_SF, weight_manager, common_shifts,
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


def select_lepton(events, channel, campaign, iso_mode="tight"):
    if channel == "mu":
        tight = (events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & mu_idiso(events, campaign)
        if iso_mode == "tight":
            mask = tight
        elif iso_mode == "sbiso":
            iso   = ak.fill_none(events.Muon.pfRelIso04_all, 999.)
            loose = (events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & ak.fill_none(events.Muon.tightId, False)
            mask  = loose & (~tight) & (iso > 0.15) & (iso < 0.40)
        return events.Muon[mask], mask

    elif channel == "el":
        tight = (events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4) & ele_mvatightid(events, campaign)
        if iso_mode == "tight":
            mask = tight
        elif iso_mode == "sbiso":
            loose = (events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)
            mva   = ak.fill_none(getattr(events.Electron, "mvaIso", ak.zeros_like(events.Electron.pt)), -99.)
            mask  = loose & (~tight) & (mva > 0.85) & (mva < 0.95)
        return events.Electron[mask], mask
    else:
        raise ValueError(channel)


def solve_nu_pz(px_l, py_l, pz_l, e_l, px_n, py_n, mW=80.4):
    muW = (mW**2)/2.0 + px_l*px_n + py_l*py_n
    a = (e_l**2 - pz_l**2)
    b = -2.0*muW*pz_l
    c = e_l**2*(px_n**2 + py_n**2) - muW**2
    disc = b**2 - 4.0*a*c
    disc_pos   = ak.where(disc > 0, disc, 0.0)
    sqrt_disc  = np.sqrt(disc_pos)
    pz1 = (-b + sqrt_disc)/(2.0*a)
    pz2 = (-b - sqrt_disc)/(2.0*a)
    e1  = np.sqrt(px_n**2 + py_n**2 + pz1**2)
    e2  = np.sqrt(px_n**2 + py_n**2 + pz2**2)
    return (pz1, e1), (pz2, e2)


def ttbar_reco(jets, lepton, met, maxjets=6, mW=80.4, mT=172.5,
               sig_w=30.0, sig_t=40.0):
    # limit jets
    jets = jets[:, :maxjets]
    idx  = ak.local_index(jets, axis=1)

    # build combos
    b_pairs = ak.cartesian({"bl": idx, "bh": idx}, axis=1, nested=False)
    b_pairs = b_pairs[b_pairs.bl != b_pairs.bh]
    w_pairs = ak.combinations(idx, 2, axis=1, replacement=False)
    comb = ak.cartesian({"b": b_pairs, "w": w_pairs}, axis=1, nested=False)
    distinct = ((comb.b.bl != comb.w["0"]) & (comb.b.bl != comb.w["1"]) &
                (comb.b.bh != comb.w["0"]) & (comb.b.bh != comb.w["1"]))
    comb = comb[distinct]
    has_cand = ak.num(comb.b.bl, axis=1) > 0
    if not ak.any(has_cand):
        return {"has_cand": has_cand}  # empty payload

    iBL, iBH, iJa, iJb = comb.b.bl, comb.b.bh, comb.w["0"], comb.w["1"]
    BL, BH, JA, JB = jets[iBL], jets[iBH], jets[iJa], jets[iJb]

    # broadcast lepton & MET
    lep = lepton
    lep_eta = ak.broadcast_arrays(lep.eta, BL.eta)[0]
    lep_phi = ak.broadcast_arrays(lep.phi, BL.phi)[0]
    lep_px  = ak.broadcast_arrays(lep.px , BL.px )[0]
    lep_py  = ak.broadcast_arrays(lep.py , BL.py )[0]
    lep_pz  = ak.broadcast_arrays(lep.pz , BL.pz )[0]
    lep_e   = ak.broadcast_arrays(lep.energy, BL.energy)[0]
    met_x   = ak.broadcast_arrays(met.x, BL.px)[0]
    met_y   = ak.broadcast_arrays(met.y, BL.px)[0]

    # lepton subtraction if ΔR<0.4
    dphi = ak.where((BL.phi - lep_phi) > np.pi, BL.phi - lep_phi - 2*np.pi,
                    ak.where((BL.phi - lep_phi) < -np.pi, BL.phi - lep_phi + 2*np.pi, BL.phi - lep_phi))
    dR = np.sqrt((BL.eta - lep_eta)**2 + dphi**2)
    overlap = dR < 0.4
    BL_px = ak.where(overlap, BL.px - lep_px, BL.px)
    BL_py = ak.where(overlap, BL.py - lep_py, BL.py)
    BL_pz = ak.where(overlap, BL.pz - lep_pz, BL.pz)
    BL_e  = ak.where(overlap, BL.energy - lep_e, BL.energy)

    # neutrino pz: choose root closer to mT
    (pz1, en1), (pz2, en2) = solve_nu_pz(lep_px, lep_py, lep_pz, lep_e, met_x, met_y, mW)
    def four(e, px, py, pz): return (e, px, py, pz)
    t1 = four(lep_e + en1 + BL_e, lep_px + met_x + BL_px, lep_py + met_y + BL_py, lep_pz + pz1 + BL_pz)
    t2 = four(lep_e + en2 + BL_e, lep_px + met_x + BL_px, lep_py + met_y + BL_py, lep_pz + pz2 + BL_pz)
    m  = lambda e,px,py,pz: np.sqrt(np.maximum(0.0, e*e - (px*px + py*py + pz*pz)))
    m_tlep1, m_tlep2 = m(*t1), m(*t2)
    choose_1 = (abs(m_tlep1 - mT) < abs(m_tlep2 - mT))
    pz_n = ak.where(choose_1, pz1, pz2)
    e_n  = ak.where(choose_1, en1, en2)
    tlep_e  = ak.where(choose_1, t1[0], t2[0])
    tlep_px = ak.where(choose_1, t1[1], t2[1])
    tlep_py = ak.where(choose_1, t1[2], t2[2])
    tlep_pz = ak.where(choose_1, t1[3], t2[3])
    m_tlep  = ak.where(choose_1, m_tlep1, m_tlep2)
    valid_dnu = (m_tlep > 100.0) & (m_tlep < 240.0)

    # had side
    W_px, W_py, W_pz = JA.px + JB.px, JA.py + JB.py, JA.pz + JB.pz
    W_e  = JA.energy + JB.energy
    mW_h = m(W_e, W_px, W_py, W_pz)
    th_px, th_py, th_pz = BH.px + W_px, BH.py + W_py, BH.pz + W_pz
    th_e  = BH.energy + W_e
    mT_h  = m(th_e, th_px, th_py, th_pz)

    chi2 = ((mW_h - mW)/sig_w)**2 + ((mT_h - mT)/sig_t)**2 + ((m_tlep - mT)/sig_t)**2
    chi2 = ak.where(valid_dnu, chi2, np.inf)
    best_idx = ak.where(has_cand, ak.argmin(chi2, axis=1), -1)
    li = ak.local_index(chi2, axis=1)
    best_mask = (li == best_idx) & (best_idx >= 0)

    # pick best objects
    best = {
        "has_cand": has_cand,
        "best_mask": best_mask,
        "BL": ak.firsts(BL[best_mask]),
        "BH": ak.firsts(BH[best_mask]),
        "JA": ak.firsts(JA[best_mask]),
        "JB": ak.firsts(JB[best_mask]),
        "nu": ak.zip({"x": ak.firsts(met_x[best_mask]), "y": ak.firsts(met_y[best_mask]),
                      "z": ak.firsts(pz_n[best_mask]), "t": ak.firsts(e_n[best_mask])},
                     with_name="FourVector"),
        "tlep": ak.zip({"x": ak.firsts(tlep_px[best_mask]), "y": ak.firsts(tlep_py[best_mask]),
                        "z": ak.firsts(tlep_pz[best_mask]), "t": ak.firsts(tlep_e[best_mask])},
                       with_name="FourVector"),
        "thad": ak.zip({"x": ak.firsts(th_px[best_mask]), "y": ak.firsts(th_py[best_mask]),
                        "z": ak.firsts(th_pz[best_mask]), "t": ak.firsts(th_e[best_mask])},
                       with_name="FourVector"),
        "chi2": ak.firsts(chi2[best_mask]),
    }
    return best


def tt_truth_category(best_BL, best_BH, is_mc):
    """Return an Awkward string array: 'sig' / 'bkg' (MC) or 'other' (data)."""
    n = len(best_BL)
    if not is_mc:
        return ak.Array(np.full(n, "other", dtype=np.str_))

    bl_b = ak.fill_none(ak.values_astype(getattr(best_BL, "hadronFlavour", 0) == 5, bool), False)
    bh_b = ak.fill_none(ak.values_astype(getattr(best_BH, "hadronFlavour", 0) == 5, bool), False)
    is_sig = ak.to_numpy(bl_b & bh_b)

    lab = np.where(is_sig, "sig", "bkg").astype(np.str_)
    return ak.Array(lab)



class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=25000,
        selectionModifier="mu",  # "mu" or "el"
        tag_tagger="UParTAK4",
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

        # -------------------- Common preselection --------------------
        req_lumi = np.ones(len(events), dtype=bool)
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        # Triggers
        triggers_mu = ["IsoMu24", "Mu50"]
        triggers_el = ["Ele32_WPTight_Gsf", "Ele50_CaloIdVT_GsfTrkIdT_PFJet165"]
        triggers = triggers_mu if self.channel == "mu" else triggers_el
        req_trig = HLT_helper(events, triggers)

        # MET filters
        req_metf = MET_filters(events, self._campaign)

        # Loose veto objects (for exactly-one-lepton selection and jet cleaning)
        mu_loose = (
            (events.Muon.pt > 15) & (abs(events.Muon.eta) < 2.4) &
            ak.fill_none(events.Muon.looseId, False)
        )
        el_loose = (
            (events.Electron.pt > 15) & (abs(events.Electron.eta) < 2.4) &
            ~((abs(events.Electron.eta) > 1.4442) & (abs(events.Electron.eta) < 1.566)) &
            (ak.fill_none(events.Electron.cutBased >= 2, False) |
             ak.fill_none(getattr(events.Electron, "mvaIso_WP90", False), False))
        )

        # Jet cleaning vs extra loose leptons
        def _clean_jets(ev, other_mu_mask, other_el_mask):
            dr_mu = ev.Jet.metric_table(ev.Muon[other_mu_mask])
            dr_el = ev.Jet.metric_table(ev.Electron[other_el_mask])
            all_true = ak.ones_like(ev.Jet.pt, dtype=bool)
            has_mu = ak.num(ev.Muon[other_mu_mask], axis=1) > 0
            has_el = ak.num(ev.Electron[other_el_mask], axis=1) > 0
            clean_mu = ak.where(has_mu, ak.all(dr_mu > 0.4, axis=-1, mask_identity=True), all_true)
            clean_el = ak.where(has_el, ak.all(dr_el > 0.4, axis=-1, mask_identity=True), all_true)
            base_jet_mask = jet_id(ev, self._campaign)
            return ak.fill_none(base_jet_mask & clean_mu & clean_el, False, axis=-1)

        # cutflow helper
        if not self.noHist and "cutflow" not in output:
            cf_axis = Hist.axis.StrCategory([], name="step", growth=True)
            output["cutflow"] = Hist.Hist(cf_axis, Hist.storage.Weight())
        def _cf(step, mask):
            if self.noHist: return
            if isinstance(mask, ak.Array): mask = ak.to_numpy(mask)
            output["cutflow"].fill(step, weight=float(np.asarray(mask, dtype=bool).sum()))
        _cf("all", np.ones(len(events), dtype=bool))
        _cf("lumi", req_lumi)
        _cf("trig", req_lumi & req_trig)
        _cf("metf", req_lumi & req_trig & req_metf)

        # -------------------- Build regions --------------------
        # There are five regions:
        # - central: tag b-jet working point "M", tight isolated lepton
        # - isolation side-band: tag b-jet working point "M", lepton fail tight isolation but pass loose isolation
        # - btagM sideband: tag b-jet working point fail "M" but pass "L", tight isolated lepton
        # - btagL sideband: tag b-jet working point fail "M" and fail "L", tight isolated lepton
        # - btagM + iso sideband: tag b-jet working point fail "M" but pass "L", lepton fail tight isolation but pass loose isolation

        region_specs = [
            # name, iso_mode, tag_btagWP
            ("central", "tight", "M"),  
            ("sb_iso", "sbiso", "M"),
            ("sb_btagM", "tight", "L"),
            ("sb_btagL", "tight", "No"),
            ("sb_iso_btagM", "sbiso", "L")
        ]

        for iso_mode in ["tight", "sbiso"]:
            # Select 1 lepton of requested isolation working point
            sel_leps, sel_mask = select_lepton(events, self.channel, self._campaign, iso_mode=iso_mode)


            # exactly one of those leptons, zero of the other flavor (veto)
            if self.channel == "mu":
                req_lepveto = (ak.num(events.Muon[mu_loose], axis=1) == 1) & (ak.num(events.Electron[el_loose], axis=1) == 0)
                other_mu = events.Muon[mu_loose & ~ak.fill_none(sel_mask, False)]
                other_el = events.Electron[el_loose]
            else:
                req_lepveto = (ak.num(events.Muon[mu_loose], axis=1) == 0) & (ak.num(events.Electron[el_loose], axis=1) == 1)
                other_mu = events.Muon[mu_loose]
                other_el = events.Electron[el_loose & ~ak.fill_none(sel_mask, False)]

            req_lep = ak.num(sel_leps, axis=1) == 1
            _cf(f"{iso_mode}:lepveto", req_lumi & req_trig & req_metf & req_lepveto)
            _cf(f"{iso_mode}:tightlep", req_lumi & req_trig & req_metf & req_lepveto & req_lep)

            # Jets (cleaned)
            jet_mask = _clean_jets(events, other_mu_mask=mu_loose & ~ak.fill_none((events.Muon.pt > 30), False),
                                            other_el_mask=el_loose & ~ak.fill_none((events.Electron.pt > 30), False))
            jets_all = events.Jet[jet_mask]
            req_jets = ak.num(jets_all, axis=1) == 4
            _cf(f"{iso_mode}:jets==4", req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets)

            # Base mask for this family (no btag condition yet)
            evmask_base = req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets
            if not ak.any(evmask_base):
                continue

            # Slice to base
            ev_base   = events[evmask_base]
            jets_base = jets_all[evmask_base]
            lep_base  = ak.firsts(sel_leps[evmask_base])

            # B-tag counts at L and M for event-level sideband logic
            bmask_L = btag_wp(jets_base, self._year, self._campaign, tagger=self.tag_tagger, borc="b", wp="L")
            bmask_M = btag_wp(jets_base, self._year, self._campaign, tagger=self.tag_tagger, borc="b", wp="M")
            nb_L = ak.sum(ak.fill_none(bmask_L, False), axis=1)
            nb_M = ak.sum(ak.fill_none(bmask_M, False), axis=1)

            # Central requires ≥1 M tag
            mask_central = (nb_M >= 1)

            # Sidebands in b-tag:  sb_btagM: (0 M) & (≥1 L), sb_btagL: (0 L)
            mask_sb_btagM = (nb_M == 0) & (nb_L >= 1)
            mask_sb_btagL = (nb_L == 0)

            # MET 4-vector
            met_b = ak.zip(
                {"x": ev_base.MET.pt * np.cos(ev_base.MET.phi),
                 "y": ev_base.MET.pt * np.sin(ev_base.MET.phi),
                 "z": ak.zeros_like(ev_base.MET.pt),
                 "t": ev_base.MET.pt},
                with_name="FourVector",
            )

            # ttbar reco on this family once
            best = ttbar_reco(jets_base, lep_base, met_b, maxjets=4)
            has_cand = best.get("has_cand", ak.Array([]))
            if not ak.any(has_cand):
                continue

            BH = best["BH"]
            BL = best["BL"]
            JA = best["JA"]
            JB = best["JB"]
            nu = best["nu"]
            tlep = best["tlep"]
            thad = best["thad"]
            chi2 = best["chi2"]

            had_tag = ak.fill_none(btag_wp(BH, self._year, self._campaign, tagger=self.tag_tagger, borc="b", wp="M"), False)
            lep_tag = ak.fill_none(btag_wp(BL, self._year, self._campaign, tagger=self.tag_tagger, borc="b", wp="M"), False)
            had_pt_ok = ak.fill_none(BH.pt >= 30.0, False)
            lep_pt_ok = ak.fill_none(BL.pt >= 30.0, False)

            def _four_mass(v):
                arg = v.t*v.t - (v.x*v.x + v.y*v.y + v.z*v.z)
                return np.sqrt(ak.where(arg > 0, arg, 0.0))

            def _four_pt(v):
                arg = v.x*v.x + v.y*v.y
                return np.sqrt(ak.where(arg > 0, arg, 0.0))

            # build pruned for each region
            for rname, riso_mode, rtag_btagwp in region_specs:
                if iso_mode != riso_mode:
                    continue

                if rtag_btagwp == "M":
                    rmask = has_cand & mask_central
                elif rtag_btagwp == "L":
                    rmask = has_cand & mask_sb_btagM
                else:
                    rmask = has_cand & mask_sb_btagL

                rmask_np = ak.to_numpy(ak.fill_none(rmask, False))

                # TnP fills (gate on opposite-side tag only for central)
                require_tag = (rname == "central")
                ones = ak.ones_like(had_pt_ok, dtype=bool)
                tnp_had_fill = ak.to_numpy(had_pt_ok & (lep_tag if require_tag else ones))[rmask_np]
                tnp_lep_fill = ak.to_numpy(lep_pt_ok & (had_tag if require_tag else ones))[rmask_np]
                tnp_had_pt   = ak.to_numpy(BH.pt)[rmask_np]
                tnp_lep_pt   = ak.to_numpy(BL.pt)[rmask_np]

                # per-region slices
                ev_r   = ev_base[rmask_np]
                jets_r = jets_base[rmask_np]
                lep_r  = lep_base[rmask_np]
                BL_r   = BL[rmask_np]
                BH_r = BH[rmask_np]
                JA_r   = JA[rmask_np]
                JB_r = JB[rmask_np]
                nu_r   = nu[rmask_np]
                tl_r = tlep[rmask_np]
                th_r = thad[rmask_np]
                chi2_r = ak.to_numpy(chi2)[rmask_np]

                # compute metric for quality of kinematic fit
                met_r = ak.to_numpy(ev_r.MET.pt)
                met_edges = np.array([0., 40., 80., 200.], dtype=float)
                prob_edges = np.arange(11., 21., 1.0)  # 11,12,...,20  => 9 bins
                # proper implementation:
                #    from scipy.special import gammaincc
                #    p = np.clip(gammaincc(1.5, ak.to_numpy(chi2_r)/2.0), 1e-300, 1.0)
                #    q = -np.log10(p)
                # fast workaround:
                q = ak.to_numpy(chi2_r) / (2.0*np.log(10.0))
                met_bin = np.clip(np.digitize(met_r, met_edges) - 1, 0, len(met_edges)-2)
                prob_bin = np.clip(np.digitize(q,      prob_edges) - 1, 0, len(prob_edges) - 2)
                NQ = len(prob_edges) - 1
                kinbin = (met_bin * NQ + prob_bin).astype(np.int32)

                # masses, pts, ΔR
                W_r = ak.zip({"x": JA_r.x + JB_r.x, "y": JA_r.y + JB_r.y,
                              "z": JA_r.z + JB_r.z, "t": JA_r.t + JB_r.t},
                             with_name="FourVector")
                tlep_mass = _four_mass(tl_r)
                thad_mass = _four_mass(th_r)
                whad_mass = _four_mass(W_r)
                tlep_pt   = _four_pt(tl_r)
                thad_pt   = _four_pt(th_r)

                def _dR(a, b):
                    dphi = ak.where((a.phi - b.phi) >  np.pi, a.phi - b.phi - 2*np.pi,
                        ak.where((a.phi - b.phi) < -np.pi, a.phi - b.phi + 2*np.pi, a.phi - b.phi))
                    arg = (a.eta - b.eta)**2 + dphi**2
                    return np.sqrt(ak.where(arg > 0, arg, 0.0))

                dr_lep_blep = _dR(lep_r, BL_r)
                dr_lep_bhad = _dR(lep_r, BH_r)
                dr_ja_jb    = _dR(JA_r,  JB_r)

                # TT truth category (fast proxy): 'sig' if both chosen jets are b-hadron jets
                tt_cat = tt_truth_category(BL_r, BH_r, is_mc=not isRealData)

                # assemble pruned view
                pr = ev_r
                pr = ak.with_field(pr, BL_r,              "blep")
                pr = ak.with_field(pr, BH_r,              "bhad")
                pr = ak.with_field(pr, JA_r,              "ja")
                pr = ak.with_field(pr, JB_r,              "jb")
                pr = ak.with_field(pr, nu_r,              "nu")
                pr = ak.with_field(pr, tl_r,              "tlep")
                pr = ak.with_field(pr, th_r,              "thad")
                pr = ak.with_field(pr, chi2_r,            "chi2")
                pr = ak.with_field(pr, self._year,        "run_year")
                pr = ak.with_field(pr, self._campaign,    "run_campaign")
                pr = ak.with_field(pr, jets_r[:, :4],"SelJet")
                if self.channel == "mu":
                    pr = ak.with_field(pr, lep_r,     "SelMuon")
                else:
                    pr = ak.with_field(pr, lep_r,     "SelElectron")
                pr = ak.with_field(pr, ak.num(pr.SelJet, axis=1), "njet")

                pr = ak.with_field(pr, tlep_mass,         "tlep_mass")
                pr = ak.with_field(pr, thad_mass,         "thad_mass")
                pr = ak.with_field(pr, whad_mass,         "whad_mass")
                pr = ak.with_field(pr, tlep_pt,           "tlep_pt")
                pr = ak.with_field(pr, thad_pt,           "thad_pt")
                pr = ak.with_field(pr, dr_lep_blep,       "dr_lep_blep")
                pr = ak.with_field(pr, dr_lep_bhad,       "dr_lep_bhad")
                pr = ak.with_field(pr, dr_ja_jb,          "dr_ja_jb")

                # tag-and-probe payload
                pr = ak.with_field(pr, tnp_had_fill, "tnp_had_fill")
                pr = ak.with_field(pr, tnp_lep_fill, "tnp_lep_fill")
                pr = ak.with_field(pr, tnp_had_pt,   "tnp_had_pt")
                pr = ak.with_field(pr, tnp_lep_pt,   "tnp_lep_pt")

                # region & tt category labels for histogramming
                pr = ak.with_field(pr, np.full(len(pr), rname, dtype=np.str_), "tnp_region")
                pr = ak.with_field(pr, tt_cat,                                    "tt_cat")
                pr = ak.with_field(pr, ak.Array(kinbin), "kinbin")

                # -------------------- Output --------------------
                weights = weight_manager(pr, self.SF_map, self.isSyst)
                systematics = [shift_name] if shift_name is not None else ["nominal"] + list(weights.variations)

                if not self.noHist:
                    output = histo_writter(pr, output, weights, systematics, self.isSyst, self.SF_map)

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
