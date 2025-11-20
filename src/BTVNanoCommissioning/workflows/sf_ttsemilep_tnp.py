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
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    MET_filters,
    mu_idiso,
    ele_mvatightid,
    btag_wp,
    wp_dict,
)

from BTVNanoCommissioning.helpers.definitions import get_discriminators, get_definitions


def select_lepton(events, channel, campaign, iso_mode="tight"):
    if channel == "mu":
        tight = (
            (events.Muon.pt > 30)
            & (abs(events.Muon.eta) < 2.4)
            & mu_idiso(events, campaign)
        )
        if iso_mode == "tight":
            mask = tight
        elif iso_mode == "sbiso":
            iso = ak.fill_none(events.Muon.pfRelIso04_all, 999.0)
            loose = (
                (events.Muon.pt > 30)
                & (abs(events.Muon.eta) < 2.4)
                & ak.fill_none(events.Muon.tightId, False)
            )
            mask = loose & (~tight) & (iso > 0.15) & (iso < 0.40)
        return events.Muon[mask], mask

    elif channel == "el":
        tight = (
            (events.Electron.pt > 30)
            & (abs(events.Electron.eta) < 2.4)
            & ele_mvatightid(events, campaign)
        )
        if iso_mode == "tight":
            mask = tight
        elif iso_mode == "sbiso":
            loose = (events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)
            mva = ak.fill_none(
                getattr(events.Electron, "mvaIso", ak.zeros_like(events.Electron.pt)),
                -99.0,
            )
            mask = loose & (~tight) & (mva > 0.85) & (mva < 0.95)
        return events.Electron[mask], mask
    else:
        raise ValueError(channel)


def solve_nu_pz(px_l, py_l, pz_l, e_l, px_n, py_n, mW=80.4):
    muW = (mW**2) / 2.0 + px_l * px_n + py_l * py_n
    a = e_l**2 - pz_l**2
    b = -2.0 * muW * pz_l
    c = e_l**2 * (px_n**2 + py_n**2) - muW**2
    disc = b**2 - 4.0 * a * c
    disc_pos = ak.where(disc > 0, disc, 0.0)
    sqrt_disc = np.sqrt(disc_pos)
    pz1 = (-b + sqrt_disc) / (2.0 * a)
    pz2 = (-b - sqrt_disc) / (2.0 * a)
    e1 = np.sqrt(px_n**2 + py_n**2 + pz1**2)
    e2 = np.sqrt(px_n**2 + py_n**2 + pz2**2)
    return (pz1, e1), (pz2, e2)
def gaussian(m, mu = 0, sigma = 1):

    return (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5 * ((m-mu)/sigma)**2)

def calculate_mass_probability(m_type, masses):
    
    #PLACEHOLDER FIXME assume probabilities as gaussians
    #Where are these values from?
    mW = 80.4
    mT = 172.5
    sig_w = 30.0
    sig_t = 40.0
    if m_type == "W, T":
        #FIXME these are not independent masses
        return gaussian(masses[0],mu = mW, sigma = sig_w) * gaussian(masses[1], mu = mT, sigma = sig_t)
    elif m_type == "T":
        return gaussian(masses[0], mu = mT, sigma = sig_t)
    else:
        raise ExceptionType("Not a valid mass distribution")
def ttbar_reco(jets, lepton, met, maxjets=6, mW=80.4, mT=172.5, sig_w=30.0, sig_t=40.0):
    # limit jets
    jets = jets[:, :maxjets]
    idx = ak.local_index(jets, axis=1)

    # build combos
    b_pairs = ak.cartesian({"bl": idx, "bh": idx}, axis=1, nested=False)
    b_pairs = b_pairs[b_pairs.bl != b_pairs.bh]
    w_pairs = ak.combinations(idx, 2, axis=1, replacement=False)
    comb = ak.cartesian({"b": b_pairs, "w": w_pairs}, axis=1, nested=False)
    distinct = (
        (comb.b.bl != comb.w["0"])
        & (comb.b.bl != comb.w["1"])
        & (comb.b.bh != comb.w["0"])
        & (comb.b.bh != comb.w["1"])
    )
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
    lep_px = ak.broadcast_arrays(lep.px, BL.px)[0]
    lep_py = ak.broadcast_arrays(lep.py, BL.py)[0]
    lep_pz = ak.broadcast_arrays(lep.pz, BL.pz)[0]
    lep_pt = np.sqrt(lep_px**2 + lep_py**2)
    lep_e = ak.broadcast_arrays(lep.energy, BL.energy)[0]
    met_x = ak.broadcast_arrays(met.x, BL.px)[0]
    met_y = ak.broadcast_arrays(met.y, BL.px)[0]
    met_phi = ak.broadcast_arrays(np.arctan(met.y/met.x), BL.phi)[0]
    met_pt = np.sqrt(met_x**2 + met_y**2)
    # lepton subtraction if ΔR<0.4
    dphi = ak.where(
        (BL.phi - lep_phi) > np.pi,
        BL.phi - lep_phi - 2 * np.pi,
        ak.where(
            (BL.phi - lep_phi) < -np.pi, BL.phi - lep_phi + 2 * np.pi, BL.phi - lep_phi
        ),
    )
    dR = np.sqrt((BL.eta - lep_eta) ** 2 + dphi**2)
    overlap = dR < 0.4
    BL_px = ak.where(overlap, BL.px - lep_px, BL.px)
    BL_py = ak.where(overlap, BL.py - lep_py, BL.py)
    BL_pz = ak.where(overlap, BL.pz - lep_pz, BL.pz)
    BL_e = ak.where(overlap, BL.energy - lep_e, BL.energy)

    # neutrino pz: choose root closer to mT
    (pz1, en1), (pz2, en2) = solve_nu_pz(
        lep_px, lep_py, lep_pz, lep_e, met_x, met_y, mW
    )

    def four(e, px, py, pz):
        return (e, px, py, pz)

    t1 = four(
        lep_e + en1 + BL_e,
        lep_px + met_x + BL_px,
        lep_py + met_y + BL_py,
        lep_pz + pz1 + BL_pz,
    )
    t2 = four(
        lep_e + en2 + BL_e,
        lep_px + met_x + BL_px,
        lep_py + met_y + BL_py,
        lep_pz + pz2 + BL_pz,
    )
    m = lambda e, px, py, pz: np.sqrt(
        np.maximum(0.0, e * e - (px * px + py * py + pz * pz))
    )
    m_tlep1, m_tlep2 = m(*t1), m(*t2)
    choose_1 = abs(m_tlep1 - mT) < abs(m_tlep2 - mT)
    pz_n = ak.where(choose_1, pz1, pz2)
    e_n = ak.where(choose_1, en1, en2)
    tlep_e = ak.where(choose_1, t1[0], t2[0])
    tlep_px = ak.where(choose_1, t1[1], t2[1])
    tlep_py = ak.where(choose_1, t1[2], t2[2])
    tlep_pz = ak.where(choose_1, t1[3], t2[3])
    m_tlep = ak.where(choose_1, m_tlep1, m_tlep2)
    valid_dnu = (m_tlep > 100.0) & (m_tlep < 240.0)
    # had side
    W_px, W_py, W_pz = JA.px + JB.px, JA.py + JB.py, JA.pz + JB.pz
    W_e = JA.energy + JB.energy
    mW_h = m(W_e, W_px, W_py, W_pz)

    #calculate and store transverse W mass for fit variable bins
    delta_phi = (lep_phi - met_phi) % (2 * np.pi) - np.pi
    mTW_l = np.sqrt((2 * lep_pt * met_pt) * (1 - np.cos(delta_phi)))

    th_px, th_py, th_pz = BH.px + W_px, BH.py + W_py, BH.pz + W_pz
    th_e = BH.energy + W_e
    mT_h = m(th_e, th_px, th_py, th_pz)

    #FIXME replace chi2 with lambda calculation
    P_m2_m3 = calculate_mass_probability("W, T", [mW_h, mT_h])
    P_m_t1 = calculate_mass_probability("T", [m_tlep])

    #chi2 = (
    #    ((mW_h - mW) / sig_w) ** 2
    #    + ((mT_h - mT) / sig_t) ** 2
    #    + ((m_tlep - mT) / sig_t) ** 2
    #)
    #chi2 = ak.where(valid_dnu, chi2, np.inf)


    neg_log_lambda = ak.where(valid_dnu, -np.log(P_m2_m3 * P_m_t1), np.inf)
    #best_idx = ak.where(has_cand, ak.argmin(chi2, axis=1), -1)
    best_idx = ak.where(has_cand, ak.argmin(neg_log_lambda, axis=1), -1)
    li = ak.local_index(neg_log_lambda, axis=1)
    best_mask = (li == best_idx) & (best_idx >= 0)

    # pick best objects
    best = {
        "has_cand": has_cand,
        "best_mask": best_mask,
        "BL": ak.firsts(BL[best_mask]),
        "BH": ak.firsts(BH[best_mask]),
        "JA": ak.firsts(JA[best_mask]),
        "JB": ak.firsts(JB[best_mask]),
        "nu": ak.zip(
            {
                "x": ak.firsts(met_x[best_mask]),
                "y": ak.firsts(met_y[best_mask]),
                "z": ak.firsts(pz_n[best_mask]),
                "t": ak.firsts(e_n[best_mask]),
            },
            with_name="FourVector",
        ),
        "tlep": ak.zip(
            {
                "x": ak.firsts(tlep_px[best_mask]),
                "y": ak.firsts(tlep_py[best_mask]),
                "z": ak.firsts(tlep_pz[best_mask]),
                "t": ak.firsts(tlep_e[best_mask]),
            },
            with_name="FourVector",
        ),
        "thad": ak.zip(
            {
                "x": ak.firsts(th_px[best_mask]),
                "y": ak.firsts(th_py[best_mask]),
                "z": ak.firsts(th_pz[best_mask]),
                "t": ak.firsts(th_e[best_mask]),
            },
            with_name="FourVector",
        ),
        #"chi2": ak.firsts(chi2[best_mask]),
        "neg_log_lambda": ak.firsts(neg_log_lambda[best_mask]),
        "mTW_l" : ak.firsts(mTW_l[best_mask]),
    }
    return best


def tt_truth_category(best_BL, best_BH, is_mc):
    """Return an Awkward string array: 'sig' / 'bkg' (MC) or 'other' (data)."""
    n = len(best_BL)
    if not is_mc:
        return ak.Array(np.full(n, "other", dtype="U5"))

    bl_b = ak.fill_none(
        ak.values_astype(getattr(best_BL, "hadronFlavour", 0) == 5, bool), False
    )
    bh_b = ak.fill_none(
        ak.values_astype(getattr(best_BH, "hadronFlavour", 0) == 5, bool), False
    )
    is_sig = ak.to_numpy(bl_b & bh_b)

    lab = np.where(is_sig, "sig", "bkg").astype("U5")
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
        chunksize=10000,
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
        self._regions = ["central", "sbiso", "sbbtagM", "sbbtagL", "sbisobtagM"]

    @property
    def accumulator(self):
        return self._accumulator

    def define_histograms(self, events):
        """
        Define histograms to be written out by workflow
        """
        _hist_dict = {}

        ## Common axes
        flav_axis = Hist.axis.IntCategory(
            [0, 1, 4, 5, 6], name="flav", label="Genflavour"
        )
        syst_axis = Hist.axis.StrCategory([], name="syst", growth=True)
        pt_axis = Hist.axis.Regular(60, 0, 300, name="pt", label=" $p_{T}$ [GeV]")
        mass_axis = Hist.axis.Regular(50, 0, 300, name="mass", label=" $p_{T}$ [GeV]")
        eta_axis = Hist.axis.Regular(25, -2.5, 2.5, name="eta", label=" $\eta$")
        phi_axis = Hist.axis.Regular(30, -3, 3, name="phi", label="$\phi$")
        mt_axis = Hist.axis.Regular(30, 0, 300, name="mt", label=" $m_{T}$ [GeV]")

        mTW_l_axis = Hist.axis.Regular(40, 0, 200, name="mTW_l", label=r"$m_T(W_l)$")
        neg_log_lambda_axis = Hist.axis.Regular(40, 9, 30, name="neg_log_lambda", label=r"$-log(\\lambda)$")

        #chi2_axis = Hist.axis.Regular(60, 0, 60, name="chi2", label=r"$\chi^2$")
        tpt_axis = Hist.axis.Regular(60, 0, 600, name="pt", label=r"$p_T$ [GeV]")
        dr_axis = Hist.axis.Regular(20, 0, 8, name="dr", label="$\Delta$R")
        n_axis = Hist.axis.Integer(0, 10, name="n", label="N obj")

        # TnP-specific axes
        cat_axis = Hist.axis.StrCategory(["had", "lep"], name="cat")
        wp_axis = Hist.axis.StrCategory(["L", "M", "T", "XT", "XXT"], name="wp")
        kin_axis = Hist.axis.Regular(24, 0, 24, name="kinbin", label="Bin: $M_T(W_h)$ v. $-log(\\lambda)$")
        ttcat_axis = Hist.axis.StrCategory(["sig", "bkg", "other"], name="tt_cat")
        result_axis = Hist.axis.StrCategory(["pass", "fail"], name="result")
        ptb_edges = [30.0, 50.0, 80.0, 120.0, 200.0, 400.0, 1000.0]
        ptb_axis = Hist.axis.Variable(ptb_edges, name="ptb", label="$p_{T}(b)$ [GeV]")

        # objects for common kinematics
        obj_list = ["MET", "mu"]
        for i in range(4):
            obj_list.append(f"jet{i}")

        # Create histograms for each region
        for region in self._regions:
            # Basic kinematics
            _hist_dict[f"{region}_njet"] = Hist.Hist(
                syst_axis, n_axis, Hist.storage.Weight()
            )

            # Bookkeeping / categories
            #_hist_dict[f"{region}_chi2"] = Hist.Hist(
            #    syst_axis, chi2_axis, Hist.storage.Weight()
            #)

            # ttbar reconstruction summaries
            _hist_dict[f"{region}_tlep_mass"] = Hist.Hist(
                syst_axis, mt_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_neg_log_lambda"] = Hist.Hist(
                syst_axis, neg_log_lambda_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_mTW_l"] = Hist.Hist(
                syst_axis, mTW_l_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_kinbin"] = Hist.Hist(
                syst_axis, kin_axis, Hist.storage.Weight()
            )

            _hist_dict[f"{region}_thad_mass"] = Hist.Hist(
                syst_axis, mt_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_whad_mass"] = Hist.Hist(
                syst_axis, mt_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_tlep_pt"] = Hist.Hist(
                syst_axis, tpt_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_thad_pt"] = Hist.Hist(
                syst_axis, tpt_axis, Hist.storage.Weight()
            )

            # Angular variables
            _hist_dict[f"{region}_dr_lep_blep"] = Hist.Hist(
                syst_axis, dr_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_dr_lep_bhad"] = Hist.Hist(
                syst_axis, dr_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{region}_dr_ja_jb"] = Hist.Hist(
                syst_axis, dr_axis, Hist.storage.Weight()
            )

            # Basic kinematics for objects
            for obj in obj_list:
                if "jet" in obj:
                    _hist_dict[f"{region}_{obj}_pt"] = Hist.Hist(
                        syst_axis, flav_axis, pt_axis, Hist.storage.Weight()
                    )
                    _hist_dict[f"{region}_{obj}_eta"] = Hist.Hist(
                        syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
                    )
                    _hist_dict[f"{region}_{obj}_phi"] = Hist.Hist(
                        syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
                    )
                    _hist_dict[f"{region}_{obj}_mass"] = Hist.Hist(
                        syst_axis, flav_axis, mass_axis, Hist.storage.Weight()
                    )
                else:
                    _hist_dict[f"{region}_{obj}_pt"] = Hist.Hist(
                        syst_axis, pt_axis, Hist.storage.Weight()
                    )
                    _hist_dict[f"{region}_{obj}_phi"] = Hist.Hist(
                        syst_axis, phi_axis, Hist.storage.Weight()
                    )
                    if obj != "MET":
                        _hist_dict[f"{region}_{obj}_eta"] = Hist.Hist(
                            syst_axis, eta_axis, Hist.storage.Weight()
                        )

            # B-tagger scores
            for disc in get_discriminators():
                if disc not in events.Jet.fields:
                    continue
                njet = 1
                for i in range(njet):
                    if "btag" in disc or "ProbaN" == disc:
                        _hist_dict[f"{region}_{disc}_{i}"] = Hist.Hist(
                            syst_axis,
                            flav_axis,
                            Hist.axis.Regular(50, 0.0, 1, name="discr", label=disc),
                            Hist.storage.Weight(),
                        )

            # TnP yields per region per tagger
            for tagger in self._all_taggers:
                _hist_dict[f"{region}_{tagger}_tnp_yields"] = Hist.Hist(
                    syst_axis,
                    cat_axis,
                    wp_axis,
                    result_axis,
                    ttcat_axis,
                    kin_axis,
                    ptb_axis,
                    Hist.storage.Weight(),
                )

        return _hist_dict

    def write_histograms(self, pruned_ev, output, weights, systematics, isSyst, SF_map):
        """
        Write histograms to the output dictionary based on pruned events and other parameters.
        Fills each region's histograms using a boolean mask built from pruned_ev.tnp_region.
        """

        # exclude b tag SFs for btag inputs
        exclude_btv = [
            v
            for v in weights.variations
            if any(
                k in v.upper()
                for k in ("DEEP", "PNET", "ROBUST", "UPART", "BTV", "BTAG", "CTAG")
            )
        ]

        nj = 4
        pruned_ev.SelJet = pruned_ev.SelJet[:, :nj]

        # data or MC
        if "hadronFlavour" in pruned_ev.SelJet.fields:
            isRealData = False
            genflavor = ak.values_astype(
                pruned_ev.SelJet.hadronFlavour
                + 1
                * (
                    (pruned_ev.SelJet.partonFlavour == 0)
                    & (pruned_ev.SelJet.hadronFlavour == 0)
                ),
                int,
            )
        else:
            isRealData = True
            genflavor = ak.zeros_like(pruned_ev.SelJet.pt, dtype=int)

        # region labels from the pruned view
        # every event in pruned_ev belongs to exactly one region
        if "tnp_region" in pruned_ev.fields:
            region_labels = np.asarray(
                ak.to_numpy(pruned_ev.tnp_region), dtype="U20"
            )  # Up to 20 chars
        else:
            # if you ever call this before assigning the region field, treat all as "central"
            region_labels = np.full(
                len(pruned_ev), "central", dtype="U20"
            )  # Up to 20 chars

        # Loop over the systematic variations
        for syst in systematics:
            if isSyst is False and syst != "nominal":
                continue

            # pick weights for this syst
            evt_w = (
                weights.weight()
                if syst == "nominal" or syst not in list(weights.variations)
                else weights.weight(modifier=syst)
            )
            # a version that excludes BTV like systematics for b tag score plots
            exclude_list = [k for k in exclude_btv if k in weights.variations]
            if exclude_list:
                evt_w_excl_btv = weights.partial_weight(exclude=exclude_list)
            else:
                evt_w_excl_btv = evt_w  # Nothing to exclude, use full weight

            # fill each histogram, but only with events that match its region prefix
            for histname, h in output.items():
                region_prefix = histname.split("_", 1)[0]
                if region_prefix not in self._regions:
                    continue

                # build the mask for this region, slice all inputs with it
                rmask = region_labels == region_prefix
                if not np.any(rmask):
                    continue

                # light aliases for sliced things
                ev = pruned_ev[rmask]
                w = evt_w[rmask]
                w_excl_btv = evt_w_excl_btv[rmask]

                # Selected electron histograms
                if (
                    "ele" in ev.fields
                    and ("ele_" in histname)
                    and (histname.replace(f"{region_prefix}_ele_", "") in ev.ele.fields)
                ):
                    fld = histname.replace(f"{region_prefix}_ele_", "")
                    h.fill(syst, ak.to_numpy(ev.ele[fld]), weight=w)
                    continue

                # Selected muon histograms
                if (
                    "mu" in ev.fields
                    and ("mu_" in histname)
                    and (histname.replace(f"{region_prefix}_mu_", "") in ev.mu.fields)
                ):
                    fld = histname.replace(f"{region_prefix}_mu_", "")
                    h.fill(syst, ak.to_numpy(ev.mu[fld]), weight=w)
                    continue

                # njet
                if histname == f"{region_prefix}_njet":
                    h.fill(syst, ak.to_numpy(ev.njet), weight=w)
                    continue

                # Jet kinematics and flavours
                if "jet" in histname:
                    for i in range(nj):
                        if f"jet{i}_" not in histname:
                            continue
                        # e.g. "<region>_jet0_pt" gives fld "pt"
                        fld = histname.replace(f"{region_prefix}_jet{i}_", "")
                        if fld not in ev.SelJet.fields:
                            continue
                        sel_jet = ev.SelJet[:, i]
                        flav = (
                            genflavor[rmask][:, i]
                            if genflavor.ndim == 2
                            else genflavor[rmask]
                        )
                        h.fill(
                            syst,
                            ak.to_numpy(flav),
                            ak.to_numpy(sel_jet[fld]),
                            weight=w,
                        )
                    continue

                # b tag discriminants, filled with BTV excluded weights
                if any(k in histname for k in ("btag", "PNet", "ProbaN")):
                    # hist names like "<region>_btagDeepFlavB_0"
                    # grab the jet index suffix
                    idx_str = histname.rsplit("_", 1)[-1]
                    if not idx_str.isdigit():
                        continue
                    i = int(idx_str)
                    if i >= nj:
                        continue
                    disc_name = histname.replace(f"{region_prefix}_", "").rsplit(
                        "_", 1
                    )[0]
                    if disc_name not in ev.SelJet.fields:
                        continue
                    seljet = ev.SelJet[:, i]
                    flav = (
                        genflavor[rmask][:, i]
                        if genflavor.ndim == 2
                        else genflavor[rmask]
                    )
                    h.fill(
                        syst=syst,
                        flav=ak.to_numpy(flav),
                        discr=ak.to_numpy(seljet[disc_name]),
                        weight=w_excl_btv,
                    )
                    continue

                # TnP yields, per region and per tagger

                # keys look like "<region>_<TAGGER>_tnp_yields"
                if histname.endswith("_tnp_yields"):
                    # parse the tagger for sanity if you need it
                    # tagger = histname.replace(f"{region_prefix}_", "").replace("_tnp_yields", "")
                    needed = {
                        "tnp_had_fill",
                        "tnp_lep_fill",
                        "tnp_had_pt",
                        "tnp_lep_pt",
                        "bhad",
                        "blep",
                        "tt_cat",
                        "kinbin",
                    }
                    if not needed.issubset(set(ev.fields)):
                        continue

                    # build available (tagger, wp, thr) from current run info
                    WPD = wp_dict(self._year, self._campaign)

                    def _has_score(obj, tagger):
                        return hasattr(obj, f"btag{tagger}B")

                    available = []
                    for tagger, sub in WPD.items():
                        if not (
                            _has_score(ev.bhad, tagger) and _has_score(ev.blep, tagger)
                        ):
                            continue
                        for wp_name, thr in sub.get("b", {}).items():
                            if wp_name == "No":
                                continue
                            available.append((tagger, wp_name, float(thr)))

                    ttcat = np.asarray(ak.to_numpy(ev.tt_cat), dtype="U5")
                    kinbin = np.asarray(ak.to_numpy(ev.kinbin), dtype=np.int32)

                    # had side
                    had_fill = np.asarray(
                        ak.to_numpy(ak.fill_none(ev.tnp_had_fill, False)), dtype=bool
                    )
                    had_pt = np.asarray(ak.to_numpy(ev.tnp_had_pt), dtype=float)
                    bhad = ev.bhad

                    if had_fill.any():
                        sel = had_fill
                        nsel = int(sel.sum())
                        for tagger, wp_name, thr in available:
                            scores = getattr(bhad, f"btag{tagger}B")
                            tagbit = np.asarray(ak.to_numpy(scores > thr), dtype=bool)
                            if nsel == 0:
                                continue
                            h.fill(
                                syst=syst,
                                cat=np.full(nsel, "had", dtype="U3"),
                                wp=np.full(nsel, wp_name, dtype="U3"),
                                result=np.where(tagbit[sel], "pass", "fail").astype(
                                    "U4"
                                ),
                                tt_cat=ttcat[sel],
                                kinbin=kinbin[sel],
                                ptb=had_pt[sel],
                                weight=evt_w[rmask][sel],
                            )

                    # lep side
                    lep_fill = np.asarray(
                        ak.to_numpy(ak.fill_none(ev.tnp_lep_fill, False)), dtype=bool
                    )
                    lep_pt = np.asarray(ak.to_numpy(ev.tnp_lep_pt), dtype=float)
                    blep = ev.blep

                    if lep_fill.any():
                        sel = lep_fill
                        nsel = int(sel.sum())
                        for tagger, wp_name, thr in available:
                            scores = getattr(blep, f"btag{tagger}B")
                            tagbit = np.asarray(ak.to_numpy(scores > thr), dtype=bool)
                            if nsel == 0:
                                continue
                            h.fill(
                                syst=syst,
                                cat=np.full(nsel, "lep", dtype="U3"),
                                wp=np.full(nsel, wp_name, dtype="U3"),
                                result=np.where(tagbit[sel], "pass", "fail").astype(
                                    "U4"
                                ),
                                tt_cat=ttcat[sel],
                                kinbin=kinbin[sel],
                                ptb=lep_pt[sel],
                                weight=evt_w[rmask][sel],
                            )
                    continue

                # chi2 and ttbar reco summaries
                base_name = histname.replace(f"{region_prefix}_", "")
                if base_name in (
                    #"chi2",
                    "neg_log_lambda",
                    "mTW_l",
                    "kinbin",
                    "tlep_mass",
                    "thad_mass",
                    "whad_mass",
                    "tlep_pt",
                    "thad_pt",
                    "dr_lep_blep",
                    "dr_lep_bhad",
                    "dr_ja_jb",
                ):
                    if base_name in ev.fields:
                        h.fill(syst, ak.to_numpy(ev[base_name]), weight=w)
                    continue

                # MET, use the region scoped keys
                if histname == f"{region_prefix}_MET_pt":
                    h.fill(syst, ak.to_numpy(ev.MET.pt), weight=w)
                    continue
                if histname == f"{region_prefix}_MET_phi":
                    h.fill(syst, ak.to_numpy(ev.MET.phi), weight=w)
                    continue

        return output

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
        output = {} if self.noHist else self.define_histograms(events)
        # print(f"=== process_shift: {dataset}, shift={shift_name}, isData={isRealData}, n={len(events)} ===")

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

        # Loose veto objects
        mu_loose = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < 2.4)
            & mu_idiso(events, self._campaign)
        )
        el_loose = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.4)
            & ele_mvatightid(events, self._campaign)
        )

        # Jet cleaning
        def _clean_jets(ev, other_mu_mask, other_el_mask):
            dr_mu = ev.Jet.metric_table(ev.Muon[other_mu_mask])
            dr_el = ev.Jet.metric_table(ev.Electron[other_el_mask])
            all_true = ak.ones_like(ev.Jet.pt, dtype=bool)
            has_mu = ak.num(ev.Muon[other_mu_mask], axis=1) > 0
            has_el = ak.num(ev.Electron[other_el_mask], axis=1) > 0
            clean_mu = ak.where(
                has_mu, ak.all(dr_mu > 0.4, axis=-1, mask_identity=True), all_true
            )
            clean_el = ak.where(
                has_el, ak.all(dr_el > 0.4, axis=-1, mask_identity=True), all_true
            )
            base_jet_mask = jet_id(ev, self._campaign, max_eta=2.4, min_pt=25)
            return ak.fill_none(base_jet_mask & clean_mu & clean_el, False, axis=-1)

        # Cutflow helper
        if not self.noHist and "cutflow" not in output:
            cf_axis = Hist.axis.StrCategory([], name="step", growth=True)
            output["cutflow"] = Hist.Hist(cf_axis, Hist.storage.Weight())

        def _cf(step, mask):
            if self.noHist:
                return
            if isinstance(mask, ak.Array):
                mask = ak.to_numpy(mask)
            output["cutflow"].fill(
                step, weight=float(np.asarray(mask, dtype=bool).sum())
            )

        _cf("all", np.ones(len(events), dtype=bool))
        _cf("lumi", req_lumi)
        _cf("trig", req_lumi & req_trig)
        _cf("metf", req_lumi & req_trig & req_metf)

        # Region specs
        region_specs = [
            ("central", "tight", "M"),
            ("sbiso", "sbiso", "M"),
            ("sbbtagM", "tight", "L"),
            ("sbbtagL", "tight", "No"),
            ("sbisobtagM", "sbiso", "L"),
        ]

        # Helper functions for 4-vector operations
        def _four_mass(v):
            arg = v.t * v.t - (v.x * v.x + v.y * v.y + v.z * v.z)
            return np.sqrt(ak.where(arg > 0, arg, 0.0))

        def _four_pt(v):
            arg = v.x * v.x + v.y * v.y
            return np.sqrt(ak.where(arg > 0, arg, 0.0))

        def _dR(a, b):
            dphi = ak.where(
                (a.phi - b.phi) > np.pi,
                a.phi - b.phi - 2 * np.pi,
                ak.where(
                    (a.phi - b.phi) < -np.pi, a.phi - b.phi + 2 * np.pi, a.phi - b.phi
                ),
            )
            arg = (a.eta - b.eta) ** 2 + dphi**2
            return np.sqrt(ak.where(arg > 0, arg, 0.0))

        # Loop over isolation modes
        for iso_mode in ["tight", "sbiso"]:
            sel_leps, sel_mask = select_lepton(
                events, self.channel, self._campaign, iso_mode=iso_mode
            )

            # Lepton veto
            if self.channel == "mu":
                req_lepveto = (ak.num(events.Muon[mu_loose], axis=1) == 1) & (
                    ak.num(events.Electron[el_loose], axis=1) == 0
                )
                other_mu = events.Muon[mu_loose & ~ak.fill_none(sel_mask, False)]
                other_el = events.Electron[el_loose]
            else:
                req_lepveto = (ak.num(events.Muon[mu_loose], axis=1) == 0) & (
                    ak.num(events.Electron[el_loose], axis=1) == 1
                )
                other_mu = events.Muon[mu_loose]
                other_el = events.Electron[el_loose & ~ak.fill_none(sel_mask, False)]

            req_lep = ak.num(sel_leps, axis=1) == 1
            _cf(f"{iso_mode}:lepveto", req_lumi & req_trig & req_metf & req_lepveto)
            _cf(
                f"{iso_mode}:tightlep",
                req_lumi & req_trig & req_metf & req_lepveto & req_lep,
            )

            # Jets
            jet_mask = _clean_jets(
                events,
                other_mu_mask=mu_loose,
                other_el_mask=el_loose,
            )
            jets_all = events.Jet[jet_mask]
            req_jets = ak.num(jets_all, axis=1) == 4
            _cf(
                f"{iso_mode}:jets==4",
                req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets,
            )

            # Base mask
            evmask_base = (
                req_lumi & req_trig & req_metf & req_lepveto & req_lep & req_jets
            )

            if not ak.any(evmask_base):
                continue

            # Slice to base once
            ev_base = events[evmask_base]
            jets_base = jets_all[evmask_base]
            lep_base = ak.firsts(sel_leps[evmask_base])

            # B-tag counts
            bmask_L = btag_wp(
                jets_base,
                self._year,
                self._campaign,
                tagger=self.tag_tagger,
                borc="b",
                wp="L",
            )
            bmask_M = btag_wp(
                jets_base,
                self._year,
                self._campaign,
                tagger=self.tag_tagger,
                borc="b",
                wp="M",
            )
            nb_L = ak.sum(ak.fill_none(bmask_L, False), axis=1)
            nb_M = ak.sum(ak.fill_none(bmask_M, False), axis=1)

            mask_central = nb_M >= 1
            mask_sb_btagM = (nb_M == 0) & (nb_L >= 1)
            mask_sb_btagL = nb_L == 0

            # MET 4-vector
            #FIXME why isn't this just ptetaphim ?
            met_b = ak.zip(
                {
                    "x": ev_base.MET.pt * np.cos(ev_base.MET.phi),
                    "y": ev_base.MET.pt * np.sin(ev_base.MET.phi),
                    "z": ak.zeros_like(ev_base.MET.pt),
                    "t": ev_base.MET.pt,
                },
                with_name="FourVector",
            )

            # ttbar reco ONCE per iso family
            best = ttbar_reco(jets_base, lep_base, met_b, maxjets=4)
            has_cand = best.get("has_cand", ak.Array([]))
            if not ak.any(has_cand):
                continue

            # Extract results once
            BH, BL, JA, JB = best["BH"], best["BL"], best["JA"], best["JB"]
#            nu, tlep, thad, chi2 = best["nu"], best["tlep"], best["thad"], best["chi2"]
            nu, tlep, thad = best["nu"], best["tlep"], best["thad"]
            neg_log_lambda, mTW_l = best["neg_log_lambda"],best["mTW_l"]

            # Compute derived quantities once
            had_tag = ak.fill_none(
                btag_wp(
                    BH,
                    self._year,
                    self._campaign,
                    tagger=self.tag_tagger,
                    borc="b",
                    wp="M",
                ),
                False,
            )
            lep_tag = ak.fill_none(
                btag_wp(
                    BL,
                    self._year,
                    self._campaign,
                    tagger=self.tag_tagger,
                    borc="b",
                    wp="M",
                ),
                False,
            )
            had_pt_ok = ak.fill_none(BH.pt >= 30.0, False)
            lep_pt_ok = ak.fill_none(BL.pt >= 30.0, False)

            # Compute masses and kinematics once
            W_full = ak.zip(
                {
                    "x": JA.x + JB.x,
                    "y": JA.y + JB.y,
                    "z": JA.z + JB.z,
                    "t": JA.t + JB.t,
                },
                with_name="FourVector",
            )
            tlep_mass_full = _four_mass(tlep)
            thad_mass_full = _four_mass(thad)
            whad_mass_full = _four_mass(W_full)
            tlep_pt_full = _four_pt(tlep)
            thad_pt_full = _four_pt(thad)
            dr_lep_blep_full = _dR(lep_base, BL)
            dr_lep_bhad_full = _dR(lep_base, BH)
            dr_ja_jb_full = _dR(JA, JB)

            # Compute kinbin once
            #met_vals = ak.to_numpy(ev_base.MET.pt)
            #chi2_np = ak.to_numpy(chi2)
            #met_edges = np.array([0.0, 40.0, 80.0, 200.0], dtype=float)
            #prob_edges = np.arange(11.0, 21.0, 1.0)
            #q = chi2_np / (2.0 * np.log(10.0))
            #met_bin = np.clip(
            #    np.digitize(met_vals, met_edges) - 1, 0, len(met_edges) - 2
            #)
            #prob_bin = np.clip(np.digitize(q, prob_edges) - 1, 0, len(prob_edges) - 2)
            #NQ = len(prob_edges) - 1
            

            mTW_l_bins = np.linspace(0,200,8)
            neg_log_lambda_bins = np.linspace(9, 30,3)
            mTW_l_bin = np.digitize(ak.to_numpy(mTW_l_bin), mTW_l_bins)
            neg_log_lambda_bin = np.digitize(ak.to_numpy(neg_log_lambda), neg_log_lambda_bins)

            #build kinbin by unfolding 2D binning into 1D for fitting
            kinbin_full = (mTW_l_bin * 3 + neg_log_lambda_bin).astype(np.int32)

            # Truth category once
            tt_cat_full = tt_truth_category(BL, BH, is_mc=not isRealData)

            # Now loop over regions and write immediately
            for rname, riso_mode, rtag_btagwp in region_specs:
                if iso_mode != riso_mode:
                    continue

                # Select region mask
                if rtag_btagwp == "M":
                    rmask = has_cand & mask_central
                elif rtag_btagwp == "L":
                    rmask = has_cand & mask_sb_btagM
                else:
                    rmask = has_cand & mask_sb_btagL

                rmask_np = ak.to_numpy(ak.fill_none(rmask, False))
                if not np.any(rmask_np):
                    continue

                # TnP fills
                require_tag = rname == "central"
                ones = ak.ones_like(had_pt_ok, dtype=bool)
                tnp_had_fill = ak.to_numpy(
                    had_pt_ok & (lep_tag if require_tag else ones)
                )[rmask_np]
                tnp_lep_fill = ak.to_numpy(
                    lep_pt_ok & (had_tag if require_tag else ones)
                )[rmask_np]
                tnp_had_pt = ak.to_numpy(BH.pt)[rmask_np]
                tnp_lep_pt = ak.to_numpy(BL.pt)[rmask_np]

                # Slice everything for this region
                ev_r = ev_base[rmask_np]
                jets_r = jets_base[rmask_np]
                lep_r = lep_base[rmask_np]
                BL_r = BL[rmask_np]
                BH_r = BH[rmask_np]
                JA_r = JA[rmask_np]
                JB_r = JB[rmask_np]
                nu_r = nu[rmask_np]
                tl_r = tlep[rmask_np]
                th_r = thad[rmask_np]

                # Assemble pruned view for this region only
                pr = ev_r
                pr = ak.with_field(pr, BL_r, "blep")
                pr = ak.with_field(pr, BH_r, "bhad")
                pr = ak.with_field(pr, JA_r, "ja")
                pr = ak.with_field(pr, JB_r, "jb")
                pr = ak.with_field(pr, nu_r, "nu")
                pr = ak.with_field(pr, tl_r, "tlep")
                pr = ak.with_field(pr, th_r, "thad")
                #pr = ak.with_field(pr, chi2_np[rmask_np], "chi2")
                pr = ak.with_field(pr, jets_r[:, :4], "SelJet")

                if self.channel == "mu":
                    pr = ak.with_field(pr, lep_r, "mu")
                else:
                    pr = ak.with_field(pr, lep_r, "ele")

                pr = ak.with_field(pr, ak.num(pr.SelJet, axis=1), "njet")
                pr = ak.with_field(pr, tlep_mass_full[rmask_np], "tlep_mass")
                pr = ak.with_field(pr, thad_mass_full[rmask_np], "thad_mass")
                pr = ak.with_field(pr, whad_mass_full[rmask_np], "whad_mass")
                pr = ak.with_field(pr, tlep_pt_full[rmask_np], "tlep_pt")
                pr = ak.with_field(pr, thad_pt_full[rmask_np], "thad_pt")
                pr = ak.with_field(pr, mTW_l[rmask_np], "mTW_l")
                pr = ak.with_field(pr, neg_log_lambda[rmask_np], "neg_log_lambda")
                pr = ak.with_field(pr, kinbin_full[rmask_np], "kinbin")
                pr = ak.with_field(pr, dr_lep_blep_full[rmask_np], "dr_lep_blep")
                pr = ak.with_field(pr, dr_lep_bhad_full[rmask_np], "dr_lep_bhad")
                pr = ak.with_field(pr, dr_ja_jb_full[rmask_np], "dr_ja_jb")
                pr = ak.with_field(pr, tnp_had_fill, "tnp_had_fill")
                pr = ak.with_field(pr, tnp_lep_fill, "tnp_lep_fill")
                pr = ak.with_field(pr, tnp_had_pt, "tnp_had_pt")
                pr = ak.with_field(pr, tnp_lep_pt, "tnp_lep_pt")
                pr = ak.with_field(
                    pr, np.full(len(pr), rname, dtype="U12"), "tnp_region"
                )
                pr = ak.with_field(pr, tt_cat_full[rmask_np], "tt_cat")
                pr = ak.with_field(pr, ak.Array(kinbin_full[rmask_np]), "kinbin")

                # Write immediately for this region
                if not self.noHist:
                    weights = weight_manager(pr, self.SF_map, self.isSyst)
                    systematics = (
                        [shift_name]
                        if shift_name is not None
                        else ["nominal"] + list(weights.variations)
                    )
                    output = self.write_histograms(
                        pr, output, weights, systematics, self.isSyst, self.SF_map
                    )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
