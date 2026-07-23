import collections
import uproot
import numpy as np
import awkward as ak
import hist as Hist
from coffea import processor
import re

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
    reweighting,
)
from BTVNanoCommissioning.helpers.func import update, dump_lumi
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_promptmvaid,
    ele_promptmvaid,
    MET_filters,
    mu_idiso,
    ele_mvatightid,
    btag_wp,
    wp_dict,
)

from BTVNanoCommissioning.helpers.definitions import get_discriminators, get_definitions


def select_lepton(events, channel, campaign, iso_mode="tight"):
    eta_cut = 2.5 if (("24" in campaign) or ("25" in campaign)) else 2.4
    use_promptmva = False
    if "24" in campaign or "25" in campaign:
        use_promptmva = True
    if channel == "mu":
        tight = (
            (events.Muon.pt > 30)
            & (abs(events.Muon.eta) < eta_cut)
            & (
                mu_promptmvaid(events, campaign)
                if use_promptmva
                else mu_idiso(events, campaign)
            )
        )
        if iso_mode == "tight":
            mask = tight
        elif iso_mode == "sbiso":
            iso = ak.fill_none(events.Muon.pfRelIso04_all, 999.0)
            loose = (
                (events.Muon.pt > 30)
                & (abs(events.Muon.eta) < eta_cut)
                & ak.fill_none(events.Muon.tightId, False)
            )
            mask = loose & (~tight) & (iso > 0.15) & (iso < 0.40)
        return events.Muon[mask], mask

    elif channel == "el":
        tight = (
            (events.Electron.pt > 30)
            & (abs(events.Electron.eta) < eta_cut)
            & (
                ele_promptmvaid(events, campaign)
                if use_promptmva
                else ele_mvatightid(events, campaign)
            )
        )
        if iso_mode == "tight":
            mask = tight
        elif iso_mode == "sbiso":
            loose = (events.Electron.pt > 30) & (abs(events.Electron.eta) < eta_cut)
            mva = ak.fill_none(
                getattr(events.Electron, "mvaIso", ak.zeros_like(events.Electron.pt)),
                -99.0,
            )
            mask = loose & (~tight) & (mva > 0.85) & (mva < 0.95)
        return events.Electron[mask], mask
    else:
        raise ValueError(channel)


def solve_nu_pz(px_l, py_l, pz_l, e_l, px_n, py_n, mW=80.4):
    """All inputs/outputs are flat 1D numpy arrays (caller flattens/unflattens)."""
    muW = (mW**2) / 2.0 + px_l * px_n + py_l * py_n
    a = e_l**2 - pz_l**2
    b = -2.0 * muW * pz_l
    c = e_l**2 * (px_n**2 + py_n**2) - muW**2
    disc = b**2 - 4.0 * a * c

    sqrt_disc = np.sqrt(np.maximum(disc, 0.0))
    inv2a = 0.5 / a
    pz1 = (-b + sqrt_disc) * inv2a
    pz2 = (-b - sqrt_disc) * inv2a

    # fast path: all combos have real solutions — skip gradient descent entirely
    ii = np.where(disc <= 0)[0]
    if len(ii) == 0:
        e1 = np.sqrt(px_n**2 + py_n**2 + pz1**2)
        e2 = np.sqrt(px_n**2 + py_n**2 + pz2**2)
        return (pz1, e1), (pz2, e2)

    # imaginary-discriminant subset only: gradient descent to rescale MET
    pzl_i = pz_l[ii]
    el_i = e_l[ii]
    el2_i = el_i**2
    pxn_i = px_n[ii]
    pyn_i = py_n[ii]
    plpnx_i = px_l[ii] * pxn_i
    plpny_i = py_l[ii] * pyn_i

    C_i = (pzl_i / el_i) ** 2 - 1.0
    mW2 = mW**2
    pzl_el2 = pzl_i / el2_i  # pz_l / e_l^2
    pzl_el2_sq = pzl_el2**2  # (pz_l / e_l^2)^2
    mW_el_sq = mW2 / el2_i  # (mW / e_l)^2

    # quadratic-form coefficients — constant across gradient-descent iterations
    As_i = 0.25 * (mW2**2 * pzl_el2_sq / C_i - mW_el_sq) / C_i
    bcoef = mW2 * pzl_el2_sq / C_i - mW_el_sq
    Bsx_i = bcoef * plpnx_i / C_i
    Bsy_i = bcoef * plpny_i / C_i
    ccoef = pzl_el2_sq / C_i - 1.0 / el2_i
    Csxx_i = (ccoef * plpnx_i**2 + pxn_i**2) / C_i
    Csxy_i = ccoef * 2.0 * plpnx_i * plpny_i / C_i
    Csyy_i = (ccoef * plpny_i**2 + pyn_i**2) / C_i

    x = np.ones(len(ii))
    y = np.ones(len(ii))
    step = np.full(len(ii), 0.1)
    for _ in range(3):
        U = (
            As_i
            + x * Bsx_i
            + y * Bsy_i
            + x**2 * Csxx_i
            + x * y * Csxy_i
            + y**2 * Csyy_i
        )
        dx = Bsx_i + 2.0 * Csxx_i * x + Csxy_i * y
        dy = Bsy_i + 2.0 * Csyy_i * y + Csxy_i * x
        norm = np.sqrt(dx**2 + dy**2 + 1e-10)
        x += step * dx / norm
        y += step * dy / norm
        U_new = (
            As_i
            + x * Bsx_i
            + y * Bsy_i
            + x**2 * Csxx_i
            + x * y * Csxy_i
            + y**2 * Csyy_i
        )
        step = np.where(U * U_new < 0, -0.5 * step, step)

    pxn_im = pxn_i * x
    pyn_im = pyn_i * y
    cross = plpnx_i * x + plpny_i * y
    pz_im = -0.5 / C_i * (mW_el_sq * pzl_i + 2.0 * pzl_el2 * cross)
    e_im = np.sqrt(pxn_im**2 + pyn_im**2 + pz_im**2)

    pz1_out = pz1.copy()
    pz2_out = pz2.copy()
    pz1_out[ii] = pz_im
    pz2_out[ii] = pz_im
    e1_out = np.sqrt(px_n**2 + py_n**2 + pz1**2)
    e2_out = np.sqrt(px_n**2 + py_n**2 + pz2**2)
    e1_out[ii] = e_im
    e2_out[ii] = e_im
    return (pz1_out, e1_out), (pz2_out, e2_out)


# open histograms of top/W mass distribution used for likelihood calculation
tf = uproot.open("src/BTVNanoCommissioning/helpers/sf_ttsemilep_likelihoods_pas.root")
# print axis ranges and peak locations
h2 = tf["had_tmwm"]
x_edges = h2.axis(0).edges()
y_edges = h2.axis(1).edges()
z_2d = h2.values().T
h1 = tf["had_mnu"]
nu_edges = h1.axis().edges()
z_1d = h1.values()


def centers(edges):
    return (edges[:-1] + edges[1:]) / 2.0


cw = centers(x_edges)
ct = centers(y_edges)
cnu = centers(nu_edges)


def interp1d(xaxis, zval, x):
    x = np.asarray(x)
    i = np.clip(np.searchsorted(xaxis, x, side="right") - 1, 0, len(xaxis) - 2)
    x1, x2 = xaxis[i], xaxis[i + 1]
    z1, z2 = zval[i], zval[i + 1]
    return z1 + (z2 - z1) * (x - x1) / (x2 - x1)


def interp2d(xaxis, yaxis, zgrid, x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    ix = np.clip(np.searchsorted(xaxis, x, side="right") - 1, 0, len(xaxis) - 2)
    iy = np.clip(np.searchsorted(yaxis, y, side="right") - 1, 0, len(yaxis) - 2)

    x1, x2 = xaxis[ix], xaxis[ix + 1]
    y1, y2 = yaxis[iy], yaxis[iy + 1]
    z11 = zgrid[ix, iy]
    z12 = zgrid[ix, iy + 1]
    z21 = zgrid[ix + 1, iy]
    z22 = zgrid[ix + 1, iy + 1]

    t = (x - x1) / (x2 - x1)
    u = (y - y1) / (y2 - y1)
    return (1 - t) * (1 - u) * z11 + t * (1 - u) * z21 + (1 - t) * u * z12 + t * u * z22


def calculate_mass_probability(m_type, masses):
    if m_type == "W, T":
        mW_val, mT_val = masses[0], masses[1]
        mW_val = np.clip(mW_val, cw[0], cw[-1])
        mT_val = np.clip(mT_val, ct[0], ct[-1])
        prob = interp2d(cw, ct, z_2d, mW_val, mT_val)
        return np.maximum(prob, 1e-24)
    elif m_type == "T":
        mTlep_val = np.clip(masses[0], cnu[0], cnu[-1])
        prob = interp1d(cnu, z_1d, mTlep_val)
        return np.maximum(prob, 1e-6)
    else:
        raise ValueError("Not a valid mass distribution")


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

    # broadcast lepton & MET via np.repeat — avoids awkward ragged overhead
    lep = lepton
    nc = ak.to_numpy(ak.num(BL.px, axis=1))  # combo count per event

    _lep_px = np.repeat(ak.to_numpy(lep.px), nc)
    _lep_py = np.repeat(ak.to_numpy(lep.py), nc)
    _lep_pz = np.repeat(ak.to_numpy(lep.pz), nc)
    _lep_e = np.repeat(ak.to_numpy(lep.energy), nc)
    _met_x = np.repeat(ak.to_numpy(met.x), nc)
    _met_y = np.repeat(ak.to_numpy(met.y), nc)

    (pz1_f, en1_f), (pz2_f, en2_f) = solve_nu_pz(
        _lep_px, _lep_py, _lep_pz, _lep_e, _met_x, _met_y, mW
    )

    pz1 = ak.unflatten(pz1_f, nc)
    en1 = ak.unflatten(en1_f, nc)
    pz2 = ak.unflatten(pz2_f, nc)
    en2 = ak.unflatten(en2_f, nc)

    lep_px = ak.unflatten(_lep_px, nc)
    lep_py = ak.unflatten(_lep_py, nc)
    lep_pz = ak.unflatten(_lep_pz, nc)
    lep_e = ak.unflatten(_lep_e, nc)
    lep_phi = ak.unflatten(np.repeat(ak.to_numpy(lep.phi), nc), nc)
    lep_pt = ak.unflatten(np.sqrt(_lep_px**2 + _lep_py**2), nc)
    met_x = ak.unflatten(_met_x, nc)
    met_y = ak.unflatten(_met_y, nc)
    met_phi = ak.unflatten(
        np.repeat(np.arctan2(ak.to_numpy(met.y), ak.to_numpy(met.x)), nc), nc
    )
    met_pt = ak.unflatten(np.sqrt(_met_x**2 + _met_y**2), nc)

    BL_px = BL.px
    BL_py = BL.py
    BL_pz = BL.pz
    BL_e = BL.energy

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

    def m(e, px, py, pz):
        return np.sqrt(np.maximum(0.0, e * e - (px * px + py * py + pz * pz)))

    m_tlep1, m_tlep2 = m(*t1), m(*t2)
    choose_1 = abs(m_tlep1 - mT) < abs(m_tlep2 - mT)
    pz_n = ak.where(choose_1, pz1, pz2)
    e_n = ak.where(choose_1, en1, en2)
    tlep_e = ak.where(choose_1, t1[0], t2[0])
    tlep_px = ak.where(choose_1, t1[1], t2[1])
    tlep_py = ak.where(choose_1, t1[2], t2[2])
    tlep_pz = ak.where(choose_1, t1[3], t2[3])
    m_tlep = ak.where(choose_1, m_tlep1, m_tlep2)
    valid_dnu = (m_tlep > 100) & (m_tlep < 240)
    # had side
    W_px, W_py, W_pz = JA.px + JB.px, JA.py + JB.py, JA.pz + JB.pz
    W_e = JA.energy + JB.energy
    mW_h = m(W_e, W_px, W_py, W_pz)

    # calculate and store transverse W mass for fit variable bins
    delta_phi = (lep_phi - met_phi + np.pi) % (2 * np.pi) - np.pi
    mTW_l = np.sqrt((2 * lep_pt * met_pt) * (1 - np.cos(delta_phi)))

    th_px, th_py, th_pz = BH.px + W_px, BH.py + W_py, BH.pz + W_pz
    th_e = BH.energy + W_e
    mT_h = m(th_e, th_px, th_py, th_pz)
    P_m2_m3 = calculate_mass_probability("W, T", [mW_h, mT_h])
    P_m_t1 = calculate_mass_probability("T", [m_tlep])
    valid_combo = P_m2_m3 > 1e-24

    neg_log_lambda = ak.where(valid_dnu, -np.log(P_m2_m3 * P_m_t1), np.inf)
    neg_log_lambda = ak.where(valid_combo, neg_log_lambda, np.inf)
    best_idx = ak.where(has_cand, ak.argmin(neg_log_lambda, axis=1), -1)
    li = ak.local_index(neg_log_lambda, axis=1)
    best_mask = (li == best_idx) & (best_idx >= 0)

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
        "neg_log_lambda": ak.firsts(neg_log_lambda[best_mask]),
        "mTW_l": ak.firsts(mTW_l[best_mask]),
    }
    return best


def sort_category(is_sig, cat_number: int):
    n = len(is_sig)

    if cat_number == 0:
        return ak.Array(np.full(n, "data", dtype="U5"))

    if cat_number == 1:
        gen_complete_np = ak.to_numpy(ak.fill_none(is_sig, False))
        return ak.Array(np.where(gen_complete_np, "sig", "ttbkg").astype("U5"))

    if cat_number == 2:
        return ak.Array(np.full(n, "st", dtype="U5"))
    if cat_number == 3:
        return ak.Array(np.full(n, "ew", dtype="U5"))
    if cat_number == 4:
        return ak.Array(np.full(n, "qcd", dtype="U5"))


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
        selectionModifier="",  # "tt_semilep_mu" or "tt_semilep_el"
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
        if selectionModifier not in ["tt_semilep_el", "tt_semilep_mu"]:
            raise ValueError(f"Invalid selectionModifier: {selectionModifier}")
        self.channel = "el" if (selectionModifier == "tt_semilep_el") else "mu"
        self.SF_map = load_SF(self._year, self._campaign, selectionModifier, isSyst)
        self._wp_table = wp_dict(year, campaign)
        self._all_taggers = sorted(self._wp_table.keys())
        self.tag_tagger = tag_tagger
        #        self._regions = ["central", "sbiso", "sbbtagM", "sbbtagL", "sbisobtagM"]
        self._regions = ["central", "sbbtagL", "sbbtagM"]

    @property
    def accumulator(self):
        return self._accumulator

    def define_histograms(self, events):
        """
        Define histograms to be written out by workflow
        """
        _hist_dict = {}

        # Common axes
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
        neg_log_lambda_axis = Hist.axis.Regular(
            48, 10, 22, name="neg_log_lambda", label=r"$-log(\\lambda)$"
        )

        tpt_axis = Hist.axis.Regular(60, 0, 600, name="pt", label=r"$p_T$ [GeV]")
        dr_axis = Hist.axis.Regular(20, 0, 8, name="dr", label="$\Delta$R")
        n_axis = Hist.axis.Integer(0, 10, name="n", label="N obj")

        # TnP-specific axes
        cat_axis = Hist.axis.StrCategory(["had", "lep"], name="cat")
        wp_axis = Hist.axis.StrCategory(["L", "M", "T", "XT", "XXT"], name="wp")
        kin_axis = Hist.axis.Regular(
            24, 0, 24, name="kinbin", label="Bin: $M_T(W_h)$ v. $-log(\\lambda)$"
        )
        ttcat_axis = Hist.axis.StrCategory(
            ["data", "sig", "ttbkg", "st", "ew", "qcd"], name="tt_cat"
        )
        result_axis = Hist.axis.StrCategory(["pass", "fail"], name="result")
        ptb_edges = [30.0, 50.0, 70.0, 100.0, 140.0, 200.0, 300.0]  # Update ptbin
        ptb_axis = Hist.axis.Variable(ptb_edges, name="ptb", label="$p_{T}(b)$ [GeV]")

        # objects for common kinematics
        obj_list = ["MET", "mu", "ele"]
        for i in range(4):
            obj_list.append(f"jet{i}")

        # Create histograms for each region
        for region in self._regions:
            # Basic kinematics
            _hist_dict[f"{region}_njet"] = Hist.Hist(
                syst_axis, n_axis, Hist.storage.Weight()
            )

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
        exclude_btv = [
            v
            for v in weights.variations
            if any(
                k in v.upper()
                for k in ("DEEP", "PNET", "ROBUST", "UPART", "BTV", "BTAG", "CTAG")
            )
        ]

        nj = 4
        pruned_ev["SelJet"] = pruned_ev.SelJet[:, :nj]

        if "hadronFlavour" in pruned_ev.SelJet.fields:
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
            genflavor = ak.zeros_like(pruned_ev.SelJet.pt, dtype=int)

        if "tnp_region" in pruned_ev.fields:
            region_labels = np.asarray(ak.to_numpy(pruned_ev.tnp_region), dtype="U20")
        else:
            region_labels = np.full(len(pruned_ev), "central", dtype="U20")

        # Pre-group output histograms by region prefix (one-time O(nhist) pass,
        # amortised across all systematics).
        region_hists = {r: {} for r in self._regions}
        for histname, h in output.items():
            prefix = histname.split("_", 1)[0]
            if prefix in self._regions:
                region_hists[prefix][histname] = h

        # Pre-compute available (tagger, wp, thr) triples once — fields are
        # constant across systs and regions.
        available = []
        if "bhad" in pruned_ev.fields and "blep" in pruned_ev.fields:
            for tagger, sub in self._wp_table.items():
                if not (
                    hasattr(pruned_ev.bhad, f"btag{tagger}B")
                    and hasattr(pruned_ev.blep, f"btag{tagger}B")
                ):
                    continue
                for wp_name, thr in sub.get("b", {}).items():
                    if wp_name == "No":
                        continue
                    available.append((tagger, wp_name, float(thr)))

        # Pre-compute per-region event slices and tnp arrays (constant across systs).
        _tnp_needed = {
            "tnp_had_fill",
            "tnp_lep_fill",
            "tnp_had_pt",
            "tnp_lep_pt",
            "bhad",
            "blep",
            "tt_cat",
            "kinbin",
        }
        region_cache = {}
        for region_prefix, hists in region_hists.items():
            rmask = region_labels == region_prefix
            if not np.any(rmask):
                continue
            ev = pruned_ev[rmask]
            gf_r = genflavor[rmask]
            tnp_data = None
            if _tnp_needed.issubset(set(ev.fields)):
                tnp_data = {
                    "ttcat": np.asarray(ak.to_numpy(ev.tt_cat), dtype="U5"),
                    "kinbin": np.asarray(ak.to_numpy(ev.kinbin), dtype=np.int32),
                    "had_fill": np.asarray(
                        ak.to_numpy(ak.fill_none(ev.tnp_had_fill, False)), dtype=bool
                    ),
                    "had_pt": np.asarray(ak.to_numpy(ev.tnp_had_pt), dtype=float),
                    "lep_fill": np.asarray(
                        ak.to_numpy(ak.fill_none(ev.tnp_lep_fill, False)), dtype=bool
                    ),
                    "lep_pt": np.asarray(ak.to_numpy(ev.tnp_lep_pt), dtype=float),
                    "bhad": ev.bhad,
                    "blep": ev.blep,
                }
            region_cache[region_prefix] = (rmask, ev, gf_r, tnp_data)

        # Loop over systematic variations.
        for syst in systematics:
            if isSyst is False and syst != "nominal":
                continue

            evt_w = (
                weights.weight()
                if syst == "nominal" or syst not in list(weights.variations)
                else weights.weight(modifier=syst)
            )
            exclude_list = [k for k in exclude_btv if k in weights.variations]
            evt_w_excl_btv = (
                weights.partial_weight(exclude=exclude_list) if exclude_list else evt_w
            )

            # One pass per region (not per histogram): compute weight slices once.
            for region_prefix, hists in region_hists.items():
                if region_prefix not in region_cache:
                    continue
                rmask, ev, gf_r, tnp_data = region_cache[region_prefix]
                w = evt_w[rmask]
                w_excl_btv = evt_w_excl_btv[rmask]

                for histname, h in hists.items():
                    # Selected electron histograms
                    if (
                        "SelElectron" in ev.fields
                        and "ele_" in histname
                        and histname.replace(f"{region_prefix}_ele_", "")
                        in ev.SelElectron.fields
                    ):
                        fld = histname.replace(f"{region_prefix}_ele_", "")
                        h.fill(syst, ak.to_numpy(ev.SelElectron[fld]), weight=w)
                        continue

                    # Selected muon histograms
                    if (
                        "SelMuon" in ev.fields
                        and "mu_" in histname
                        and histname.replace(f"{region_prefix}_mu_", "")
                        in ev.SelMuon.fields
                    ):
                        fld = histname.replace(f"{region_prefix}_mu_", "")
                        h.fill(syst, ak.to_numpy(ev.SelMuon[fld]), weight=w)
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
                            fld = histname.replace(f"{region_prefix}_jet{i}_", "")
                            if fld not in ev.SelJet.fields:
                                continue
                            h.fill(
                                syst,
                                ak.to_numpy(gf_r[:, i]),
                                ak.to_numpy(ev.SelJet[:, i][fld]),
                                weight=w,
                            )
                        continue

                    # b tag discriminants, filled with BTV-excluded weights
                    if any(
                        k in histname.replace(f"{region_prefix}_", " ")
                        for k in ("btag", "PNet", "ProbaN")
                    ):
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
                        h.fill(
                            syst=syst,
                            flav=ak.to_numpy(gf_r[:, i]),
                            discr=ak.to_numpy(ev.SelJet[:, i][disc_name]),
                            weight=w_excl_btv,
                        )
                        continue

                    # TnP yields — keys look like "<region>_<TAGGER>_tnp_yields"
                    if histname.endswith("_tnp_yields"):
                        if tnp_data is None:
                            continue

                        ttcat = tnp_data["ttcat"]
                        kinbin = tnp_data["kinbin"]

                        # had side
                        if tnp_data["had_fill"].any():
                            sel = tnp_data["had_fill"]
                            nsel = int(sel.sum())
                            bhad = tnp_data["bhad"]
                            for tagger, wp_name, thr in available:
                                scores = getattr(bhad, f"btag{tagger}B")
                                tagbit = np.asarray(
                                    ak.to_numpy(scores > thr), dtype=bool
                                )
                                h.fill(
                                    syst=syst,
                                    cat=np.full(nsel, "had", dtype="U3"),
                                    wp=np.full(nsel, wp_name, dtype="U3"),
                                    result=np.where(tagbit[sel], "pass", "fail").astype(
                                        "U4"
                                    ),
                                    tt_cat=ttcat[sel],
                                    kinbin=kinbin[sel],
                                    ptb=tnp_data["had_pt"][sel],
                                    weight=w[sel],
                                )

                        # lep side
                        if tnp_data["lep_fill"].any():
                            sel = tnp_data["lep_fill"]
                            nsel = int(sel.sum())
                            blep = tnp_data["blep"]
                            for tagger, wp_name, thr in available:
                                scores = getattr(blep, f"btag{tagger}B")
                                tagbit = np.asarray(
                                    ak.to_numpy(scores > thr), dtype=bool
                                )
                                h.fill(
                                    syst=syst,
                                    cat=np.full(nsel, "lep", dtype="U3"),
                                    wp=np.full(nsel, wp_name, dtype="U3"),
                                    result=np.where(tagbit[sel], "pass", "fail").astype(
                                        "U4"
                                    ),
                                    tt_cat=ttcat[sel],
                                    kinbin=kinbin[sel],
                                    ptb=tnp_data["lep_pt"][sel],
                                    weight=w[sel],
                                )
                        continue

                    # ttbar reco summaries
                    base_name = histname.replace(f"{region_prefix}_", "")
                    if base_name in (
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

                    # MET
                    if histname == f"{region_prefix}_MET_pt":
                        h.fill(syst, ak.to_numpy(ev.MET.pt), weight=w)
                        continue
                    if histname == f"{region_prefix}_MET_phi":
                        h.fill(syst, ak.to_numpy(ev.MET.phi), weight=w)
                        continue

        return output

    def process(self, events):
        events = missing_branch(events)
        sumws = reweighting(events, self.isSyst)
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), sumws, name)
            for collections, name in shifts
        )

    def process_shift(self, events, sumws, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        if isRealData:
            cat_number = 0
        elif re.search(r"QCD", dataset):
            cat_number = 4
        elif re.search(r"DY|Wto|WW|WZ|ZZ|WJet", dataset):
            cat_number = 3
        elif re.search(r"TTto|TT_|TTTo", dataset):
            cat_number = 1
        elif re.search(r"T[W-]|Tbar|TBbar|ST", dataset):
            cat_number = 2
        else:
            raise RuntimeError(f"Unknown MC category for dataset '{dataset}'. ")

        output = {} if self.noHist else self.define_histograms(events)
        # print(f"=== process_shift: {dataset}, shift={shift_name}, isData={isRealData}, n={len(events)} ===")

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
                    print("I AM HERE")
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

        # -------------------- Common preselection --------------------
        req_lumi = np.ones(len(events), dtype=bool)
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        # Triggers
        triggers_mu = ["IsoMu24"]
        triggers_el = []
        if "2016" in self._campaign:
            triggers_mu.append("IsoTkMu24")
            triggers_el = ["Ele27_WPTight_Gsf"]
        elif "17" in self._campaign:
            triggers_mu = ["IsoMu27"]
            triggers_el = ["Ele27_WPTight_Gsf", "Ele32_WPTight_Gsf"]
        elif "18" in self._campaign:
            triggers_mu = ["IsoMu24"]
            triggers_el = ["Ele32_WPTight_Gsf"]
        else:
            triggers_el = ["Ele30_WPTight_Gsf"]
        triggers = triggers_mu if self.channel == "mu" else triggers_el
        req_trig = HLT_helper(events, triggers)

        # MET filters
        req_metf = MET_filters(events, self._campaign)
        eta_cut = 2.5 if (("24" in self._campaign) or ("25" in self._campaign)) else 2.4
        use_promptmva = "24" in self._campaign or "25" in self._campaign
        # Loose veto objects
        mu_loose = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < eta_cut)
            & (
                mu_promptmvaid(events, self._campaign)
                if use_promptmva
                else mu_idiso(events, self._campaign)
            )
        )
        el_loose = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < eta_cut)
            & (
                ele_promptmvaid(events, self._campaign)
                if use_promptmva
                else ele_mvatightid(events, self._campaign)
            )
        )

        # Jet cleaning
        def _clean_jets(ev, veto_mu_mask, veto_el_mask):
            dr_mu = ev.Jet.metric_table(ev.Muon[veto_mu_mask])
            dr_el = ev.Jet.metric_table(ev.Electron[veto_el_mask])
            all_true = ak.ones_like(ev.Jet.pt, dtype=bool)
            has_mu = ak.num(ev.Muon[veto_mu_mask], axis=1) > 0
            has_el = ak.num(ev.Electron[veto_el_mask], axis=1) > 0
            clean_mu = ak.where(
                has_mu, ak.all(dr_mu > 0.4, axis=-1, mask_identity=True), all_true
            )
            clean_el = ak.where(
                has_el, ak.all(dr_el > 0.4, axis=-1, mask_identity=True), all_true
            )
            base_jet_mask = jet_id(ev, self._campaign, max_eta=eta_cut, min_pt=30)
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
        #        for iso_mode in ["tight", "sbiso"]:
        for iso_mode in ["tight"]:
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
                veto_mu_mask=mu_loose,
                veto_el_mask=el_loose,
            )
            jets_all = events.Jet[jet_mask]
            # require exactly 4 jets AND all jet pairs separated by DeltaR > 0.8
            pairs = ak.combinations(jets_all, 2, axis=1, fields=["0", "1"])
            dr_pairs = _dR(pairs["0"], pairs["1"])
            pairwise_ok = ak.all(dr_pairs > 0.8, axis=1)
            req_jets = (ak.num(jets_all, axis=1) == 4) & pairwise_ok
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
            mask_sb_btagM = nb_L >= 1
            mask_sb_btagL = nb_L == 0

            # MET 4-vector
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
            nu, tlep, thad = best["nu"], best["tlep"], best["thad"]
            neg_log_lambda, mTW_l = best["neg_log_lambda"], best["mTW_l"]
            log_lambda_mask = ak.to_numpy(ak.fill_none(neg_log_lambda < 50.0, False))

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
            had_tag_L = ak.fill_none(
                btag_wp(
                    BH,
                    self._year,
                    self._campaign,
                    tagger=self.tag_tagger,
                    borc="b",
                    wp="L",
                ),
                False,
            )
            lep_tag_L = ak.fill_none(
                btag_wp(
                    BL,
                    self._year,
                    self._campaign,
                    tagger=self.tag_tagger,
                    borc="b",
                    wp="L",
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

            # Update kinbin definition
            mTW_l_bins = np.array([0.0, 40.0, 80.0, 200.0], dtype=np.double)
            neg_log_lambda_bins = np.array(
                [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0], dtype=np.double
            )

            mTW_l_bin = (
                np.digitize(
                    ak.to_numpy(np.clip(ak.to_numpy(mTW_l), 0.1, 199.9)), mTW_l_bins
                )
                - 1
            )
            neg_log_lambda_bin = (
                np.digitize(
                    ak.to_numpy(np.clip(ak.to_numpy(neg_log_lambda), 11.01, 19.99)),
                    neg_log_lambda_bins,
                )
                - 1
            )

            kinbin_full = (mTW_l_bin * 8 + neg_log_lambda_bin).astype(np.int32)

            # GEN MATCHING:
            if cat_number == 1:
                genpart_all = ev_base.GenPart
                genpart = genpart_all[genpart_all.hasFlags("isHardProcess")]

                b_from_top = genpart[
                    (abs(genpart.pdgId) == 5)
                    & (abs(ak.fill_none(genpart.parent.pdgId, 0)) == 6)
                ]
                gen_charged_leptons = genpart[
                    ((abs(genpart.pdgId) == 11) | (abs(genpart.pdgId) == 13))
                    & (abs(ak.fill_none(genpart.parent.pdgId, 0)) == 24)
                ]
                gen_quarks_from_w = genpart[
                    (abs(genpart.pdgId) <= 4)
                    & (abs(genpart.pdgId) >= 1)
                    & (abs(ak.fill_none(genpart.parent.pdgId, 0)) == 24)
                ]

                has_2b = ak.num(b_from_top, axis=1) == 2
                n_clep = ak.num(gen_charged_leptons, axis=1)
                n_q = ak.num(gen_quarks_from_w, axis=1)
                gen_complete = has_2b & (n_clep == 1) & (n_q == 2)

                gen_charged_lep = ak.pad_none(gen_charged_leptons, 1, axis=1)
                lepton_pdgId = gen_charged_lep[:, 0].pdgId
                b_lep_pdgId = -np.sign(lepton_pdgId) * 5
                b_had_pdgId = np.sign(lepton_pdgId) * 5

                gen_b_lep = b_from_top[b_from_top.pdgId == b_lep_pdgId]
                gen_b_had = b_from_top[b_from_top.pdgId == b_had_pdgId]
                gen_b_lep = ak.pad_none(gen_b_lep, 1, axis=1)
                gen_b_had = ak.pad_none(gen_b_had, 1, axis=1)
                gen_qs = ak.pad_none(gen_quarks_from_w, 2, axis=1)

                dr = 0.4

                def get_dr(reco, gen):
                    return ak.fill_none(reco.delta_r(gen), 99.0)

                match_bhad = get_dr(BH, gen_b_had[:, 0]) < dr
                match_blep = get_dr(BL, gen_b_lep[:, 0]) < dr
                match_llep = get_dr(lep_base, gen_charged_lep[:, 0]) < dr

                dr_ja_q0 = get_dr(JA, gen_qs[:, 0])
                dr_ja_q1 = get_dr(JA, gen_qs[:, 1])
                dr_jb_q0 = get_dr(JB, gen_qs[:, 0])
                dr_jb_q1 = get_dr(JB, gen_qs[:, 1])
                w_match_1 = (dr_ja_q0 < dr) & (dr_jb_q1 < dr)
                w_match_2 = (dr_ja_q1 < dr) & (dr_jb_q0 < dr)
                match_whad = w_match_1 | w_match_2

                is_correct_kin = match_bhad & match_blep & match_llep & match_whad
                is_sig = gen_complete & is_correct_kin

                bmask_M_base = btag_wp(
                    jets_base,
                    self._year,
                    self._campaign,
                    tagger=self.tag_tagger,
                    borc="b",
                    wp="M",
                )
                has_2M = ak.sum(ak.fill_none(bmask_M_base, False), axis=1) >= 2
                gc2m = gen_complete & has_2M
            else:
                is_sig = ak.Array(np.zeros(len(ev_base), dtype=bool))

            tt_cat_full = sort_category(is_sig, cat_number)

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

                rmask_np = ak.to_numpy(ak.fill_none(rmask, False)) & log_lambda_mask
                if not np.any(rmask_np):
                    continue
                require_tag = rname == "central"
                ones = ak.ones_like(had_pt_ok, dtype=bool)
                if rname == "central":
                    tnp_had_fill = ak.to_numpy(had_pt_ok & (lep_tag))[rmask_np]
                    tnp_lep_fill = ak.to_numpy(lep_pt_ok & (had_tag))[rmask_np]
                elif rname == "sbbtagM":
                    tnp_had_fill = ak.to_numpy(had_pt_ok & (lep_tag_L & ~lep_tag))[
                        rmask_np
                    ]
                    tnp_lep_fill = ak.to_numpy(lep_pt_ok & (had_tag_L & ~had_tag))[
                        rmask_np
                    ]
                else:  # rname = sbbtagL
                    tnp_had_fill = ak.to_numpy(had_pt_ok & (~lep_tag & ~lep_tag_L))[
                        rmask_np
                    ]
                    tnp_lep_fill = ak.to_numpy(lep_pt_ok & (~had_tag & ~had_tag_L))[
                        rmask_np
                    ]

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
                pr = ak.with_field(pr, jets_r[:, :4], "SelJet")

                if self.channel == "mu":
                    pr = ak.with_field(pr, lep_r, "SelMuon")
                else:
                    pr = ak.with_field(pr, lep_r, "SelElectron")

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

                # Write immediately for this region
                if not self.noHist:
                    weights = weight_manager(
                        pr,
                        self.SF_map,
                        self.isSyst,
                        ttbar_reweights=getattr(self, "ttbar_reweights", "none"),
                        campaign=self._campaign,
                    )
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
