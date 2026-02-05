import awkward as ak
import numpy as np
from BTVNanoCommissioning.helpers.func import flatten
from .hist_helpers import get_hist_collections, get_axes_collections


def histogrammer(
    jet_fields=None,
    obj_list=[],
    hist_collections=["fourvector"],
    axes_collections=["common"],
    custom_axes=None,
    **kwargs,
):
    """
    Get a dictionary of histograms from the given histogram collections with the specified axes.

    Parameters:
    jet_fields (list): List of jet fields to be included in the histograms, e.g. events.Jet.fields
    obj_list (list): List of objects to be included in the histograms, e.g. "jet", "mu", "jet0". Accessed by the collections.
    hist_collections (list): List of histogram collections to be included in the histograms. Defined in the utils/histogramming/histograms directory.
    axes_collections (list): List of axes collections to be included in the histograms. Defined in the utils/histogramming/axes directory.
    custom_axes (dict): Dictionary of custom axes to be included in the histograms. Overwrites the axes defined in the axes_collections.
    **kwargs: Additional keyword arguments to be passed to the histogram collection functions.

    Example:
    ```python
    output = histogrammer(
        jet_fields=events.Jet.fields,
        obj_list=["jet", "mu"],
        hist_collections=["example", "common"],
        axes_collections=["common"],
    )
    ```

    Returns:
    dict: A dictionary containing the defined histograms.
    """

    ## Common axes
    if custom_axes:
        axes = custom_axes
    else:
        axes = get_axes_collections(axes_collections)

    ## Histograms
    _hist_dict = get_hist_collections(
        axes,
        hist_collections,
        obj_list=obj_list,
        jet_fields=jet_fields,
        **kwargs,
    )

    return _hist_dict


def histo_writer(pruned_ev, output, weights, systematics, isSyst, SF_map):
    return histo_writter(pruned_ev, output, weights, systematics, isSyst, SF_map)


# Filled common histogram
def histo_writter(pruned_ev, output, weights, systematics, isSyst, SF_map):
    """
    Write histograms to the output dictionary based on pruned events and other parameters.

    This function processes the pruned events and writes the histograms to the `output` dictionary. It takes into account the weights, systematics, and scale factors.

    Parameters:
    pruned_ev (coffea.nanoaodevents): The pruned events data to be histogrammed.
    output (dict): The output dictionary where histograms will be stored.
    weights (coffea.analysis_tools.Weights): The weights object for the events.
    systematics (list): A list of systematic variations to be considered.
    isSyst (str,bool): Indicating whether systematic variations are to be applied.
    SF_map (dict): A dictionary containing scale factors for different variables.

    Example:
    ```python
    histo_writter(pruned_ev, output, weights, systematics, isSyst, SF_map)
    ```

    Returns:
    None
    """
    exclude_btv = [
        "DeepCSVC",
        "DeepCSVB",
        "DeepJetB",
        "DeepJetC",
    ]  # exclude b-tag SFs for btag inputs
    # define Jet flavor

    # Reduce the jet to the correct dimension in the plot
    found4jets = False
    found2jets = False
    for key in output.keys():
        if "jet3" in key:  # Because 0-indexing
            found4jets = True
            break
        if "jet1" in key:  # Because 0-indexing
            found2jets = True
    nj = 4 if found4jets else 2 if found2jets else 1
    if nj != 1:
        pruned_ev["SelJet"] = pruned_ev.SelJet[:, :nj]
    if "var" in str(ak.type(pruned_ev.SelJet.pt)) and nj == 1:
        pruned_ev.SelJet = pruned_ev.SelJet[:, 0]

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
        if "MuonJet" in pruned_ev.fields:
            smflav = ak.values_astype(
                1
                * (
                    (pruned_ev.MuonJet.partonFlavour == 0)
                    & (pruned_ev.MuonJet.hadronFlavour == 0)
                )
                + pruned_ev.MuonJet.hadronFlavour,
                int,
            )
    else:
        isRealData = True
        genflavor = ak.zeros_like(pruned_ev.SelJet.pt, dtype=int)
        if "MuonJet" in pruned_ev.fields:
            smflav = ak.zeros_like(pruned_ev.MuonJet.pt, dtype=int)

    # Loop over the systematic variations
    for syst in systematics:
        if isSyst == False and syst != "nominal":
            break
        # Weight modifications for systematics
        weight = (
            weights.weight()
            if syst == "nominal" or syst not in list(weights.variations)
            else weights.weight(modifier=syst)
        )
        syst = np.full(len(weight), syst)
        # Loop over the histograms
        for histname, h in output.items():
            # Tagger score histograms
            if (
                "Deep" in histname
                and "btag" not in histname
                and histname in pruned_ev.SelJet.fields
            ):
                temp_weights = flatten(
                    ak.broadcast_arrays(
                        weights.partial_weight(exclude=exclude_btv),
                        pruned_ev.SelJet["pt"],
                    )[0]
                )
                temp_syst = np.full(len(temp_weights), syst[0])
                h.fill(
                    temp_syst,
                    flatten(genflavor),
                    flatten(pruned_ev.SelJet[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(
                            weights.partial_weight(exclude=exclude_btv),
                            pruned_ev.SelJet["pt"],
                        )[0]
                    ),
                )
            # PFcands histograms
            elif (
                "PFCands" in pruned_ev.fields
                and "PFCands" in histname
                and histname.split("_")[1] in pruned_ev.PFCands.fields
            ):
                if "MuonJet" in pruned_ev.fields:
                    h.fill(
                        syst,
                        flatten(
                            ak.broadcast_arrays(smflav, pruned_ev.PFCands["pt"])[0]
                        ),
                        flatten(pruned_ev.PFCands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                pruned_ev.PFCands["pt"],
                            )[0]
                        ),
                    )
                else:
                    h.fill(
                        syst,
                        flatten(
                            ak.broadcast_arrays(
                                genflavor[:, 0], pruned_ev.PFCands["pt"]
                            )[0]
                        ),
                        flatten(pruned_ev.PFCands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                pruned_ev.PFCands["pt"],
                            )[0]
                        ),
                    )
            # Leading lepton histograms
            elif (
                "hl_" in histname
                and "hl" in pruned_ev.fields
                and histname.replace("hl_", "") in pruned_ev.hl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.hl[histname.replace("hl_", "")]),
                    weight=weight,
                )
            # Subleading lepton histograms
            elif (
                "sl_" in histname
                and "sl" in pruned_ev.fields
                and histname.replace("sl_", "") in pruned_ev.sl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.sl[histname.replace("sl_", "")]),
                    weight=weight,
                )
            # Selected electron histograms
            elif (
                "ele_" in histname
                and histname.replace("ele_", "") in pruned_ev.SelElectron.fields
                and "Plus" not in pruned_ev.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.SelElectron[histname.replace("ele_", "")]),
                    weight=weight,
                )
            # Selected muon histograms
            elif (
                "mu_" in histname
                and histname.replace("mu_", "") in pruned_ev.SelMuon.fields
                and "Plus" not in pruned_ev.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.SelMuon[histname.replace("mu_", "")]),
                    weight=weight,
                )
            # Negatively charged lepton histograms-in DY workflow
            elif (
                "negl_" in histname
                and histname.replace("negl_", "") in pruned_ev.negl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.negl[histname.replace("negl_", "")]),
                    weight=weight,
                )
            # Posively charged lepton histograms-in DY workflow
            elif (
                "posl_" in histname
                and histname.replace("posl_", "") in pruned_ev.posl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.posl[histname.replace("posl_", "")]),
                    weight=weight,
                )
            # Soft muon histograms
            elif "soft_l" in histname and "ptratio" not in histname:
                h.fill(
                    syst,
                    smflav,
                    flatten(pruned_ev.SoftMuon[histname.replace("soft_l_", "")]),
                    weight=weight,
                )
            elif "njet" == histname:
                output["njet"].fill(syst, pruned_ev.njet, weight=weight)
            elif "npv" == histname:
                output["npv"].fill(syst, pruned_ev.PV.npvsGood, weight=weight)
            # Jet kinematics & deltaR between jet and lepton
            elif (
                "jet" in histname and "posl" not in histname and "negl" not in histname
            ):
                for i in range(nj):
                    if f"jet{i}" not in histname:
                        continue
                    if nj == 1:
                        sel_jet, flav = pruned_ev.SelJet, genflavor
                    else:
                        sel_jet, flav = pruned_ev.SelJet[:, i], genflavor[:, i]
                    if str(i) in histname:
                        if "dr_mujet" in histname:
                            h.fill(
                                syst,
                                flatten(flav),
                                flatten(sel_jet.delta_r(pruned_ev.SelMuon)),
                                weight=weight,
                            )
                        else:
                            h.fill(
                                syst,
                                flatten(flav),
                                flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                                weight=weight,
                            )
                # fill positively tagged jets, negatively tagged jets, and inclusive jets, binned in pt
                if histname.endswith("_postag_jet_pt") or histname.endswith(
                    "_negtag_jet_pt"
                ):
                    h.fill(
                        syst,
                        flatten(pruned_ev[histname.replace("_pt", "")].flavor),
                        flatten(pruned_ev[histname.replace("_pt", "")].pt),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weight, pruned_ev[histname.replace("_pt", "")].pt
                            )[0]
                        ),
                    )
                elif histname.endswith("jet_pt") and "AllSelJet" in pruned_ev.fields:
                    h.fill(
                        syst,
                        flatten(pruned_ev["AllSelJet"].flavor),
                        flatten(pruned_ev["AllSelJet"].pt),
                        weight=flatten(
                            ak.broadcast_arrays(weight, pruned_ev["AllSelJet"].pt)[0]
                        ),
                    )
            # Mu-jets distribution
            elif "lmujet_" in histname:
                h.fill(
                    syst,
                    smflav,
                    flatten(pruned_ev.MuonJet[histname.replace("lmujet_", "")]),
                    weight=weight,
                )
            # Filled discriminants
            elif "btag" in histname or "PNet" in histname:
                # Events with muon jet
                if "MuonJet" in pruned_ev.fields:
                    flavs, seljets = smflav, pruned_ev.MuonJet
                    nj = 1
                else:
                    flavs, seljets = genflavor, pruned_ev.SelJet
                for i in range(nj):
                    if not histname.endswith(str(i)):
                        continue
                    if nj > 1:
                        flav, seljet = flavs[:, i], seljets[:, i]
                    else:
                        flav, seljet = flavs, seljets
                    h.fill(
                        syst=syst,
                        flav=flav,
                        discr=seljet[histname.replace(f"_{i}", "")],
                        weight=weights.partial_weight(exclude=exclude_btv),
                    )

        if "dr_poslnegl" in output.keys():
            # DY histograms
            output["dr_poslnegl"].fill(
                syst, pruned_ev.posl.delta_r(pruned_ev.negl), weight=weight
            )
            output["dr_posljet"].fill(
                syst,
                genflavor,
                pruned_ev.posl.delta_r(pruned_ev.SelJet),
                weight=weight,
            )
            output["dr_negljet"].fill(
                syst,
                genflavor,
                pruned_ev.negl.delta_r(pruned_ev.SelJet),
                weight=weight,
            )

        # Muon enriched jet histograms
        if "MuonJet" in pruned_ev.fields:
            if (
                "hl" not in pruned_ev.fields
                and "SelElectron" in pruned_ev.fields
                and "SelMuon" in pruned_ev.fields
            ):
                pruned_ev["hl"] = ak.zip(
                    {
                        "pt": (
                            pruned_ev.SelMuon.pt
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.pt
                        ),
                        "eta": (
                            pruned_ev.SelMuon.eta
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.eta
                        ),
                        "phi": (
                            pruned_ev.SelMuon.phi
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.phi
                        ),
                        "energy": (
                            pruned_ev.SelMuon.energy
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.energy
                        ),
                    },
                    with_name="PtEtaPhiECandidate",
                )
            output["soft_l_ptratio"].fill(
                syst,
                flav=smflav,
                ratio=pruned_ev.SoftMuon.pt / pruned_ev.MuonJet.pt,
                weight=weight,
            )
            output["dr_lmujetsmu"].fill(
                syst,
                flav=smflav,
                dr=pruned_ev.MuonJet.delta_r(pruned_ev.SoftMuon),
                weight=weight,
            )
            if "hl" in pruned_ev.fields:
                output["hl_ptratio"].fill(
                    syst,
                    smflav,
                    ratio=pruned_ev.hl.pt / pruned_ev.MuonJet.pt,
                    weight=weight,
                )
            if "sl" in pruned_ev.fields:
                output["sl_ptratio"].fill(
                    syst,
                    smflav,
                    ratio=pruned_ev.sl.pt / pruned_ev.MuonJet.pt,
                    weight=weight,
                )
            if "SelMuon" in pruned_ev.fields and "hl" not in pruned_ev.fields:
                output["dr_lmujethmu"].fill(
                    syst,
                    flav=smflav,
                    dr=pruned_ev.MuonJet.delta_r(pruned_ev.SelMuon),
                    weight=weight,
                )
                output["dr_hmusmu"].fill(
                    syst, pruned_ev.SelMuon.delta_r(pruned_ev.SoftMuon), weight=weight
                )

        # dilepton system histograms: DY workflow
        if "dilep" in pruned_ev.fields:
            output["dilep_pt"].fill(syst, flatten(pruned_ev.dilep.pt), weight=weight)
            output["dilep_eta"].fill(syst, flatten(pruned_ev.dilep.eta), weight=weight)
            output["dilep_phi"].fill(syst, flatten(pruned_ev.dilep.phi), weight=weight)
            output["dilep_mass"].fill(
                syst, flatten(pruned_ev.dilep.mass), weight=weight
            )

        if "MET_pt" in output.keys():
            output["MET_pt"].fill(syst, flatten(pruned_ev.MET.pt), weight=weight)
            output["MET_phi"].fill(syst, flatten(pruned_ev.MET.phi), weight=weight)

        # ttbar dilepton kin workflow
        if "kindisc" in output.keys():
            output["kindisc"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.kindisc),
                weight=flatten(
                    ak.broadcast_arrays(
                        weights.partial_weight(exclude=exclude_btv), pruned_ev.kindisc
                    )[0]
                ),
            )
            output["close_mlj"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_mlj),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_mlj)[0]),
            )
            output["close_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_deta),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_deta)[0]),
            )
            output["close_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_dphi),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_dphi)[0]),
            )
            output["close_ptrel"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_ptrel),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_ptrel)[0]),
            )
            output["close_lj2ll_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_lj2ll_deta),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.close_lj2ll_deta)[0]
                ),
            )
            output["close_lj2ll_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_lj2ll_dphi),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.close_lj2ll_dphi)[0]
                ),
            )
            output["far_mlj"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_mlj),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_mlj)[0]),
            )
            output["far_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_deta),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_deta)[0]),
            )
            output["far_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_dphi),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_dphi)[0]),
            )
            output["far_ptrel"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_ptrel),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_ptrel)[0]),
            )
            output["far_lj2ll_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_lj2ll_deta),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.far_lj2ll_deta)[0]
                ),
            )
            output["far_lj2ll_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_lj2ll_dphi),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.far_lj2ll_dphi)[0]
                ),
            )
            output["j2ll_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.j2ll_deta),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.j2ll_deta)[0]),
            )
            output["j2ll_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.j2ll_dphi),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.j2ll_dphi)[0]),
            )

    return output
