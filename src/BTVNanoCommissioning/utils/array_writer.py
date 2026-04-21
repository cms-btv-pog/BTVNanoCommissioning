from BTVNanoCommissioning.helpers.func import uproot_writeable
import numpy as np
import awkward as ak
import os, uproot

arraySchema = {
    "CFM": [
        "SelJet_btag",
        "SelJet_pt",
        "SelJet_muEF",
        "SelJet_chEmEF",
        "SelJet_chHEF",
        "SelJet_chMultiplicity",
        "SelJet_electronIdx1",
        "SelJet_hfEmEF",
        "SelJet_hfHEF",
        "SelJet_muonIdx1",
        "SelJet_nConstituents",
        "SelJet_nMuons",
        "SelJet_nElectrons",
        "SelJet_nSVs",
        "SelJet_neEmEF",
        "SelJet_neHEF",
        "SelJet_neMultiplicity",
        "SelJet_puIdDisc",
        "SelJet_rawFactor",
        "SelJet_eta",
        "SelJet_phi",
        "SelJet_mass",
        "SelJet_hadronFlavour",
        "SelJet_partonFlavour",
        "SelJet_isMuonJet",
        "njet",
        "PuppiMET_pt",
        "PuppiMET_phi",
        "dijet_pt",
        "dijet_eta",
        "dijet_phi",
        "dijet_mass",
        "top_pt",
        "top_eta",
        "top_phi",
        "top_mass",
        "PV_npvs",
        "PV_npvsGood",
        "dilep_mass",
        "dilep_pt",
        "dilep_eta",
        "dilep_phi",
        "SoftMuon_dxySig",
        "MuonJet_muneuEF",
        "soft_l_ptratio",
        "osss",
        "W_transmass",
        "W_pt",
        "W_eta",
        "W_phi",
        "W_mass",
        "Pileup_nTrueInt",
        "Pileup_nPU",
    ]
}


def array_writer(
    processor_class,  # the NanoProcessor class ("self")
    pruned_event,  # the event with specific calculated variables stored
    nano_event,  # entire NanoAOD/PFNano event with many variables
    weights,  # weight for the event
    systname,  # name of systematic shift
    dataset,  # dataset name
    isRealData,  # boolean
    out_dir_base="",  # string
    remove=[
        "SoftMuon",
        "MuonJet",
        "dilep",
        "OtherJets",
        "Jet",
    ],  # remove from variable list
    kinOnly=[
        "Muon",
        "Jet",
        "SoftMuon",
        "dilep",
        "charge",
        "MET",
    ],  # variables for which only kinematic properties are kept
    kins=[
        "pt",
        "eta",
        "phi",
        "mass",
        "pfRelIso04_all",
        "pfRelIso03_all",
        "dxy",
        "dz",
    ],  # kinematic propoerties for the above variables
    othersData=[
        "PFCands_*",
        "MuonJet_*",
        "SV_*",
        "PV_npvs",
        "PV_npvsGood",
        "Rho_*",
        "SoftMuon_dxySig",
        "Muon_sip3d",
    ],  # other fields, for Data and MC
    doOnly=None,
    schema=None,
    othersMC=["Pileup_nTrueInt", "Pileup_nPU"],  # other fields, for MC only
    empty=False,
):
    if weights is not None:
        pruned_event["weight"] = weights.weight()
        for ind_wei in weights.weightStatistics.keys():
            pruned_event[f"{ind_wei}_weight"] = weights.partial_weight(
                include=[ind_wei]
            )
        if len(systname) > 1:
            for syst in systname:
                if syst == "nominal":
                    continue
                pruned_event[f"weight_syst_{syst}"] = weights.weight(modifier=syst)

    if empty:
        print("WARNING: No events selected. Writing blank file.")
        out_branch = []
    elif doOnly is not None:
        if "weight" not in doOnly:
            doOnly.extend([b for b in pruned_event.fields if "weight" in b])
        out_branch = np.array(doOnly)
        if not isRealData:
            out_branch = np.append(out_branch, othersMC)
    elif schema is not None:
        netout = arraySchema[schema] + [b for b in pruned_event.fields if "weight" in b]
        out_branch = np.array(netout)
    else:
        # Get only the variables that were added newly
        out_branch = np.setdiff1d(
            np.array(pruned_event.fields), np.array(nano_event.fields)
        )

        # Handle kinOnly vars
        remove = remove + ["PFCands", "hl", "sl", "posl", "negl"]
        for v in remove:
            out_branch = np.delete(out_branch, np.where((out_branch == v)))

        for kin in kins:
            for obj in kinOnly:
                if "MET" in obj and ("pt" != kin or "phi" != kin):
                    continue
                if (obj != "SelMuon" and obj != "SoftMuon") and (
                    "pfRelIso04_all" == kin or "d" in kin
                ):
                    continue
                out_branch = np.append(out_branch, [f"{obj}_{kin}"])

        # Handle data vars
        out_branch = np.append(out_branch, othersData)

        if not isRealData:
            out_branch = np.append(out_branch, othersMC)

    # Write to root files
    print("Branches to write:", out_branch)
    outdir = f"{out_dir_base}{processor_class.name}/{systname[0]}/{dataset}/"
    os.system(f"mkdir -p {outdir}")

    with uproot.recreate(
        f"{outdir}/{nano_event.metadata['filename'].split('/')[-1].replace('.root','')}_{int(nano_event.metadata['entrystop']/processor_class.chunksize)}.root"
    ) as fout:
        if not empty:
            fout["Events"] = uproot_writeable(pruned_event, include=out_branch)
        fout["TotalEventCount"] = ak.Array(
            [nano_event.metadata["entrystop"] - nano_event.metadata["entrystart"]]
        )
        if not isRealData:
            fout["TotalEventWeight"] = ak.Array([ak.sum(nano_event.genWeight)])
