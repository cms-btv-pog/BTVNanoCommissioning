import awkward as ak
import numpy as np


def HLT_helper(events, triggers):
    checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
    if ak.all(checkHLT == False):
        raise ValueError(
            "HLT paths:", triggers, " are all invalid in", events.metadata["dataset"]
        )
    elif ak.any(checkHLT == False):
        print(
            np.array(triggers)[~checkHLT], " not exist in", events.metadata["dataset"]
        )
    trig_arrs = [events.HLT[_trig] for _trig in triggers if hasattr(events.HLT, _trig)]
    req_trig = np.zeros(len(events), dtype="bool")
    for t in trig_arrs:
        req_trig = req_trig | t
    return req_trig


def jet_id(events, campaign, max_eta=2.5, min_pt=20):
    # Run 3 NanoAODs have a bug in jetId,
    # Implement fix from:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#nanoAOD_Flags
    # Note: this is only the jetId==6, ie. passJetIdTightLepVeto. Looser selection is not implemented.
    if campaign in ["Summer22", "Summer22EE", "Summer23", "Summer23BPix"]:
        # NanoV12
        jetid = ak.where(
            abs(events.Jet.eta) <= 2.7,
            (events.Jet.jetId >= 2)
            & (events.Jet.muEF < 0.8)
            & (events.Jet.chEmEF < 0.8),
            ak.where(
                (abs(events.Jet.eta) > 2.7) & (abs(events.Jet.eta) <= 3.0),
                (events.Jet.jetId >= 2) & (events.Jet.neHEF < 0.99),
                ak.where(
                    (abs(events.Jet.eta) > 3.0),
                    (events.Jet.jetId & (1 << 1)) & (events.Jet.neEmEF < 0.4),
                    ak.zeros_like(events.Jet.pt, dtype=bool),
                ),
            ),
        )
    elif campaign in ["Winter24", "Summer24"]:
        # NanoV13 & NanoV14
        barrel = (
            (events.Jet.neHEF < 0.99)
            & (events.Jet.neEmEF < 0.9)
            & (events.Jet.chMultiplicity + events.Jet.neMultiplicity > 1)
            & (events.Jet.chHEF > 0.01)
            & (events.Jet.chMultiplicity > 0)
        )
        t1 = (events.Jet.neHEF < 0.9) & (events.Jet.neEmEF < 0.99)
        t2 = events.Jet.neHEF < 0.99
        endcap = (events.Jet.neMultiplicity >= 2) & (events.Jet.neEmEF < 0.4)

        jetid = ak.where(
            abs(events.Jet.eta) <= 2.6,
            barrel,
            ak.where(
                (abs(events.Jet.eta) > 2.6) & (abs(events.Jet.eta) <= 2.7),
                t1,
                ak.where(
                    (abs(events.Jet.eta) > 2.7) & (abs(events.Jet.eta) <= 3.0),
                    t2,
                    ak.where(
                        (abs(events.Jet.eta) > 3.0),
                        endcap,
                        ak.zeros_like(events.Jet.pt, dtype=bool),
                    ),
                ),
            ),
        )
    else:
        jetid = events.Jet.jetId >= 5

    jetid = ak.values_astype(jetid, np.bool)

    if campaign == "Rereco17_94X":
        # Use puId for Run2
        jetmask = (
            (events.Jet.pt > min_pt)
            & (abs(events.Jet.eta) <= max_eta)
            & (jetid)
            & ((events.Jet.pt > 50) | (events.Jet.puId >= 7))
        )
    else:
        jetmask = (events.Jet.pt > min_pt) & (abs(events.Jet.eta) <= max_eta) & (jetid)
    return jetmask


## FIXME: Electron cutbased Id & MVA ID not exist in Winter22Run3 sample
def ele_cuttightid(events, campaign):
    elemask = (
        (abs(events.Electron.eta) < 1.4442)
        | ((abs(events.Electron.eta) < 2.5) & (abs(events.Electron.eta) > 1.566))
    ) & (events.Electron.cutBased > 3)
    return elemask


def ele_mvatightid(events, campaign):
    elemask = (
        (abs(events.Electron.eta) < 1.4442)
        | ((abs(events.Electron.eta) < 2.5) & (abs(events.Electron.eta) > 1.566))
    ) & (events.Electron.mvaIso_WP80 > 0.5)
    return elemask


def softmu_mask(events, campaign, dxySigCut=0):
    softmumask = (
        (events.Muon.pt < 25)
        & (abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId > 0.5)
        & (events.Muon.pfRelIso04_all > 0.2)
        & (abs(events.Muon.dxy / events.Muon.dxyErr) > dxySigCut)
        & (events.Muon.jetIdx != -1)
    )
    return softmumask


def mu_idiso(events, campaign):
    mumask = (
        (abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId > 0.5)
        & (events.Muon.pfRelIso04_all <= 0.15)
    )
    return mumask


def btag_mu_idiso(events, campaign):
    mumask = (
        (abs(events.Muon.eta) < 2.4)
        & (events.Muon.tightId > 0.5)
        & (events.Muon.pfRelIso04_all < 0.12)
    )
    return mumask


def jet_cut(events, campaign):
    multijetmask = (
        (abs(events.Jet.eta) < 2.4) & (events.Jet.pt > 180) & (events.Jet.jetId >= 5)
    )
    return multijetmask


def MET_filters(events, campaign):
    # apply MET filter
    metfilter = ak.ones_like(events.run, dtype=bool)
    for flag in met_filters[campaign]["data" if "Run" else "mc"]:
        metfilter = events.Flag[flag] & metfilter
    ## Flag_ecalBadCalibFilter
    badjet = (
        (events.Jet.pt > 50)
        & (events.Jet.eta >= -0.5)
        & (events.Jet.eta <= -0.1)
        & (events.Jet.phi >= -2.1)
        & (events.Jet.phi <= -1.8)
        & ((events.Jet.neEmEF > 0.9) | (events.Jet.chEmEF > 0.9))
        & (events.Jet.delta_phi(events.PuppiMET) > 2.9)
    )
    ecalBadCalibFilter = (
        (ak.sum(badjet, axis=-1) >= 1)
        & (events.PuppiMET.pt > 100)
        & (events.run >= 362433)
        & (events.run <= 367144)
    )
    metfilter = metfilter & ~ecalBadCalibFilter
    return metfilter


def btag_wp(jets, year, campaign, tagger, borc, wp):
    WP = wp_dict(year, campaign)
    if borc == "b":
        jet_mask = jets[f"btag{tagger}B"] > WP[tagger]["b"][wp]
    elif borc == "c":
        jet_mask = (jets[f"btag{tagger}CvB"] > WP[tagger]["c"][wp][1]) & (jets[f"btag{tagger}CvL"] > WP[tagger]["c"][wp][0])
    return jet_mask


def calculate_new_discriminators(ith_jets):
    probudg = ith_jets.btagUParTAK4UDG
    SvUDG = ith_jets.btagUParTAK4SvUDG
    probs = ak.Array(np.where(
        (SvUDG > 0.0) & (probudg > 0.0),
        SvUDG * probudg / (1.0 - SvUDG), -1.0
    ))
    CvL = ith_jets.btagUParTAK4CvL
    probc = ak.Array(np.where(
        (CvL > 0.0) & (probs > 0.0) & (probudg > 0.0),
        CvL * (probs + probudg) / (1.0 - CvL), -1.0
    ))
    CvB = ith_jets.btagUParTAK4CvB
    probbbblepb = ak.Array(np.where(
        (CvB > 0.0) & (probc > 0.0),
        (1.0 - CvB) * probc / CvB, -1.0
    ))
    BvC = ak.Array(np.where(
        CvB > 0.0,
        1.0 - CvB, -1.0
    ))
    HFvLF = ak.Array(np.where(
        (probbbblepb > 0.0) & (probc > 0.0) & (probs > 0.0) & (probudg > 0.0),
        (probbbblepb + probc) / (probbbblepb + probc + probs + probudg), -1.0
    ))
    return HFvLF, BvC


def get_wp_2D(ith_jet_HFvLF_val, ith_jet_BvC_val, year, campaign, tagger):
    WP = wp_dict(year, campaign)
    for wp in WP[tagger]["2D"]["mapping"].keys():
        HFvLF_low = WP[tagger]["2D"][wp][0]
        HFvLF_high = WP[tagger]["2D"][wp][1]
        BvC_low = WP[tagger]["2D"][wp][2]
        BvC_high = WP[tagger]["2D"][wp][3]
        if ith_jet_HFvLF_val is None or ith_jet_HFvLF_val == -1:
            return None
        if ith_jet_BvC_val is None or ith_jet_BvC_val == -1:
            return None
        if HFvLF_low == 0:
            passWP = (ith_jet_HFvLF_val >= HFvLF_low)
        else:
            passWP = (ith_jet_HFvLF_val > HFvLF_low)
        passWP = passWP & (ith_jet_HFvLF_val <= HFvLF_high)
        if BvC_low == 0:
            passWP = passWP & (ith_jet_BvC_val >= BvC_low)
        else:
            passWP = passWP & (ith_jet_BvC_val > BvC_low)
        passWP = passWP & (ith_jet_BvC_val <= BvC_high)
        if passWP:
            return WP[tagger]["2D"]["mapping"][wp]
    raise ValueError(
        "Somehow did not find the working point for values HFvLF=",
        ith_jet_HFvLF_val, "and BvC=", ith_jet_BvC_val
    )


def btag_wp_2D(ith_jet_HFvLF, ith_jet_BvC, year, campaign, tagger, wp_low, wp_high):
    WP = wp_dict(year, campaign)
    jet_wps = ak.Array([get_wp_2D(ith_jet_HFvLF[i], ith_jet_BvC[i], year, campaign, tagger) for i in range(len(ith_jet_HFvLF))])
    return (jet_wps >= WP[tagger]["2D"]["mapping"][wp_low]) & (jet_wps >= WP[tagger]["2D"]["mapping"][wp_high])


btag_wp_dict = {
    "2022_Summer22": {
        "DeepFlav": {
            "b": {
                "No": 0.0,
                "L": 0.0583,
                "M": 0.3086,
                "T": 0.7183,
                "XT": 0.8111,
                "XXT": 0.9512,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.042, 0.208], # CvL, then CvB
                "M": [0.108, 0.299],
                "T": [0.303, 0.243],
            },
        },
        "RobustParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0849,
                "M": 0.4319,
                "T": 0.8482,
                "XT": 0.9151,
                "XXT": 0.9874,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.039, 0.068],
                "M": [0.117, 0.130],
                "T": [0.360, 0.095],
            },
        },
        "PNet": {
            "b": {
                "No": 0.0,
                "L": 0.047,
                "M": 0.245,
                "T": 0.6734,
                "XT": 0.7862,
                "XXT": 0.961,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.054, 0.181],
                "M": [0.160, 0.306],
                "T": [0.492, 0.259],
            },
        },
    },
    "2022_Summer22EE": {
        "DeepFlav": {
            "b": {
                "No": 0.0,
                "L": 0.0614,
                "M": 0.3196,
                "T": 0.73,
                "XT": 0.8184,
                "XXT": 0.9542,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.042, 0.206], # CvL, then CvB
                "M": [0.108, 0.298],
                "T": [0.305, 0.241],
            },
        },
        "RobustParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0897,
                "M": 0.451,
                "T": 0.8604,
                "XT": 0.9234,
                "XXT": 0.9893,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.039, 0.067],
                "M": [0.117, 0.128],
                "T": [0.358, 0.095],
            },
        },
        "PNet": {
            "b": {
                "No": 0.0,
                "L": 0.0499,
                "M": 0.2605,
                "T": 0.6915,
                "XT": 0.8033,
                "XXT": 0.9664,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.054, 0.182],
                "M": [0.160, 0.304],
                "T": [0.491, 0.258],
            },
        },
    },
    "2023_Summer23": {
        "DeepFlav": {
            "b": {
                "No": 0.0,
                "L": 0.0479,
                "M": 0.2431,
                "T": 0.6553,
                "XT": 0.7667,
                "XXT": 0.9459,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.042, 0.234], # CvL, then CvB
                "M": [0.102, 0.322],
                "T": [0.250, 0.262],
                "XT": [0.371, 0.440],
            },
        },
        "RobustParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0681,
                "M": 0.3487,
                "T": 0.7969,
                "XT": 0.8882,
                "XXT": 0.9883,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.038, 0.086],
                "M": [0.109, 0.153],
                "T": [0.308, 0.113],
                "XT": [0.469, 0.275],
            },
        },
        "PNet": {
            "b": {
                "No": 0.0,
                "L": 0.0358,
                "M": 0.1917,
                "T": 0.6172,
                "XT": 0.7515,
                "XXT": 0.9659,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.052, 0.220],
                "M": [0.148, 0.353],
                "T": [0.434, 0.300],
                "XT": [0.634, 0.549],
            },
        },
    },
    "2023_Summer23BPix": {
        "DeepFlav": {
            "b": {
                "No": 0.0,
                "L": 0.048,
                "M": 0.2435,
                "T": 0.6563,
                "XT": 0.7671,
                "XXT": 0.9483,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.042, 0.242], # CvL, then CvB
                "M": [0.102, 0.328],
                "T": [0.250, 0.267],
                "XT": [0.371, 0.444],
            },
        },
        "RobustParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0683,
                "M": 0.3494,
                "T": 0.7994,
                "XT": 0.8877,
                "XXT": 0.9883,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.038, 0.091],
                "M": [0.109, 0.157],
                "T": [0.308, 0.116],
                "XT": [0.469, 0.281],
            },
        },
        "PNet": {
            "b": {
                "No": 0.0,
                "L": 0.0359,
                "M": 0.1919,
                "T": 0.6133,
                "XT": 0.7544,
                "XXT": 0.9688,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.052, 0.228],
                "M": [0.149, 0.358],
                "T": [0.436, 0.303],
                "XT": [0.634, 0.5552],
            },
        },
    },
    "2024_Summer24" : {
        "UParTAK4": {
            "b": {
                "No": 0.0,
                "L": 0.0246,
                "M": 0.1272,
                "T": 0.4648,
                "XT": 0.6298,
                "XXT": 0.9739,
            },
            "c": {
                "No": [0.0, 0.0],
                "L": [0.086, 0.233], # CvL, then CvB
                "M": [0.291, 0.457],
                "T": [0.650, 0.421],
                "XT": [0.810, 0.736],
            },
            "2D": {
                "No": [0.0,   1.0,   0.0,   1.0  ],
                "L0": [0.0,   0.264, 0.0,   1.0  ], # [HFvLF low, HFvLF high, BvC low, BvC high]
                "C0": [0.264, 0.448, 0.0,   1.0  ],
                "C1": [0.448, 0.767, 0.0,   1.0  ],
                "C2": [0.767, 1.0,   0.028, 0.094],
                "C3": [0.767, 1.0,   0.01,  0.028],
                "C4": [0.767, 1.0,   0.0,   0.01 ],
                "B0": [0.767, 1.0,   0.094, 0.69 ],
                "B1": [0.767, 1.0,   0.69,  0.918],
                "B2": [0.767, 1.0,   0.918, 0.978],
                "B3": [0.767, 1.0,   0.978, 0.994],
                "B4": [0.767, 1.0,   0.994, 1.0  ],
                "mapping": {
                    "L0": 0,
                    "C0": 1,
                    "C1": 2,
                    "C2": 3,
                    "C3": 4,
                    "C4": 5,
                    "B0": 6,
                    "B1": 7,
                    "B2": 8,
                    "B3": 9,
                    "B4": 10,
                }
            },
        },
    },
}


import os, correctionlib

def wp_dict(year, campaign):
    """
    year :
    """
    # btag_wp_dicts={}
    cache_key = f"{year}_{campaign}"

    if cache_key in btag_wp_dict:
        return btag_wp_dict[cache_key]

    name_map = {
        "deepJet": "DeepFlav",
        "robustParticleTransformer": "RobustParTAK4",
        "particleNet": "PNet",
        "unifiedParticleTransformer": "UParTAK4",
    }

    wps_dict = {}
    if os.path.exists(
        f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{year}_{campaign}"
    ):
        btag = correctionlib.CorrectionSet.from_file(
            f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{year}_{campaign}/btagging.json.gz"
        )
        ctag = correctionlib.CorrectionSet.from_file(
            f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{year}_{campaign}/ctagging.json.gz"
        )
        tagger_list = [i for i in list(btag.keys()) if "wp_values" in i]

        if len(tagger_list) == 0:
            btag_wp_dict[cache_key] = wps_dict
            return wps_dict

        for tagger in tagger_list:
            wps_dict[name_map[tagger.replace("_wp_values", "")]] = {"b": {}, "c": {}}
            # Get b WPs
            bwp = btag[tagger].inputs[0].description.split("/")
            wps_dict[name_map[tagger.replace("_wp_values", "")]]["b"] = {
                wp: btag[tagger].evaluate(wp) for wp in bwp
            }
            # Get c WPs in [CvL, CvB]
            cwp = ctag[tagger].inputs[0].description.split("/")
            wps_dict[name_map[tagger.replace("_wp_values", "")]]["c"] = {
                wp: [ctag[tagger].evaluate(wp, "CvL"), ctag[tagger].evaluate(wp, "CvB")]
                for wp in cwp
            }
        btag_wp_dict[cache_key] = wps_dict
        return wps_dict

    else:
        btag_wp_dict[cache_key] = wps_dict
        return wps_dict


met_filters = {
    "2016preVFP_UL": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "eeBadScFilter",
        ],
    },
    "2016postVFP_UL": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "eeBadScFilter",
        ],
    },
    "2017_UL": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
            "ecalBadCalibFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
            "ecalBadCalibFilter",
        ],
    },
    "2018_UL": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
            "ecalBadCalibFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "HBHENoiseFilter",
            "HBHENoiseIsoFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
            "ecalBadCalibFilter",
        ],
    },
    "Summer22": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
    },
    "Summer22EE": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
    },
    "Summer23": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
    },
    "Summer23BPix": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
    },
    "Summer24": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
    },
    "prompt_dataMC": {
        "data": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
        "mc": [
            "goodVertices",
            "globalSuperTightHalo2016Filter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "BadPFMuonFilter",
            "BadPFMuonDzFilter",
            "hfNoisyHitsFilter",
            "eeBadScFilter",
        ],
    },
}
