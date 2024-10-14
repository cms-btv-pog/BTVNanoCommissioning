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


## Jet pu ID not exist in Winter22Run3 sample
def jet_id(events, campaign):
    if campaign == "Rereco17_94X":
        jetmask = (
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (events.Jet.jetId >= 5)
            & ((events.Jet.pt > 50) | (events.Jet.puId >= 7))
        )
    else:
        jetmask = (
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (events.Jet.jetId >= 5)
        )
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
    else:
        jet_mask = (jets[f"btag{tagger}CvB"] > WP[tagger]["c"][wp][1]) & (
            jets[f"btag{tagger}CvL"] > WP[tagger]["c"][wp][0]
        )
    return jet_mask


import os, correctionlib

btag_wp_dict = {}


def wp_dict(year, campaign):
    """
    year :
    """
    global btag_wp_dict
    cache_key = f"{year}_{campaign}"

    if cache_key in btag_wp_dict:
        return btag_wp_dict[cache_key]
    name_map = {
        "deepJet": "DeepFlav",
        "robustParticleTransformer": "RobustParTAK4",
        "particleNet": "PNet",
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
            bwp = btag[tagger].inputs[0].description.split("/")
            # Get b WPs
            wps_dict[name_map[tagger.replace("_wp_values", "")]]["b"] = {
                wp: btag[tagger].evaluate(wp) for wp in bwp
            }
            cwp = ctag[tagger].inputs[0].description.split("/")
            wps_dict[name_map[tagger.replace("_wp_values", "")]]["c"] = {
                wp: [ctag[tagger].evaluate(wp, "CvL"), ctag[tagger].evaluate(wp, "CvB")]
                for wp in cwp
            }  # [CvL, CvB]
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
}
