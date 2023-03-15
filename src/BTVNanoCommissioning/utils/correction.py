import importlib.resources
import gzip
import pickle
import contextlib
import cloudpickle
import os
import numpy as np
import awkward as ak
from coffea.lookup_tools import extractor

from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
import correctionlib
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.cTagSFReader import getSF


def load_SF(campaign):
    correction_map = {}
    for SF in correction_config[campaign].keys():
        if SF == "JME" or SF == "lumiMask":
            continue
        if SF == "PU":
            ## Check whether files in jsonpog-integration exist
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/LUM/{campaign}"
            ):
                correction_map["PU"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/LUM/{campaign}/puWeights.json.gz"
                )
            ## Otherwise custom files
            else:
                _pu_path = f"BTVNanoCommissioning.data.PU.{campaign}"
                with importlib.resources.path(
                    _pu_path, correction_config[campaign]["PU"]
                ) as filename:
                    if str(filename).endswith(".pkl.gz"):
                        with gzip.open(filename) as fin:
                            correction_map["PU"] = cloudpickle.load(fin)[
                                "2017_pileupweight"
                            ]
                    elif str(filename).endswith(".histo.root"):
                        ext = extractor()
                        ext.add_weight_sets([f"* * {filename}"])
                        ext.finalize()
                        correction_map["PU"] = ext.make_evaluator()["PU"]

        elif SF == "BTV":
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/BTV/{campaign}"
            ):
                correction_map["btag"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/BTV/{campaign}/btagging.json.gz"
                )
                correction_map["ctag"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/BTV/{campaign}/ctagging.json.gz"
                )
            else:
                correction_map["btag"] = {}
                correction_map["ctag"] = {}
                _btag_path = f"BTVNanoCommissioning.data.BTV.{campaign}"
                for tagger in correction_config[campaign]["BTV"]:
                    with importlib.resources.path(
                        _btag_path, correction_config[campaign]["BTV"][tagger]
                    ) as filename:
                        if "B" in tagger:
                            correction_map["btag"][tagger] = BTagScaleFactor(
                                filename,
                                BTagScaleFactor.RESHAPE,
                                methods="iterativefit,iterativefit,iterativefit",
                            )
                        else:
                            if campaign == "Rereco17_94X":
                                correction_map["ctag"][tagger] = (
                                    "BTVNanoCommissioning/data/BTV/"
                                    + campaign
                                    + "/"
                                    + correction_config[campaign]["BTV"][tagger]
                                )
                            else:
                                correction_map["ctag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )

        elif SF == "LSF":
            correction_map["MUO_cfg"] = {
                mu: f
                for mu, f in correction_config[campaign]["LSF"].items()
                if "mu" in mu
            }
            correction_map["EGM_cfg"] = {
                e: f
                for e, f in correction_config[campaign]["LSF"].items()
                if "ele" in e
            }
            ## Muon
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/MUO/{campaign}"
            ) and os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/EGM/{campaign}"
            ):
                correction_map["MUO"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/MUO/{campaign}/muon_Z.json.gz"
                )
                correction_map["EGM"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/EGM/{campaign}/electron.json.gz"
                )
            ### Check if any custom corrections needed
            # FIXME: (some low pT muons not supported in jsonpog-integration at the moment)

            if (
                ".json" in "\t".join(list(correction_config[campaign]["LSF"].values()))
                or ".txt"
                in "\t".join(list(correction_config[campaign]["LSF"].values()))
                or ".root"
                in "\t".join(list(correction_config[campaign]["LSF"].values()))
            ):
                _mu_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    inputs, real_paths = [
                        k
                        for k in correction_map["MUO_cfg"].keys()
                        if ".json" in correction_map["MUO_cfg"][k]
                        or ".txt" in correction_map["MUO_cfg"][k]
                        or ".root" in correction_map["MUO_cfg"][k]
                    ], [
                        stack.enter_context(importlib.resources.path(_mu_path, f))
                        for f in correction_map["MUO_cfg"].values()
                        if ".json" in f or ".txt" in f or ".root" in f
                    ]

                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(inputs, real_paths)
                            if ".json" in str(file)
                            or ".txt" in str(file)
                            or ".root" in str(file)
                        ]
                    )
                ext.finalize()
                correction_map["MUO_custom"] = ext.make_evaluator()

                _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    inputs, real_paths = [
                        k
                        for k in correction_map["EGM_cfg"].keys()
                        if ".json" in correction_map["EGM_cfg"][k]
                        or ".txt" in correction_map["EGM_cfg"][k]
                        or ".root" in correction_map["EGM_cfg"][k]
                    ], [
                        stack.enter_context(importlib.resources.path(_ele_path, f))
                        for f in correction_map["EGM_cfg"].values()
                        if ".json" in f or ".txt" in f or ".root" in f
                    ]
                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(inputs, real_paths)
                            if ".json" in str(file)
                            or ".txt" in str(file)
                            or ".root" in str(file)
                        ]
                    )
                ext.finalize()
                correction_map["EGM_custom"] = ext.make_evaluator()

    return correction_map


def load_lumi(campaign):
    _lumi_path = "BTVNanoCommissioning.data.lumiMasks"
    with importlib.resources.path(
        _lumi_path, correction_config[campaign]["lumiMask"]
    ) as filename:
        return LumiMask(filename)


## MET filters
met_filters = {
    "2016preVFP_UL": {
        "data": [
            "goodVertices",
            "globalSuperTightHaloUL16Filter",
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
    "2016postVFP_UL": {
        "data": [
            "goodVertices",
            "globalSuperTightHaloUL16Filter",
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
}


##JEC
# FIXME: would be nicer if we can move to correctionlib in the future together with factory and workable
def load_jmefactory(campaign):
    _jet_path = f"BTVNanoCommissioning.data.JME.{campaign}"
    with importlib.resources.path(
        _jet_path, correction_config[campaign]["JME"]
    ) as filename:
        with gzip.open(filename) as fin:
            jmestuff = cloudpickle.load(fin)
    return jmestuff


def add_jec_variables(jets, event_rho):
    jets["pt_raw"] = (1 - jets.rawFactor) * jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor) * jets.mass
    if hasattr(jets, "genJetIdxG"):
        jets["pt_gen"] = ak.values_astype(
            ak.fill_none(jets.matched_gen.pt, 0), np.float32
        )
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets


## PU weight
def puwei(correct_map, nPU):
    if "correctionlib" in str(type(correct_map["PU"])):
        return correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
            nPU, "nominal"
        )
    else:
        return correct_map["PU"](nPU)


def btagSFs(jet, correct_map, SFtype, syst="central"):
    if "correctionlib" in str(type(correct_map["btag"])):
        if SFtype.endswith("C"):
            ## hard coded
            if SFtype == "DeepJetC":
                return correct_map["ctag"]["deepJet_shape"].evaluate(
                    syst, jet.hadronFlavour, jet.btagDeepFlavCvL, jet.btagDeepFlavCvB
                )
            if SFtype == "DeepCSVC":
                return correct_map["ctag"]["deepCSV_shape"].evaluate(
                    syst, jet.hadronFlavour, jet.btagDeepCvL, jet.btagDeepCvB
                )
        else:
            if SFtype == "DeepJetB":
                return correct_map["btag"]["deepJet_shape"].evaluate(
                    syst, jet.hadronFlavour, abs(jet.eta), jet.pt, jet.btagDeepFlavB
                )
            if SFtype == "DeepCSVB":
                return correct_map["btag"]["deepCSV_shape"].evaluate(
                    syst, jet.hadronFlavour, abs(jet.eta), jet.pt, jet.btagDeepB
                )
    else:
        if SFtype.endswith("C"):
            if SFtype == "DeepJetC":
                return getSF(
                    jet.hadronFlavour,
                    jet.btagDeepFlavCvL,
                    jet.btagDeepFlavCvB,
                    correct_map["ctag"][SFtype],
                    syst,
                )
            else:
                return getSF(
                    jet.hadronFlavour,
                    jet.btagDeepCvL,
                    jet.btagDeepCvB,
                    correct_map["ctag"][SFtype],
                    syst,
                )

        else:
            if SFtype == "DeepJetB":
                return correct_map["btag"][SFtype].eval(
                    syst,
                    jet.hadronFlavour,
                    abs(jet.eta),
                    jet.pt,
                    jet.btagDeepFlavB,
                )
            else:
                return correct_map["btag"][SFtype].eval(
                    syst,
                    jet.hadronFlavour,
                    abs(jet.eta),
                    jet.pt,
                    jet.btagDeepB,
                )


### Lepton SFs
def eleSFs(ele, correct_map, isHLT=False):
    ele_eta = ele.eta
    ele_pt = np.where(ele.pt < 20, 20.0, ele.pt)
    weight = 1.0
    for sf in correct_map["EGM_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        sf_type = sf[: sf.find(" ")]
        if "correctionlib" in str(type(correct_map["EGM"])):
            if "Reco" in sf:
                weight = weight * correct_map["EGM"][
                    list(correct_map["EGM"].keys())[0]
                ].evaluate(sf[sf.find(" ") + 1 :], "sf", "RecoAbove20", ele_eta, ele_pt)
            else:
                weight = weight * correct_map["EGM"][
                    list(correct_map["EGM"].keys())[0]
                ].evaluate(
                    sf[sf.find(" ") + 1 :],
                    "sf",
                    correct_map["EGM_cfg"][sf],
                    ele_eta,
                    ele_pt,
                )
        else:
            if "ele_Trig" in sf:
                weight = weight * correct_map["EGM_custom"][sf_type](ele_pt)
            elif "ele" in sf:
                weight = weight * correct_map["EGM_custom"][sf_type](ele_eta, ele_pt)
    return weight


def muSFs(mu, correct_map, isHLT=False):
    mu_eta = np.abs(mu.eta)
    mu_pt = mu.pt
    weight = 1.0
    sfs = 1.0
    for sf in correct_map["MUO_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        mask = mu_pt > 15.0
        if "low" not in sf:
            mu_pt = np.where(mu.pt < 15.0, 15.0, mu.pt)
        else:
            mu_pt = np.where(mu.pt >= 15.0, 15.0, mu.pt)
        sf_type = sf[: sf.find(" ")]
        if (
            "correctionlib" in str(type(correct_map["MUO"]))
            and "MUO_custom" in correct_map
        ):
            if "ID" in sf:
                if "low" in sf:
                    sfs = np.where(
                        ~mask,
                        correct_map["MUO_custom"][
                            "mu_ID_lowNUM_TightID_DEN_TrackerMuons/abseta_pt_value"
                        ](mu_eta, mu_pt),
                        1.0,
                    )
                else:
                    sfs = np.where(
                        mask,
                        correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                            sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "sf"
                        ),
                        1.0,
                    )
            elif "Reco" in sf:
                if "low" in sf:
                    sfs = np.where(
                        ~mask,
                        correct_map["MUO_custom"][
                            "mu_Reco_lowNUM_TrackerMuons_DEN_genTracks/abseta_pt_value"
                        ](mu_eta, mu_pt),
                        1.0,
                    )
                else:
                    sfs = np.where(
                        mu_pt > 15.0,
                        correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                            sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "sf"
                        ),
                        1.0,
                    )

            weight = weight * sfs
        elif "correctionlib" in str(type(correct_map["MUO"])):
            weight = weight * correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "sf"
            )
        else:
            if "mu" in sf:
                weight = weight * correct_map["MUO_custom"][sf_type](mu_eta, mu_pt)

    return weight
