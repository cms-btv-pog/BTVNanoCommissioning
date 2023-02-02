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
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor
import correctionlib
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.cTagSFReader import getSF


def load_SF(campaign):
    correction_map = {}
    for SF in correction_config[campaign].keys():
        if SF == "JME" or SF == "lumiMask":
            continue
        if SF == "PU":
            ## load from correction config
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/LUM/{campaign}"
            ):
                correction_map["PU"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/LUM/{campaign}/puWeights.json.gz"
                )
            ##
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
                if campaign == "Rereco17_94X":
                    _btag_path = f"BTVNanoCommissioning.data.BTV.{campaign}"
                    for tagger in correction_config[campaign]["BTV"]:
                        if tagger == "DeepCSVB" or tagger == "DeepJetB":
                            with importlib.resources.path(
                                _btag_path, correction_config[campaign]["BTV"][tagger]
                            ) as filename:
                                correction_map["btag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )
                        elif tagger == "DeepCSVC" or tagger == "DeepJetC":
                            correction_map["ctag"][tagger] = (
                                "BTVNanoCommissioning/data/BTV/"
                                + campaign
                                + "/"
                                + correction_config[campaign]["BTV"][tagger]
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
            else:
                _mu_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    real_paths = [
                        stack.enter_context(importlib.resources.path(_mu_path, f))
                        for f in correction_config[campaign]["LSF"].values()
                    ]
                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(
                                correction_config[campaign]["LSF"].keys(), real_paths
                            )
                            if "mu" in paths
                        ]
                    )
                ext.finalize()
                correction_map["MUO"] = ext.make_evaluator()

                _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    real_paths = [
                        stack.enter_context(importlib.resources.path(_ele_path, f))
                        for f in correction_config[campaign]["LSF"].values()
                    ]
                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(
                                correction_config[campaign]["LSF"].keys(), real_paths
                            )
                            if "ele" in paths
                        ]
                    )
                ext.finalize()
                correction_map["EGM"] = ext.make_evaluator()
    return correction_map


def load_lumi(campaign):
    _lumi_path = "BTVNanoCommissioning.data.lumiMasks"
    with importlib.resources.path(
        _lumi_path, correction_config[campaign]["lumiMask"]
    ) as filename:
        return LumiMask(filename)


## MET filters
met_filters = {
    "UL16": {
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
    "UL17": {
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
    "UL18": {
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
def load_pu(campaign, path):
    _pu_path = f"BTVNanoCommissioning.data.PU.{campaign}"
    with importlib.resources.path(_pu_path, path) as filename:
        if str(filename).endswith(".pkl.gz"):
            with gzip.open(filename) as fin:
                compiled = cloudpickle.load(fin)
        elif str(filename).endswith(".histo.root"):
            ext = extractor()
            ext.add_weight_sets([f"* * {filename}"])
            ext.finalize()
            compiled = ext.make_evaluator()

    return compiled


## BTag SFs
def load_BTV(campaign, path, tagger):
    if campaign == "Rereco17_94X":
        _btag_path = f"BTVNanoCommissioning.data.BTV.{campaign}"
        if tagger == "DeepCSVB" or tagger == "DeepJetB":
            with importlib.resources.path(_btag_path, path[tagger]) as filename:
                deepsf = BTagScaleFactor(
                    filename,
                    BTagScaleFactor.RESHAPE,
                    methods="iterativefit,iterativefit,iterativefit",
                )
        elif tagger == "DeepCSVC" or tagger == "DeepJetC":
            deepsf = "BTVNanoCommissioning/data/BTV/" + campaign + "/" + path[tagger]
        return deepsf


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
            if "DeepJetC":
                return correct_map["ctag"]["deepJet_shape"].evaluate(
                    syst, jet.hadronFlavour, jet.btagDeepFlavCvL, jet.btagDeepFlavCvB
                )
            else:
                return correct_map["ctag"]["deepCSV_shape"].evaluate(
                    syst, jet.hadronFlavour, jet.btagDeepCvL, jet.btagDeepCvB
                )
        else:
            if "DeepJetB":
                return correct_map["btag"]["deepJet_shape"].evaluate(
                    syst, jet.hadronFlavour, abs(jet.eta), jet.pt, jet.btagDeepFlavB
                )
            if "DeepCSVB":
                return correct_map["btag"]["deepJet_shape"].evaluate(
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
            if SFtype == "DeepJetC":
                return correct_map["btag"][SFtype].eval(
                    syst, jet.hadronFlavour, abs(jet.eta), jet.pt, jet.btagDeepFlavB
                )
            else:
                return correct_map["btag"][SFtype].eval(
                    syst, jet.hadronFlavour, abs(jet.eta), jet.pt, jet.btagDeepB
                )


### Lepton SFs
def eleSFs(ele, correct_map):
    ele_eta = ak.fill_none(ele.eta, 0.0)
    ele_pt = ak.fill_none(ele.pt, 20.0)
    weight = 1.0
    for sf in correct_map["EGM_cfg"].keys():
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
                weight = weight * correct_map["EGM"][sf_type](ele_pt)
            elif "ele" in sf:
                weight = weight * correct_map["EGM"][sf_type](ele_eta, ele_pt)
    return weight


def muSFs(mu, correct_map):
    mu_eta = np.abs(ak.fill_none(mu.eta, 0.0))

    weight = 1.0
    for sf in correct_map["MUO_cfg"].keys():
        min_pt = 15.0
        if "HLT" in sf:
            min_pt = 29.0
        mu_pt = ak.fill_none(mu.pt, min_pt)
        mu_pt = np.where(mu_pt < min_pt, min_pt, mu_pt)

        sf_type = sf[: sf.find(" ")]
        if "correctionlib" in str(type(correct_map["MUO"])):
            weight = weight * correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "sf"
            )
        else:
            if "mu" in sf:
                weight = weight * correct_map["MUO"][sf_type](mu_eta, mu_pt)

    return weight
