import importlib.resources
import gzip
import pickle
import contextlib
import cloudpickle
import os
import numpy as np
import awkward as ak
from coffea.lookup_tools import extractor, txt_converters, rochester_lookup

from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
import correctionlib

from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.AK4_parameters import correction_config as config


def load_SF(campaign, syst=False):
    correction_map = {}
    for SF in config[campaign].keys():
        if SF == "lumiMask":
            continue
        ## pileup weight
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
                    _pu_path, config[campaign]["PU"]
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
        ## btag weight
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
                for tagger in config[campaign]["BTV"]:
                    with importlib.resources.path(
                        _btag_path, config[campaign]["BTV"][tagger]
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
                                    + config[campaign]["BTV"][tagger]
                                )
                            else:
                                correction_map["ctag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )
        ## lepton SFs
        elif SF == "LSF":
            correction_map["MUO_cfg"] = {
                mu: f for mu, f in config[campaign]["LSF"].items() if "mu" in mu
            }
            correction_map["EGM_cfg"] = {
                e: f for e, f in config[campaign]["LSF"].items() if "ele" in e
            }
            ## Muon
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/MUO/{campaign}"
            ):
                correction_map["MUO"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/MUO/{campaign}/muon_Z.json.gz"
                )
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/EGM/{campaign}"
            ):
                correction_map["EGM"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/EGM/{campaign}/electron.json.gz"
                )
            ### Check if any custom corrections needed
            # FIXME: (some low pT muons not supported in jsonpog-integration at the moment)

            if (
                ".json" in "\t".join(list(config[campaign]["LSF"].values()))
                or ".txt" in "\t".join(list(config[campaign]["LSF"].values()))
                or ".root" in "\t".join(list(config[campaign]["LSF"].values()))
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

                    inputs = [
                        i.split(" ")[0] + " *" if "_low" in i else i for i in inputs
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
                    if syst:
                        ext.add_weight_sets(
                            paths.split(" ")[0]
                            + "_error "
                            + paths.split(" ")[1]
                            + "_error "
                            + file
                            for paths, file in zip(inputs, real_paths)
                            if ".root" in str(file)
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
                    if syst:
                        ext.add_weight_sets(
                            paths.split(" ")[0]
                            + "_error "
                            + paths.split(" ")[1]
                            + "_error "
                            + file
                            for paths, file in zip(inputs, real_paths)
                            if ".root" in str(file)
                        )
                ext.finalize()
                correction_map["EGM_custom"] = ext.make_evaluator()
        ## rochester muon momentum correction
        elif SF == "roccor":
            if "2016postVFP_UL" == campaign:
                filename = "RoccoR2016bUL.txt"
            elif "2016preVFP_UL" in campaign:
                filename = "RoccoR2016aUL.txt"
            elif "2017_UL" in campaign:
                filename = "RoccoR2017UL.txt"
            if "2018_UL" in campaign:
                filename = "RoccoR2018UL.txt"

            full_path = "src/BTVNanoCommissioning/data/LSF/roccor/" + filename
            rochester_data = txt_converters.convert_rochester_file(
                full_path, loaduncs=True
            )
            correction_map["roccor"] = rochester_lookup.rochester_lookup(rochester_data)
        elif SF == "JME":
            correction_map["JME"] = load_jmefactory(campaign)
        elif SF == "JMAR":
            if os.path.exists(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/JME/{campaign}/jmar.json.gz"
            ):
                correction_map["JMAR_cfg"] = {
                    j: f for j, f in config[campaign]["JMAR"].items()
                }
                correction_map["JMAR"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/JME/{campaign}/jmar.json.gz"
                )

    return correction_map


def load_lumi(campaign):
    _lumi_path = "BTVNanoCommissioning.data.lumiMasks"
    with importlib.resources.path(_lumi_path, config[campaign]["lumiMask"]) as filename:
        return LumiMask(filename)


## MET filters
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
}

import matplotlib.pyplot as plt

ext_jetvetomap = extractor()
ext_jetvetomap.add_weight_sets(
    [
        "RunCD jetvetomap src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunCD_v1.histo.root",
        "RunE jetvetomap src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunE_v1.histo.root",
    ]
)

ext_jetvetomap.finalize()
jetvetomap = ext_jetvetomap.make_evaluator()


def jetveto(events):
    if (
        "Run2022C" in events.metadata["dataset"]
        or "Run2022D" in events.metadata["dataset"]
    ):
        return ak.where(
            jetvetomap["RunCD"](events.Jet.phi, events.Jet.eta) > 0,
            ak.ones_like(events.Jet.eta),
            ak.zeros_like(events.Jet.eta),
        )
    elif "Run2022E" in events.metadata["dataset"]:
        return ak.where(
            jetvetomap["RunE"](events.Jet.phi, events.Jet.eta) > 0,
            ak.ones_like(events.Jet.eta),
            ak.zeros_like(events.Jet.eta),
        )


##JEC
# FIXME: would be nicer if we can move to correctionlib in the future together with factory and workable
def load_jmefactory(campaign):
    _jet_path = f"BTVNanoCommissioning.data.JME.{campaign}"
    with importlib.resources.path(_jet_path, config[campaign]["JME"]) as filename:
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
    else:
        jets["pt_gen"] = ak.zeros_like(jets.pt)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets


## JERC
def JME_shifts(
    shifts,
    correct_map,
    events,
    campaign,
    isRealData,
    systematic=False,
    exclude_jetveto=False,
):
    dataset = events.metadata["dataset"]
    if isRealData:
        if "2016preVFP_UL" == campaign:
            if "2016B" in dataset or "2016C" in dataset or "2016D" in dataset:
                jecname = "BCD"
            elif "2016E" in dataset or "2016F" in dataset:
                jecname = "EF"
        elif "2016postVFP_UL" == campaign:
            jecname = "FGH"
        elif campaign == "Rereco17_94X":
            jecname = ""
        elif "un" in dataset:
            jecname = dataset[dataset.find("un") + 6]
        else:
            print("No valid jec name")
            raise NameError
        jecname = "data" + jecname
    else:
        jecname = "mc"
    jets = correct_map["JME"]["jet_factory"][jecname].build(
        add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
        lazy_cache=events.caches[0],
    )
    met = correct_map["JME"]["met_factory"].build(events.MET, jets, {})
    ## HEM 18 issue
    if isRealData:
        if "2018" in events.metadata["dataset"]:
            _runid = events.run >= 319077
            j_mask = ak.where(
                _runid
                & (jets.phi > -1.57)
                & (jets.phi < -0.87)
                & (jets.eta > -2.5)
                & (jets.eta < -1.3),
                0.8,
                1,
            )
            j_high_eta_mask = ak.where(
                _runid
                & (jets.phi > -1.57)
                & (jets.phi < -0.87)
                & (jets.eta > -3.0)
                & (jets.eta < -2.5),
                0.65,
                1,
            )

            for var in ["mass", "pt"]:
                jets[var] = j_mask * j_high_eta_mask * jets[var]
        if (
            "Run2022C" in events.metadata["dataset"]
            or "Run2022D" in events.metadata["dataset"]
            or "Run2022E" in events.metadata["dataset"]
        ) and not exclude_jetveto:
            jets["pt"] = ak.where(jets.veto == 0, jets.pt, 0.0)
    shifts += [({"Jet": jets, "MET": met}, None)]

    ## systematics
    if not isRealData:
        if systematic != False:
            if systematic == "split":
                for jes in met.fields:
                    if "JES" not in jes or "Total" in jes:
                        continue

                    shifts += [
                        (
                            {
                                "Jet": jets[jes]["up"],
                                "MET": met[jes]["up"],
                            },
                            f"{jes}Up",
                        ),
                        (
                            {
                                "Jet": jets[jes]["down"],
                                "MET": met[jes]["down"],
                            },
                            f"{jes}Down",
                        ),
                    ]

            else:
                shifts += [
                    (
                        {
                            "Jet": jets.JES_Total.up,
                            "MET": met.JES_Total.up,
                        },
                        "JESUp",
                    ),
                    (
                        {
                            "Jet": jets.JES_Total.down,
                            "MET": met.JES_Total.down,
                        },
                        "JESDown",
                    ),
                ]
            shifts += [
                (
                    {
                        "Jet": jets,
                        "MET": met.MET_UnclusteredEnergy.up,
                    },
                    "UESUp",
                ),
                (
                    {
                        "Jet": jets,
                        "MET": met.MET_UnclusteredEnergy.down,
                    },
                    "UESDown",
                ),
                (
                    {
                        "Jet": jets.JER.up,
                        "MET": met.JER.up,
                    },
                    "JERUp",
                ),
                (
                    {
                        "Jet": jets.JER.down,
                        "MET": met.JER.down,
                    },
                    "JERDown",
                ),
            ]
    return shifts


## Muon Rochester correction
def Roccor_shifts(shifts, correct_map, events, isRealData, systematic=False):
    mu = events.Muon
    if isRealData:
        SF = correct_map["roccor"].kScaleDT(
            events.Muon.charge, events.Muon.pt, events.Muon.eta, events.Muon.phi
        )

    else:
        hasgen = ~np.isnan(ak.fill_none(events.Muon.matched_gen.pt, np.nan))
        mc_kspread = correct_map["roccor"].kSpreadMC(
            events.Muon.charge[hasgen],
            events.Muon.pt[hasgen],
            events.Muon.eta[hasgen],
            events.Muon.phi[hasgen],
            events.Muon.matched_gen.pt[hasgen],
        )
        mc_rand = np.random.rand(len(ak.flatten(events.Muon.pt, axis=1)))
        mc_rand = ak.unflatten(mc_rand, ak.num(events.Muon.pt))
        mc_ksmear = correct_map["roccor"].kSmearMC(
            events.Muon.charge[~hasgen],
            events.Muon.pt[~hasgen],
            events.Muon.eta[~hasgen],
            events.Muon.phi[~hasgen],
            events.Muon.nTrackerLayers[~hasgen],
            mc_rand[~hasgen],
        )
        SF = np.array(ak.flatten(ak.ones_like(events.Muon.pt)))
        hasgen_flat = np.array(ak.flatten(hasgen))
        SF[hasgen_flat] = np.array(ak.flatten(mc_kspread))
        SF[~hasgen_flat] = np.array(ak.flatten(mc_ksmear))
        SF = ak.unflatten(SF, ak.num(events.Muon.pt))

    mu["pt"] = SF * events.Muon.pt
    # add rochester correction to shift
    for i in range(len(shifts)):
        shifts[i][0]["Muon"] = mu

    if systematic:
        if isRealData:
            err = correct_map["roccor"].kScaleDTerror(
                events.Muon.charge, events.Muon.pt, events.Muon.eta, events.Muon.phi
            )
        else:
            mc_errspread = correct_map["roccor"].kSpreadMCerror(
                events.Muon.charge[hasgen],
                events.Muon.pt[hasgen],
                events.Muon.eta[hasgen],
                events.Muon.phi[hasgen],
                events.Muon.matched_gen.pt[hasgen],
            )
            mc_errsmear = correct_map["roccor"].kSmearMCerror(
                events.Muon.charge[~hasgen],
                events.Muon.pt[~hasgen],
                events.Muon.eta[~hasgen],
                events.Muon.phi[~hasgen],
                events.Muon.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
            )
            err = np.array(ak.flatten(ak.ones_like(events.Muon.pt)))
            err[hasgen_flat] = np.array(ak.flatten(mc_errspread))
            err[~hasgen_flat] = np.array(ak.flatten(mc_errsmear))
            err = ak.unflatten(err, ak.num(events.Muon.pt))
        muup, mudown = events.Muon, events.Muon
        muup["pt"] = (SF + err) * events.Muon.pt
        mudown["pt"] = (SF - err) * events.Muon.pt
        shifts += [
            (
                {"Jet": shifts[0][0]["Jet"], "MET": shifts[0][0]["MET"], "Muon": muup},
                "RoccorUp",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": mudown,
                },
                "RoccorDown",
            )
        ]
    return shifts


## PU weight
def puwei(nPU, correct_map, weights, syst=False):
    if "correctionlib" in str(type(correct_map["PU"])):
        if syst:
            return weights.add(
                "puweight",
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "nominal"
                ),
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "up"
                ),
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "down"
                ),
            )
        else:
            return weights.add(
                "puweight",
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "nominal"
                ),
            )
    else:
        if syst:
            weights.add(
                "puweight",
                correct_map["PU"](nPU),
                correct_map["PUup"](nPU),
                correct_map["PUdown"](nPU),
            )
        else:
            weights.add("puweight", correct_map["PU"](nPU))


def btagSFs(jet, correct_map, weights, SFtype, syst=False):
    if SFtype == "DeepJetC" or SFtype == "DeepCSVC":
        systlist = [
            "Extrap",
            "Interp",
            "LHEScaleWeight_muF",
            "LHEScaleWeight_muR",
            "PSWeightFSR",
            "PSWeightISR",
            "PUWeight",
            "Stat",
            "XSec_BRUnc_DYJets_b",
            "XSec_BRUnc_DYJets_c",
            "XSec_BRUnc_WJets_c",
            "jer",
            "jesTotal",
        ]
    elif SFtype == "DeepJetB" or SFtype == "DeepCSVB":
        systlist = [
            "hf",
            "lf",
            "cferr1",
            "cferr2",
            "hfstat1",
            "hfstat2",
            "lfstats1",
            "lfstats2",
        ]
    sfs_up_all, sfs_down_all = {}, {}
    alljet = jet if jet.ndim > 1 else ak.singletons(jet)
    for i, sys in enumerate(systlist):
        sfs, sfs_down, sfs_up = (
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
        )
        for nj in range(ak.num(alljet.pt)[0]):
            jet = alljet[:, nj]
            masknone = ak.is_none(jet.pt)
            jet.btagDeepFlavCvL = ak.fill_none(jet.btagDeepFlavCvL, 0.0)
            jet.btagDeepFlavCvB = ak.fill_none(jet.btagDeepFlavCvB, 0.0)
            jet.btagDeepCvL = ak.fill_none(jet.btagDeepCvL, 0.0)
            jet.btagDeepCvB = ak.fill_none(jet.btagDeepCvB, 0.0)
            jet.hadronFlavour = ak.fill_none(jet.hadronFlavour, 0)
            if "correctionlib" in str(type(correct_map["ctag"])):
                if SFtype == "DeepJetC":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepJet_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepFlavCvL,
                            jet.btagDeepFlavCvB,
                        ),
                    )
                    if syst:
                        tmp_sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["ctag"]["deepJet_shape"].evaluate(
                                f"up_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                        tmp_sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["ctag"]["deepJet_shape"].evaluate(
                                f"down_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                if SFtype == "DeepCSVC":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepCSV_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_up = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepCSV_shape"].evaluate(
                            f"up_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_down = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepCSV_shape"].evaluate(
                            f"down_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
            if "correctionlib" in str(type(correct_map["btag"])):
                if SFtype == "DeepJetB":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepJet_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepFlavCvL,
                            jet.btagDeepFlavCvB,
                        ),
                    )
                    if syst:
                        tmp_sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["btag"]["deepJet_shape"].evaluate(
                                f"up_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                        tmp_sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["btag"]["deepJet_shape"].evaluate(
                                f"down_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                if SFtype == "DeepCSVB":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepCSV_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_up = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepCSV_shape"].evaluate(
                            f"up_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_down = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepCSV_shape"].evaluate(
                            f"down_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )

            sfs = sfs * tmp_sfs
            if syst:
                sfs_up = sfs_up * tmp_sfs_up
                sfs_down = sfs_down * tmp_sfs_down

        if i == 0 and syst == False:
            weights.add(SFtype, sfs)
            break
        else:
            sfs_up_all[sys] = sfs_up
            sfs_down_all[sys] = sfs_down
    if syst == True:
        weights.add_multivariation(
            SFtype,
            sfs,
            systlist,
            np.array(list(sfs_up_all.values())),
            np.array(list(sfs_down_all.values())),
        )
    return weights


### Lepton SFs
def eleSFs(ele, correct_map, weights, syst=True, isHLT=False):
    allele = ele if ele.ndim > 1 else ak.singletons(ele)
    for sf in correct_map["EGM_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        sf_type = sf[: sf.find(" ")]
        if "low" in sf:
            continue
        for nele in range(ak.num(allele.pt)[0]):
            ele = allele[:, nele]
            ele_eta = ak.fill_none(ele.eta, -2.5)
            ele_pt = ak.fill_none(ele.pt, 20)
            mask = ele.pt > 20.0
            masknone = ak.is_none(ele.pt)
            sfs_alle, sfs_alle_up, sfs_alle_down = (
                np.ones_like(allele[:, 0].pt),
                np.ones_like(allele[:, 0].pt),
                np.ones_like(allele[:, 0].pt),
            )

            if "correctionlib" in str(type(correct_map["EGM"])):
                ## Reco SF -splitted pT
                if "Reco" in sf:
                    ele_pt = np.where(ele.pt < 20.0, 20.0, ele.pt)
                    ele_pt_low = np.where(ele.pt >= 20.0, 19.9, ele.pt)

                    sfs_low = np.where(
                        (~mask) & (~masknone),
                        correct_map["EGM"][list(correct_map["EGM"].keys())[0]].evaluate(
                            sf[sf.find(" ") + 1 :],
                            "sf",
                            "RecoBelow20",
                            ele_eta,
                            ele_pt_low,
                        ),
                        1.0,
                    )
                    sfs = np.where(
                        mask & (~masknone),
                        correct_map["EGM"][list(correct_map["EGM"].keys())[0]].evaluate(
                            sf[sf.find(" ") + 1 :], "sf", "RecoAbove20", ele_eta, ele_pt
                        ),
                        sfs_low,
                    )
                    sfs = np.where(masknone, 1.0, sfs)

                    if syst:
                        sfs_up_low = np.where(
                            ~mask & ~masknone,
                            correct_map["EGM"][
                                list(correct_map["EGM"].keys())[0]
                            ].evaluate(
                                sf[sf.find(" ") + 1 :],
                                "sfup",
                                "RecoBelow20",
                                ele_eta,
                                ele_pt_low,
                            ),
                            0.0,
                        )
                        sfs_down_low = np.where(
                            ~mask & ~masknone,
                            correct_map["EGM"][
                                list(correct_map["EGM"].keys())[0]
                            ].evaluate(
                                sf[sf.find(" ") + 1 :],
                                "sfdown",
                                "RecoBelow20",
                                ele_eta,
                                ele_pt_low,
                            ),
                            0.0,
                        )
                        sfs_up = np.where(
                            mask & ~masknone,
                            correct_map["EGM"][
                                list(correct_map["EGM"].keys())[0]
                            ].evaluate(
                                sf[sf.find(" ") + 1 :],
                                "sfup",
                                "RecoAbove20",
                                ele_eta,
                                ele_pt,
                            ),
                            sfs_up_low,
                        )
                        sfs_down = np.where(
                            mask & ~masknone,
                            correct_map["EGM"][
                                list(correct_map["EGM"].keys())[0]
                            ].evaluate(
                                sf[sf.find(" ") + 1 :],
                                "sfdown",
                                "RecoAbove20",
                                ele_eta,
                                ele_pt,
                            ),
                            sfs_down_low,
                        )
                        sfs_up, sfs_down = np.where(masknone, 1.0, sfs_up), np.where(
                            masknone, 1.0, sfs_down
                        )
                ## Other files
                else:
                    sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["EGM"][list(correct_map["EGM"].keys())[0]].evaluate(
                            sf[sf.find(" ") + 1 :],
                            "sf",
                            correct_map["EGM_cfg"][sf],
                            ele_eta,
                            ele_pt,
                        ),
                    )

                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM"][
                                list(correct_map["EGM"].keys())[0]
                            ].evaluate(
                                sf[sf.find(" ") + 1 :],
                                "sfup",
                                correct_map["EGM_cfg"][sf],
                                ele_eta,
                                ele_pt,
                            ),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM"][
                                list(correct_map["EGM"].keys())[0]
                            ].evaluate(
                                sf[sf.find(" ") + 1 :],
                                "sfdown",
                                correct_map["EGM_cfg"][sf],
                                ele_eta,
                                ele_pt,
                            ),
                        )
            else:
                if "ele_Trig" in sf:
                    sfs = np.where(
                        masknone, 1.0, correct_map["EGM_custom"][sf_type](ele_pt)
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_pt)
                            + correct_map["EGM_custom"][f"{sf_type}_error"](ele_pt),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_pt)
                            - correct_map["EGM_custom"][f"{sf_type}_error"](ele_pt),
                        )
                elif "ele" in sf:
                    sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["EGM_custom"][sf_type](ele_eta, ele_pt),
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_eta, ele_pt)
                            + correct_map["EGM_custom"][f"{sf_type}_error"](
                                ele_eta, ele_pt
                            ),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_eta, ele_pt)
                            - correct_map["EGM_custom"][f"{sf_type}_error"](
                                ele_eta, ele_pt
                            ),
                        )
            sfs_alle = sfs_alle * sfs
            if syst:
                sfs_alle_down = sfs_alle_down * sfs_down
                sfs_alle_up = sfs_alle_up * sfs_up

        if syst:
            weights.add(sf.split(" ")[0], sfs_alle, sfs_alle_up, sfs_alle_down)
        else:
            weights.add(sf.split(" ")[0], sfs_alle)
    return weights


def muSFs(mu, correct_map, weights, syst=False, isHLT=False):
    allmu = mu if mu.ndim > 1 else ak.singletons(mu)
    for sf in correct_map["MUO_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        if "low" in sf:
            continue
        sfs_allmu, sfs_allmu_up, sfs_allmu_down = (
            np.ones_like(allmu[:, 0].pt),
            np.ones_like(allmu[:, 0].pt),
            np.ones_like(allmu[:, 0].pt),
        )
        sf_type = sf[: sf.find(" ")]
        for nmu in range(ak.num(allmu.pt)[0]):
            mu = allmu[:, nmu]
            masknone = ak.is_none(mu.pt)

            mu_pt = np.where(mu.pt < 15, 15, mu.pt)
            mu_eta = np.where(np.abs(mu.eta) >= 2.4, 2.39, np.abs(mu.eta))
            mask = mu_pt > 30
            sfs = 1.0
            if "correctionlib" in str(type(correct_map["MUO"])):
                if "ID" in sf or "Reco" in sf:
                    mu_pt = ak.fill_none(np.where(mu.pt < 30, 30, mu.pt), 30)
                    mu_pt_low = ak.fill_none(np.where(mu.pt >= 30, 30, mu.pt), 30)
                    sfs_low = np.where(
                        ~mask & ~masknone,
                        correct_map["MUO_custom"][
                            f'{sf.split(" ")[0]}_low{correct_map["MUO_cfg"][sf]}/abseta_pt_value'
                        ](mu_eta, mu_pt_low),
                        1.0,
                    )

                    sfs = np.where(
                        mask & ~masknone,
                        correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                            sf.split(" ")[1], mu_eta, mu_pt, "sf"
                        ),
                        sfs_low,
                    )
                    sfs = np.where(masknone, 1.0, sfs)
                    sfs_forerr = sfs

                    if syst:
                        sfs_err_low = np.where(
                            ~mask & ~masknone,
                            correct_map["MUO_custom"][
                                f'{sf.split(" ")[0]}_low{correct_map["MUO_cfg"][sf]}/abseta_pt_error'
                            ](mu_eta, mu_pt_low),
                            0.0,
                        )
                        sfs_up = np.where(
                            mask & ~masknone,
                            correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                                sf.split(" ")[1], mu_eta, mu_pt, "systup"
                            ),
                            sfs_forerr + sfs_err_low,
                        )
                        sfs_down = np.where(
                            mask & ~masknone,
                            correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                                sf.split(" ")[1], mu_eta, mu_pt, "systdown"
                            ),
                            sfs_forerr - sfs_err_low,
                        )
                        sfs_up, sfs_down = np.where(masknone, 1.0, sfs_up), np.where(
                            masknone, 1.0, sfs_down
                        )

                else:
                    sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                            sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "sf"
                        ),
                    )

                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                                sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "systup"
                            ),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                                sf[sf.find(" ") + 1 :], mu_eta, mu_pt, "systdown"
                            ),
                        )
            else:
                if "mu" in sf:
                    sfs = np.where(
                        masknone, 1.0, correct_map["MUO_custom"][sf_type](mu_eta, mu_pt)
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            sfs
                            + correct_map["MUO_custom"][f"{sf_type}_error"](
                                mu_eta, mu_pt
                            ),
                        )
                        sf_down = np.where(
                            masknone,
                            1.0,
                            sfs
                            - correct_map["MUO_custom"][f"{sf_type}_error"](
                                mu_eta, mu_pt
                            ),
                        )

            sfs_allmu = sfs_allmu * sfs
            if syst:
                sfs_allmu_down = sfs_allmu_down * sfs_down
                sfs_allmu_up = sfs_allmu_up * sfs_up
        if syst:
            weights.add(sf.split(" ")[0], sfs_allmu, sfs_allmu_up, sfs_allmu_down)
        else:
            weights.add(sf.split(" ")[0], sfs_allmu)
    return weights


def jmar_sf(jet, correct_map, weights, syst=False):
    alljet = jet if jet.ndim > 1 else ak.singletons(jet)
    for sf in correct_map["JMAR_cfg"].keys():
        sfs_all, sfs_all_up, sfs_all_down = (
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
        )
        for njet in range(ak.num(alljet.pt)[0]):
            jet = alljet[:, njet]
            jet_eta = ak.fill_none(jet.eta, -2.5)
            jet_pt = ak.fill_none(jet.pt, 20)
            masknone = ak.is_none(jet.pt)
            # PU Jet ID applied only jet pT<50GeV
            if sf == "PUJetID_eff":
                jet_pt = np.where(jet.pt > 50, 20, jet.pt)
                masknone = (ak.is_none(jet.pt)) & (jet.pt > 50)
            if "correctionlib" in str(type(correct_map["JMAR"])):
                sfs = np.where(
                    masknone,
                    1.0,
                    correct_map["JMAR"][sf].evaluate(
                        jet_eta,
                        jet_pt,
                        "nom",
                        correct_map["JMAR_cfg"][sf],
                    ),
                )
                if syst:
                    sfs_up = np.where(
                        masknone,
                        1.0,
                        correct_map["JMAR"][sf].evaluate(
                            jet_eta,
                            jet_pt,
                            "up",
                            correct_map["JMAR_cfg"][sf],
                        ),
                    )
                    sfs_down = np.where(
                        masknone,
                        1.0,
                        correct_map["JMAR"][sf].evaluate(
                            jet_eta,
                            jet_pt,
                            "down",
                            correct_map["JMAR_cfg"][sf],
                        ),
                    )
                    sfs_all_up = sfs_up * sfs_all_up
                    sfs_all_down = sfs_down * sfs_all_down
            sfs_all = sfs * sfs_all
        if syst:
            weights.add(sf, sfs_all, sfs_all_up, sfs_all_down)
        else:
            weights.add(sf, sfs_all)


def add_pdf_weight(weights, pdf_weights):
    nom = np.ones(len(weights.weight()))
    up = np.ones(len(weights.weight()))
    down = np.ones(len(weights.weight()))

    # NNPDF31_nnlo_hessian_pdfas
    # https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas/NNPDF31_nnlo_hessian_pdfas.info
    if pdf_weights is not None and "306000 - 306102" in pdf_weights.__doc__:
        # Hessian PDF weights
        # Eq. 21 of https://arxiv.org/pdf/1510.03865v1.pdf
        arg = pdf_weights[:, 1:-2] - np.ones((len(weights.weight()), 100))
        summed = ak.sum(np.square(arg), axis=1)
        pdf_unc = np.sqrt((1.0 / 99.0) * summed)
        weights.add("PDF_weight", nom, pdf_unc + nom)

        # alpha_S weights
        # Eq. 27 of same ref
        as_unc = 0.5 * (pdf_weights[:, 102] - pdf_weights[:, 101])
        weights.add("aS_weight", nom, as_unc + nom)

        # PDF + alpha_S weights
        # Eq. 28 of same ref
        pdfas_unc = np.sqrt(np.square(pdf_unc) + np.square(as_unc))
        weights.add("PDFaS_weight", nom, pdfas_unc + nom)

    else:
        weights.add("aS_weight", nom, up, down)
        weights.add("PDF_weight", nom, up, down)
        weights.add("PDFaS_weight", nom, up, down)


# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
def top_pT_sf_formula(pt):
    return 0.103 * np.exp(-0.0118 * pt) - 0.000134 * pt + 0.973


def top_pT_reweighting(gen):
    #     """
    #     Apply this SF only to TTbar datasets! Updated to latest suggestion
    #     Documentation:
    #         - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
    #         - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
    #     """
    top = gen[(gen.pdgId == 6) & gen.hasFlags(["isLastCopy"])]
    anti_top = gen[(gen.pdgId == -6) & gen.hasFlags(["isLastCopy"])]
    return np.sqrt(
        top_pT_sf_formula(ak.flatten(top.pt, axis=-1))
        * top_pT_sf_formula(ak.flatten(anti_top.pt, axis=-1))
    )


# Jennet adds PS weights
# https://github.com/andrzejnovak/boostedhiggs/blob/master/boostedhiggs/corrections.py#L88-L108
def add_ps_weight(weights, ps_weights):
    nom = np.ones(len(weights.weight()))
    up_isr = np.ones(len(weights.weight()))
    down_isr = np.ones(len(weights.weight()))
    up_fsr = np.ones(len(weights.weight()))
    down_fsr = np.ones(len(weights.weight()))

    if ps_weights is not None:
        if len(ps_weights[0]) == 4:
            up_isr = ps_weights[:, 0]
            down_isr = ps_weights[:, 2]
            up_fsr = ps_weights[:, 1]
            down_fsr = ps_weights[:, 3]
        # else:
        #   warnings.warn(f"PS weight vector has length {len(ps_weights[0])}")

    weights.add("UEPS_ISR", nom, up_isr, down_isr)
    weights.add("UEPS_FSR", nom, up_fsr, down_fsr)


def add_scalevar_7pt(weights, lhe_weights):
    nom = np.ones(len(weights.weight()))

    if len(lhe_weights) > 0:
        if len(lhe_weights[0]) == 9:
            up = np.maximum.reduce(
                [
                    lhe_weights[:, 0],
                    lhe_weights[:, 1],
                    lhe_weights[:, 3],
                    lhe_weights[:, 5],
                    lhe_weights[:, 7],
                    lhe_weights[:, 8],
                ]
            )
            down = np.minimum.reduce(
                [
                    lhe_weights[:, 0],
                    lhe_weights[:, 1],
                    lhe_weights[:, 3],
                    lhe_weights[:, 5],
                    lhe_weights[:, 7],
                    lhe_weights[:, 8],
                ]
            )
        elif len(lhe_weights[0]) > 1:
            print("Scale variation vector has length ", len(lhe_weights[0]))
    else:
        up = np.ones(len(weights.weight()))
        down = np.ones(len(weights.weight()))

    weights.add("scalevar_7pt", nom, up, down)


def add_scalevar_3pt(weights, lhe_weights):
    nom = np.ones(len(weights.weight()))

    if len(lhe_weights) > 0:
        if len(lhe_weights[0]) == 9:
            up = np.maximum(lhe_weights[:, 0], lhe_weights[:, 8])
            down = np.minimum(lhe_weights[:, 0], lhe_weights[:, 8])
        elif len(lhe_weights[0]) > 1:
            print("Scale variation vector has length ", len(lhe_weights[0]))
    else:
        up = np.ones(len(weights.weight()))
        down = np.ones(len(weights.weight()))

    weights.add("scalevar_3pt", nom, up, down)
