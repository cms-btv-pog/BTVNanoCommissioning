import importlib.resources
import gzip
import pickle
import contextlib
import cloudpickle

import numpy as np
import awkward as ak
from coffea.lookup_tools.lookup_base import lookup_base
from coffea.lookup_tools import extractor

from coffea.analysis_tools import Weights
from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor

# FIXME make jsons campaing configurable
lumiMasks = {}
_lumi_path = "BTVNanoCommissioning.data.lumiMasks"
with importlib.resources.path(
    _lumi_path, "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
) as filename:
    lumiMasks["2016"] = LumiMask(filename)
with importlib.resources.path(
    _lumi_path, "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
) as filename:
    lumiMasks["2017"] = LumiMask(filename)
with importlib.resources.path(
    _lumi_path,
    "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt",
) as filename:
    lumiMasks["2018"] = LumiMask(filename)
with importlib.resources.path(
    _lumi_path,
    "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
) as filename:
    lumiMasks["UL16"] = LumiMask(filename)
with importlib.resources.path(
    _lumi_path,
    "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
) as filename:
    lumiMasks["UL17"] = LumiMask(filename)
with importlib.resources.path(
    _lumi_path,
    "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
) as filename:
    lumiMasks["UL18"] = LumiMask(filename)
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
def load_jetfactory(campaign, path):
    _jet_path = f"BTVNanoCommissioning.data.JME.{campaign}"
    with importlib.resources.path(_jet_path, path) as filename:
        with gzip.open(filename) as fin:
            jmestuff = cloudpickle.load(fin)

    jet_factory = jmestuff["jet_factory"]
    return jet_factory

def load_jmefactory(campaign, path):
    _jet_path = f"BTVNanoCommissioning.data.JME.{campaign}"
    with importlib.resources.path(_jet_path, path) as filename:
        with gzip.open(filename) as fin:
            jmestuff = cloudpickle.load(fin)

        #jet_factory = jmestuff["jet_factory"]
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


def load_metfactory(campaign, path):
    _jet_path = f"BTVNanoCommissioning.data.JME.{campaign}"
    with importlib.resources.path(_jet_path, path) as filename:
        with gzip.open(filename) as fin:
            jmestuff = cloudpickle.load(fin)

    met_factory = jmestuff["met_factory"]
    return met_factory


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


### Lepton SFs


# with contextlib.ExitStack() as stack:
#     # this would work even in zipballs but since extractor keys on file extension and
#     # importlib make a random tempfile, it won't work. coffea needs to enable specifying the type manually
#     # for now we run this whole module as $ python -m boostedhiggs.build_jec boostedhiggs/data/jec_compiled.pkl.gz
#     # so the compiled value can be loaded using the importlib tool in corrections.py
#     _ele_path = "BTVNanoCommissioning.data.LSF.Rereco17_94X"
#     real_paths = [stack.enter_context(importlib.resources.path(_ele_path, f)) for f in ele_sf_mapping.values()]
#     ext.add_weight_sets([f"{path} {file}" for path, file in zip(ele_sf_mapping.keys(), real_paths)])
#     ext.finalize()
# evaluator = ext.make_evaluator()
### Lepton SFs
def eleSFs(ele, campaign, path):
    _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
    ext = extractor()
    with contextlib.ExitStack() as stack:
        real_paths = [
            stack.enter_context(importlib.resources.path(_ele_path, f))
            for f in path.values()
        ]
        ext.add_weight_sets(
            [f"{path} {file}" for path, file in zip(path.keys(), real_paths)]
        )

    ext.finalize()
    evaluator = ext.make_evaluator()
    ele_eta = ak.fill_none(ele.eta, 0.0)
    ele_pt = ak.fill_none(ele.pt, 0.0)
    weight = 1.0
    for paths in path.keys():
        if "ele" in paths:
            if "above20" in paths:
                weight = weight * np.where(
                    ele_pt < 20.0,
                    1.0,
                    evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt),
                )
            elif "below20" in paths:
                weight = weight * np.where(
                    ele_pt > 20.0,
                    1.0,
                    evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt),
                )
            else:
                weight = weight * evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt)
    return weight


def add_eleSFs(ele, campaign, path, weights, cut):
    _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
    ext = extractor()
    with contextlib.ExitStack() as stack:
        real_paths = [
            stack.enter_context(importlib.resources.path(_ele_path, f))
            for f in path.values()
        ]
        ext.add_weight_sets(
            [f"{path} {file}" for path, file in zip(path.keys(), real_paths)]
        )

    ext.finalize()
    evaluator = ext.make_evaluator()
    ele_eta = ak.fill_none(ele.eta, 0.0)
    ele_pt = ak.fill_none(ele.pt, 0.0)
    weight = 1.0
    weight_err = 0.0
    for paths in path.keys():
        if "ele" in paths:
            if "above20" in paths:
                if "error" not in paths:
                    weight = weight * np.where(
                        ele_pt < 20.0,
                        1.0,
                        evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt),
                    )
                else:
                    weight_err = weight_err + np.power(
                        np.where(
                            ele_pt < 20.0,
                            0.0,
                            evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt),
                        ),
                        2,
                    )
            elif "below20" in paths:
                if "error" not in paths:
                    weight = weight * np.where(
                        ele_pt > 20.0,
                        1.0,
                        evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt),
                    )
                else:
                    weight_err = weight_err + np.power(
                        np.where(
                            ele_pt > 20.0,
                            0.0,
                            evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt),
                        ),
                        2,
                    )
            else:
                if "error" not in paths:
                    weight = weight * evaluator[paths[: paths.find(" ")]](
                        ele_eta, ele_pt
                    )
                else:
                    weight_err = weight_err + np.power(
                        evaluator[paths[: paths.find(" ")]](ele_eta, ele_pt), 2
                    )
    weight = np.where(cut & (ele.lep_flav == 11), weight, 1.0)
    weight_err = np.where(cut & (ele.lep_flav == 11), weight_err, 0.0)
    weights.add(
        "eleSFs", weight, weight + np.sqrt(weight_err), weight - np.sqrt(weight_err)
    )
    return weights


def muSFs(mu, campaign, path):
    _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
    ext = extractor()
    with contextlib.ExitStack() as stack:
        real_paths = [
            stack.enter_context(importlib.resources.path(_ele_path, f))
            for f in path.values()
        ]
        ext.add_weight_sets(
            [f"{path} {file}" for path, file in zip(path.keys(), real_paths)]
        )

    ext.finalize()
    evaluator = ext.make_evaluator()
    mu_eta = ak.fill_none(mu.eta, 0.0)
    mu_pt = ak.fill_none(mu.pt, 0.0)
    weight = 1.0
    for paths in path.keys():
        if "mu" in paths:
            weight = weight * evaluator[paths[: paths.find(" ")]](mu_eta, mu_pt)
    return weight


def add_muSFs(mu, campaign, path, weights, cut):
    _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
    ext = extractor()
    with contextlib.ExitStack() as stack:
        real_paths = [
            stack.enter_context(importlib.resources.path(_ele_path, f))
            for f in path.values()
        ]
        ext.add_weight_sets(
            [f"{path} {file}" for path, file in zip(path.keys(), real_paths)]
        )

    ext.finalize()
    evaluator = ext.make_evaluator()
    mu_eta = ak.fill_none(mu.eta, 0.0)
    mu_pt = ak.fill_none(mu.pt, 0.0)
    weight = 1.0
    weight_err = 0.0
    for paths in path.keys():
        if "mu" in paths:
            if "error" in paths:
                weight_err = weight_err + np.power(
                    evaluator[paths[: paths.find(" ")]](mu_eta, mu_pt), 2
                )
            else:
                weight = weight * evaluator[paths[: paths.find(" ")]](mu_eta, mu_pt)
    weight = np.where(cut & (mu.lep_flav == 13), weight, 1.0)
    weight_err = np.where(cut & (mu.lep_flav == 13), weight_err, 0.0)
    weights.add(
        "muSFs", weight, weight + np.sqrt(weight_err), weight - np.sqrt(weight_err)
    )
    return weights


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


# Jennet adds PS weights
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
        #else:
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
