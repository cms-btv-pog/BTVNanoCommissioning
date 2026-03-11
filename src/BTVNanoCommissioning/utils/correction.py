import importlib.resources
import contextlib
import copy
import os
import re
import warnings

import numpy as np
import awkward as ak
import uproot
import correctionlib

from coffea.lookup_tools import extractor, txt_converters, rochester_lookup
from coffea.lumi_tools import LumiMask
from coffea.jetmet_tools.CorrectedMETFactory import corrected_polar_met
from coffea.analysis_tools import Weights
from coffea.btag_tools import BTagScaleFactor

from BTVNanoCommissioning.helpers.MuonScaRe import (
    pt_resol,
    pt_scale,
    pt_resol_var,
    pt_scale_var,
)
from BTVNanoCommissioning.helpers.func import (
    _compile_jec_,
    _load_jmefactory,
    campaign_map,
)
from BTVNanoCommissioning.utils.AK4_parameters import correction_config as config


def load_SF(year, campaign, syst=False):
    """
    Load scale factors (SF) for a given year and campaign.

    This function reads scale factors from the specified campaign configuration and returns them in a suitable format.
    It handles different types of scale factors, such as pileup weights, and checks for the existence of files in
    the jsonpog-integration directory or custom files.

    Example:
    ```python
    ## Initialization, add EGM map from correctionlib
    correction_map["EGM"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/jsonpog-integration/POG/EGM/{campaign}/electron.json.gz"
                )
    ## Initialization, add EGM map from custom file by extractor
    ext = extractor()
    ext.add_weight_sets(["eleID EGamma2D {filename}.root"])
    ext.finalize()
    correction_map["EGM"] = ext.make_evaluator()
    ```

    Parameters:
    year (str): The year for which to load the scale factors.
    campaign (str): The name of the campaign for which to load the scale factors.
    syst (bool, optional): A flag to indicate whether to load systematic variations. Default is False.

    Returns:
    dict: A dictionary containing the scale factors, where keys are the relevant identifiers and values are the scale factors.

    Raises:
    FileNotFoundError: If the specified file does not exist.
    ValueError: If the file content is not in the expected format.
    KeyError: If the specified campaign or year is not found in the configuration.
    """
    # read the configuration file to get the correct SFs
    correct_map = {"campaign": campaign}

    for SF in config[campaign].keys():
        if SF == "DC":
            continue
        ## pileup weight
        if SF == "LUM":
            ## Check whether files in jsonpog-integration exist
            if os.path.exists(
                f"/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/{campaign_map()[campaign]}/latest/"
            ):
                try:
                    correct_map["LUM"] = correctionlib.CorrectionSet.from_file(
                        f"/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/{campaign_map()[campaign]}/latest/puWeights.json.gz"
                    )
                except FileNotFoundError:
                    correct_map["LUM"] = correctionlib.CorrectionSet.from_file(
                        f"/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/{campaign_map()[campaign]}/latest/puWeights_BCDEFGHI.json.gz"
                    )
            ## Otherwise custom files
            else:
                _pu_path = f"BTVNanoCommissioning.data.LUM.{campaign}"
                with importlib.resources.path(
                    _pu_path, config[campaign]["LUM"]
                ) as filename:
                    if str(filename).endswith(".json.gz"):
                        correct_map["LUM"] = correctionlib.CorrectionSet.from_file(
                            str(filename)
                        )
                    elif str(filename).endswith(".histo.root"):
                        ext = extractor()
                        ext.add_weight_sets([f"* * {filename}"])
                        ext.finalize()
                        correct_map["LUM"] = ext.make_evaluator()

        ## btag weight
        elif SF == "BTV":
            if "btag" in config[campaign]["BTV"].keys() and config[campaign]["BTV"][
                "btag"
            ].endswith(".json.gz"):
                correct_map["btag"] = correctionlib.CorrectionSet.from_file(
                    importlib.resources.path(
                        f"BTVNanoCommissioning.data.BTV.{campaign}", filename
                    )
                )
            if "ctag" in config[campaign]["BTV"].keys() and config[campaign]["BTV"][
                "ctag"
            ].endswith(".json.gz"):
                correct_map["btag"] = correctionlib.CorrectionSet.from_file(
                    importlib.resources.path(
                        f"BTVNanoCommissioning.data.BTV.{campaign}", filename
                    )
                )
            if os.path.exists(
                f"/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/{campaign_map()[campaign]}/latest/"
            ):
                correct_map["btag"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/{campaign_map()[campaign]}/latest/btagging.json.gz"
                )
                correct_map["ctag"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/{campaign_map()[campaign]}/latest/ctagging.json.gz"
                )
            else:
                correct_map["btag"] = {}
                correct_map["ctag"] = {}
                correct_map["btv_cfg"] = config[campaign]["BTV"]
                _btag_path = f"BTVNanoCommissioning.data.BTV.{campaign}"
                for tagger in config[campaign]["BTV"]:
                    with importlib.resources.path(
                        _btag_path, config[campaign]["BTV"][tagger]
                    ) as filename:
                        if "B" in tagger:
                            if filename.endswith(".json.gz"):
                                correct_map["btag"] = (
                                    correctionlib.CorrectionSet.from_file(filename)
                                )
                            else:
                                correct_map["btag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )
                        else:
                            if filename.endswith(".json.gz"):
                                correct_map["ctag"] = (
                                    correctionlib.CorrectionSet.from_file(filename)
                                )
                            else:
                                correct_map["ctag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )

        ## lepton SFs
        elif SF == "MUO" or SF == "EGM":
            correct_map["MUO_cfg"] = {
                mu: f
                for mu, f in config[campaign]["MUO"].items()
                if "mu" in mu and "_json" not in mu
            }
            correct_map["EGM_cfg"] = {
                e: f
                for e, f in config[campaign]["EGM"].items()
                if "ele" in e and "_json" not in e
            }
            ## muon
            _mu_path = f"/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/{campaign_map()[campaign]}/latest//muon_Z.json.gz"
            if not os.path.exists(_mu_path):
                _mu_path = f"src/BTVNanoCommissioning/data/MUO/{campaign_map()[campaign]}/latest//muon_Z.json.gz"
            if os.path.exists(_mu_path):
                correct_map["MUO"] = correctionlib.CorrectionSet.from_file(_mu_path)
            ## electron
            for _ele_file, _ele_map in {
                "electron": "EGM",
                "electronHlt": "EGM_HLT",
            }.items():
                _ele_path = f"/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/{campaign_map()[campaign]}/latest//{_ele_file}.json.gz"
                if not os.path.exists(_ele_path):
                    _ele_path = f"src/BTVNanoCommissioning/data/EGM/{campaign_map()[campaign]}/latest//{_ele_file}.json.gz"
                if os.path.exists(_ele_path):
                    correct_map[_ele_map] = correctionlib.CorrectionSet.from_file(
                        _ele_path
                    )
            ## json
            if any(
                np.char.find(np.array(list(config[campaign]["MUO"].keys())), "mu_json")
                != -1
            ):
                correct_map["MUO"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/data/MUO/{campaign_map()[campaign]}/latest/{config[campaign]['MUO']['mu_json']}"
                )
            if any(
                np.char.find(np.array(list(config[campaign]["EGM"].keys())), "ele_json")
                != -1
            ):
                correct_map["EGM"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/data/EGM/{campaign_map()[campaign]}/latest/{config[campaign]['EGM']['ele_json']}"
                )

            ## check if any custom corrections needed
            # FIXME: (some low pT muons not supported in jsonpog-integration at the moment)
            if (
                "histo.json" in "\t".join(list(config[campaign]["MUO"].values()))
                or "histo.txt" in "\t".join(list(config[campaign]["MUO"].values()))
                or "histo.root" in "\t".join(list(config[campaign]["MUO"].values()))
                or "histo.json" in "\t".join(list(config[campaign]["EGM"].values()))
                or "histo.txt" in "\t".join(list(config[campaign]["EGM"].values()))
                or "histo.root" in "\t".join(list(config[campaign]["EGM"].values()))
            ):
                _mu_path = f"BTVNanoCommissioning.data.MUO.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    inputs, real_paths = (
                        [
                            k
                            for k in correct_map["MUO_cfg"].keys()
                            if "histo.json" in correct_map["MUO_cfg"][k]
                            or "histo.txt" in correct_map["MUO_cfg"][k]
                            or "histo.root" in correct_map["MUO_cfg"][k]
                        ],
                        [
                            stack.enter_context(importlib.resources.path(_mu_path, f))
                            for f in correct_map["MUO_cfg"].values()
                            if ".json" in f or ".txt" in f or ".root" in f
                        ],
                    )

                    inputs = [
                        i.split(" ")[0] + " *" if "_low" in i else i for i in inputs
                    ]

                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(inputs, real_paths)
                            if "histo.json" in str(file)
                            or "histo.txt" in str(file)
                            or "histo.root" in str(file)
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
                correct_map["MUO_custom"] = ext.make_evaluator()

                _ele_path = f"BTVNanoCommissioning.data.EGM.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    inputs, real_paths = (
                        [
                            k
                            for k in correct_map["EGM_cfg"].keys()
                            if "histo.json" in correct_map["EGM_cfg"][k]
                            or "histo.txt" in correct_map["EGM_cfg"][k]
                            or "histo.root" in correct_map["EGM_cfg"][k]
                        ],
                        [
                            stack.enter_context(importlib.resources.path(_ele_path, f))
                            for f in correct_map["EGM_cfg"].values()
                            if "histo.json" in f or ".txt" in f or ".root" in f
                        ],
                    )
                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(inputs, real_paths)
                            if "histo.json" in str(file)
                            or "histo.txt" in str(file)
                            or "histo.root" in str(file)
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
                correct_map["EGM_custom"] = ext.make_evaluator()

        ## lepton scale & smearing
        elif SF == "muonSS":
            _mu_path = f"/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/{campaign_map()[campaign]}/latest/muon_scalesmearing.json.gz"
            if not os.path.exists(_mu_path):
                _mu_path = f"src/BTVNanoCommissioning/data/MUO/{campaign_map()[campaign]}/latest/muon_scalesmearing.json.gz"
            if os.path.exists(_mu_path):
                correct_map["muonSS"] = correctionlib.CorrectionSet.from_file(_mu_path)
        elif SF == "electronSS":
            _ele_path = f"/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/{campaign_map()[campaign]}/latest/electronSS_EtDependent{'_v1' if year == '2024' else ''}.json.gz"
            if not os.path.exists(_ele_path):
                _ele_path = f"src/BTVNanoCommissioning/data/EGM/{campaign_map()[campaign]}/latest/electronSS_EtDependent.json.gz"
            if os.path.exists(_ele_path):
                correct_map["electronSS"] = correctionlib.CorrectionSet.from_file(
                    _ele_path
                )
            correct_map["electronSS_cfg"] = config[campaign]["electronSS"]

        ## Rochester muon momentum correction (Run 2)
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
            correct_map["roccor"] = rochester_lookup.rochester_lookup(rochester_data)

        ## JME corrections
        elif SF == "JME":
            if "name" in config[campaign]["JME"].keys():
                if not os.path.exists(
                    f"src/BTVNanoCommissioning/data/JME/{campaign_map()[campaign]}/latest/jec_compiled_{config[campaign]['JME']['name']}.pkl.gz"
                ):
                    _compile_jec_(
                        year,
                        campaign,
                        config[campaign]["JME"],
                        f"jec_compiled_{config[campaign]['JME']['name']}",
                    )

                correct_map["JME"] = _load_jmefactory(
                    year,
                    campaign,
                    f"jec_compiled_{config[campaign]['JME']['name']}.pkl.gz",
                )
            elif os.path.exists(
                f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{campaign_map()[campaign]}/latest/jet_jerc.json.gz"
            ):
                correct_map["JME"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{campaign_map()[campaign]}/latest/jet_jerc.json.gz"
                )
                correct_map["JME_cfg"] = config[campaign]["JME"]
                for dataset in correct_map["JME_cfg"].keys():
                    if (
                        np.all(
                            np.char.find(
                                np.array(list(correct_map["JME"].keys())),
                                correct_map["JME_cfg"][dataset],
                            )
                        )
                        == -1
                    ):
                        raise (
                            f"{dataset} has no JEC map : {correct_map['JME_cfg'][dataset]} available"
                        )
        elif SF == "JMAR":
            if os.path.exists(
                f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{campaign_map()[campaign]}/latest/jmar.json.gz"
            ):
                correct_map["JMAR_cfg"] = {
                    j: f for j, f in config[campaign]["JMAR"].items()
                }
                correct_map["JMAR"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{campaign_map()[campaign]}/latest/jmar.json.gz"
                )

        elif SF == "jetveto":
            correct_map["jetveto_cfg"] = {
                j: f for j, f in config[campaign]["jetveto"].items()
            }

            isRootFile = False
            for val in correct_map["jetveto_cfg"].values():
                if ".root" in val:
                    isRootFile = True

            if not isRootFile and os.path.exists(
                f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{campaign_map()[campaign]}/latest/jetvetomaps.json.gz"
            ):
                correct_map["jetveto"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/{campaign_map()[campaign]}/latest/jetvetomaps.json.gz"
                )
            else:
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    ext.add_weight_sets(
                        [
                            f"{run} {stack.enter_context(importlib.resources.path(f'BTVNanoCommissioning.data.JME.{campaign}', file))}"
                            for run, file in config[campaign]["jetveto"].items()
                        ]
                    )
                    ext.finalize()
                    correct_map["jetveto"] = ext.make_evaluator()

    return correct_map


def load_lumi(campaign):
    """
    Load luminosity mask for a given campaign.

    This function reads the luminosity mask file for the specified campaign and returns a `LumiMask` object.

    Parameters:
    campaign (str): The name of the campaign for which to load the luminosity mask.

    Returns:
    LumiMask: An object representing the luminosity mask for the specified campaign.

    Raises:
    KeyError: If the specified campaign is not found in the configuration.
    FileNotFoundError: If the luminosity mask file does not exist.
    """

    _lumi_path = "BTVNanoCommissioning.data.DC"
    if os.path.exists(
        f'/cvmfs/cms-griddata.cern.ch/cat/metadata/DC/Collisions{campaign[-2:]}/config[campaign]["DC"]'
    ):
        return LumiMask(
            f"/cvmfs/cms-griddata.cern.ch/cat/metadata/DC/Collisions{campaign[-2:]}/{config[campaign]['DC']}"
        )
    else:
        with importlib.resources.path(_lumi_path, config[campaign]["DC"]) as filename:
            return LumiMask(filename)


## JEC
# FIXME: would be nicer if we can move to correctionlib in the future together with factory and workable


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


## Jet Veto
def jetveto(jets, correct_map):
    """
    Apply a veto to jets based on predefined transverse momentum (pt) and pseudorapidity (eta) thresholds.

    This function filters out jets that do not meet the predefined pt and eta criteria. It also utilizes a correction map
    to apply additional corrections or selections to the jets.

    Parameters:
    jets (iterable): A collection of jet objects or dictionaries containing jet properties.
    correct_map (dict): A dictionary containing correction factors or additional selection criteria for the jets.

    Returns:
    jets: A collection of jets that pass the predefined pt and eta criteria and any additional criteria from the correction map.

    Raises:
    TypeError: If the jets parameter is not an iterable.
    KeyError: If the jet objects do not contain the required 'pt' or 'eta' properties.
    """

    if "correctionlib" in str(
        type(correct_map["jetveto"][list(correct_map["jetveto"].keys())[0]])
    ):
        j, nj = ak.flatten(jets), ak.num(jets)
        return ak.unflatten(
            correct_map["jetveto"][list(correct_map["jetveto"].keys())[0]].evaluate(
                correct_map["jetveto_cfg"][list(correct_map["jetveto"].keys())[0]],
                np.clip(j.eta, -5.191, 5.191),
                np.clip(j.phi, -3.141592653589793, 3.141592653589793),
            ),
            nj,
        )
    else:
        return ak.where(
            correct_map["jetveto"][list(correct_map["jetveto"].keys())[0]](
                jets.phi, jets.eta
            )
            > 0,
            ak.broadcast_arrays(jets.eta, 100),
            ak.zeros_like(jets.eta),
        )


# from https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/jercExample.py
def get_corr_inputs(input_dict, corr_obj, jersyst="nom"):
    """
    Helper function for getting values of input variables
    given a dictionary and a correction object.
    """
    input_values = []
    for inputs in corr_obj.inputs:
        if "systematic" in inputs.name:
            input_values.append(jersyst)
        else:
            input_values.append(
                np.array(
                    input_dict[
                        inputs.name.replace("Jet", "")
                        .replace("Pt", "pt")
                        .replace("Phi", "phi")
                        .replace("Eta", "eta")
                        .replace("Mass", "mass")
                        .replace("Rho", "rho")
                        .replace("A", "area")
                    ]
                )
            )
    return input_values


cset_jersmear_paths = [
    "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/JER-Smearing/latest/jer_smear.json.gz",
    "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/jer_smear.json.gz",
]
for _jersmear_path in cset_jersmear_paths:
    if os.path.exists(_jersmear_path):
        cset_jersmear = correctionlib.CorrectionSet.from_file(_jersmear_path)
        break
else:
    warnings.warn(
        "JER smearing JSON not found in CVMFS. JER smearing corrections will be disabled.",
        RuntimeWarning,
        stacklevel=2,
    )
    cset_jersmear = {"JERSmear": None}
sf_jersmear = cset_jersmear["JERSmear"]


## JERC
def JME_shifts(
    shifts,
    correct_map,
    events,
    year,
    campaign,
    isRealData,
    systematic=False,
):
    """
    Apply Jet Energy Corrections (JEC) and Jet Energy Resolutions (JER) shifts to events.

    This function applies JEC and JER shifts to the jets in the events based on the provided correction map and campaign.
    It handles both real data and simulated data, and can optionally apply systematic variations and exclude jet vetoes.

    Parameters:
    shifts (list): A list of shift types to apply (e.g., 'up', 'down').
    correct_map (dict): A dictionary containing correction factors and settings for JEC and JER.
    events (awkward.Array): An array of events containing jet information.
    year (str): The year for which to apply the corrections.
    campaign (str): The name of the campaign for which to apply the corrections.
    isRealData (bool): A flag indicating whether the data is real or simulated.
    systematic (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

    Returns:
    awkward.Array: The events array with applied JEC and JER shifts.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the campaign is not recognized or supported.
    """
    dataset = events.metadata["dataset"]
    jecname = ""
    # https://cms-jerc.web.cern.ch/JECUncertaintySources/, currently no recommendation of reduced/full split sources
    syst_list = [
        i.split("_")[3]
        for i in correct_map["JME"].keys()
        if "MC" in i and "L1" not in i and "L2" not in i and "L3" not in i
    ]
    if "JME" in correct_map.keys():
        ## correctionlib
        if "JME_cfg" in correct_map.keys():
            if isRealData:
                jecname = [
                    v
                    for k, v in correct_map["JME_cfg"].items()
                    if k in events.metadata["dataset"]
                ]
                if len(jecname) > 1:
                    raise ValueError("Multiple uncertainties match to this era")
                elif len(jecname) == 0:
                    raise ValueError(
                        "Available JEC variations in this era are not compatible with this file. Did you choose the correct dataset-era combination?"
                    )
                else:
                    jecname = jecname[0] + "_DATA"
            else:
                jecname = correct_map["JME_cfg"]["MC"].split(" ")[0] + "_MC"
                jrname = correct_map["JME_cfg"]["MC"].split(" ")[1] + "_MC"

            # store the original jet info
            nocorrjet = events.Jet
            nocorrjet["pt_raw"] = (1 - nocorrjet["rawFactor"]) * nocorrjet["pt"]
            nocorrjet["mass_raw"] = (1 - nocorrjet["rawFactor"]) * nocorrjet["mass"]
            nocorrjet["rho"] = ak.broadcast_arrays(
                events.fixedGridRhoFastjetAll, nocorrjet.pt
            )[0]
            nocorrjet["EventID"] = ak.broadcast_arrays(events.event, nocorrjet.pt)[0]
            nocorrjet["run"] = ak.broadcast_arrays(events.run, nocorrjet.pt)[0]
            if not isRealData:
                genjetidx = ak.where(nocorrjet.genJetIdx == -1, 0, nocorrjet.genJetIdx)
                nocorrjet["Genpt"] = ak.where(
                    nocorrjet.genJetIdx == -1, -1, events.GenJet[genjetidx].pt
                )
            jets = copy.copy(nocorrjet)
            jets["orig_pt"] = ak.values_astype(nocorrjet["pt"], np.float32)

            ## flatten jets
            j, nj = ak.flatten(nocorrjet), ak.num(nocorrjet)

            # JEC
            JECcorr = correct_map["JME"].compound[f"{jecname}_L1L2L3Res_AK4PFPuppi"]
            JEC_input = get_corr_inputs(j, JECcorr)
            JECflatCorrFactor = JECcorr.evaluate(*JEC_input)

            ## JER
            if isRealData:
                # in data only the JEC is applied
                corrFactor = JECflatCorrFactor
            else:
                JERSF = correct_map["JME"][f"{jrname}_ScaleFactor_AK4PFPuppi"]
                JERptres = correct_map["JME"][f"{jrname}_PtResolution_AK4PFPuppi"]
                # for MC, correct the jet pT with JEC first
                j["pt"] = j["pt_raw"] * JECflatCorrFactor
                j["mass"] = j["mass_raw"] * JECflatCorrFactor
                JERSF_input = get_corr_inputs(j, JERSF)
                JERptres_input = get_corr_inputs(j, JERptres)
                j["JER"] = JERptres.evaluate(*JERptres_input)
                j["JERSF"] = JERSF.evaluate(*JERSF_input)
                if sf_jersmear is None:
                    warnings.warn(
                        "JER smearing coefficients unavailable. "
                        "Proceeding without applying JER smearing.",
                        RuntimeWarning,
                        stacklevel=2,
                    )
                    corrFactor = JECflatCorrFactor
                else:
                    JERsmear_input = get_corr_inputs(j, sf_jersmear)
                    corrFactor = JECflatCorrFactor * sf_jersmear.evaluate(
                        *JERsmear_input
                    )
            corrFactor = ak.unflatten(corrFactor, nj)

            jets["pt"] = ak.values_astype(nocorrjet["pt_raw"] * corrFactor, np.float32)
            jets["mass"] = ak.values_astype(
                nocorrjet["mass_raw"] * corrFactor, np.float32
            )

            # Type-I MET correction, from corrected MET factory
            # https://github.com/scikit-hep/coffea/blob/master/src/coffea/jetmet_tools/CorrectedMETFactory.py
            nocorrmet = events.PuppiMET if int(year) > 2020 else events.MET
            met = copy.copy(nocorrmet)
            metinfo = [nocorrmet.pt, nocorrmet.phi, jets.pt, jets.phi, jets.pt_raw]
            met["pt"], met["phi"] = (
                ak.values_astype(corrected_polar_met(*metinfo).pt, np.float32),
                ak.values_astype(corrected_polar_met(*metinfo).phi, np.float32),
            )
            met["orig_pt"], met["orig_phi"] = nocorrmet["pt"], nocorrmet["phi"]

            ## JEC variations
            if not isRealData and systematic != False:
                if systematic != "JERC_split":
                    jesuncmap = correct_map["JME"][f"{jecname}_Total_AK4PFPuppi"]
                    jesunc = ak.unflatten(jesuncmap.evaluate(j.eta, j.pt), nj)
                    unc_jets, unc_met = {}, {}

                    for var in ["up", "down"]:
                        fac = 1.0 if var == "up" else -1.0
                        ## JES total
                        unc_jets[f"JES_Total{var}"] = copy.copy(nocorrjet)
                        unc_met[f"JES_Total{var}"] = copy.copy(nocorrmet)

                        unc_jets[f"JES_Total{var}"]["pt"] = ak.values_astype(
                            jets["pt"] * (1 + fac * jesunc),
                            np.float32,
                        )
                        unc_jets[f"JES_Total{var}"]["mass"] = ak.values_astype(
                            jets["mass"] * (1 + fac * jesunc),
                            np.float32,
                        )
                        unc_met[f"JES_Total{var}"]["pt"] = corrected_polar_met(
                            nocorrmet.pt,
                            nocorrmet.phi,
                            unc_jets[f"JES_Total{var}"]["pt"],
                            jets.phi,
                            jets.pt_raw,
                        ).pt
                        unc_met[f"JES_Total{var}"]["phi"] = corrected_polar_met(
                            nocorrmet.pt,
                            nocorrmet.phi,
                            unc_jets[f"JES_Total{var}"]["pt"],
                            jets.phi,
                            jets.pt_raw,
                        ).phi

                        JERSF_input_var = get_corr_inputs(j, JERSF, var)

                        ## JER variations
                        if sf_jersmear is None:
                            warnings.warn(
                                "Skipping JER smearing variations because the "
                                "correction file is unavailable.",
                                RuntimeWarning,
                                stacklevel=2,
                            )
                            unc_jets[f"JER{var}"] = copy.copy(jets)
                            unc_met[f"JER{var}"] = copy.copy(met)
                        else:
                            unc_jets[f"JER{var}"] = copy.copy(nocorrjet)
                            unc_met[f"JER{var}"] = copy.copy(nocorrmet)
                            j["JERSF"] = JERSF.evaluate(*JERSF_input_var)
                            JERsmear_input_var = get_corr_inputs(j, sf_jersmear)
                            corrFactor_var = JECflatCorrFactor * sf_jersmear.evaluate(
                                *JERsmear_input_var
                            )
                            corrFactor_var = ak.unflatten(corrFactor_var, nj)

                            unc_jets[f"JER{var}"]["pt"] = ak.values_astype(
                                nocorrjet["pt_raw"] * corrFactor_var,
                                np.float32,
                            )
                            unc_jets[f"JER{var}"]["mass"] = ak.values_astype(
                                nocorrjet["mass_raw"] * corrFactor_var,
                                np.float32,
                            )
                            unc_met[f"JER{var}"]["pt"] = corrected_polar_met(
                                nocorrmet.pt,
                                nocorrmet.phi,
                                unc_jets[f"JER{var}"]["pt"],
                                jets.phi,
                                jets.pt_raw,
                            ).pt
                            unc_met[f"JER{var}"]["phi"] = corrected_polar_met(
                                nocorrmet.pt,
                                nocorrmet.phi,
                                unc_jets[f"JER{var}"]["pt"],
                                jets.phi,
                                jets.pt_raw,
                            ).phi
                    jets["JES_Total"] = ak.zip(
                        {
                            "up": unc_jets["JES_Totalup"],
                            "down": unc_jets["JES_Totaldown"],
                        }
                    )
                    jets["JER"] = ak.zip(
                        {
                            "up": unc_jets["JERup"],
                            "down": unc_jets["JERdown"],
                        }
                    )
                    met["JES_Total"] = ak.zip(
                        {
                            "up": unc_met["JES_Totalup"],
                            "down": unc_met["JES_Totaldown"],
                        }
                    )
                    met["JER"] = ak.zip(
                        {
                            "up": unc_met["JERup"],
                            "down": unc_met["JERdown"],
                        }
                    )
                else:
                    raise NotImplementedError

        else:
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
                elif campaign == "Summer23":
                    if "v4" in dataset:
                        jecname = "Cv4"
                    else:
                        jecname = "Cv123"
                elif re.search(r"[Rr]un20\d{2}([A-Z])", dataset):
                    jecname = re.search(r"[Rr]un20\d{2}([A-Z])", dataset).group(1)
                else:
                    print("No valid jec name")
                    raise NameError
                jecname = "data" + jecname
            else:
                jecname = "MC"

            jets = correct_map["JME"]["jet_factory"][jecname].build(
                add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
                lazy_cache=events.caches[0],
            )
            met = correct_map["JME"]["met_factory"].build(events.PuppiMET, jets, {})

        # systematics
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
                    if "JES_Total" in jets.fields:
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
                    if "MET_UnclusteredEnergy" in met.fields:
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
                        ]
                    if "JER" in jets.fields:
                        shifts += [
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

    else:
        met = events.PuppiMET
        jets = events.Jet

    shifts.insert(0, ({"Jet": jets, "MET": met}, None))
    return shifts


## Muon Rochester correction (Run 2)
def Roccor_shifts(shifts, correct_map, events, isRealData, systematic=False):
    """
    Apply Rochester corrections (Roccor) shifts to muons in events.

    This function applies Rochester corrections to the muons in the events based on the provided correction map and campaign.
    It handles both real data and simulated data, and can optionally apply systematic variations.

    Parameters:
    shifts (list): A list of shift types to apply (e.g., 'up', 'down').
    correct_map (dict): A dictionary containing correction factors and settings for Rochester corrections.
    events (awkward.Array): An array of events containing muon information.
    campaign (str): The name of the campaign for which to apply the corrections.
    isRealData (bool): A flag indicating whether the data is real or simulated.
    systematic (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

    Returns:
    awkward.Array: The events array with applied Rochester corrections.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the campaign is not recognized or supported.
    """

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


def MUO_shifts(shifts, correct_map, events, isRealData, systematic=False):
    """
    Applies the Run 3 recommended muon scale and smearing corrections.
    Returns the corrected muon objects, including systematics.
    Adapted from this example of muon SS correction usage:
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/muoScaleAndSmearingCoffeaExample.py
    """

    mu = events.Muon

    if isRealData:
        mu_pt_corr = pt_scale(
            1,  # 1 for data, 0 for mc
            mu.pt,
            mu.eta,
            mu.phi,
            mu.charge,
            correct_map["muonSS"],
            nested=True,  # for awkward arrays, set False for 1d arrays
        )

    else:
        mu_pt_scalecorr = pt_scale(
            0, mu.pt, mu.eta, mu.phi, mu.charge, correct_map["muonSS"], nested=True
        )
        mu_pt_corr = pt_resol(
            mu_pt_scalecorr,
            mu.eta,
            mu.phi,
            mu.nTrackerLayers,
            events.event,
            events.luminosityBlock,
            correct_map["muonSS"],
            nested=True,
        )

    # scale and smearing uncertainties should be evaluated and applied on MC only
    if systematic and not isRealData:
        mu_pt_corr_scaleup = pt_scale_var(
            mu_pt_corr,
            mu.eta,
            mu.phi,
            mu.charge,
            "up",
            correct_map["muonSS"],
            nested=True,
        )
        mu_pt_corr_scaledown = pt_scale_var(
            mu_pt_corr,
            mu.eta,
            mu.phi,
            mu.charge,
            "dn",
            correct_map["muonSS"],
            nested=True,
        )
        mu_pt_corr_resolup = pt_resol_var(
            mu_pt_scalecorr,
            mu_pt_corr,
            mu.eta,
            "up",
            correct_map["muonSS"],
            nested=True,
        )
        mu_pt_corr_resoldown = pt_resol_var(
            mu_pt_scalecorr,
            mu_pt_corr,
            mu.eta,
            "dn",
            correct_map["muonSS"],
            nested=True,
        )

    mu["pt"] = mu_pt_corr
    # add nominal scale & smearing correction to shifts
    for i in range(len(shifts)):
        shifts[i][0]["Muon"] = mu

    if systematic:
        mu_up, mu_down = events.Muon, events.Muon
        mu_resol_up, mu_resol_down = events.Muon, events.Muon

        if not isRealData:
            mu_up["pt"] = mu_pt_corr_scaleup
            mu_down["pt"] = mu_pt_corr_scaledown
            mu_resol_up["pt"] = mu_pt_corr_resolup
            mu_resol_down["pt"] = mu_pt_corr_resoldown

        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": mu_up,
                },
                "MuonScaleUp",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": mu_down,
                },
                "MuonScaleDown",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": mu_resol_up,
                },
                "MuonResolUp",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": mu_resol_down,
                },
                "MuonResolDown",
            )
        ]

    return shifts


def EGM_shifts(shifts, correct_map, events, isRealData, systematic=False):
    """
    Applies the Run 3 recommended electron scale and smearing corrections.
    Returns the corrected electron objects, including systematics.
    Adapted from this example of electron SS correction usage:
    https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/egmScaleAndSmearingExample.py
    """

    ele = events.Electron
    n_ele = ak.num(ele)
    events_run = ak.flatten(ak.broadcast_arrays(events.run, ele.eta)[0])
    ele_etaSC = (
        ak.flatten(ele.superclusterEta)
        if "Summer24" in correct_map["campaign"]
        else ak.flatten(ele.eta + ele.deltaEtaSC)
    )
    ele_r9 = ak.flatten(ele.r9)
    ele_pt = ak.flatten(ele.pt)
    ele_seedGain = ak.flatten(ele.seedGain)

    if isRealData:  # scale correction is only applied to data
        scale_evaluator = correct_map["electronSS"].compound[
            correct_map["electronSS_cfg"][0]
        ]
        if "Summer24" in correct_map["campaign"]:
            scale = scale_evaluator.evaluate(
                "scale", events_run, ele_etaSC, ele_r9, ele_pt, ele_seedGain
            )
        else:
            scale = scale_evaluator.evaluate(
                "scale",
                events_run,
                ele_etaSC,
                ele_r9,
                ele_pt,
                ele_seedGain,
            )
        scale = ak.unflatten(scale, n_ele)
        ele_pt_corr = scale * ele.pt
    else:  # smear correction is only applied to MC
        smear_and_syst_evaluator = correct_map["electronSS"][
            correct_map["electronSS_cfg"][1]
        ]
        smear = smear_and_syst_evaluator.evaluate(
            "smear", ele_pt, ele_r9, np.abs(ele_etaSC)
        )
        smear = ak.unflatten(smear, n_ele)
        # since the smearing is stochastic, a random number is needed for each event
        rng = np.random.default_rng(seed=125)
        random_numbers = rng.normal(loc=0.0, scale=1.0, size=len(ele.pt))
        ele_pt_corr = ele.pt * (1 + smear * random_numbers)

    # scale and smearing uncertainties should be evaluated on the original MC only
    if systematic and not isRealData:
        unc_scale = smear_and_syst_evaluator.evaluate(
            "escale", ele_pt, ele_r9, np.abs(ele_etaSC)
        )
        unc_scale = ak.unflatten(unc_scale, n_ele)
        unc_smear = smear_and_syst_evaluator.evaluate(
            "esmear", ele_pt, ele_r9, np.abs(ele_etaSC)
        )
        unc_smear = ak.unflatten(unc_smear, n_ele)

    ele["pt"] = ele_pt_corr
    # add nominal scale & smearing correction to shifts
    for i in range(len(shifts)):
        shifts[i][0]["Electron"] = ele

    if systematic:
        ele_scale_up, ele_scale_down = events.Electron, events.Electron
        ele_smear_up, ele_smear_down = events.Electron, events.Electron

        if not isRealData:
            ele_scale_up["pt"] = (1 + unc_scale) * ele_pt_corr
            ele_scale_down["pt"] = (1 - unc_scale) * ele_pt_corr
            ele_smear_up["pt"] = events.Electron.pt * (
                1 + (smear + unc_smear) * random_numbers
            )
            ele_smear_down["pt"] = events.Electron.pt * (
                1 + np.maximum(0.0, (smear - unc_smear)) * random_numbers
            )

        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": shifts[0][0]["Muon"],
                    "Electron": ele_scale_up,
                },
                "ElectronScaleUp",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": shifts[0][0]["Muon"],
                    "Electron": ele_scale_down,
                },
                "ElectronScaleDown",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": shifts[0][0]["Muon"],
                    "Electron": ele_smear_up,
                },
                "ElectronSmearUp",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": shifts[0][0]["Muon"],
                    "Electron": ele_smear_down,
                },
                "ElectronSmearDown",
            )
        ]

    return shifts


# wrapped up common shifts


def puwei(nPU, correct_map, weights, syst=False):
    """
    Return pileup weight
    Parameters
    ----------
    nPU: ak.Array
    correct_map : dict
    weights : coffea.analysis_tool.weights
    syst: "split", "weight_only"
    Apply pileup weights to events based on the number of primary vertices (nPU).

    This function applies pileup weights to the events using the provided correction map and weights.
    It can optionally apply systematic variations.

    Parameters:
    nPU (awkward.Array(int)): The number of primary vertices in the event.
    correct_map (dict): A dictionary containing correction factors and settings for pileup weights.
    weights (): A dictionary to store the calculated weights.
    syst (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

    Returns:
    None: The function modifies the weights dictionary in place.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the nPU value is not recognized or supported.
    """
    if "correctionlib" in str(type(correct_map["LUM"])):
        if syst:
            return weights.add(
                "puweight",
                correct_map["LUM"][list(correct_map["LUM"].keys())[0]].evaluate(
                    nPU, "nominal"
                ),
                correct_map["LUM"][list(correct_map["LUM"].keys())[0]].evaluate(
                    nPU, "up"
                ),
                correct_map["LUM"][list(correct_map["LUM"].keys())[0]].evaluate(
                    nPU, "down"
                ),
            )
        else:
            return weights.add(
                "puweight",
                correct_map["LUM"][list(correct_map["LUM"].keys())[0]].evaluate(
                    nPU, "nominal"
                ),
            )
    else:
        # Legacy ROOT histos may be keyed as PU/PUup/PUdown instead of LUM/PUup/PUdown
        if "LUM" in correct_map["LUM"]:
            central_key = "LUM"
        elif "PU" in correct_map["LUM"]:
            central_key = "PU"
        else:
            raise KeyError("Pileup central value not found in correct_map['LUM']")
        if "PUup" in correct_map["LUM"]:
            up_key = "PUup"
        elif "LUMup" in correct_map["LUM"]:
            up_key = "LUMup"
        else:
            raise KeyError("Pileup up-variation not found in correct_map['LUM']")

        if "PUdown" in correct_map["LUM"]:
            down_key = "PUdown"
        elif "LUMdown" in correct_map["LUM"]:
            down_key = "LUMdown"
        else:
            raise KeyError("Pileup down-variation not found in correct_map['LUM']")
        if syst:
            weights.add(
                "puweight",
                correct_map["LUM"][central_key](nPU),
                correct_map["LUM"][up_key](nPU),
                correct_map["LUM"][down_key](nPU),
            )
        else:
            weights.add("puweight", correct_map["LUM"][central_key](nPU))


def btagSFs(jet, correct_map, weights, SFtype, syst=False):
    """
    Apply b-tagging scale factors (SFs) to a single jet.

    This function applies b-tagging scale factors to the given jet based on the provided correction map, weights, and scale factor type.
    It can optionally apply systematic variations.

    Parameters:
    jet (dict): A dictionary containing the properties of the jet.
    correct_map (dict): A dictionary containing correction factors and settings for b-tagging scale factors.x
    weights (coffea.weight.Weight): An instance of coffea's Weight class to store the calculated weights.
    SFtype (str): The type of scale factor to apply , only shape-based C, B are supported.
    syst (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

    Returns:
    None: The function modifies the weights instance in place.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the SFtype is not recognized or supported.
    """
    if SFtype.endswith("C"):
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
    elif SFtype.endswith("B"):
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


def eleSFs(ele, correct_map, weights, syst=True, isHLT=False):
    allele = ele if ele.ndim > 1 else ak.singletons(ele)

    for sf in correct_map["EGM_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        sf_type = sf[: sf.find(" ")]
        for nele in range(ak.num(allele.pt)[0]):
            ele = allele[:, nele]
            ele_etaSC = (
                ak.fill_none(ele.eta + ele.deltaEtaSC, -2.5)
                if "Summer24" not in correct_map["campaign"]
                else ak.fill_none(ele.superclusterEta, -2.5)
            )
            ele_pt = ak.fill_none(ele.pt, 20)
            ele_pt = np.clip(ele_pt, 20, 999)
            masknone = ak.is_none(ele.pt)
            sfs_alle, sfs_alle_up, sfs_alle_down = (
                np.ones_like(allele[:, 0].pt),
                np.ones_like(allele[:, 0].pt),
                np.ones_like(allele[:, 0].pt),
            )

            if "correctionlib" in str(type(correct_map["EGM"])):
                ## reco SFs, split by pT
                if "Reco" in sf:
                    ## phi is used in Summer23
                    ele_pt = np.clip(ele.pt, 20.1, 74.9)
                    ele_pt_low = np.where(ele.pt >= 20.0, 19.9, ele.pt)
                    ele_pt_high = np.clip(ele.pt, 75.0, 500.0)
                    if "Summer23" in correct_map["campaign"]:
                        sfs_low = np.where(
                            (ele.pt <= 20.0) & ~masknone,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoBelow20",
                                ele_etaSC,
                                ele_pt_low,
                                ele.phi,
                            ),
                            1.0,
                        )
                        sfs_high = np.where(
                            (ele.pt > 75.0) & ~masknone,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoAbove75",
                                ele_etaSC,
                                ele_pt_high,
                                ele.phi,
                            ),
                            sfs_low,
                        )
                        sfs = np.where(
                            (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "Reco20to75",
                                ele_etaSC,
                                ele_pt,
                                ele.phi,
                            ),
                            sfs_high,
                        )
                        sfs = np.where(masknone, 1.0, sfs)

                        if syst != False:
                            sfs_up_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoBelow20",
                                    ele_etaSC,
                                    ele_pt_low,
                                    ele.phi,
                                ),
                                0.0,
                            )
                            sfs_down_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoBelow20",
                                    ele_etaSC,
                                    ele_pt_low,
                                    ele.phi,
                                ),
                                0.0,
                            )
                            sfs_up_high = np.where(
                                (ele.pt > 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoAbove75",
                                    ele_etaSC,
                                    ele_pt_high,
                                    ele.phi,
                                ),
                                sfs_up_low,
                            )
                            sfs_down_high = np.where(
                                (ele.pt > 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoAbove75",
                                    ele_etaSC,
                                    ele_pt_high,
                                    ele.phi,
                                ),
                                sfs_down_low,
                            )
                            sfs_up = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "Reco20to75",
                                    ele_etaSC,
                                    ele_pt,
                                    ele.phi,
                                ),
                                sfs_up_high,
                            )
                            sfs_down = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "Reco20to75",
                                    ele_etaSC,
                                    ele_pt,
                                    ele.phi,
                                ),
                                sfs_down_high,
                            )
                            sfs_up = np.where(masknone, 1.0, sfs_up)
                            sfs_down = np.where(masknone, 1.0, sfs_down)

                    else:
                        sfs_low = np.where(
                            (ele.pt <= 20.0) & ~masknone,
                            (
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sf",
                                    "RecoBelow20",
                                    ele_etaSC,
                                    ele_pt_low,
                                )
                                if "Summer24" not in correct_map["campaign"]
                                else 1.0
                            ),  # TODO: temporary until RecoBelow20 is released for 2024
                            1.0,
                        )
                        sfs_high = np.where(
                            (ele.pt > 75.0) & ~masknone,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoAbove75",
                                ele_etaSC,
                                ele_pt_high,
                            ),
                            sfs_low,
                        )
                        sfs = np.where(
                            (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1], "sf", "Reco20to75", ele_etaSC, ele_pt
                            ),
                            sfs_high,
                        )
                        sfs = np.where(masknone, 1.0, sfs)

                        if syst:
                            sfs_up_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                (
                                    correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                        sf.split(" ")[1],
                                        "sfup",
                                        "RecoBelow20",
                                        ele_etaSC,
                                        ele_pt_low,
                                    )
                                    if "Summer24" not in correct_map["campaign"]
                                    else 0.0
                                ),  # TODO: temporary until RecoBelow20 is released for 2024
                                0.0,
                            )
                            sfs_down_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                (
                                    correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                        sf.split(" ")[1],
                                        "sfdown",
                                        "RecoBelow20",
                                        ele_etaSC,
                                        ele_pt_low,
                                    )
                                    if "Summer24" not in correct_map["campaign"]
                                    else 0.0
                                ),  # TODO: temporary until RecoBelow20 is released for 2024
                                0.0,
                            )
                            sfs_up_high = np.where(
                                (ele.pt > 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoAbove75",
                                    ele_etaSC,
                                    ele_pt_high,
                                ),
                                sfs_up_low,
                            )
                            sfs_down_high = np.where(
                                (ele.pt > 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoAbove75",
                                    ele_etaSC,
                                    ele_pt_high,
                                ),
                                sfs_down_low,
                            )
                            sfs_up = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "Reco20to75",
                                    ele_etaSC,
                                    ele_pt,
                                ),
                                sfs_up_high,
                            )
                            sfs_down = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "Reco20to75",
                                    ele_etaSC,
                                    ele_pt,
                                ),
                                sfs_down_high,
                            )
                            sfs_up = np.where(masknone, 1.0, sfs_up)
                            sfs_down = np.where(masknone, 1.0, sfs_down)

                else:
                    # trigger SFs
                    if "Trig" in sf and "correctionlib" in str(
                        type(correct_map["EGM_HLT"])
                    ):
                        _ele_map = "EGM_HLT"
                    # ID SFs
                    else:
                        _ele_map = "EGM"

                    if "Summer23" in correct_map["campaign"]:
                        sfs = np.where(
                            masknone | (ele.pt > 1000.0),
                            1.0,
                            correct_map[_ele_map][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                correct_map["EGM_cfg"][sf],
                                ele_etaSC,
                                ele_pt,
                                ele.phi,
                            ),
                        )

                        if syst:
                            sfs_up = np.where(
                                masknone | (ele.pt > 1000.0),
                                1.0,
                                correct_map[_ele_map][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    correct_map["EGM_cfg"][sf],
                                    ele_etaSC,
                                    ele_pt,
                                    ele.phi,
                                ),
                            )
                            sfs_down = np.where(
                                masknone | (ele.pt > 1000.0),
                                1.0,
                                correct_map[_ele_map][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    correct_map["EGM_cfg"][sf],
                                    ele_etaSC,
                                    ele_pt,
                                    ele.phi,
                                ),
                            )
                    else:
                        sfs = np.where(
                            masknone | (ele.pt > 1000.0),
                            1.0,
                            correct_map[_ele_map][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                correct_map["EGM_cfg"][sf],
                                ele_etaSC,
                                ele_pt,
                            ),
                        )

                        if syst:
                            sfs_up = np.where(
                                masknone | (ele.pt > 1000.0),
                                1.0,
                                correct_map[_ele_map][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    correct_map["EGM_cfg"][sf],
                                    ele_etaSC,
                                    ele_pt,
                                ),
                            )
                            sfs_down = np.where(
                                masknone | (ele.pt > 1000.0),
                                1.0,
                                correct_map[_ele_map][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    correct_map["EGM_cfg"][sf],
                                    ele_etaSC,
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
                        correct_map["EGM_custom"][sf_type](ele_etaSC, ele_pt),
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_etaSC, ele_pt)
                            + correct_map["EGM_custom"][f"{sf_type}_error"](
                                ele_etaSC, ele_pt
                            ),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_etaSC, ele_pt)
                            - correct_map["EGM_custom"][f"{sf_type}_error"](
                                ele_etaSC, ele_pt
                            ),
                        )

            sfs_alle = sfs_alle * sfs
            if syst:
                sfs_alle_down = sfs_alle_down * sfs_down
                sfs_alle_up = sfs_alle_up * sfs_up

        sfname = sf.split(" ")[0]
        if syst:
            weights.add(sfname, sfs_alle, sfs_alle_up, sfs_alle_down)
        else:
            weights.add(sfname, sfs_alle)

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
            mu_pt = np.clip(mu.pt, 15.0, 199.9)
            mu_eta = np.clip(np.abs(mu.eta), 0.0, 2.4)
            sfs = 1.0
            if "correctionlib" in str(type(correct_map["MUO"])):
                sfs = np.where(
                    masknone,
                    1.0,
                    correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                        mu_eta, mu_pt, "nominal"
                    ),
                )

                if syst:
                    sf_unc = np.where(
                        masknone,
                        0.0,
                        correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                            mu_eta, mu_pt, "syst"
                        ),
                    )
                    sfs_up, sfs_down = 1.0 + sf_unc, 1.0 - sf_unc
            else:
                if "mu" in sf:
                    sfs = np.where(
                        masknone, 1.0, correct_map["MUO_cfg"][sf_type](mu_eta, mu_pt)
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
                        sfs_down = np.where(
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


def add_pdf_weight(weights, pdf_weights, isSyst=False):
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

        # alpha_S weights
        # Eq. 27 of same ref
        as_unc = 0.5 * (pdf_weights[:, 102] - pdf_weights[:, 101])

        # PDF + alpha_S weights
        # Eq. 28 of same ref
        pdfas_unc = np.sqrt(np.square(pdf_unc) + np.square(as_unc))
        if isSyst != False:
            weights.add("PDF_weight", nom, pdf_unc + nom)
            weights.add("aS_weight", nom, as_unc + nom)
            weights.add("PDFaS_weight", nom, pdfas_unc + nom)

        else:
            weights.add("PDF_weight", nom)
            weights.add("aS_weight", nom)
            weights.add("PDFaS_weight", nom)
    else:
        warnings.warn("PDF weights are not available")
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
def add_ps_weight(weights, ps_weights, isSyst=False):
    nom = np.ones(len(weights.weight()))
    up_isr = np.ones(len(weights.weight()))
    down_isr = np.ones(len(weights.weight()))
    up_fsr = np.ones(len(weights.weight()))
    down_fsr = np.ones(len(weights.weight()))

    if ps_weights is not None and isSyst != False:
        if len(ps_weights[0]) == 4:
            up_isr = ps_weights[:, 0]
            down_isr = ps_weights[:, 2]
            up_fsr = ps_weights[:, 1]
            down_fsr = ps_weights[:, 3]
            weights.add("UEPS_ISR", nom, up_isr, down_isr)
            weights.add("UEPS_FSR", nom, up_fsr, down_fsr)

        else:
            warnings.warn(f"PS weight vector has length {len(ps_weights[0])}")
            weights.add("UEPS_FSR", nom, nom, nom)


def add_scalevar_weight(weights, lhe_weights, isSyst=False):
    """
    Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Factorization_and_renormalizatio

    __doc__:
    ['LHE scale variation weights (w_var / w_nominal)',
    ' [0] is renscfact=0.5d0 facscfact=0.5d0 ',
    ' [1] is renscfact=0.5d0 facscfact=1d0 ',
    ' [2] is renscfact=0.5d0 facscfact=2d0 ',
    ' [3] is renscfact=1d0 facscfact=0.5d0 ',
    ' [4] is renscfact=1d0 facscfact=1d0 ',
    ' [5] is renscfact=1d0 facscfact=2d0 ',
    ' [6] is renscfact=2d0 facscfact=0.5d0 ',
    ' [7] is renscfact=2d0 facscfact=1d0 ',
    ' [8] is renscfact=2d0 facscfact=2d0 ']
    """

    nom = np.ones(len(weights.weight()))
    if isSyst != False:
        if len(lhe_weights) > 0:
            if len(lhe_weights[0]) == 9:
                nom = lhe_weights[:, 4]
                weights.add(
                    "scalevar_muR",
                    nom,
                    lhe_weights[:, 1] / nom,
                    lhe_weights[:, 7] / nom,
                )
                weights.add(
                    "scalevar_muF",
                    nom,
                    lhe_weights[:, 3] / nom,
                    lhe_weights[:, 5] / nom,
                )
                weights.add(
                    "scalevar_muR_muF", nom, lhe_weights[:, 0], lhe_weights[:, 8]
                )
            elif len(lhe_weights[0]) > 1:
                print("Scale variation vector has length ", len(lhe_weights[0]))
        else:
            warnings.warn(
                "LHE scale variation weights are not available, put nominal weights"
            )
            weights.add("scalevar_muR", nom, nom, nom)
            weights.add("scalevar_muF", nom, nom, nom)
            weights.add("scalevar_muR_muF", nom, nom, nom)

    else:
        weights.add("scalevar_3pt", nom)


# JP calibration utility
class JPCalibHandler(object):
    def __init__(self, year, campaign, isRealData, dataset, isSyst=False):
        """
        A tool for calculating the track probability and jet probability
            campaign: campaign name
            isRealData: whether the dataset is real data
            dataset: dataset name from events.metadata["dataset"]
        """
        if "JPCalib" not in config[campaign].keys():
            templates = uproot.open(
                "src/BTVNanoCommissioning/data/JPCalib/Summer22Run3/calibeHistoWrite_MC2022_NANO130X_v2.root"
            )
        else:
            if isRealData:
                if isSyst is not False:
                    filename = config[campaign]["JPCalib"]["MC"]
                else:
                    filename = "default"
                    for key in config[campaign]["JPCalib"]:
                        if key in dataset:
                            filename = config[campaign]["JPCalib"][key]
                            break
                    if filename == "default":
                        raise ValueError(f"No JPCalib file found for dataset {dataset}")
            else:
                filename = config[campaign]["JPCalib"]["MC"]

            templates = uproot.open(
                f"src/BTVNanoCommissioning/data/JPCalib/{campaign}/{filename}"
            )
        self.ipsig_histo_val = np.array(
            [templates[f"histoCat{i}"].values() for i in range(10)]
        )
        self.ipsig_histo_tot = np.sum(self.ipsig_histo_val, axis=1)
        self.values_cumsum = np.cumsum(self.ipsig_histo_val[:, ::-1], axis=1)[:, ::-1]
        self.edges = templates["histoCat0"].axes[0].edges()

    def flatten(self, array):
        """
        Get the fully flattened array and its layout for each layer
        """
        layouts = []
        array_fl = array
        while str(ak.type(array_fl)).count("*") > 1:
            layouts.append(ak.num(array_fl))
            array_fl = ak.flatten(array_fl)
        return array_fl, layouts

    def unflatten(self, array_fl, layouts):
        """
        Recover a flattened array using the original layouts
        """

        array = array_fl
        for layout in layouts[::-1]:
            array = ak.unflatten(array, layout)
        return array

    def calc_track_proba(self, ipsig: ak.Array, cat: ak.Array):
        """
        Calculate the track probability from the integral of the track IPsig templates, given the IPsig and category.
        Reference code: https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/RecoBTag/TrackProbability/src/HistogramProbabilityEstimator.cc
            ipsig: IP significance array
            cat: category array (0-9)
        """

        if ak.any(cat < 0) or ak.any(cat > 9):
            raise ValueError("Category out of range [0, 9]")

        # get the fully flattened array of the input while storing its layouts for later recovery
        ipsig_fl, layouts = self.flatten(ipsig)
        cat_fl = ak.flatten(cat, axis=None)

        # index of the IPsig bins
        ipsig_fl = abs(ipsig_fl)
        ipsig_fl_index = np.minimum(
            np.searchsorted(self.edges, ipsig_fl), self.ipsig_histo_val.shape[1] - 1
        )

        # retrieve the cumsum value (\int_{ipsig}^{inf} p(ipsig') d(ipsig')) from the correct template
        ipsig_cumsum_fl = self.values_cumsum[cat_fl, ipsig_fl_index]

        # calculate the track probability as (\int_{ipsig}^{inf} ..) / (\int_{0}^{inf} ..) * sign(IPsig)
        proba_fl = (ipsig_cumsum_fl / self.ipsig_histo_tot[cat_fl]) * np.sign(ipsig_fl)

        # recover the original layout
        proba = self.unflatten(proba_fl, layouts)
        return proba

    def calc_jet_proba(self, proba):
        # Calculate jet probability (JP)
        # according to jetProbability func in https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/RecoBTag/ImpactParameter/interface/TemplatedJetProbabilityComputer.h

        # minium proba = 0.5%
        proba = np.maximum(proba, 0.005)  # dim: (evt, jet, trk)

        ntrk = ak.num(proba, axis=-1)  # dim: (evt, jet), the number of tracks in a jet
        prodproba_log = ak.sum(
            np.log(proba), axis=-1
        )  # dim: (evt, jet), the log(Π(proba)) of all tracks in a jet
        prodproba_log_m_log = ak.where(
            (ak.num(proba, axis=-1) >= 2) & (prodproba_log < 0),
            np.log(-prodproba_log),
            0,
        )  # log(-logΠ), if >=2 tracks in a jet

        # now calculating Σ_tr{0..N-1} ((-logΠ)^tr / tr!)
        trk_index = ak.local_index(proba)
        fact_array = ak.concatenate(
            [
                [1.0],
                np.arange(1, max(5, ak.max(trk_index) + 1), dtype=np.float64).cumprod(),
            ]
        )  # construct a factorial array
        trk_index_fl, _layouts = self.flatten(trk_index)
        lfact = self.unflatten(
            fact_array[trk_index_fl], _layouts
        )  # dim: (evt, jet, trk), nested factorial array given the track index

        prob = ak.sum(
            np.exp(trk_index * prodproba_log_m_log - np.log(lfact)), axis=-1
        )  # dim: (evt, jet), Σ_tr{0..N-1} ((-logΠ)^tr / tr!)

        prob_jet = np.minimum(
            np.exp(np.maximum(np.log(np.maximum(prob, 1e-30)) + prodproba_log, -30.0)),
            1.0,
        )  # dim: (evt, jet), calculating Π * Σ_tr{0..N-1} ((-logΠ)^tr / tr!)

        prob_jet = ak.where(prodproba_log < 0, prob_jet, 1.0)
        prob_jet = np.maximum(prob_jet, 1e-30)

        return prob_jet


def common_shifts(self, events):
    """
    Apply common shifts to a events(mostly affect energy resolution/scale of objects).

    This function applies common shifts to the input DataFrame based on the specified shift type.
    It modifies the DataFrame in place to reflect the systematic variations/dedicated corrections.
    This includes JERC, Rochester (Run 2), muon scale & smearing (SS), and electron SS corrections.

    Scale/Resolution Corrections:
    Construct a shift list with tuples of (obj_dict, shift_name).
    These corrections are applied independently on all objects by updating the contents of the branch.
    Normally done before selection to apply updated objects.
    Uncertainties are handled by updating object collections with up/down variations.

    Example for Shift List:
    ```python
    # nominal correction
    shift = [({"Jet": jets, "MET": met, "Muon": muon, "Electron": ele}, None)]
    # add variations
    shifts += [
                    (
                        {
                            "Jet": jets.JES_Total.up,
                            "MET": met.JES_Total.up,
                        },
                        "JESUp",
                    )]
    shifts += [
                    (
                        {
                            "Jet": jets.JES_Total.down,
                            "MET": met.JES_Total.down,
                        },
                        "JESDown",
                    )]
    ```

    Different treatment for weights and scale/resolution shifts is necessary
    to ensure accurate corrections and uncertainties are applied to the data.

    Parameters:
    self (dict): The configuration dictionary from SF_map containing the scale factors and other settings.
    events (events): The input events containing the data to be shifted.

    Returns:
    pandas.DataFrame: The DataFrame with the applied systematic shifts.
    """

    isRealData = not hasattr(events, "genWeight")
    dataset = events.metadata["dataset"]

    shifts = []

    if "JME" in self.SF_map.keys():
        syst_JERC = self.isSyst
        if self.isSyst == "JERC_split":
            syst_JERC = "split"
        shifts = JME_shifts(
            shifts,
            self.SF_map,
            events,
            self._year,
            self._campaign,
            isRealData,
            syst_JERC,
        )
    else:
        ## Using PFMET
        if int(self._year) < 2020:
            shifts = [
                (
                    {
                        "Jet": events.Jet,
                        "MET": events.MET,
                    },
                    None,
                )
            ]
        ## Using PuppiMET
        else:
            shifts = [
                (
                    {
                        "Jet": events.Jet,
                        "MET": events.PuppiMET,
                    },
                    None,
                )
            ]

    if "roccor" in self.SF_map.keys():
        shifts = Roccor_shifts(shifts, self.SF_map, events, isRealData, False)
    elif "muonSS" in self.SF_map.keys():
        shifts = MUO_shifts(shifts, self.SF_map, events, isRealData, False)
    else:
        for shift in shifts:
            shift[0]["Muon"] = events.Muon

    if "electronSS" in self.SF_map.keys():
        shifts = EGM_shifts(shifts, self.SF_map, events, isRealData, False)
    else:
        for shift in shifts:
            shift[0]["Electron"] = events.Electron

    # Apply jet veto
    if "jetveto" in self.SF_map.keys():
        jet_veto = jetveto(events.Jet, self.SF_map)
        event_veto = ak.any(jet_veto > 0, axis=1)
        vetoed_events = events[~event_veto]
        for collections, _ in shifts:
            for key in collections:
                collections[key] = collections[key][~event_veto]
    else:
        vetoed_events = events

    return vetoed_events, shifts


# common weights
def weight_manager(pruned_ev, SF_map, isSyst):
    """
    Example for Scaling Factors (SFs):
    ```python
    # evaluation depends on file types...
    ## add SFs & uncertainties to weight function
    weights.add(sf.split(" ")[0], sfs_alle, sfs_alle_up, sfs_alle_down)
    """
    weights = Weights(len(pruned_ev), storeIndividual=True)
    # Gen info
    if "genWeight" in pruned_ev.fields:
        weights.add("genweight", pruned_ev.genWeight)
    if "PSWeight" in pruned_ev.fields:
        # PS ISR/FSR weights
        add_ps_weight(weights, pruned_ev.PSWeight, isSyst)
    if "LHEPdfWeight" in pruned_ev.fields:
        add_pdf_weight(weights, pruned_ev.LHEPdfWeight, isSyst)
    if "LHEScaleWeight" in pruned_ev.fields:
        add_scalevar_weight(weights, pruned_ev.LHEScaleWeight, isSyst)
    if "TT" in pruned_ev.metadata["dataset"]:
        weights.add(
            "ttbar_weight",
            top_pT_reweighting(pruned_ev.GenPart),
            (
                top_pT_reweighting(pruned_ev.GenPart)
                - ak.ones_like(top_pT_reweighting(pruned_ev.GenPart))
            )
            * 2.0
            + ak.ones_like(top_pT_reweighting(pruned_ev.GenPart)),
        )
        if isSyst != False:
            weights.add(
                "ttbar_weight",
                top_pT_reweighting(pruned_ev.GenPart),
                (
                    top_pT_reweighting(pruned_ev.GenPart)
                    - ak.ones_like(top_pT_reweighting(pruned_ev.GenPart))
                )
                * 2.0,
                ak.ones_like(top_pT_reweighting(pruned_ev.GenPart)),
            )

    if "hadronFlavour" in pruned_ev.Jet.fields:
        syst_wei = True if isSyst != False else False
        if "LUM" in SF_map.keys():
            puwei(
                pruned_ev.Pileup.nTrueInt,
                SF_map,
                weights,
                syst_wei,
            )
        if "MUO" in SF_map.keys() and "SelMuon" in pruned_ev.fields:
            muSFs(pruned_ev.SelMuon, SF_map, weights, syst_wei, False)
        if "EGM" in SF_map.keys() and "SelElectron" in pruned_ev.fields:
            eleSFs(pruned_ev.SelElectron, SF_map, weights, syst_wei, False)
        if "BTV" in SF_map.keys() and "SelJet" in pruned_ev.fields:
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepJetC", syst_wei)
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepJetB", syst_wei)
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepCSVB", syst_wei)
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepCSVC", syst_wei)

    return weights
