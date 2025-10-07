## code source: Master Thesis Annika Stein
## https://git.rwth-aachen.de/annika-stein/ai_safety_2021/-/blob/master/june_21/evaluate/definitions.py
from .definition_vars.pfcands import get_variables as pfcands_definitions
from .definition_vars.deepjet import get_variables as deepjet_definitions
from .definition_vars.deepcsv import get_variables as deepcsv_definitions
from .definition_vars.sv import get_variables as sv_definitions


def get_definitions(
    include_definitions: list = ["DeepJet", "DeepCSV", "PFCands"],
):
    """
    Add new definitions to the definitions dictionary.

    Parameters:
    definitions_dict (dict): The dictionary to which new definitions will be added.

    Example:
    ```python
    definitions_dict["DeepCSV_jetNSelectedTracks"] = {
        "displayname": "Jet N Selected Tracks",
        "manual_ranges": [0.0, 25],
        "ylabel_text": "Jets",
        "format_unit": "2f",
        "format_unit_digits": 2,
        "bins": 25,
        "inputVar_units": None,
    }
    ```

    Returns:
    dict: with defitions of the tagger input variables added to the dictionary.
    """

    available_definitions = {
        "PFCands": pfcands_definitions(),
        "DeepJet": deepjet_definitions(),
        "DeepCSV": deepcsv_definitions(),
        "SV": sv_definitions(),
    }

    output = {}

    for d in include_definitions:
        if d not in available_definitions:
            raise ValueError(
                f"Definitions for variable {d} not found {available_definitions.keys()}. Check utils/histogramming/definitions"
            )
        output = output | available_definitions[d]

    return output


def get_discriminators():
    disc_list = [
        "btagDeepFlavB",
        "btagDeepFlavC",
        "btagDeepFlavCvL",
        "btagDeepFlavCvB",
        "btagDeepFlavB_b",
        "btagDeepFlavB_bb",
        "btagPNetB",
        "btagPNetCvB",
        "btagPNetCvL",
        "btagPNetCvNotB",
        "btagPNetProbB",
        "btagPNetProbC",
        "btagPNetProbG",
        "btagPNetProbUDS",
        "btagPNetQvG",
        "btagPNetTauVJet",
        "btagRobustParTAK4B",
        "btagRobustParTAK4B_b",
        "btagRobustParTAK4B_bb",
        "btagRobustParTAK4B_lepb",
        "btagRobustParTAK4C",
        "btagRobustParTAK4G",
        "btagRobustParTAK4UDS",
        "btagRobustParTAK4CvB",
        "btagRobustParTAK4CvL",
        "btagRobustParTAK4QG",
        "btagUParTAK4B",
        "btagUParTAK4CvL",
        "btagUParTAK4CvB",
        "btagUParTAK4CvNotB",
        "btagUParTAK4QvG",
        "btagUParTAK4TauVJet",
        ## Negative tagger
        "btagNegDeepFlavB",
        "btagNegDeepFlavB_b",
        "btagNegDeepFlavB_bb",
        "btagNegDeepFlavB_lepb",
        "btagNegDeepFlavC",
        "btagNegDeepFlavCvB",
        "btagNegDeepFlavCvL",
        "btagNegDeepFlavG",
        "btagNegDeepFlavQG",
        "btagNegDeepFlavUDS",
        "btagNegPNetB",
        "btagNegPNetC",
        "btagNegPNetCvB",
        "btagNegPNetCvL",
        "btagNegPNetProbB",
        "btagNegPNetProbC",
        "btagNegPNetProbG",
        "btagNegPNetProbUDS",
        "btagNegRobustParTAK4B",
        "btagNegRobustParTAK4B_b",
        "btagNegRobustParTAK4B_bb",
        "btagNegRobustParTAK4B_lepb",
        "btagNegRobustParTAK4C",
        "btagNegRobustParTAK4CvB",
        "btagNegRobustParTAK4CvL",
        "btagNegRobustParTAK4G",
        "btagNegRobustParTAK4QG",
        "btagNegRobustParTAK4UDS",
        # other prob info
        "PNetRegPtRawCorr",
        "PNetRegPtRawCorrNeutrino",
        "PNetRegPtRawRes",
        "UParTAK4RegPtRawRes",
        "UParTAK4RegPtRawCorrNeutrino",
        "UParTAK4RegPtRawCorr",
        "Bprob",
        "BprobN",
        "ProbaN",
    ]

    return disc_list


def axes_name(var):

    output = {}

    unit = ""
    obj = ""
    if "dr" in var:
        obj = "$\\Delta$R"
        if "lmujethmu" in var:
            unit = "($\\mu$-Jet,$\\mu$)"
        elif "lmujetsmu" in var:
            unit = "($\\mu$-Jet,soft-$\\mu$)"
        elif "lmusmu" in var:
            unit = "($\\mu$,soft-$\\mu$)"
        elif "mujet" in var:
            unit = "($\\mu$,Jet)"
        elif "SVjet0" in var:
            unit = "(SV,Jet)"
    elif "MET_" in var:
        obj = "MET"
    elif "ele_" in var:
        obj = "e"
    elif "posl_" in var:
        obj = "$\\ell^+$"
    elif "negl_" in var:
        obj = "$\\ell^-$"
    elif "mu_" in var:
        obj = "$\\mu$"
    elif "hl_" in var:
        obj = "$\\ell_1$"
    elif "sl_" in var:
        obj = "$\\ell_2$"
    elif "soft_l" in var:
        obj = "soft-$\\mu$"
    elif "mujet" in var:
        obj = "$\\mu$-Jet"
    elif "jet" in var:
        obj = "Jet"
    elif "w_" in var:
        obj = "W"
    elif "z_" in var:
        obj = "Z"
    if "pt" in var:
        unit = " $p_T$ [GeV]"
    elif "mass" in var:
        unit = " mass [GeV]"
    elif "eta" in var:
        unit = " $\\eta$"
    elif "phi" in var:
        unit = " $\\phi$"
    elif "dxy" in var:
        unit = " $d_{xy}$ [cm]"
    elif "dz" in var:
        unit = " $d_{z}$ [cm]"
    elif "pfRelIso" in var:
        unit = " Rel. Iso"
    elif "btag" in var:
        unit = "Jet"
        ## Negative tagger
        if "Neg" in var:
            unit = "Neg. " + unit
        if "tanh" in var:
            unit = "Trans. " + unit
        ## Different taggers
        if "DeepFlav" in var:
            unit = unit + " DeepJet"
        elif "RobustParTAK4" in var:
            unit = unit + " RobustParTAK4"
        elif "PNet" in var:
            unit = unit + " PNet"
        elif "UParT" in var:
            unit = unit + " UParTAK4"
        else:
            unit = unit
        # output node
        if "CvL" in var:
            unit = unit + " CvL"
        elif "CvB" in var:
            unit = unit + " CvB"
        elif "CvNotB" in var:
            unit = unit + " CvNotB"
        elif "B_b" in var or "ProbB" in var:
            unit = unit + " Prob(b)"
        elif "B_bb" in var:
            unit = unit + " Prob(bb)"
        elif "B_lepb" in var:
            unit = unit + " Prob(lepb)"
        elif "QvG" in var or "QG" in var:
            unit = unit + " QvG"
        elif "G" in var:
            unit = unit + " Prob(g)"
        elif "UDS" in var:
            unit = unit + " Prob(uds)"
        elif "TauVJet" in var:
            unit = unit + " $\\tau$vJ"
        elif "B" in var:
            unit = unit + " BvAll"
        elif "C" in var:
            unit = unit + " Prob(c)"
        elif "PNetReg" in var:
            unit = "Jet PNet Reg."
            if "Corr" in var:
                unit = " Corr"
            elif "CorrNeutrino" in var:
                unit = " Corr (w/$\\nu$)"
            elif "RegPtRawRes" in var:
                unit = " Resoultion"
        elif "ProbaN" in var:
            unit = "Neg. Jet Probability"
        elif "Proba" in var:
            unit = "Jet Probability"
        elif "BprobN" in var:
            unit = "Neg. b-Jet Probability"
        elif "Bprob" in var:
            unit = "b-Jet Probability"
    label = obj + unit
    if var.endswith("0") or "jet0" in var:
        if "btagDeepFlav" not in var:
            label = label.replace("Jet", "$1^{st}$-Jet")
        else:
            label = label.replace("Jet", "$1^{st}$-Jet", 1)
    elif var.endswith("1") or "jet1" in var:
        if "btagDeepFlav" not in var:
            label = label.replace("Jet", "$2^{nd}$-Jet")
        else:
            label = label.replace("Jet", "$2^{nd}$-Jet", 1)
    elif var.endswith("2") or "jet2" in var:
        if "btagDeepFlav" not in var:
            label = label.replace("Jet", "$3^{rd}$-Jet")
        else:
            label = label.replace("Jet", "$3^{rd}$-Jet", 1)
    elif var.endswith("3") or "jet3" in var:
        if "btagDeepFlav" not in var:
            label = label.replace("Jet", "$4^{th}$-Jet")
        else:
            label = label.replace("Jet", "$4^{th}$-Jet", 1)
    if "ptratio" in var:
        if var == "hl_ptratio":
            label = "$p_T^{\\ell_1}/p_T^{Jet}$"
        if var == "sl_ptratio":
            label = "$p_T^{\\ell_2}/p_T^{Jet}$"
        if var == "soft_l_ptratio":
            label = "$p_T^{soft-\\mu}/p_T^{Jet}$"
    return label
