import copy
import hist
from coffea.processor import accumulate
import os
from BTVNanoCommissioning.helpers.xsection import xsection

"""
Scale histograms to corresponding cross-section. Merge mutiple `.coffea` and collate the MC samples into sub-class in this function.
"""
from BTVNanoCommissioning.helpers.xsection_13TeV import xsection_13TeV
import numpy as np


def scale_xs(hist, lumi, events):
    xs_dict = {}
    for obj in xsection + xsection_13TeV:
        xs_dict[obj["process_name"]] = float(obj["cross_section"])
    scales = {}
    for key in events:
        if type(key) != str or "Run" in key:
            continue
        scales[key] = xs_dict[key] * lumi / events[key]
    hist.scale(scales, axis="dataset")
    return hist


def scaleSumW(output, lumi):
    scaled = {}
    xs_dict = {}
    for obj in xsection + xsection_13TeV:
        xs_dict[obj["process_name"]] = float(obj["cross_section"])
        # if k-factor in the xsection: multiply by the k-factor
        if "kFactor" in obj.keys() and obj["kFactor"] != "":
            xs_dict[obj["process_name"]] = xs_dict[obj["process_name"]] * float(
                obj["kFactor"]
            )
    merged_output = merge_output(output)

    for sample, accu in merged_output.items():
        scaled[sample] = {}
        for key, h_obj in accu.items():
            scaled[sample]["sumw"] = merged_output[sample]["sumw"]
            if isinstance(h_obj, hist.Hist):
                h = copy.deepcopy(h_obj)
                h_PDF_weightUp = copy.deepcopy(h_obj)
                h_PDF_weightDown = copy.deepcopy(h_obj)
                h_aS_weightUp = copy.deepcopy(h_obj)
                h_aS_weightDown = copy.deepcopy(h_obj)
                h_scalevar_muRUp = copy.deepcopy(h_obj)
                h_scalevar_muRDown = copy.deepcopy(h_obj)
                h_scalevar_muFUp = copy.deepcopy(h_obj)
                h_scalevar_muFDown = copy.deepcopy(h_obj)
                h_PS_ISRUp = copy.deepcopy(h_obj)
                h_PS_ISRDown = copy.deepcopy(h_obj)
                h_PS_FSRUp = copy.deepcopy(h_obj)
                h_PS_FSRDown = copy.deepcopy(h_obj)

                if sample in xs_dict.keys():
                    for syst in ["PDF", "aS", "PDFaS", "muR", "muF", "ISR", "FSR"]:
                        for var in ["Up", "Down"]:
                            key_reweight = f"{syst}_sumw{var}"
                            if key_reweight in merged_output[sample].keys():
                                scaled[sample][key_reweight] = merged_output[sample][
                                    key_reweight
                                ]
                            else:
                                scaled[sample][key_reweight] = merged_output[sample][
                                    "sumw"
                                ]
                                print(f"WARNING: {key_reweight} not found!")

                    h = h * xs_dict[sample] * lumi / merged_output[sample]["sumw"]
                    h_PDF_weightUp = (
                        h_PDF_weightUp
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["PDF_sumwUp"]
                    )
                    h_PDF_weightDown = (
                        h_PDF_weightDown
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["PDF_sumwDown"]
                    )
                    h_aS_weightUp = (
                        h_aS_weightUp
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["aS_sumwUp"]
                    )
                    h_aS_weightDown = (
                        h_aS_weightDown
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["aS_sumwDown"]
                    )
                    h_scalevar_muRUp = (
                        h_scalevar_muRUp
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["muR_sumwUp"]
                    )
                    h_scalevar_muRDown = (
                        h_scalevar_muRDown
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["muR_sumwDown"]
                    )
                    h_scalevar_muFUp = (
                        h_scalevar_muFUp
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["muF_sumwUp"]
                    )
                    h_scalevar_muFDown = (
                        h_scalevar_muFDown
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["muF_sumwDown"]
                    )
                    h_PS_ISRUp = (
                        h_PS_ISRUp
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["ISR_sumwUp"]
                    )
                    h_PS_ISRDown = (
                        h_PS_ISRDown
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["ISR_sumwDown"]
                    )
                    h_PS_FSRUp = (
                        h_PS_FSRUp
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["FSR_sumwUp"]
                    )
                    h_PS_FSRDown = (
                        h_PS_FSRDown
                        * xs_dict[sample]
                        * lumi
                        / scaled[sample]["FSR_sumwDown"]
                    )
                else:
                    if ("data" in sample) or ("Run" in sample) or ("Double" in sample):
                        h = h
                    else:
                        raise KeyError(
                            sample,
                            "is not found in xsection.py. If you're using 13TeV samples, please use xsection_13TeV.py",
                        )

                scaled[sample][key] = h
                if sample in xs_dict.keys():
                    scaled[sample][f"{key}_PDF_weightUp"] = h_PDF_weightUp
                    scaled[sample][f"{key}_PDF_weightDown"] = h_PDF_weightDown
                    scaled[sample][f"{key}_aS_weightUp"] = h_aS_weightUp
                    scaled[sample][f"{key}_aS_weightDown"] = h_aS_weightDown
                    scaled[sample][f"{key}_scalevar_muRUp"] = h_scalevar_muRUp
                    scaled[sample][f"{key}_scalevar_muRDown"] = h_scalevar_muRDown
                    scaled[sample][f"{key}_scalevar_muFUp"] = h_scalevar_muFUp
                    scaled[sample][f"{key}_scalevar_muFDown"] = h_scalevar_muFDown
                    scaled[sample][f"{key}_UEPS_ISRUp"] = h_PS_ISRUp
                    scaled[sample][f"{key}_UEPS_ISRDown"] = h_PS_ISRDown
                    scaled[sample][f"{key}_UEPS_FSRUp"] = h_PS_FSRUp
                    scaled[sample][f"{key}_UEPS_FSRDown"] = h_PS_FSRDown
    return scaled


## Additional rescale for MC
def additional_scale(output, scale, sample_to_scale):
    scaled = {}
    for files in output.keys():
        scaled[files] = {}
        if "sumw" not in output[files].keys():
            for sample, accu in output[files].items():
                scaled[files][sample] = {}
                for key, h_obj in accu.items():
                    if isinstance(h_obj, hist.Hist):
                        h = copy.deepcopy(h_obj)
                        if sample in sample_to_scale:
                            h = h * scale
                        else:
                            h = h
                        scaled[files][sample][key] = h
        else:
            for sample, accu in output.items():
                scaled[sample] = {}
                for key, h_obj in accu.items():
                    if isinstance(h_obj, hist.Hist):
                        h = copy.deepcopy(h_obj)
                        if sample in sample_to_scale:
                            h = h * scale
                        else:
                            h = h
                        scaled[sample][key] = h
    return scaled


def dict_depth(d):
    if isinstance(d, dict):
        return 1 + (max(map(dict_depth, d.values())) if d else 0)
    return 0


def merge_output(output):
    if dict_depth(output) <= 2:
        return output
    else:
        return accumulate([output[f] for f in output.keys()])


def collate(merged_output, mergemap):
    out = {}
    merged_output = merge_output(merged_output)
    for group, names in mergemap.items():
        out[group] = accumulate(
            [
                copy.deepcopy(v)
                for k, v in merged_output.items()
                if k.split("_FNAME_")[0] in names
            ]
        )
    return out


def getSumW(accumulator):
    sumw = {}
    for key, accus in accumulator.items():
        sumw[key] = accus["sumw"]
    return sumw
