import copy
import hist
from coffea.processor import accumulate
import os
from BTVNanoCommissioning.helpers.xsection import xsection

"""
Scale histograms to corresponding cross-section. Merge mutiple `.coffea` and collate the MC samples into sub-class in this function.
"""
# from BTVNanoCommissioning.helpers.xsection_13TeV import xsection_13TeV
import numpy as np


def scale_xs(hist, lumi, events):
    xs_dict = {}
    for obj in xsection:
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
    for obj in xsection:
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
                h_scalevar_muFUp = copy.deepcopy(h_obj)
                h_scalevar_muFDown = copy.deepcopy(h_obj)
                h_scalevar_muRUp = copy.deepcopy(h_obj)
                h_scalevar_muRDown = copy.deepcopy(h_obj)
                h_PS_ISRUp = copy.deepcopy(h_obj)
                h_PS_ISRDown = copy.deepcopy(h_obj)
                h_PS_FSRUp = copy.deepcopy(h_obj)
                h_PS_FSRDown = copy.deepcopy(h_obj)
                h_PDF_weightUp = copy.deepcopy(h_obj)
                h_PDF_weightDown = copy.deepcopy(h_obj)
                h_aS_weightUp = copy.deepcopy(h_obj)
                h_aS_weightDown = copy.deepcopy(h_obj)

                if sample in xs_dict.keys():

                    nLHEScaleSumw = 9
                    temp_reshaped_LHEScaleSumw = np.reshape(
                        merged_output[sample]["LHEScaleSumw"],
                        (
                            int(len(merged_output[sample]["LHEScaleSumw"]) / nLHEScaleSumw),
                            nLHEScaleSumw
                        )
                    )
                    reshaped_LHEScaleSumw = []
                    for i in range(nLHEScaleSumw):
                        reshaped_LHEScaleSumw.append(
                            np.sum(temp_reshaped_LHEScaleSumw[:, i])
                        )

                    nLHEPdfSumw = 103
                    temp_reshaped_LHEPdfSumw = np.reshape(
                        merged_output[sample]["LHEPdfSumw"],
                        (
                            int(len(merged_output[sample]["LHEPdfSumw"]) / nLHEPdfSumw),
                            nLHEPdfSumw
                        )
                    )
                    reshaped_LHEPdfSumw = []
                    for i in range(nLHEPdfSumw):
                        reshaped_LHEPdfSumw.append(
                            np.sum(temp_reshaped_LHEPdfSumw[:, i])
                        )

                    nPSSumw = 4
                    temp_reshaped_PSSumw = np.reshape(
                        merged_output[sample]["PSSumw"],
                        (
                            int(len(merged_output[sample]["PSSumw"]) / nPSSumw),
                            nPSSumw
                        )
                    )
                    reshaped_PSSumw = []
                    for i in range(nPSSumw):
                        reshaped_PSSumw.append(
                            np.sum(temp_reshaped_PSSumw[:, i])
                        )

                    scaled[sample]["LHEScaleSumw"] = reshaped_LHEScaleSumw
                    scaled[sample]["LHEPdfSumw"] = reshaped_LHEPdfSumw
                    scaled[sample]["PSSumw"] = reshaped_PSSumw

                    ############################################################################################
                    # WARNING: THIS IS WRONG
                    delta = scaled[sample]["LHEPdfSumw"][1:-2] - scaled[sample]["LHEPdfSumw"][0]
                    pdfWeightSumw = np.sqrt(np.sum(np.square(delta)))
                    aSWeightSumw = 0.5 * (scaled[sample]["LHEPdfSumw"][102] - scaled[sample]["LHEPdfSumw"][101])
                    ############################################################################################

                    h = h * xs_dict[sample] * lumi / merged_output[sample]["sumw"]
                    h_scalevar_muFUp = h_scalevar_muFUp * xs_dict[sample] * lumi / scaled[sample]["LHEScaleSumw"][5]
                    h_scalevar_muFDown = h_scalevar_muFDown * xs_dict[sample] * lumi / scaled[sample]["LHEScaleSumw"][3]
                    h_scalevar_muRUp = h_scalevar_muRUp * xs_dict[sample] * lumi / scaled[sample]["LHEScaleSumw"][7]
                    h_scalevar_muRDown = h_scalevar_muRDown * xs_dict[sample] * lumi / scaled[sample]["LHEScaleSumw"][1]
                    h_PS_ISRUp = h_PS_ISRUp * xs_dict[sample] * lumi / scaled[sample]["PSSumw"][0]
                    h_PS_ISRDown = h_PS_ISRDown * xs_dict[sample] * lumi / scaled[sample]["PSSumw"][2]
                    h_PS_FSRUp = h_PS_FSRUp * xs_dict[sample] * lumi / scaled[sample]["PSSumw"][1]
                    h_PS_FSRDown = h_PS_FSRDown * xs_dict[sample] * lumi / scaled[sample]["PSSumw"][3]
                    h_PDF_weightUp = h_PDF_weightUp * xs_dict[sample] * lumi / (scaled[sample]["sumw"] + pdfWeightSumw)
                    h_PDF_weightDown = h_PDF_weightDown * xs_dict[sample] * lumi / (scaled[sample]["sumw"] - pdfWeightSumw)
                    h_aS_weightUp = h_aS_weightUp * xs_dict[sample] * lumi / (scaled[sample]["sumw"] + aSWeightSumw)
                    h_aS_weightDown = h_aS_weightDown * xs_dict[sample] * lumi / (scaled[sample]["sumw"] - aSWeightSumw)
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
                    scaled[sample][f"{key}_scalevar_muFUp"] = h_scalevar_muFUp
                    scaled[sample][f"{key}_scalevar_muFDown"] = h_scalevar_muFDown
                    scaled[sample][f"{key}_scalevar_muRUp"] = h_scalevar_muRUp
                    scaled[sample][f"{key}_scalevar_muRDown"] = h_scalevar_muRDown
                    scaled[sample][f"{key}_UEPS_ISRUp"] = h_PS_ISRUp
                    scaled[sample][f"{key}_UEPS_ISRDown"] = h_PS_ISRDown
                    scaled[sample][f"{key}_UEPS_FSRUp"] = h_PS_FSRUp
                    scaled[sample][f"{key}_UEPS_FSRDown"] = h_PS_FSRDown
                    scaled[sample][f"{key}_PDF_weightUp"] = h_PDF_weightUp
                    scaled[sample][f"{key}_PDF_weightDown"] = h_PDF_weightDown
                    scaled[sample][f"{key}_aS_weightUp"] = h_aS_weightUp
                    scaled[sample][f"{key}_aS_weightDown"] = h_aS_weightDown
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
