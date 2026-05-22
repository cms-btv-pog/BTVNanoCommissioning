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


def scaleSumW(output, lumi, syst_shapes=True):
    scaled = {}
    xs_dict = {}
    for obj in xsection + xsection_13TeV:
        xs_dict[obj["process_name"]] = float(obj["cross_section"])
        # if k-factor in the xsection: multiply by the k-factor
        if "kFactor" in obj.keys() and obj["kFactor"] != "":
            xs_dict[obj["process_name"]] = xs_dict[obj["process_name"]] * float(
                obj["kFactor"]
            )
    # Maps output histogram suffix -> sumw accumulator key used for normalisation.
    # When syst_shapes=False every variation is normalised with the nominal sumw instead.
    syst_sumw_map = {
        "PDF_weightUp": "PDF_sumwUp",
        "PDF_weightDown": "PDF_sumwDown",
        "aS_weightUp": "aS_sumwUp",
        "aS_weightDown": "aS_sumwDown",
        "scalevar_muRUp": "muR_sumwUp",
        "scalevar_muRDown": "muR_sumwDown",
        "scalevar_muFUp": "muF_sumwUp",
        "scalevar_muFDown": "muF_sumwDown",
        "UEPS_ISRUp": "ISR_sumwUp",
        "UEPS_ISRDown": "ISR_sumwDown",
        "UEPS_FSRUp": "FSR_sumwUp",
        "UEPS_FSRDown": "FSR_sumwDown",
    }

    merged_output = merge_output(output)

    for sample, accu in merged_output.items():
        scaled[sample] = {}
        if "sumw" not in accu.keys():
            continue
        for key, h_obj in accu.items():
            scaled[sample]["sumw"] = merged_output[sample]["sumw"]
            if not isinstance(h_obj, hist.Hist):
                continue

            if sample in xs_dict.keys():
                xs = xs_dict[sample] * lumi
                nominal_sumw = merged_output[sample]["sumw"]

                # Collect per-syst sumw accumulators (fall back to nominal if missing)
                for syst in ["PDF", "aS", "PDFaS", "muR", "muF", "ISR", "FSR"]:
                    for var in ["Up", "Down"]:
                        key_reweight = f"{syst}_sumw{var}"
                        if key_reweight in merged_output[sample].keys():
                            scaled[sample][key_reweight] = merged_output[sample][
                                key_reweight
                            ]
                        else:
                            scaled[sample][key_reweight] = nominal_sumw
                            print(f"WARNING: {key_reweight} not found!")

                scaled[sample][key] = copy.deepcopy(h_obj) * xs / nominal_sumw

                for suffix, sumw_key in syst_sumw_map.items():
                    sumw = scaled[sample][sumw_key] if syst_shapes else nominal_sumw
                    scaled[sample][f"{key}_{suffix}"] = copy.deepcopy(h_obj) * xs / sumw
            else:
                if ("data" in sample) or ("Run" in sample) or ("Double" in sample):
                    scaled[sample][key] = copy.deepcopy(h_obj)
                else:
                    raise KeyError(
                        sample,
                        "is not found in xsection.py. If you're using 13TeV samples, please use xsection_13TeV.py",
                    )
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
            [v for k, v in merged_output.items() if k.split("_FNAME_")[0] in names]
        )
    return out


def getSumW(accumulator):
    sumw = {}
    for key, accus in accumulator.items():
        sumw[key] = accus["sumw"]
    return sumw
