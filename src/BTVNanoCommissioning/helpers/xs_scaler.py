import copy
import hist
from coffea.processor import accumulate
import os
from BTVNanoCommissioning.helpers.xsection import xsection

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
    merged = {}
    merged_output = accumulate([output[f] for f in output.keys()])

    for sample, accu in merged_output.items():
        scaled[sample] = {}
        for key, h_obj in accu.items():
            scaled[sample]["sumw"] = merged_output[sample]["sumw"]
            if isinstance(h_obj, hist.Hist):
                h = copy.deepcopy(h_obj)
                if sample in xs_dict.keys():
                    h = h * xs_dict[sample] * lumi / merged_output[sample]["sumw"]
                else:
                    if ("data" in sample) or ("Run" in sample) or ("Double" in sample):
                        h = h
                    else:
                        raise KeyError(
                            sample,
                            "is not found in xsection.py. If you're using 13TeV samples, please use xsection_13TeV.py",
                        )

                scaled[sample][key] = h
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


def collate(merged_output, mergemap):
    out = {}
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
