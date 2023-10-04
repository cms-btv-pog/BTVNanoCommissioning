import copy
import hist
from coffea.processor import accumulate
import os
from BTVNanoCommissioning.helpers.xsection import xsection
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
    duplicated_name = False
    sumw = {}
    flist = []
    for f in output.keys():
        flist.extend([m for m in output[f].keys() if "Run" not in m])
    for files in output.keys():
        if "sumw" not in output[files].keys() and len(flist) != len(set(flist)):
            duplicated_name = True
            for sample in output[files].keys():
                if "Run" in str(output[files][sample]):
                    continue
                if sample in sumw.keys():
                    sumw[sample] = sumw[sample] + float(output[files][sample]["sumw"])
                else:
                    sumw[sample] = float(output[files][sample]["sumw"])
    for files in output.keys():
        if "sumw" not in output[files].keys():
            scaled[files] = {}
            for sample, accu in output[files].items():
                scaled[files][sample] = {}
                scaled[files][sample]["sumw"] = output[files][sample]["sumw"]
                if duplicated_name:
                    scaled[files][sample]["sumw"] = sumw[sample]
                for key, h_obj in accu.items():
                    if isinstance(h_obj, hist.Hist):
                        h = copy.deepcopy(h_obj)
                        if sample in xs_dict.keys():
                            h = (
                                h
                                * xs_dict[sample]
                                * lumi
                                / scaled[files][sample]["sumw"]
                            )
                        else:
                            if not (("data" in sample) or ("Run" in sample)):
                                raise KeyError(sample, "is not found in xsection.py")
                            else:
                                h = h
                        scaled[files][sample][key] = h
        else:
            for sample, accu in output[files].items():
                scaled[sample] = {}
                for key, h_obj in accu.items():
                    scaled[sample]["sumw"] = output[files]["sumw"]
                    if isinstance(h_obj, hist.Hist):
                        h = copy.deepcopy(h_obj)
                        if sample in xs_dict.keys():
                            h = h * xs_dict[sample] * lumi / output[files]["sumw"]
                        else:
                            if not (("data" in sample) or ("Run" in sample)) or (
                                "Double" in sample
                            ):
                                raise KeyError(sample, "is not found in xsection.py")
                            else:
                                h = h
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


def collate(output, mergemap):
    out = {}
    merged = {}
    merged_output = accumulate([output[f] for f in output.keys()])
    for files in merged_output.keys():
        if "sumw" not in merged_output[files].keys():
            for m in output[files].keys():
                merged[m] = dict(merged_output[files][m].items())
        else:
            merged[files] = dict(merged_output[files].items())
    for group, names in mergemap.items():
        out[group] = accumulate(
            [v for k, v in merged.items() if k.split("_FNAME_")[0] in names]
        )
    return out


def getSumW(accumulator):
    sumw = {}
    for key, accus in accumulator.items():
        sumw[key] = accus["sumw"]
    return sumw
