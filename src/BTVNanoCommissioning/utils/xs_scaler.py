import copy
import hist 
from coffea import processor
import os
from BTVNanoCommissioning.helpers.xsection import xsection



def scale_xs(hist, lumi, events):
    xs_dict = {}
    for obj in xsection:
        xs_dict[obj["process_name"]] = float(obj["cross_section"])
    scales = {}
    for key in events:
        if type(key) != str or key == "Data" or "Run" in key:
            continue
        scales[key] = xs_dict[key] * lumi / events[key]
    hist.scale(scales, axis="dataset")
    return hist


def scale_xs_arr(events, lumi):
    xs_dict = {}
    for obj in xsection:
        xs_dict[obj["process_name"]] = float(obj["cross_section"])
    scales = {}
    wei_array = {}
    for key in events:
        if type(key) != str or key == "Data" or "Run" in key:
            continue
        scales[key] = xs_dict[key] * lumi / events[key]

        wei_array[key] = scales[key]
    return wei_array
def scaleSumW(accumulator, lumi, sumw, dyscale=1.,xsfile="xsection.json"):
    scaled = {}
    xs_dict = {}
    for obj in xsection:
        xs_dict[obj["process_name"]] = float(obj["cross_section"])
    for sample, accu in accumulator.items():
        scaled[sample] = {}
        for key, h_obj in accu.items():
            scaled[sample]['sumw']=sumw[sample]
            if isinstance(h_obj, hist.Hist):
                h = copy.deepcopy(h_obj)
                if sample in xs_dict.keys():h = h * xs_dict[sample] *lumi/sumw[sample]
                else:
                    if not (("data" in sample) or ("Run" in sample)):
                        continue
                    else: h = h
                scaled[sample][key] = h
    
    return scaled
def additional_scale(accumulator,scale, target):
    scaled = {}
    for sample, accu in accumulator.items():
        scaled[sample] = {}
        for key, h_obj in accu.items():
            if isinstance(h_obj, hist.Hist):
                h = copy.deepcopy(h_obj)
                if sample in target:h = h * scale
                
                else: h = h
                scaled[sample][key] = h
    return scaled
def collate(output, mergemap):
    out = {}    
    merged = {}
    for files in output.keys():
        if 'sumw' not in output[files].keys():
            for m in output[files].keys():
                merged[m]=dict(output[files][m].items())
        else:
            
            merged[files]=dict(output[files].items())
    for group, names in mergemap.items():
        print(group,names)
        out[group] = processor.accumulate([v for k, v in merged.items() if k in names])
    
    return out
def getSumW(accumulator):
    sumw = {}
    for key, accus in accumulator.items():
        sumw[key] = accus['sumw']
    return sumw
