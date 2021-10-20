def read_xs(file):
    import json

    f = open("../metadata/%s" %file)
    data = json.load(f)
    xs_dict={}
    for obj in data:
        xs_dict[obj['process_name']]=float(obj['cross_section'])
    return xs_dict
def scale_xs(hist,lumi,events,xsfile="../metadata/xsection.json"):
    xs_dict = read_xs(xsfile)
    scales={}
    for key in events:
        if type(key) != str or key=="Data": continue
        scales[key]=xs_dict[key]*lumi/events[key]
    hist.scale(scales, axis="dataset")
    return hist
