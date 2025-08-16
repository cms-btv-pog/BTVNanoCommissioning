import os
from collections import defaultdict
import time, json


# Based on https://github.com/PocketCoffea/PocketCoffea/blob/main/pocket_coffea/utils/rucio.py
def get_xrootd_sites_map():

    in_ci_env = "CI" in os.environ or "GITLAB_CI" in os.environ

    sites_xrootd_access = defaultdict(dict)
    # Check if the cache file has been modified in the last 10 minutes
    cache_valid = False
    if os.path.exists(".sites_map.json"):
        file_time = os.path.getmtime(".sites_map.json")
        current_time = time.time()
        # ten_minutes_ago = current_time - 600
        twenty_minutes_ago = current_time - 1200
        sixty_minutes_ago = current_time - 3600
        if file_time > sixty_minutes_ago:
            cache_valid = True

    if not os.path.exists(".sites_map.json") or not cache_valid or in_ci_env:
        print("Loading SITECONF info")
        sites = [
            (s, "/cvmfs/cms.cern.ch/SITECONF/" + s + "/storage.json")
            for s in os.listdir("/cvmfs/cms.cern.ch/SITECONF/")
            if s.startswith("T")
        ]
        for site_name, conf in sites:
            if not os.path.exists(conf):
                continue
            try:
                data = json.load(open(conf))
            except:
                continue
            for site in data:
                if site["type"] != "DISK":
                    continue
                if site["rse"] == None:
                    continue
                for proc in site["protocols"]:
                    if proc["protocol"] == "XRootD":
                        if proc["access"] not in ["global-ro", "global-rw"]:
                            continue
                        if "prefix" not in proc:
                            if "rules" in proc:
                                for rule in proc["rules"]:
                                    sites_xrootd_access[site["rse"]][rule["lfn"]] = (
                                        rule["pfn"]
                                    )
                        else:
                            sites_xrootd_access[site["rse"]] = proc["prefix"]
        json.dump(sites_xrootd_access, open(".sites_map.json", "w"))

    return json.load(open(".sites_map.json"))


def find_other_file(filepath, sitemap):
    if filepath.startswith("root:/"):
        rootpref = filepath.split("/store/")[0]
        file = "/store/" + filepath.split("/store/")[1]
    else:
        rootpref = None
        file = filepath

    command = f'/cvmfs/cms.cern.ch/common/dasgoclient -query="site file={file}"'
    sites = os.popen(command, "r").read().split()
    for site in sites:
        if site not in sitemap:
            continue
        sitepath = sitemap[site]
        if not isinstance(sitepath, str):
            continue
        if rootpref:
            if rootpref in sitepath or sitepath in rootpref:
                continue
        return sitepath + file

    return filepath
