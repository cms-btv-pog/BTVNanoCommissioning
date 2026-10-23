"""Instantiate and pickle every registered workflow, without reading any events.

Catches runner/processor signature mismatches (e.g. a new keyword passed by
runner.py that some processor's __init__ does not accept) and processors that
cannot be shipped to futures/dask/parsl/condor workers.

    python scripts/check_workflows.py --campaign Summer23 --year 2023
"""

import argparse
import os
import pickle
import sys
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from runner import build_processor, get_main_parser  # noqa: E402
from BTVNanoCommissioning.workflows import workflows  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--campaign", default="Summer23")
parser.add_argument("--year", default="2023")
parser.add_argument(
    "--isSyst",
    default="False",
    help="Forwarded to runner.py's --isSyst (default: %(default)s)",
)
opts = parser.parse_args()

# Workflows that refuse to build without a --selectionModifier
required_modifiers = {"sf_ttsemilep_tnp": ["tt_semilep_mu", "tt_semilep_el"]}

failed = []
ntot = 0
for wf in workflows:
    for selmod in required_modifiers.get(wf, [""]):
        ntot += 1
        label = f"{wf} --selectionModifier {selmod}" if selmod else wf
        args = get_main_parser().parse_args(
            [
                "--wf",
                wf,
                "--campaign",
                opts.campaign,
                "--year",
                opts.year,
                "--isSyst",
                opts.isSyst,
                f"--selectionModifier={selmod}",
            ]
        )
        if args.isSyst == "False":
            args.isSyst = False
        try:
            pickle.loads(pickle.dumps(build_processor(args, f"arrays_{wf}")))
            print(f"OK    {label}")
        except Exception:
            print(f"FAIL  {label}")
            traceback.print_exc()
            failed.append(label)

print(f"\n{ntot - len(failed)}/{ntot} workflows built and pickled")
if failed:
    print("Failed: " + ", ".join(failed))
    sys.exit(1)
