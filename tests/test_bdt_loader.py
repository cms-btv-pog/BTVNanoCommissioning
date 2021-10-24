import pytest
import xgboost

from hgg_coffea.tools.xgb_loader import load_bdt


@pytest.mark.parametrize(
    "fname",
    ["tests/samples/PhoID_EB.json.gz", "tests/samples/PhoID_EB.json", "DoesNotExist"],
)
def test_load_bdt(fname):
    bdt = load_bdt(fname)
    if fname == "DoesNotExist":
        assert bdt is None
    else:
        assert isinstance(bdt, xgboost.Booster)
