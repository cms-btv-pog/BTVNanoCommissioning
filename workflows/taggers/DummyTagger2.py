import awkward as ak
import numpy as np


class DummyTagger2:
    def __init__(self):
        pass

    @property
    def name(self):
        return "DummyTagger2"

    # lower priority is better
    # first decimal point is category within tag (in case for untagged)
    @property
    def priority(self):
        return 10

    def __call__(self, events):
        counts = ak.num(events.diphotons, axis=1)
        x = np.random.exponential(size=ak.sum(counts))
        return ak.unflatten(self.priority * (x > 2), counts)
