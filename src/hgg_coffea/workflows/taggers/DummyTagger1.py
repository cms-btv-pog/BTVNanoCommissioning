import awkward as ak
import numpy as np


class DummyTagger1:
    def __init__(self):
        pass

    @property
    def name(self):
        return "DummyTagger1"

    # lower priority is better
    # first decimal point is category within tag (in case for untagged)
    @property
    def priority(self):
        return 20

    def __call__(self, events):
        return self.priority * ak.ones_like(events.diphotons.pt, dtype=np.int32)
