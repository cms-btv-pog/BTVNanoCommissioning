import awkward
import numpy


class DummyTagger2:
    def __init__(self) -> None:
        pass

    @property
    def name(self) -> str:
        return "DummyTagger2"

    # lower priority is better
    # first decimal point is category within tag (in case for untagged)
    @property
    def priority(self) -> int:
        return 10

    def __call__(self, events: awkward.Array) -> awkward.Array:
        counts = awkward.num(events.diphotons, axis=1)
        x = numpy.random.exponential(size=awkward.sum(counts))
        return awkward.unflatten(self.priority * (x > 2), counts), {}
