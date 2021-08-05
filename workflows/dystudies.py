import coffea
from coffea import hist, processor
import awkward as ak
import functools as ft
import operator as op
from .base import HggBaseProcessor

class DYStudiesProcessor(HggBaseProcessor):
    def __init__(
        self,
        metaconditions,
        do_systematics=False,
        apply_trigger=False,
        output_location=None,
        taggers=None,
    ):
        super(DYStudiesProcessor, self).__init__(
            metaconditions,
            do_systematics=do_systematics,
            apply_trigger=apply_trigger,
            output_location=output_location,
            taggers=taggers,
            trigger_group=".*DoubleEG.*",
            analysis="mainAnalysis",
        )
        self.trigger_group = ".*DoubleEG.*"
        self.analysis = "mainAnalysis"

    def process_extra(self, events):
        return events

    def postprocess(self, accumulant):
        pass
