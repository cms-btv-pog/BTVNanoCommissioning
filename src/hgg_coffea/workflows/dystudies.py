from typing import Any, Dict, List, Optional

import awkward as ak

from hgg_coffea.workflows.base import HggBaseProcessor


class DYStudiesProcessor(HggBaseProcessor):
    def __init__(
        self,
        metaconditions: Dict[str, Any],
        do_systematics: bool = False,
        apply_trigger: bool = False,
        output_location: Optional[str] = None,
        taggers: Optional[List[Any]] = None,
    ) -> None:
        super().__init__(
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

    def process_extra(self, events: ak.Array) -> ak.Array:
        return events

    def postprocess(self, accumulant: Dict[Any, Any]) -> Any:
        pass
