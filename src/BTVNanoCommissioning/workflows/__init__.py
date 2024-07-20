from functools import partial

# Validation
from BTVNanoCommissioning.workflows.validation import (
    NanoProcessor as ValidationProcessor,
)

# TTbar
from BTVNanoCommissioning.workflows.ttsemilep_valid_sf import (
    NanoProcessor as TTsemilepValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ttdilep_valid_sf import (
    NanoProcessor as TTdilepValidSFProcessor,
)

# C-tag
from BTVNanoCommissioning.workflows.ctag_emdileptt_valid_sf import (
    NanoProcessor as CTAGEMDilepttValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_dileptt_valid_sf import (
    NanoProcessor as CTAGDilepttValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_Wctt_valid_sf import (
    NanoProcessor as CTAGWcTTValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_DY_valid_sf import (
    NanoProcessor as CTAGDYValidSFProcessor,
)

##QCD
from BTVNanoCommissioning.workflows.QCD_validation import (
    NanoProcessor as QCDValidProcessor,
)

## BTA - for SFs
from BTVNanoCommissioning.workflows.BTA_producer import (
    NanoProcessor as BTA_processor,
)
from BTVNanoCommissioning.workflows.BTA_ttbar_producer import (
    NanoProcessor as BTA_ttbar_processor,
)  # ttbar -kinFit

# from BTVNanoCommissioning.workflows.example import (
#     NanoProcessor as ExampleProcessor,
# )


# FIXME - make names more systematic?
workflows = {}
workflows["validation"] = ValidationProcessor

# TTBar
workflows["ttdilep_sf"] = TTdilepValidSFProcessor
workflows["ttsemilep_sf"] = TTsemilepValidSFProcessor

workflows["emctag_ttdilep_sf"] = CTAGEMDilepttValidSFProcessor
workflows["ctag_ttdilep_sf"] = partial(
    CTAGDilepttValidSFProcessor, selectionModifier="dilepttM"
)
workflows["ectag_ttdilep_sf"] = partial(
    CTAGDilepttValidSFProcessor, selectionModifier="dilepttE"
)
workflows["ctag_ttsemilep_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="semittM"
)
workflows["ectag_ttsemilep_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="semittE"
)

##QCD
workflows["QCD_sf"] = QCDValidProcessor

# W+c
workflows["ctag_Wc_sf"] = partial(CTAGWcTTValidSFProcessor, selectionModifier="WcM")
workflows["ectag_Wc_sf"] = partial(CTAGWcTTValidSFProcessor, selectionModifier="WcE")

# DY
workflows["ctag_DY_sf"] = partial(CTAGDYValidSFProcessor, selectionModifier="DYM")
workflows["ectag_DY_sf"] = partial(CTAGDYValidSFProcessor, selectionModifier="DYE")

# Tutorial
# workflows["example"] = ExampleProcessor
# BTA producers
workflows["BTA"] = BTA_processor
workflows["BTA_addPFMuons"] = partial(BTA_processor, addPFMuons=True)
workflows["BTA_addAllTracks"] = partial(BTA_processor, addAllTracks=True)
workflows["BTA_ttbar"] = BTA_ttbar_processor

__all__ = ["workflows"]
