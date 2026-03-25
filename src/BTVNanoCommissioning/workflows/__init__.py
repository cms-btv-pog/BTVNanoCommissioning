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
from BTVNanoCommissioning.workflows.sf_ttdilep_kin import (
    NanoProcessor as TTdilepKinSFProcessor,
)
from BTVNanoCommissioning.workflows.sf_ttsemilep_tnp import (
    NanoProcessor as TTsemilepTnPSFProcessor,
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
from BTVNanoCommissioning.workflows.QCD_soft_mu_validation import (
    NanoProcessor as QCDsmuValidProcessor,
)

## BTA - for SFs
from BTVNanoCommissioning.workflows.BTA_producer import (
    NanoProcessor as BTA_processor,
)
from BTVNanoCommissioning.workflows.BTA_ttbar_producer import (
    NanoProcessor as BTA_ttbar_processor,
)  # ttbar -kinFit

from BTVNanoCommissioning.workflows.example import (
    NanoProcessor as ExampleProcessor,
)

## SFl DY+Jets
from BTVNanoCommissioning.workflows.DY_sfl import (
    NanoProcessor as NegTagDYValidSFProcessor,
)

# FIXME - make names more systematic?
workflows = {}
workflows["validation"] = ValidationProcessor

# TTBar
workflows["ttdilep_sf"] = TTdilepValidSFProcessor
workflows["ttsemilep_sf"] = partial(
    TTsemilepValidSFProcessor, selectionModifier="tt_semilep"
)
workflows["sf_ttdilep_kin"] = TTdilepKinSFProcessor

workflows["c_ttsemilep_sf"] = partial(
    TTsemilepValidSFProcessor, selectionModifier="c_tt_semilep"
)

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
workflows["ctag_ttsemilep_noMuVeto_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="semittM_noMuVeto"
)
workflows["ectag_ttsemilep_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="semittE"
)

workflows["sf_ttsemilep_tnp"] = TTsemilepTnPSFProcessor

##QCD
workflows["QCD_sf"] = QCDValidProcessor
workflows["QCD_smu_sf"] = QCDsmuValidProcessor

# W+c
workflows["ctag_Wc_sf"] = partial(CTAGWcTTValidSFProcessor, selectionModifier="WcM")
workflows["ctag_Wc_noMuVeto_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="WcM_noMuVeto"
)
workflows["ectag_Wc_sf"] = partial(CTAGWcTTValidSFProcessor, selectionModifier="WcE")
workflows["ctag_Wc_WP_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="cutbased_WcM"
)
workflows["ectag_Wc_WP_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="cutbased_WcE"
)

# DY
workflows["ctag_DY_sf"] = partial(CTAGDYValidSFProcessor, selectionModifier="DYM")
workflows["ectag_DY_sf"] = partial(CTAGDYValidSFProcessor, selectionModifier="DYE")

# DY light jet SF with negative tag method
workflows["DY_sfl"] = partial(NegTagDYValidSFProcessor, selectionModifier="DYM")
workflows["eDY_sfl"] = partial(NegTagDYValidSFProcessor, selectionModifier="DYE")


# Tutorial
workflows["example"] = ExampleProcessor

# BTA producers
workflows["BTA"] = BTA_processor
workflows["BTA_addPFMuons"] = partial(BTA_processor, addPFMuons=True)
workflows["BTA_addAllTracks"] = partial(BTA_processor, addAllTracks=True)
workflows["BTA_ttbar"] = BTA_ttbar_processor

__all__ = ["workflows"]
