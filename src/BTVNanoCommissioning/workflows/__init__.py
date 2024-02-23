from functools import partial

# Validation
from BTVNanoCommissioning.workflows.validation import (
    NanoProcessor as ValidationProcessor,
)

# TTbar
from BTVNanoCommissioning.workflows.ttbar_validation import (
    NanoProcessor as TTbarValidProcessor,
)
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
from BTVNanoCommissioning.workflows.ctag_eDY_valid_sf import (
    NanoProcessor as CTAGeDYValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_eWc_valid_sf import (
    NanoProcessor as CTAGeWcValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_Wc_valid_sf import (
    NanoProcessor as CTAGWcValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_DY_valid_sf import (
    NanoProcessor as CTAGDYValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_edileptt_valid_sf import (
    NanoProcessor as CTAGEDilepttValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_ettsemilep_valid_sf import (
    NanoProcessor as CTAGETTSemilepValidSFProcessor,
)
from BTVNanoCommissioning.workflows.ctag_semileptt_valid_sf import (
    NanoProcessor as CTAGSemilepttValidSFProcessor,
)
from BTVNanoCommissioning.workflows.BTA_producer import (
    NanoProcessor as BTA_processor,
)
from BTVNanoCommissioning.workflows.BTA_ttbar_producer import (
    NanoProcessor as BTA_ttbar_processor,
)

# from BTVNanoCommissioning.workflows.example import (
#     NanoProcessor as ExampleProcessor,
# )
##QCD
from BTVNanoCommissioning.workflows.QCD_validation import (
    NanoProcessor as QCDValidProcessor,
)

# FIXME - make names more systematic?
workflows = {}
workflows["validation"] = ValidationProcessor

# TTBar
workflows["ttcom"] = TTbarValidProcessor
workflows["ttdilep_sf"] = TTdilepValidSFProcessor
workflows["ttsemilep_sf"] = TTsemilepValidSFProcessor
workflows["emctag_ttdilep_sf"] = CTAGEMDilepttValidSFProcessor
workflows["ctag_ttdilep_sf"] = CTAGDilepttValidSFProcessor
workflows["ectag_ttdilep_sf"] = CTAGEDilepttValidSFProcessor
workflows["ctag_ttsemilep_sf"] = CTAGSemilepttValidSFProcessor
workflows["ectag_ttsemilep_sf"] = CTAGETTSemilepValidSFProcessor

##QCD
workflows["QCD"] = QCDValidProcessor

# W+c
workflows["ctag_Wc_sf"] = CTAGWcValidSFProcessor
workflows["ectag_Wc_sf"] = CTAGeWcValidSFProcessor

# DY
workflows["ctag_DY_sf"] = CTAGDYValidSFProcessor
workflows["ectag_DY_sf"] = CTAGeDYValidSFProcessor

# Tutorial
# workflows["example"] = ExampleProcessor
# BTA producers
workflows["BTA"] = BTA_processor
workflows["BTA_addPFMuons"] = partial(BTA_processor, addPFMuons=True)
workflows["BTA_addAllTracks"] = partial(BTA_processor, addAllTracks=True)
workflows["BTA_ttbar"] = BTA_ttbar_processor

__all__ = ["workflows"]
