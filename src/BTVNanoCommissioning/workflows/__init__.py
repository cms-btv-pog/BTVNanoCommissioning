 
# Validation
from BTVNanoCommissioning.workflows.validation import NanoProcessor as ValidationProcessor

# TTbar
from BTVNanoCommissioning.workflows.ttbar_validation import NanoProcessor as TTbarValidProcessor
from BTVNanoCommissioning.workflows.ttsemilep_valid_jec import NanoProcessor as TTsemilepValidJECProcessor
from BTVNanoCommissioning.workflows.ttdilep_valid_jec import NanoProcessor as TTdilepValidJECProcessor
from BTVNanoCommissioning.workflows.ttsemilep_valid_sf import NanoProcessor as TTsemilepValidSFProcessor
from BTVNanoCommissioning.workflows.ttdilep_valid_sf import NanoProcessor as TTdilepValidSFProcessor

# C-tag
from BTVNanoCommissioning.workflows.ctag_emdileptt_valid_sf import NanoProcessor as CTAGEMDilepttValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_dileptt_valid_sf import NanoProcessor as CTAGDilepttValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_eDY_valid_sf import NanoProcessor as CTAGeDYValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_eWc_valid_sf import NanoProcessor as CTAGeWcValidSFProcessor 
from BTVNanoCommissioning.workflows.ctag_Wc_valid_sf import NanoProcessor as CTAGWcValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_DY_valid_sf import NanoProcessor as CTAGDYValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_edileptt_valid_sf import NanoProcessor as CTAGEDilepttValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_ettsemilep_valid_sf import NanoProcessor as CTAGETTSemilepValidSFProcessor
from BTVNanoCommissioning.workflows.ctag_valid_jec import NanoProcessor as CTAGValidJECProcessor
from BTVNanoCommissioning.workflows.ctag_semileptt_valid_sf import NanoProcessor as CTAGSemilepttValidSFProcessor

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

# W+c 
workflows["ctag_Wc_sf"] = CTAGWcValidSFProcessor
workflows["ectag_Wc_sf"] = CTAGeWcValidSFProcessor

# DY 
workflows["ctag_DY_sf"] = CTAGDYValidSFProcessor
workflows["ectag_DY_sf"] = CTAGeDYValidSFProcessor

# JECs
workflows["ctag_jec"] = CTAGValidJECProcessor
workflows["dilep_jec"] = TTdilepValidJECProcessor
workflows["semilep_jec"] = TTsemilepValidJECProcessor




__all__ = ["workflows"]
# __all__ = ["workflows", "taggers", "DYStudiesProcessor"]
