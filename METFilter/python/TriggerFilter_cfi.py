import FWCore.ParameterSet.Config as cms

TriggerFilter = cms.EDFilter(
    "FlagFilter",
    filterData = cms.InputTag("TriggerResults","","HLT"),
    filterNames = cms.vstring([
                      "HLT_Iso(Mu|TkMu)2[0-9]_v([0-9]|[1-9][0-9])$",
                      "HLT_PFMET.*",
                      "HLT_Ele([0-9]|[1-9][0-9])_WPTight_Gsf_v([0-9]|[1-9][0-9])$",
                      "HLT_PFHT([1-9][0-9][0-9]|[1-9][0-9][0-9][0-9])_v([0-9]|[1-9][0-9])$"
                  ]),
    taggingMode = cms.bool(True),
    OR_mode = cms.bool(True)
)
