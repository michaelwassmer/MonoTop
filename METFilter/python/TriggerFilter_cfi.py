import FWCore.ParameterSet.Config as cms

TriggerFilter = cms.EDFilter(
    "FlagFilter",
    filterData = cms.InputTag("TriggerResults","","HLT"),
    filterNames = cms.vstring([
                      "HLT_IsoMu2[0-9]_v([0-9]|[1-9][0-9])$",
                      "HLT_PFMET.*"
                  ]),
    taggingMode = cms.bool(True),
    OR_mode = cms.bool(True)
)
