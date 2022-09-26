import FWCore.ParameterSet.Config as cms

TriggerFilter = cms.EDFilter(
    "FlagFilter",
    filterData = cms.InputTag("TriggerResults","","HLT"),
    filterNames = cms.vstring([
                      "HLT_IsoMu.*",
                      "HLT_PFMETNoMu.*"
                  ]),
    taggingMode = cms.bool(True),
    OR_mode = cms.bool(True)
)
