import FWCore.ParameterSet.Config as cms

METFilter = cms.EDFilter(
    "FlagFilter",
    filterData = cms.InputTag("TriggerResults","","PAT"),
    filterNames = cms.vstring([
                      "Flag_goodVertices",
                      "Flag_globalSuperTightHalo2016Filter",
                      "Flag_HBHENoiseFilter",
                      "Flag_HBHENoiseIsoFilter",
                      "Flag_EcalDeadCellTriggerPrimitiveFilter",
                      "Flag_BadPFMuonFilter",
                      "Flag_BadPFMuonDzFilter",
                      "Flag_eeBadScFilter",
                      #"Flag_ecalBadCalibFilter",
                  ]),
    taggingMode = cms.bool(False),
    OR_mode = cms.bool(False)
)
