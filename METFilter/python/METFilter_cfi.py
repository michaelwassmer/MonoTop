import FWCore.ParameterSet.Config as cms

METFilter = cms.EDFilter(
    "METFilter",
    filterData = cms.InputTag("TriggerResults","","PAT"),
    filterNames = cms.vstring([
                      "Flag_goodVertices",
                      "Flag_globalSuperTightHalo2016Filter",
                      "Flag_HBHENoiseFilter",
                      "Flag_HBHENoiseIsoFilter",
                      "Flag_EcalDeadCellTriggerPrimitiveFilter",
                      "Flag_BadPFMuonFilter",
                  ])
)
