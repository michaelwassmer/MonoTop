import FWCore.ParameterSet.Config as cms

EventCounter = cms.EDAnalyzer(
    "EventCounter",
    lheInfo=cms.InputTag("externalLHEProducer"),
    genInfo=cms.InputTag("generator")
)
