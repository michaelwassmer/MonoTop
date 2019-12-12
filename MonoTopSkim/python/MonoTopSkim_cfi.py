import FWCore.ParameterSet.Config as cms

MonoTopSkim = cms.EDFilter(
    "MonoTopSkim",
    era=cms.string("NA"),
    isData=cms.bool(False),
    electrons=cms.InputTag("slimmedElectrons"),
    muons=cms.InputTag("slimmedMuons"),
    photons=cms.InputTag("slimmedPhotons"),
    AK4jets=cms.InputTag("selectedPatJetsAK4PFPuppi"),
    AK8jets=cms.InputTag("selectedPatJetsAK8PFPuppi"),
    AK15jets=cms.InputTag("selectedPatJetsAK15PFPuppi"),
    vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho=cms.InputTag("fixedGridRhoFastjetAll"),
    met=cms.InputTag("slimmedMETs"),
    met_puppi=cms.InputTag("slimmedMETsPuppi"),
    muonPtMin=cms.double(10),
    muonEtaMax=cms.double(2.5),
    electronPtMin=cms.double(10),
    electronEtaMax=cms.double(2.5),
    AK4jetPtMin=cms.double(20),
    AK4jetEtaMax=cms.double(2.5),
    AK8jetPtMin=cms.double(170),
    AK8jetEtaMax=cms.double(2.5),
    minJetsAK4=cms.int32(1),
    minJetsAK8=cms.int32(1),
    minJetsAK15=cms.int32(1),
    maxJetsAK4=cms.int32(99),
    maxJetsAK8=cms.int32(4),
    maxJetsAK15=cms.int32(3),
    photonPtMin=cms.double(15),
    photonEtaMax=cms.double(2.5),
    metPtMin=cms.double(200.0)
)
