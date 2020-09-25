import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register(
    "isData",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is it data or MC?",
)
options.parseArguments()

if options.maxEvents is -1:  # maxEvents is set in VarParsing class by default to -1
    options.maxEvents = 1001  # reset for testing

if not options.inputFiles:
    options.inputFiles = [
        "root://xrootd-cms.infn.it///store/mc/RunIIAutumn18MiniAOD/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/40000/E26C264E-6FAE-BA4B-B729-1FDD04B0B2AC.root"
    ]

process = cms.Process("EventCount")

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(0)
#)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# messages
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.load("MonoTop.EventCounter.EventCounter_cfi")

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.count = cms.EndPath(process.EventCounter)
# print(process.dumpPython())
