from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName = "THEREQUESTNAME"
config.General.workArea = "crab_skims"
config.General.transferLogs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "/nfs/dust/cms/user/mwassmer/MonoTop/Skimming/CMSSW_10_2_18/src/MonoTop/MonoTopSkim/test/skim.py"
config.JobType.outputFiles = ["Skim.root"]
# config.JobType.maxJobRuntimeMin = 2800
config.JobType.maxMemoryMB = 2000
config.JobType.pyCfgParams = [
    "isData=ISDATA",
    "maxEvents=999999999",
    "dataEra=DATAERA",
    "globalTag=GLOBALTAG",
]
config.JobType.sendPythonFolder = True

config.Data.inputDataset = "THEINPUTDATASET"
config.Data.inputDBS = "global"
config.Data.splitting = "Automatic"
config.Data.unitsPerJob = 720
config.Data.publication = False
config.Data.publishDBS = "phys03"
config.Data.outputDatasetTag = "KIT_Monotop_skims_2016_102X"

config.User.voGroup = "dcms"

config.Site.storageSite = "T2_DE_DESY"
# config.Site.blacklist = ['T2_US_*']
# config.Site.whitelist = ['T2_DE_DESY']
