import FWCore.ParameterSet.Config as cms

TriggerFilter = cms.EDFilter(
    "FlagFilter",
    filterData = cms.InputTag("TriggerResults","","HLT"),
    filterNames = cms.vstring([
                      # isolated single muon triggers
                      "HLT_Iso(Mu|TkMu)2[0-9]_v([0-9]|[1-9][0-9])$",
                      # non-isolated single muon triggers
                      "HLT_(Mu|TkMu|OldMu)([4-9][0-9]|[1-9][0-9][0-9])_v([0-9]|[1-9][0-9])$",
                      # pfmet triggers
                      "HLT_PFMET.*",
                      # isolated single electron triggers
                      "HLT_Ele[2-9][0-9]_WPTight_Gsf_v([0-9]|[1-9][0-9])$",
                      # non-isolated single electron triggers
                      "HLT_Ele([4-9][0-9]|[1-9][0-9][0-9])_CaloIdVT_GsfTrkIdT_v([0-9]|[1-9][0-9])$",
                      # high pt photon triggers
                      "HLT_Photon[1-9][0-9][0-9]_v([0-9]|[1-9][0-9])$",
                      # pfht triggers
                      "HLT_PFHT([1-9][0-9][0-9]|[1-9][0-9][0-9][0-9])_v([0-9]|[1-9][0-9])$"
                  ]),
    taggingMode = cms.bool(True),
    OR_mode = cms.bool(True)
)
