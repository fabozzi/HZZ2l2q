import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2bedmNtuples_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service("TFileService", fileName = cms.string('h2l2b450_histo.root') )

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring(
    'file:h2l2bData.root'
#    'file:h2l2bData_FJNew_9_1_F2E.root',
#    'file:h2l2bData_FJNew_99_1_2pl.root'
)

process.PUInfoNtuple = cms.EDProducer(
    "GenPUNtupleDump",
    isData = cms.bool(True)
)

process.HLTPassInfo = cms.EDProducer(
    "HLTPassInfoProducer",
    triggerEvent = cms.InputTag("patTriggerEvent"),
    # here the 1st run with a new trigger table
    runLimits = cms.vint32(160410,#5e32
                           165121,#1e33
                           167039, #1.4e33
                           170249, #2e33
                           173236, #3e33
                           178421  #5e33
                           ),
    # here insert the HLT path (without _v[n] suffix) you want to check
    # Summer11 MC path
    triggerNamesSingleMu_MC = cms.vstring(),
    triggerNamesDoubleMu_MC = cms.vstring(),
    triggerNamesSingleEl_MC = cms.vstring(),
    triggerNamesDoubleEl_MC = cms.vstring(),
    # Data: requested path(s) in the PD
    # 5e32 paths
    triggerNamesSingleMu_5e32 = cms.vstring('HLT_IsoMu17'),
    triggerNamesDoubleMu_5e32 = cms.vstring('HLT_DoubleMu7'),
    triggerNamesSingleEl_5e32 = cms.vstring(),
    triggerNamesDoubleEl_5e32 = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL'),
    # 1e33 paths
    triggerNamesSingleMu_1e33 = cms.vstring('HLT_IsoMu17'),
    triggerNamesDoubleMu_1e33 = cms.vstring('HLT_Mu13_Mu8'),
    triggerNamesSingleEl_1e33 = cms.vstring(),
    triggerNamesDoubleEl_1e33 = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL'),
##### 1.4e33 paths
    triggerNamesSingleMu_1p4e33 = cms.vstring('HLT_IsoMu17'),
    triggerNamesDoubleMu_1p4e33 = cms.vstring('HLT_Mu13_Mu8'),
    triggerNamesSingleEl_1p4e33 = cms.vstring(),
    triggerNamesDoubleEl_1p4e33 = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL'),
##### 2e33 paths
    triggerNamesSingleMu_2e33 = cms.vstring('HLT_IsoMu17'),
    triggerNamesDoubleMu_2e33 = cms.vstring('HLT_Mu13_Mu8'),
    triggerNamesSingleEl_2e33 = cms.vstring(),
    triggerNamesDoubleEl_2e33 = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL'),
##### 3e33 paths
    triggerNamesSingleMu_3e33 = cms.vstring('HLT_IsoMu20'),
    triggerNamesDoubleMu_3e33 = cms.vstring('HLT_Mu13_Mu8'),
    triggerNamesSingleEl_3e33 = cms.vstring(),
    triggerNamesDoubleEl_3e33 = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL'),
##### 5e33 paths
    triggerNamesSingleMu_5e33 = cms.vstring('HLT_IsoMu24_eta2p1'),
    triggerNamesDoubleMu_5e33 = cms.vstring('HLT_Mu17_Mu8'),
    triggerNamesSingleEl_5e33 = cms.vstring(),
    triggerNamesDoubleEl_5e33 = cms.vstring('HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL')
    )


process.edmNtuplesOut.fileName = cms.untracked.string('h2l2b_ntuple.root')
process.edmNtuplesOut.outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_HLTPassInfo_*_*",
    "keep *_eventVtxInfoNtuple_*_*",
    "keep *_PUInfoNtuple_*_*",
    "keep *_metInfoProducer_*_*",
    "keep *_kt6PFJets_rho_PAT",
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
    "keep *_Higgsemu2bEdmNtuple_*_*",
    "keep *_jetinfos_*_*"
)
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.analysisPath = cms.Path(
    process.HLTPassInfo+
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple+
    process.Higgsemu2bEdmNtuple+
    process.jetinfos
)

process.endPath = cms.EndPath(process.edmNtuplesOut)

process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
