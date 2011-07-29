import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2bedmNtuplesSF_cff")

process.load("CondCore.DBCommon.CondDBCommon_cfi")

#Data measurements from Summer11
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service("TFileService", fileName = cms.string('h2l2b450_histo.root') )

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring(
    'file:h2l2bData.root'
)

process.PUInfoNtuple = cms.EDProducer(
    "GenPUNtupleDump",
    isData = cms.bool(False)
)

process.HLTPassInfo = cms.EDProducer(
    "HLTPassInfoProducer",
    triggerEvent = cms.InputTag("patTriggerEvent"),
    # here the 1st run with a new trigger table
    # leave empty for MC
    runLimits = cms.vint32(),
    # here insert the HLT path (without _v[n] suffix) you want to check
    # Summer11 MC path
    triggerNamesSingleMu_MC = cms.vstring('HLT_IsoMu17'),
    triggerNamesDoubleMu_MC = cms.vstring('HLT_DoubleMu7'),
    triggerNamesSingleEl_MC = cms.vstring(),
    triggerNamesDoubleEl_MC = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL'),
    # Data: here all the paths making the PDs are listed
    # 5e32 paths
    triggerNamesSingleMu_5e32 = cms.vstring(),
    triggerNamesDoubleMu_5e32 = cms.vstring(),
    triggerNamesSingleEl_5e32 = cms.vstring(),
    triggerNamesDoubleEl_5e32 = cms.vstring(),
    # 1e33 paths
    triggerNamesSingleMu_1e33 = cms.vstring(),
    triggerNamesDoubleMu_1e33 = cms.vstring(),
    triggerNamesSingleEl_1e33 = cms.vstring(),
    triggerNamesDoubleEl_1e33 = cms.vstring(),
##### 1.4e33 paths
    triggerNamesSingleMu_1p4e33 = cms.vstring(),
    triggerNamesDoubleMu_1p4e33 = cms.vstring(),
    triggerNamesSingleEl_1p4e33 = cms.vstring(),
    triggerNamesDoubleEl_1p4e33 = cms.vstring(),
##### 2e33 paths
    triggerNamesSingleMu_2e33 = cms.vstring(),
    triggerNamesDoubleMu_2e33 = cms.vstring(),
    triggerNamesSingleEl_2e33 = cms.vstring(),
    triggerNamesDoubleEl_2e33 = cms.vstring()
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
#    "keep *_Higgsemu2bEdmNtuple_*_*"
)
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.hzzeejjSF = cms.EDProducer("Higgs2l2bUserDataSF",
                                     higgs = cms.InputTag("hzzeejj:h")
                                     )

process.hzzmmjjSF = cms.EDProducer("Higgs2l2bUserDataSF",
                                     higgs = cms.InputTag("hzzmmjj:h")
                                     )

process.hzzemjjSF = cms.EDProducer("Higgs2l2bUserDataSF",
                                     higgs = cms.InputTag("hzzemjj:h")
                                     )

process.analysisPath = cms.Path(
    process.HLTPassInfo+
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.hzzeejjSF +
    process.hzzmmjjSF +
    process.hzzemjjSF +
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple
#    process.Higgsemu2bEdmNtuple
)

process.endPath = cms.EndPath(process.edmNtuplesOut)


process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
