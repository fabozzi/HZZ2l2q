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
    'file:h2l2bData_FJNew_9_1_F2E.root',
    'file:h2l2bData_FJNew_99_1_2pl.root'
)

process.HLTPassInfo = cms.EDProducer(
    "HLTPassInfoProducer",
    triggerEvent = cms.InputTag("patTriggerEvent"),
    # here the 1st run with a new trigger table
    runLimits = cms.vint32(160410,#5e32
                           165121,#1e33
                           167039 #1.4e33
                           ),
    # here insert the HLT path (without _v[n] suffix) you want to check
    # Summer11 MC path
    triggerNamesSingleMu_MC = cms.vstring(),
    triggerNamesDoubleMu_MC = cms.vstring(),
    triggerNamesSingleEl_MC = cms.vstring(),
    triggerNamesDoubleEl_MC = cms.vstring(),
    # Data: here all the paths making the PDs are listed
    # 5e32 paths
    triggerNamesSingleMu_5e32 = cms.vstring('HLT_IsoMu12',
                                            'HLT_IsoMu15',
                                            'HLT_IsoMu17',
                                            'HLT_IsoMu24',
                                            'HLT_IsoMu30',
                                            'HLT_L1SingleMu10',
                                            'HLT_L1SingleMu20',
                                            'HLT_L2Mu10',
                                            'HLT_L2Mu20',
                                            'HLT_Mu12',
                                            'HLT_Mu15',
                                            'HLT_Mu20',
                                            'HLT_Mu24',
                                            'HLT_Mu30',
                                            'HLT_Mu3',
                                            'HLT_Mu5',
                                            'HLT_Mu8'
                                            ),
    triggerNamesDoubleMu_5e32 = cms.vstring('HLT_DoubleMu3',
                                            'HLT_DoubleMu4_Acoplanarity03',
                                            'HLT_DoubleMu6',
                                            'HLT_DoubleMu7',
                                            'HLT_L1DoubleMu0',
                                            'HLT_L2DoubleMu0',
                                            'HLT_L2DoubleMu23_NoVertex',
                                            'HLT_Mu8_Jet40',
                                            'HLT_TripleMu5'
                                            ),
    triggerNamesSingleEl_5e32 = cms.vstring('HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT',
                                            'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT',
                                            'HLT_Ele45_CaloIdVT_TrkIdT',
                                            'HLT_Ele90_NoSpikeFilter'
                                            ),
    triggerNamesDoubleEl_5e32 = cms.vstring('HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL',
                                            'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL',
                                            'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30',
                                            'HLT_Ele32_CaloIdL_CaloIsoVL_SC17',
                                            'HLT_Ele8_CaloIdL_CaloIsoVL_Jet40',
                                            'HLT_Ele8_CaloIdL_CaloIsoVL',
                                            'HLT_Ele8_CaloIdL_TrkIdVL',
                                            'HLT_Ele8',
                                            'HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL',
                                            'HLT_TripleEle10_CaloIdL_TrkIdVL'
                                            ),
    # 1e33 paths
    triggerNamesSingleMu_1e33 = cms.vstring('HLT_IsoMu12',
                                            'HLT_IsoMu15',
                                            'HLT_IsoMu17',
                                            'HLT_IsoMu24',
                                            'HLT_IsoMu30',
                                            'HLT_L1SingleMu10',
                                            'HLT_L1SingleMu20',
                                            'HLT_L2Mu10',
                                            'HLT_L2Mu20',
                                            'HLT_Mu100',
                                            'HLT_Mu12',
                                            'HLT_Mu15',
                                            'HLT_Mu20',
                                            'HLT_Mu24',
                                            'HLT_Mu30',
                                            'HLT_Mu3',
                                            'HLT_Mu40',
                                            'HLT_Mu5',
                                            'HLT_Mu8'
                                            ),
    triggerNamesDoubleMu_1e33 = cms.vstring('HLT_DoubleMu3',
                                            'HLT_DoubleMu45',
                                            'HLT_DoubleMu4_Acoplanarity03',
                                            'HLT_DoubleMu5_Acoplanarity03',
                                            'HLT_DoubleMu6',
                                            'HLT_DoubleMu7',
                                            'HLT_L1DoubleMu0',
                                            'HLT_L2DoubleMu0',
                                            'HLT_L2DoubleMu23_NoVertex',
                                            'HLT_Mu13_Mu8',
                                            'HLT_Mu17_Mu8',
                                            'HLT_Mu8_Jet40',
                                            'HLT_TripleMu5'
                                            ),
    triggerNamesSingleEl_1e33 = cms.vstring('HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                            'HLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                            'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT',
                                            'HLT_Ele42_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                            'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT',
                                            'HLT_Ele52_CaloIdVT_TrkIdT',
                                            'HLT_Ele65_CaloIdVT_TrkIdT'
                                            ),
    triggerNamesDoubleEl_1e33 = cms.vstring('HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10',
                                            'HLT_DoubleEle45_CaloIdL',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFT',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL',
                                            'HLT_Ele17_CaloIdL_CaloIsoVL',
                                            'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL',
                                            'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30',
                                            'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30',
                                            'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17',
                                            'HLT_Ele8_CaloIdL_CaloIsoVL_Jet40',
                                            'HLT_Ele8_CaloIdL_CaloIsoVL',
                                            'HLT_Ele8_CaloIdL_TrkIdVL',
                                            'HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL',
                                            'HLT_Ele8',
                                            'HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL',
                                            'HLT_TripleEle10_CaloIdL_TrkIdVL'
                                            ),
##### 1.4e33 paths
    triggerNamesSingleMu_1p4e33 = cms.vstring('HLT_IsoMu12',
                                              'HLT_IsoMu15',
                                              'HLT_IsoMu17_eta2p1',
                                              'HLT_IsoMu17',
                                              'HLT_IsoMu20_eta2p1',
                                              'HLT_IsoMu24',
                                              'HLT_IsoMu30',
                                              'HLT_L1SingleMu10',
                                              'HLT_L1SingleMu20',
                                              'HLT_L2Mu10',
                                              'HLT_L2Mu20',
                                              'HLT_Mu100',
                                              'HLT_Mu12',
                                              'HLT_Mu15',
                                              'HLT_Mu20',
                                              'HLT_Mu24',
                                              'HLT_Mu30',
                                              'HLT_Mu3',
                                              'HLT_Mu40',
                                              'HLT_Mu5',
                                              'HLT_Mu8'
                                              ),
    triggerNamesDoubleMu_1p4e33 = cms.vstring('HLT_DoubleMu3',
                                              'HLT_DoubleMu45',
                                              'HLT_DoubleMu4_Acoplanarity03',
                                              'HLT_DoubleMu5_Acoplanarity03',
                                              'HLT_DoubleMu5_IsoMu5',
                                              'HLT_DoubleMu6',
                                              'HLT_DoubleMu7',
                                              'HLT_L1DoubleMu0',
                                              'HLT_L2DoubleMu0',
                                              'HLT_L2DoubleMu23_NoVertex',
                                              'HLT_L2DoubleMu30_NoVertex',
                                              'HLT_Mu13_Mu8',
                                              'HLT_Mu17_Mu8',
                                              'HLT_Mu8_Jet40',
                                              'HLT_TripleMu5'
                                              ),
    triggerNamesSingleEl_1p4e33 = cms.vstring('HLT_Ele100_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                              'HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                              'HLT_Ele27_WP80_PFMT50',
                                              'HLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                              'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT',
                                              'HLT_Ele32_WP70_PFMT50',
                                              'HLT_Ele52_CaloIdVT_TrkIdT',
                                              'HLT_Ele65_CaloIdVT_TrkIdT'
                                              ),
    triggerNamesDoubleEl_1p4e33 = cms.vstring('HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_CaloIdT_TrkIdVL',
                                              'HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL',
                                              'HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFT',
                                              'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL',
                                              'HLT_Ele17_CaloIdL_CaloIsoVL',
                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                              'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30',
                                              'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30',
                                              'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17',
                                              'HLT_Ele8_CaloIdL_CaloIsoVL_Jet40',
                                              'HLT_Ele8_CaloIdL_CaloIsoVL',
                                              'HLT_Ele8_CaloIdL_TrkIdVL',
                                              'HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL',
                                              'HLT_Ele8',
                                              'HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL',
                                              'HLT_TripleEle10_CaloIdL_TrkIdVL'
                                              )
    )


process.edmNtuplesOut.fileName = cms.untracked.string('h2l2b_ntuple.root')
process.edmNtuplesOut.outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_HLTPassInfo_*_*",
    "keep *_eventVtxInfoNtuple_*_*",
    "keep *_metInfoProducer_*_*",
    "keep *_kt6PFJets_rho_PAT",
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
#    "keep *_Higgsemu2bEdmNtuple_*_*"
)
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.analysisPath = cms.Path(
    process.HLTPassInfo+
    process.eventVtxInfoNtuple+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple
#    process.Higgsemu2bEdmNtuple
)

process.endPath = cms.EndPath(process.edmNtuplesOut)


process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
