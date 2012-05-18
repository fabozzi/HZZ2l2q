import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

applyFilter = True

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2qedmNtuples_52_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring(
    'file:h2l2qSkimData.root'
)


#### Event cleaning 
process.badEventFilter = cms.EDFilter(
    "HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","PAT"),
    HLTPaths = cms.vstring('primaryVertexFilterPath',
#                           'CSCTightHaloFilterPath',
#                           'EcalDeadCellTriggerPrimitiveFilterPath',
                           'EcalDeadCellBoundaryEnergyFilterPath',
                           'noscrapingFilterPath',          
                           'hcalLaserEventFilterPath',
                           'HBHENoiseFilterPath'#,
#                           'totalKinematicsFilterPath' #only for Madgraph MC
                           ),
    eventSetupPathsKey = cms.string(''),
    # how to deal with multiple triggers: True (OR) accept if ANY is true, False
    # (AND) accept if ALL are true
    andOr = cms.bool(False),
    throw = cms.bool(True)  # throw exception on unknown path names
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
    # MC path
    triggerNamesSingleMu_MC = cms.vstring(),
    triggerNamesDoubleMu_MC = cms.vstring(),
    triggerNamesSingleEl_MC = cms.vstring(),
    triggerNamesDoubleEl_MC = cms.vstring(),
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
    triggerNamesDoubleEl_2e33 = cms.vstring(),
##### 3e33 paths
    triggerNamesSingleMu_3e33 = cms.vstring(),
    triggerNamesDoubleMu_3e33 = cms.vstring(),
    triggerNamesSingleEl_3e33 = cms.vstring(),
    triggerNamesDoubleEl_3e33 = cms.vstring(),
##### 5e33 paths
    triggerNamesSingleMu_5e33 = cms.vstring(),
    triggerNamesDoubleMu_5e33 = cms.vstring(),
    triggerNamesSingleEl_5e33 = cms.vstring(),
    triggerNamesDoubleEl_5e33 = cms.vstring()
    )

# Event rho dumper
process.rhoDumper = cms.EDProducer("EventRhoDumper",
                                    rho = cms.InputTag("kt6PFJets:rho"),
                                    restrictedRho = cms.InputTag("kt6PFJetsForIso:rho")
                                    )

# Met variables producer
process.metInfoProducer = cms.EDProducer("MetVariablesProducer",
                                    metTag = cms.InputTag("patMETsAK5"),
                                    t1CorrMetTag = cms.InputTag("patType1CorrectedPFMetAK5")
                                    )

process.edmNtuplesOut.fileName = cms.untracked.string('h2l2q_ntuple.root')
process.edmNtuplesOut.outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_HLTPassInfo_*_*",
    "keep *_eventVtxInfoNtuple_*_*",
    "keep *_PUInfoNtuple_*_*",
    "keep *_rhoDumper_*_*",
    "keep *_metInfoProducer_*_*",
# keep rho recommended for muoniso
    'keep *_kt6PFJetsCentralNeutral_rho_*',
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
    "keep *_Higgsemu2bEdmNtuple_*_*",
    "keep *_jetinfos_*_*"
)

if applyFilter:
    process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('analysisPath')
        )

process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.edmNtuplesOut.outputCommands.extend([
    'keep edmTriggerResults_TriggerResults_*_HLT',
    'keep *_TriggerResults_*_PAT',
    ])

process.analysisPath = cms.Path()

if applyFilter:
    process.analysisPath += process.badEventFilter

#process.analysisPath = cms.Path(
#    process.badEventFilter +
process.analysisSequence = cms.Sequence(
    process.HLTPassInfo+
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.rhoDumper+
    process.metInfoProducer+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple+
    process.Higgsemu2bEdmNtuple+
    process.jetinfos
)

process.analysisPath += process.analysisSequence

process.endPath = cms.EndPath(process.edmNtuplesOut)

process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
