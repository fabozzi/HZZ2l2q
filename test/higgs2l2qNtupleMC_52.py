import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

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

process.PUInfoNtuple = cms.EDProducer(
    "GenPUNtupleDump",
    isData = cms.bool(False)
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
    "keep *_eventVtxInfoNtuple_*_*",
    "keep *_PUInfoNtuple_*_*",
    "keep *_rhoDumper_*_*",
    "keep *_metInfoProducer_*_*",
    "keep *_kt6PFJets_rho_PAT",
    "keep *_kt6PFJetsForIso_rho_*",
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
    "keep *_Higgsemu2bEdmNtuple_*_*",
    "keep *_jetinfos_*_*"
)
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.analysisPath = cms.Path(
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.rhoDumper+
    process.metInfoProducer+
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
