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
)

process.PUInfoNtuple = cms.EDProducer(
    "GenPUNtupleDump"
)

process.edmNtuplesOut.fileName = cms.untracked.string('h2l2b_ntuple.root')
process.edmNtuplesOut.outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_eventVtxInfoNtuple_*_*",
    "keep *_PUInfoNtuple_*_*",
    "keep *_metInfoProducer_*_*",
    "keep *_kt6PFJets_rho_PAT",
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
#    "keep *_Higgsemu2bEdmNtuple_*_*"
)
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.analysisPath = cms.Path(
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple
#    process.Higgsemu2bEdmNtuple
)

process.endPath = cms.EndPath(process.edmNtuplesOut)


process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
