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

process.edmNtuplesOut.fileName = cms.untracked.string('h2l2b_ntuple.root')
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')


process.analysisPath = cms.Path(
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
