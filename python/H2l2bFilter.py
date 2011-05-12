import FWCore.ParameterSet.Config as cms
process = cms.Process("H2l2bFilter")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource")
process.source.fileNames=cms.untracked.vstring(
    '/store/group/local/SMHiggsToZZTo2L2Q_M-350/SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6/Skim_H2l2b_414_v2/80892fc9180aeedf1cd8e77773f43707/h2l2bData_3_1_pfD.root',
)

process.myfilter  = cms.EDFilter("H2l2bFilter",
    hzzeejjTag = cms.InputTag("hzzeejj:h"),
    hzzmmjjTag = cms.InputTag("hzzmmjj:h"),
    zLepPtCut = cms.double(20),
    zJetPtCut = cms.double(30),
    zLepMassCut = cms.double(10),
    zJetMassCut = cms.double(15),
    zLepRelIsoCut = cms.double(0.15),
    zMudBCut =  cms.double(0.02),                           
    zllPtCut = cms.double(0),
    metCut = cms.double(10000),
    jjdrCut = cms.double(10)
)

process.filterPath = cms.Path(
    process.myfilter
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("mySkimTight.root"),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring("filterPath")),
)

process.end = cms.EndPath(process.out)
