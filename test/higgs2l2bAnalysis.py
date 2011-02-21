import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2bedmNtuples_cff")
#process.load("HiggsAnalysis.Higgs2l2b.zjjEdmNtuples_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.TFileService = cms.Service("TFileService", fileName = cms.string('h2l2b350GF.root') )

process.source = cms.Source("PoolSource")

#process.source.fileNames=cms.untracked.vstring('rfio:/castor/cern.ch/user/d/decosa/Higgs/h350/skim/h2l2b300GF_1_1_Ukh.root')
process.source.fileNames=cms.untracked.vstring(
    'file:/scratch2/users/fabozzi/higgs/Skim/h350/h2l2b300GF_1_1_Ukh.root',
    'file:/scratch2/users/fabozzi/higgs/Skim/h350/h2l2b300GF_2_1_z6i.root',
    'file:/scratch2/users/fabozzi/higgs/Skim/h350/h2l2b300GF_3_1_Ar8.root'
)

####
process.elChannelAnalysis = cms.EDAnalyzer("Higgs2l2bAnalysis",
    higgsTag = cms.InputTag("hzzeejj:h"),
)

process.elChannelSelection = cms.EDAnalyzer("H2l2bSelection",
    higgsTag = cms.InputTag("hzzeejj:h"),
    zLepMassCut = cms.double(10),
    zJetMassCut = cms.double(15),
    bTagCSVCut = cms.double(0.5),
    bTagJProbCut = cms.double(0.9),
    zllPtCut = cms.double(90),
    metCut = cms.double(35),
    jjdrCut = cms.double(1.7),
    hMassMinCut = cms.double(315),
    hMassMaxCut = cms.double(385),
    lumiNormalization = cms.double(0.00115),
    output_name = cms.string("h350GF")                                            
                                            
)


process.muChannelSelection = cms.EDAnalyzer("H2l2bSelection",
    higgsTag = cms.InputTag("hzzmmjj:h"),
    zLepMassCut = cms.double(10),
    zJetMassCut = cms.double(15),
    bTagCSVCut = cms.double(0.5),
    bTagJProbCut = cms.double(0.9),
    zllPtCut = cms.double(90),
    metCut = cms.double(35),
    jjdrCut = cms.double(1.7),
    hMassMinCut = cms.double(315),
    hMassMaxCut = cms.double(385),
    lumiNormalization = cms.double(0.00115),
    output_name = cms.string("h350GF")                                            
)

process.muChannelAnalysis = cms.EDAnalyzer("Higgs2l2bAnalysis",
    higgsTag = cms.InputTag("hzzmmjj:h"),
)

process.edmNtuplesOut.fileName = cms.untracked.string('H350GF_ntuples.root')

process.analysisPath = cms.Path(
    process.elChannelAnalysis+
    process.elChannelSelection+
    process.muChannelAnalysis+
    process.muChannelSelection+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple
    )

process.endPath = cms.EndPath(process.edmNtuplesOut)


process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
