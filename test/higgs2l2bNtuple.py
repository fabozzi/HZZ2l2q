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
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_10_1_0WE.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_11_1_fgJ.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_12_1_M8d.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_1_1_Lxw.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_2_1_yjl.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_3_1_G8v.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_4_1_LRg.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_5_1_sxk.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_6_1_seh.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_7_1_l3b.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_8_1_Ukd.root',
    'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Run2011A/Skim/DoubleMu/160404-161352_v2/h2l2bData_9_1_sbR.root',
)

process.edmNtuplesOut.fileName = cms.untracked.string('ntupleMu/h2l2b_ntuple.root')

process.analysisPath = cms.Path(
    process.eventVtxInfoNtuple+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple+
    process.Higgsemu2bEdmNtuple
)

process.endPath = cms.EndPath(process.edmNtuplesOut)


process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
