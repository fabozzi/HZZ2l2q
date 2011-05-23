import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2bedmNtuplesLD_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service("TFileService", fileName = cms.string('h2l2b400_histo.root') )

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring(

    'file:/data3/scratch/users/fabozzi/Higgs/400/skim/h2l2bData_1_1_5SF.root',
    'file:/data3/scratch/users/fabozzi/Higgs/400/skim/h2l2bData_2_1_42j.root',

)

process.edmNtuplesOut.fileName = cms.untracked.string('h2l2b_ntuple_testLD.root')
process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

# select lumis interactively from a json
#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#import FWCore.ParameterSet.Types as CfgTypes
#myLumis = LumiList.LumiList(filename = 'Cert_160404-163757_7TeV_PromptReco_Collisions11_JSON.txt').getCMSSWString().split(',')
#process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#process.source.lumisToProcess.extend(myLumis)


process.hzzeejjLD = cms.EDProducer("Higgs2l2bUserDataLD",
                                     higgs = cms.InputTag("hzzeejj:h")
                                     )

process.hzzmmjjLD = cms.EDProducer("Higgs2l2bUserDataLD",
                                     higgs = cms.InputTag("hzzmmjj:h")
                                     )

process.hzzemjjLD = cms.EDProducer("Higgs2l2bUserDataLD",
                                     higgs = cms.InputTag("hzzemjj:h")
                                     )

process.analysisPath = cms.Path(
    process.eventVtxInfoNtuple+
    process.hzzeejjLD +
    process.hzzmmjjLD +
    process.hzzemjjLD +
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple
#    process.Higgsemu2bEdmNtuple
)

process.endPath = cms.EndPath(process.edmNtuplesOut)


process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
