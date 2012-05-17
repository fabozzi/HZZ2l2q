import FWCore.ParameterSet.Config as cms

process = cms.Process("SUBSK")

#### Simple cfg for subskim

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring(
    'file:h2l2qSkimData.root'
)

process.out = cms.OutputModule("PoolOutputModule",
                 fileName = cms.untracked.string('h2l2qSubSkimData.root'),
                 outputCommands =  cms.untracked.vstring(
                  'drop *_*_*_*',
                  ),
)

process.out.dropMetaData = cms.untracked.string("DROPPED")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patTriggerEventContent, patTriggerStandAloneEventContent

process.out.outputCommands += patTriggerEventContent
#process.out.outputCommands += patTriggerStandAloneEventContent

process.out.outputCommands.extend([
    'keep *_userDataSelectedElectrons_*_PAT',
    'keep *_userDataSelectedMuons_*_PAT',
    'keep *_cleanPatJetsNoPUIsoLept_*_PAT',
    # rho variables
    'keep *_*_rho_PAT',
    'keep *_kt6PFJetsCentralNeutral_rho_*',
    # ll, jj, lljj candidates
    'keep *_zee_*_PAT',
    'keep *_zmm_*_PAT',
    'keep *_zem_*_PAT',
    'keep *_zjj_*_PAT',
    'keep *_hzzeejj_*_PAT',
    'keep *_hzzmmjj_*_PAT',
    'keep *_hzzemjj_*_PAT',
    ####
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    # additional collections from AOD   
    'keep *_generalTracks_*_*',
    'keep *_electronGsfTracks_*_*',
    'keep *_muons_*_*',
    'keep recoPFCandidates_particleFlow_*_*',
    # gen Info
    'keep PileupSummaryInfos_*_*_*',
    ###### MET products
    'keep *_patMETs*_*_*',
#    'keep *_patType1CorrectedPFMet_*_*', # NOT included for the moment
    ### for HLT selection
    'keep edmTriggerResults_TriggerResults_*_HLT'])

# additional products for event cleaning
process.out.outputCommands.extend([
    'keep *_TriggerResults_*_PAT',
    ])

process.outPath = cms.EndPath(process.out)

