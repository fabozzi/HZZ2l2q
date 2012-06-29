import FWCore.ParameterSet.Config as cms

process = cms.Process("NTPFILT")

#### Simple cfg to apply event filter to ntuple

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring('file:h2l2q_ntuple.root')

#### Event cleaning 
process.badEventFilter = cms.EDFilter(
    "HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","PAT"),
    HLTPaths = cms.vstring('primaryVertexFilterPath',
#                           'CSCTightHaloFilterPath',
#                           'EcalDeadCellTriggerPrimitiveFilterPath',
                           'EcalDeadCellBoundaryEnergyFilterPath',
                           'noscrapingFilterPath',          
#### Use hcalLaserFilterFromAODPath in 5_2 ##########
#                           'hcalLaserEventFilterPath',
                           'hcalLaserFilterFromAODPath',
                           'HBHENoiseFilterPath'#,
#                           'totalKinematicsFilterPath' #only for Madgraph MC
                           ),
    eventSetupPathsKey = cms.string(''),
    # how to deal with multiple triggers: True (OR) accept if ANY is true, False
    # (AND) accept if ALL are true
    andOr = cms.bool(False),
    throw = cms.bool(True)  # throw exception on unknown path names
    )

process.cleaningPath = cms.Path(
    process.badEventFilter
    )

process.edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('h2l2q_ntuple_clean.root'),
    outputCommands = cms.untracked.vstring(
    "keep *",
    "drop *_TriggerResults_*_PAT",
    "drop *_TriggerResults_*_NTPFILT",
    )
    )

process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('cleaningPath')
    )

process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.endPath = cms.EndPath(process.edmNtuplesOut)

