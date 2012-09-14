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
    HLTPaths = cms.vstring(#'primaryVertexFilterPath',
                           'CSCTightHaloFilterPath',
                           'EcalDeadCellTriggerPrimitiveFilterPath',
#                           'noscrapingFilterPath',          
                           'hcalLaserEventFilterPath',
                           'HBHENoiseFilterPath',
                           'trackingFailureFilterPath',
                           'eeBadScFilterPath'
                           ),
    eventSetupPathsKey = cms.string(''),
    # how to deal with multiple triggers: True (OR) accept if ANY is true, False
    # (AND) accept if ALL are true
    andOr = cms.bool(False),
    throw = cms.bool(True)  # throw exception on unknown path names
    )


### To comment in case of MC samples

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.HLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.HLTFilter_2 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()

process.HLTFilter.HLTPaths  = ["HLT_Mu17_Mu8_v16",  "HLT_Mu17_TkMu8_v9" ]
process.HLTFilter_2.HLTPaths  = [ "HLT_Mu17_Mu8_v17", "HLT_Mu17_TkMu8_v10"]
process.HLTFilter.throw = cms.bool(False) 
process.HLTFilter_2.throw = cms.bool(False)


process.runFilter = cms.EDFilter(
        "FilterByRun",
        run = cms.int32(193834),
        selMode = cms.string("ge"),
        )



process.cleaningPath = cms.Path(
    process.badEventFilter * ~process.runFilter*process.HLTFilter 
    )


process.cleaningPath_2 = cms.Path(
    process.badEventFilter * process.runFilter*process.HLTFilter_2
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
    SelectEvents = cms.vstring('cleaningPath','cleaningPath_2')
    )

process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.endPath = cms.EndPath(process.edmNtuplesOut)

