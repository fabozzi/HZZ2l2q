import FWCore.ParameterSet.Config as cms

# Setup PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *




#removeMCMatching(['Muons'])

#runOnData(['All'])
# Load some generic cffs  
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

# Specify the Global Tag
process.GlobalTag.globaltag = 'GR_R_38X_V15::All'

# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# Source file : To be run on a Full RECO sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/scratch2/users/decosa/data/Run2010B/063A44C5-1BEB-DF11-87A8-E0CB4E29C4B8.root'

    )
)



# Output Module : Hopefully we keep all we need
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('h2l2bData.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("filterPath")
    ),
    outputCommands =  cms.untracked.vstring(
        'drop *_*_*_*',
        'keep *_selectedPatElectrons_*_PAT',
        'keep *_selectedPatMuons_*_PAT',
        'keep *_cleanPatJets_*_PAT',
        'keep *_zee_*_PAT',
        'keep *_zmm_*_PAT',
        'keep *_zjj_*_PAT',
        'keep *_hzzeejj_*_PAT',
        'keep *_hzzmmjj_*_PAT',
        #'keep *_TriggerResults*_*_HLT',
        #'keep *_hltTriggerSummaryAOD_*_HLT',
        #'keep *_TriggerResults*_*_REDIGI*',
        #'keep *_hltTriggerSummaryAOD_*_REDIGI*',        
        'keep *_offlinePrimaryVertices_*_*'
    ),

    #verbose = cms.untracked.bool(True)
)

# Modules for the Cut-based Electron ID in the VBTF prescription
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )

process.eidSequence = cms.Sequence(
    process.eidVBTFRel95 +
    process.eidVBTFRel80 +
    process.eidVBTFCom95 +
    process.eidVBTFCom80 
)

# Muon Selection
process.selectedPatMuons.cut = ( 
    "pt > 10 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
    "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
    "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
    "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
    "dB < 0.2 && "                                                                     +
    #"trackIso + caloIso < 0.15 * pt && "                                               +
    "numberOfMatches > 1 && abs(eta) < 2.4" 
)

# Electron Selection
process.patElectrons.electronIDSources = cms.PSet(
    eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
    eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
    eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
    eidVBTFCom80 = cms.InputTag("eidVBTFCom80")
)

process.selectedPatElectrons.cut = ( 
    "pt > 10.0 && abs(eta) < 2.5 &&"                               +
    "(isEE || isEB) && !isEBEEGap &&"                              +
    "electronID('eidVBTFCom95') == 7" 
)


# Switch to using PFJets 
switchJetCollection(process,cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])),
    doType1MET   = True,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = True
)


# Switch to using PFMET 
switchToPFMET(
    process, 
    cms.InputTag('pfMet'), 
    ""
)

# Clean the Jets from the seleted leptons
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    src = cms.InputTag("patJets"), 
    preselection = cms.string('pt > 30.0 && abs(eta) < 2.4'),
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("selectedPatMuons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False), 
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), 
        ),
        electrons = cms.PSet(
           src       = cms.InputTag("selectedPatElectrons"),
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string(""),
           deltaR              = cms.double(0.5),
           checkRecoComponents = cms.bool(False), 
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), 
        )
    ),
    finalCut = cms.string('')
)


removeMCMatching(process,['All', 'PFAll'])


# Z Candidates and Higgs Candidates
process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass >50 '),
    decay = cms.string("selectedPatElectrons@+ selectedPatElectrons@-")
)

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 50 && min(abs(daughter(0).eta), abs(daughter(1).eta)) < 2.1'),
    decay = cms.string("selectedPatMuons@+ selectedPatMuons@-")
)

process.zem = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    cut = cms.string('mass > 50 '),
    decay = cms.string("selectedPatElectrons@+ selectedPatMuons@-")
)

process.zjj = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("cleanPatJets cleanPatJets")
)   

process.hzzeejjBaseColl = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("zee zjj")
)   

process.hzzmmjjBaseColl = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("zmm zjj")
)   

process.hzzmejjBaseColl = cms.EDProducer("CandViewCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string(''),
    decay = cms.string("zem zjj")
)   


process.hzzeejj = cms.EDProducer("Higgs2l2bUserDataNoMC",
    higgs = cms.InputTag("hzzeejjBaseColl"),
    #gensTag = cms.InputTag("genParticles"),
    metTag = cms.InputTag("patMETs")                          
    )

process.hzzmmjj = cms.EDProducer("Higgs2l2bUserDataNoMC",
    higgs = cms.InputTag("hzzmmjjBaseColl"),
    #gensTag = cms.InputTag("genParticles"),
    metTag = cms.InputTag("patMETs")
    )

process.hzzmejj = cms.EDProducer("Higgs2l2bUserDataNoMC",
    higgs = cms.InputTag("hzzmejjBaseColl"),
    #gensTag = cms.InputTag("genParticles"),
    metTag = cms.InputTag("patMETs")
    )



# Define the relevant paths and schedule them
process.analysisPath = cms.Path(
    process.eidSequence + 
    process.makePatElectrons +
    process.makePatMuons +
    process.makePatJets +
    process.makePatMETs +
    process.selectedPatElectrons + 
    process.selectedPatMuons + 
    process.cleanPatJets +
    process.zee +
    process.zmm +
    process.zme + 
    Process.zjj + 
    process.hzzeejjBaseColl + 
    process.hzzmmjjBaseColl +
    process.hzzmejjBaseColl +
    process.hzzeejj + 
    process.hzzmmjj +
    Process.hzzmejj #+ 

)

# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("zee", "zmm")
) 

process.zllFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("zll"),
    minNumber = cms.uint32(1),
)

process.jetFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("cleanPatJets"),
    minNumber = cms.uint32(2),
)


process.filterPath = cms.Path(process.zll+process.zllFilter+process.jetFilter)

process.out.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("filterPath",
                                   )
        )


process.outPath = cms.EndPath(process.out)
