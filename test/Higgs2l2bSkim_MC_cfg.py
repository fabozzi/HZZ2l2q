from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

## global tag for data
process.GlobalTag.globaltag = 'START42_V12::All'


# Jet energy corrections to use: 
inputJetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
# 'L2L3Residual' not to be applied on MC  
# inputJetCorrLabel = ('AK5PF', ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

# Switch to using PFMET 
from PhysicsTools.PatAlgos.tools.pfTools import *
switchToPFMET(
    process, 
    cms.InputTag('pfMet'), 
    ""
)

# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
#process.kt6PFJets.Rho_EtaMax = cms.double(5.0)
#process.kt6PFJets.Ghost_EtaMax = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True
#process.ak5PFJets.Rho_EtaMax = cms.double(2.5)

process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')

addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5','PFOffset',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = False,
                 doL1Cleaning   = False,
                 doL1Counters   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )

switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = inputJetCorrLabel,
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )
process.patJets.addTagInfos = True
process.patJets.tagInfoSources  = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAOD"),
    )

# Apply loose PF jet ID
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                   filterParams = pfJetIDSelector.clone(),
                                   src = cms.InputTag("selectedPatJets"),
                                   filter = cms.bool(True)
                                   )

# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

# ECAL noise filtering
process.load('JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi')

# Select jets
process.selectedPatJets.cut = cms.string('pt > 25')

# Add the files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()



# Modules for the Cut-based Electron ID in the VBTF prescription
#import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
# Switch to the official Electron VBTF Selection for 2011 Data (relax H/E cut in the Endcap):
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/VbtfEleID2011
import HiggsAnalysis.Higgs2l2b.simpleCutBasedElectronIDSummer11_cfi as vbtfid
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

process.patElectrons.electronIDSources = cms.PSet(
        eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
        eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
        eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
        eidVBTFCom80 = cms.InputTag("eidVBTFCom80")
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

process.selectedPatElectrons.cut = (
        "pt > 10.0 && abs(eta) < 2.5 &&"                               +
            "(isEE || isEB) && !isEBEEGap &&"                              +
            "electronID('eidVBTFCom95') == 7"
        )

# Clean the Jets from the seleted leptons
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
                           src = cms.InputTag("selectedPatJets"),
                           preselection = cms.string('pt > 30.0 && abs(eta) < 2.4'),
                           checkOverlaps = cms.PSet(
                              muons = cms.PSet(
                                          src = cms.InputTag("selectedPatMuons"),
                                          algorithm = cms.string("byDeltaR"),
                                          preselection = cms.string(""),
                                          deltaR = cms.double(0.5),
                                          checkRecoComponents = cms.bool(False),
                                          pairCut = cms.string(""),
                                          requireNoOverlaps = cms.bool(True),
                              ),
                              electrons = cms.PSet(
                                              src = cms.InputTag("selectedPatElectrons"),
                                              algorithm = cms.string("byDeltaR"),
                                              preselection = cms.string(""),
                                              deltaR = cms.double(0.5),
                                              checkRecoComponents = cms.bool(False),
                                              pairCut = cms.string(""),
                                              requireNoOverlaps = cms.bool(True),
                             )
                           ),
                           finalCut = cms.string('')
)

process.cleanPatJetsAK5PFOffset = cms.EDProducer("PATJetCleaner",
                           src = cms.InputTag("patJetsAK5PFOffset"),
                           preselection = cms.string('pt > 30.0 && abs(eta) < 2.4'),
                           checkOverlaps = cms.PSet(
                              muons = cms.PSet(
                                          src = cms.InputTag("selectedPatMuons"),
                                          algorithm = cms.string("byDeltaR"),
                                          preselection = cms.string(""),
                                          deltaR = cms.double(0.5),
                                          checkRecoComponents = cms.bool(False),
                                          pairCut = cms.string(""),
                                          requireNoOverlaps = cms.bool(True),
                              ),
                              electrons = cms.PSet(
                                              src = cms.InputTag("selectedPatElectrons"),
                                              algorithm = cms.string("byDeltaR"),
                                              preselection = cms.string(""),
                                              deltaR = cms.double(0.5),
                                              checkRecoComponents = cms.bool(False),
                                              pairCut = cms.string(""),
                                              requireNoOverlaps = cms.bool(True),
                             )
                           ),
                           finalCut = cms.string('')
)

# Z Candidates and Higgs Candidates
process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(True),
                                 cut = cms.string('mass >50 '),
                                 decay = cms.string("selectedPatElectrons@+ selectedPatElectrons@-")
                             )

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(True),
                                 cut = cms.string('mass > 50 '),
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

process.hzzemjjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zem zjj")
                                         )


process.hzzeejj = cms.EDProducer("Higgs2l2bUserData",
                                     higgs = cms.InputTag("hzzeejjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     metTag = cms.InputTag("patMETs")
                                     )

process.hzzmmjj = cms.EDProducer("Higgs2l2bUserData",
                                     higgs = cms.InputTag("hzzmmjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     metTag = cms.InputTag("patMETs")
                                     )

process.hzzemjj = cms.EDProducer("Higgs2l2bUserData",
                                     higgs = cms.InputTag("hzzemjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     metTag = cms.InputTag("patMETs")
                                     )

# Met variables producer
process.metInfoProducer = cms.EDProducer("MetVariablesProducer",
                                    metTag = cms.InputTag("patMETs")
                                    )

readFiles.extend( [
'file:/data3/scratch/cms/mc/Summer11/GluGluToHToZZTo2L2Q_M-300/F6D80D8E-3794-E011-9800-00215E93D944.root'
 ] )

process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

# let it run

process.p = cms.Path(
    process.HBHENoiseFilter*
    process.EcalDeadCellEventFilter*
    process.recoPFJets *
    process.eidSequence *
    process.patDefaultSequence *
    process.cleanPatJetsAK5PFOffset *
    process.zee +
    process.zmm +
    process.zem +
    process.zjj +
    process.hzzeejjBaseColl +
    process.hzzmmjjBaseColl +
    process.hzzemjjBaseColl +
    process.hzzeejj +
    process.hzzmmjj +
    process.hzzemjj +
    process.metInfoProducer
)

# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
                             src = cms.VInputTag("zee", "zmm", "zem")
)

process.zllFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("zll"),
                                 minNumber = cms.uint32(1),
)

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("cleanPatJets"),
                                 minNumber = cms.uint32(2),
)


process.filterPath = cms.Path(
    process.HBHENoiseFilter *
    process.EcalDeadCellEventFilter *
    process.zll *
    process.zllFilter *
    process.jetFilter
)

# Output Module : Hopefully we keep all we need
process.out = cms.OutputModule("PoolOutputModule",
                 fileName = cms.untracked.string('h2l2bData.root'),
                 SelectEvents = cms.untracked.PSet(
                    SelectEvents = cms.vstring("filterPath")
                 ),
                 outputCommands =  cms.untracked.vstring(
                  'drop *_*_*_*',
                  'keep *_genParticles_*_*',
                  'keep *_selectedPatElectrons_*_PAT',
                  'keep *_selectedPatMuons_*_PAT',
                  'keep *_cleanPatJets_*_PAT',
                  'keep *_cleanPatJetsAK5PFOffset_*_PAT',
                  'keep *_kt6PFJets_rho_PAT',
                  'keep *_zee_*_PAT',
                  'keep *_zmm_*_PAT',
                  'keep *_zem_*_PAT',
                  'keep *_zjj_*_PAT',
                  'keep *_hzzeejj_*_PAT',
                  'keep *_hzzmmjj_*_PAT',
                  'keep *_hzzemjj_*_PAT',
                #'keep *_TriggerResults*_*_HLT',
                #'keep *_hltTriggerSummaryAOD_*_HLT',
                #'keep *_TriggerResults*_*_REDIGI*',
                #'keep *_hltTriggerSummaryAOD_*_REDIGI*',
                  'keep *_offlinePrimaryVertices_*_*',
                  'keep *_secondaryVertexTagInfos*_*_*',
                  'keep *_*_*tagInfo*_*',
                  'keep *_generalTracks_*_*',
                  'keep PileupSummaryInfos_*_*_*',
                  'keep *_metInfoProducer_*_*',
                  ),
)

process.out.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring("filterPath",
                                       )
)


# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, sequence = 'p', hltProcess = '*' )

# PAT trigger matching for muons
process.muonTriggerMatchHLTMuon = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                   src = cms.InputTag( "selectedPatMuons" ),
                                                   matched = cms.InputTag( "patTrigger" ),
                                                   matchedCuts = cms.string( 'path( "HLT_DoubleMu7" )' ),
                                                   maxDPtRel = cms.double( 1000.0 ),
                                                   maxDeltaR = cms.double( 0.2 ),
                                                   resolveAmbiguities    = cms.bool( True ),
                                                   resolveByMatchQuality = cms.bool( True )

)

# PAT trigger matching for electrons
process.electronTriggerMatchHLTElectron = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                          src = cms.InputTag( "selectedPatElectrons" ),
                                                          matched = cms.InputTag( "patTrigger" ),
                                                          matchedCuts = cms.string( 'path( "Ele17_CaloIdL_CaloIsoV_Ele8_CaloIdL_CaloIsoVL_v1" )' ),
                                                          maxDPtRel = cms.double( 1000.0 ),
                                                          maxDeltaR = cms.double( 0.2 ),
                                                          resolveAmbiguities    = cms.bool( True ),
                                                          resolveByMatchQuality = cms.bool( True )
)


switchOnTriggerMatching( process, [ 'muonTriggerMatchHLTMuon','electronTriggerMatchHLTElectron' ], sequence = 'p', hltProcess = '*' )

process.outPath = cms.EndPath(process.out)

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

# process all the events
process.maxEvents.input = 1000
process.options.wantSummary = True

