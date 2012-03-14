from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

## global tag for MC (final JEC set for 42x 2011 data)
process.GlobalTag.globaltag = 'START42_V17::All'

### event FILTERS (not applied for the moment) #########
## HB + HE noise filtering
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
#
## ECAL noise filtering
#process.load('JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi')
#
#
#########################################################

#process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi") # for electron isolation

# Switch to PFMET 
from PhysicsTools.PatAlgos.tools.pfTools import *
switchToPFMET(
    process, 
    cms.InputTag('pfMet')
)

# --> JET configuration <--
# Jet energy corrections to use for default jets: 
inputJetCorrLabel = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
# 'L2L3Residual' not to be applied on MC  
# inputJetCorrLabel = ('AK5PF', ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

##### JET energy corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True

##### PFnoPU sequence

# selection of good PV for PU subtraction
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
process.pfPileUp.Enable = True
process.pfPileUp.checkClosestZVertex = cms.bool(False)
process.pfPileUp.Vertices = 'goodOfflinePrimaryVertices'

process.pfPUSequence = cms.Sequence( process.goodOfflinePrimaryVertices * process.pfPileUp * process.pfNoPileUp )

process.kt6PFnoPUJets = process.kt6PFJets.clone()
process.kt6PFnoPUJets.src = cms.InputTag("pfNoPileUp")
process.ak5PFnoPUJets = process.ak5PFJets.clone()
process.ak5PFnoPUJets.src  = cms.InputTag("pfNoPileUp")

#### NEW rho in restricted cone
process.kt6PFJetsForIsolation = process.kt6PFJets.clone(Rho_EtaMax = 2.5)
####

# PAT Jet tools
from PhysicsTools.PatAlgos.tools.jetTools import *

#### Define additional patJets 

# PFnoPUPUFastJet
# JEC for PFJets with charged hadron subtraction used (AK5PFchs)
addJetCollection(process,cms.InputTag('ak5PFnoPUJets'),
                  'AK5','PFnoPUPUFastJet',
                  doJTA = True,
                  doBTagging = True,
                  jetCorrLabel = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                  doType1MET = True,
                  doL1Cleaning = True,
                  doL1Counters = True,
                  genJetCollection=cms.InputTag("ak5GenJets"),
                  doJetID = True
                  )

process.patJetCorrFactorsAK5PFnoPUPUFastJet.rho=cms.InputTag("kt6PFnoPUJets","rho")

#### set default patJets
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
    cms.InputTag("impactParameterTagInfosAOD")
    )

process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')

### non default embedding of AOD items for default patJets
process.patJets.embedCaloTowers = cms.bool(False)
process.patJets.embedPFCandidates = cms.bool(True)
###

# Select jets
process.selectedPatJets.cut = cms.string('pt > 25.0 && abs(eta) < 2.4')
process.selectedPatJetsAK5PFnoPUPUFastJet.cut = cms.string('pt > 25.0 && abs(eta) < 2.4')

#add user variables to PAT-jets 
process.customPFJets = cms.EDProducer('PFJetUserData',
                                      JetInputCollection=cms.untracked.InputTag("selectedPatJets"),
                                      Verbosity=cms.untracked.bool(False)
                                      )

process.customPFnoPUJetsPUFastjet = cms.EDProducer('PFJetUserData',
                                      JetInputCollection=cms.untracked.InputTag("selectedPatJetsAK5PFnoPUPUFastJet"),
                                      Verbosity=cms.untracked.bool(False)
                                      )

#### end Jet configuration

# --> electron configuration <--

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

# loosely selected pat electrons
process.selectedPatElectrons.cut = (
            "pt > 10.0 && abs(eta) < 2.5 &&"                               +
            "(isEE || isEB) && !isEBEEGap"
        )


# Electrons with UserData
process.userDataSelectedElectrons = cms.EDProducer(
    "Higgs2l2bElectronUserData",
    src = cms.InputTag("selectedPatElectrons"),
    rho = cms.InputTag("kt6PFJetsForIsolation:rho")
)

# ID-selected electrons (only ID and conversion, no isolation)
process.selectedIDElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("userDataSelectedElectrons"),
    cut = cms.string("(electronID('eidVBTFCom95') == 7) ||"               +
                     " (electronID('eidVBTFCom95') == 5) "
                     )
)

# Isolated electrons: standard isolation; can be changed to corrected one
# NOTE: isolation is different for Barrel/Endcap
process.selectedIsoElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedIDElectrons"),
    cut = cms.string("electronID('eidVBTFCom95') == 7")
#    cut = cms.string('( isEB && (userFloat("absCombIsoPUCorrected") < 0.15*pt) ) || '      +
#                     '( !isEB && isEE && (userFloat("absCombIsoPUCorrected") < 0.1*pt) ) '
)

# --> muon configuration <--

# loosely selected pat muons
process.selectedPatMuons.cut = (
            "pt > 10 && abs(eta) < 2.4"
        )

# Muons with UserData
process.userDataSelectedMuons = cms.EDProducer(
    "Higgs2l2bMuonUserData",
    src = cms.InputTag("selectedPatMuons"),
    rho = cms.InputTag("kt6PFJetsForIsolation:rho")
)

# ID-selected muons
process.selectedIDMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("userDataSelectedMuons"),
    cut = cms.string(
            "pt > 10 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
            "innerTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
            "dB < 0.2 && numberOfMatches > 1 && abs(eta) < 2.4" )
)

# Isolated muons: standard isolation; can be changed to rho-corrected one
process.selectedIsoMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedIDMuons"),
#    cut = cms.string('userFloat("absCombIsoPUCorrected") < 0.15*pt')
    cut = cms.string("trackIso + caloIso < 0.15 * pt")
)

# Jet cleaning for patJets
process.cleanPatJetsIsoLept = cms.EDProducer("PATJetCleaner",
                           src = cms.InputTag("customPFJets"),
                           preselection = cms.string(''),
                           checkOverlaps = cms.PSet(
                              muons = cms.PSet(
                                          src = cms.InputTag("selectedIsoMuons"),
                                          algorithm = cms.string("byDeltaR"),
                                          preselection = cms.string(""),
                                          deltaR = cms.double(0.5),
                                          checkRecoComponents = cms.bool(False),
                                          pairCut = cms.string(""),
                                          requireNoOverlaps = cms.bool(True),
                              ),
                              electrons = cms.PSet(
                                              src = cms.InputTag("selectedIsoElectrons"),
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


# Jet cleaning for patJets noPU
process.cleanPatJetsNoPUIsoLept = cms.EDProducer("PATJetCleaner",
                           src = cms.InputTag("customPFnoPUJetsPUFastjet"),
                           preselection = cms.string(''),
                           checkOverlaps = cms.PSet(
                              muons = cms.PSet(
                                          src = cms.InputTag("selectedIsoMuons"),
                                          algorithm = cms.string("byDeltaR"),
                                          preselection = cms.string(""),
                                          deltaR = cms.double(0.5),
                                          checkRecoComponents = cms.bool(False),
                                          pairCut = cms.string(""),
                                          requireNoOverlaps = cms.bool(True),
                              ),
                              electrons = cms.PSet(
                                              src = cms.InputTag("selectedIsoElectrons"),
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


#################################################################################
# SEQUENCE TO GET Type1 Corrected MET (use default PFJet selection and cleaning)
#--------------------------------------------------------------------------------
# NOTE: do not compute Type 1 MET corrections for |eta| > 4.7,
#       in order to work around problem with CMSSW_4_2_x JEC factors at high eta,
#       reported in
#         https://hypernews.cern.ch/HyperNews/CMS/get/jes/270.html
#         https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1259/1.html
#
process.selectedPatJetsForMETtype1p2Corr = cms.EDFilter("PATJetSelector",
    src = cms.InputTag('patJets'),
    cut = cms.string('abs(eta) < 4.7'),
    filter = cms.bool(False)
)
#--------------------------------------------------------------------------------
# produce Type 1 corrections for pat::Jets of PF-type
process.patPFJetMETtype1p2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    src = cms.InputTag('selectedPatJetsForMETtype1p2Corr'),
    offsetCorrLabel = cms.string("L1FastJet"),
# ---> NOTE: use "L3Absolute" for MC / "L2L3Residual" for Data
    jetCorrLabel = cms.string("L3Absolute"),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)
#--------------------------------------------------------------------------------
# use MET corrections to produce Type 1 corrected PFMET objects
process.patType1CorrectedPFMet = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag('patMETs'),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('patPFJetMETtype1p2Corr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)   

process.met_Type1Corr_Sequence = cms.Sequence (
    process.selectedPatJetsForMETtype1p2Corr *
    process.patPFJetMETtype1p2Corr *
    process.patType1CorrectedPFMet 
)
#################################################################################


# Z Candidates and Higgs Candidates
# relax m_ll cut and take SameSign leptons
process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedElectrons@-")
                             )

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedMuons@+ userDataSelectedMuons@-")
                             )

process.zem = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedMuons@-")
                             )

process.zjj = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("cleanPatJetsIsoLept cleanPatJetsIsoLept")
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
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                 )

process.hzzmmjj = cms.EDProducer("Higgs2l2bUserData",
                                     higgs = cms.InputTag("hzzmmjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                     )

process.hzzemjj = cms.EDProducer("Higgs2l2bUserData",
                                     higgs = cms.InputTag("hzzemjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                     )

# Met variables producer (let's comment it for the moment)
#process.metInfoProducer = cms.EDProducer("MetVariablesProducer",
#                                    metTag = cms.InputTag("patMETs"),
#                                    t1CorrMetTag = cms.InputTag("patType1CorrectedPFMet")
#                                    )

# Add the files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
#file:/data3/scratch/cms/mc/Summer11/TT/EA1C42B3-E47D-E011-BA8D-001D0967D670.root'
'file:/data3/scratch/cms/mc/Summer11/GluGluToHToZZTo2L2Q_M-300/F6D80D8E-3794-E011-9800-00215E93D944.root'
 ] )

process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

# let it run

process.p = cms.Path(
#    process.HBHENoiseFilter*
#    process.EcalDeadCellEventFilter*
    process.pfPUSequence *
    process.kt6PFJets*
    process.kt6PFJetsForIsolation*
    process.ak5PFJets*
#    process.recoPFJets *
    process.kt6PFnoPUJets*
    process.ak5PFnoPUJets*
    process.eidSequence *
    process.patDefaultSequence *    
    process.userDataSelectedMuons *
    process.userDataSelectedElectrons *
    process.selectedIDMuons *
    process.selectedIDElectrons *
    process.selectedIsoMuons *
    process.selectedIsoElectrons *
    process.customPFJets *
    process.customPFnoPUJetsPUFastjet * 
    process.cleanPatJetsIsoLept *
    process.cleanPatJetsNoPUIsoLept *
######### TypeI corrected MET    
    process.met_Type1Corr_Sequence *
#########
    process.zee +
    process.zmm +
    process.zem +
    process.zjj +
    process.hzzeejjBaseColl +
    process.hzzmmjjBaseColl +
    process.hzzemjjBaseColl +
    process.hzzeejj +
    process.hzzmmjj +
    process.hzzemjj #+
#    process.metInfoProducer
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
                                 src = cms.InputTag("customPFJets"),
                                 minNumber = cms.uint32(2),
)


process.filterPath = cms.Path(
#    process.HBHENoiseFilter *
#    process.EcalDeadCellEventFilter *
    process.zll *
    process.zllFilter *
    process.jetFilter
)

# Output Module : Hopefully we keep all we need
process.out = cms.OutputModule("PoolOutputModule",
                 fileName = cms.untracked.string('h2l2qSkimData.root'),
                 SelectEvents = cms.untracked.PSet(
                    SelectEvents = cms.vstring("filterPath")
                 ),
                 outputCommands =  cms.untracked.vstring(
                  'drop *_*_*_*',
# PAT leptons and jets
                  'keep *_userDataSelectedElectrons_*_PAT',
                  'keep *_userDataSelectedMuons_*_PAT',
                  'keep *_customPFJets_*_PAT',
                  'keep *_customPFnoPUJetsPUFastjet_*_PAT',
                  'keep *_cleanPatJetsIsoLept_*_PAT',
# rho variables
                  'keep *_kt6PFJets_rho_PAT',
                  'keep *_kt6PFnoPUJets_rho_*',
                  'keep *_kt6PFJetsForIsolation_rho_*',
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
                  'keep *_secondaryVertexTagInfos*_*_*',
                  'keep *_*_*tagInfo*_*',
# additional collections from AOD   
                  'keep *_generalTracks_*_*',
                  'keep *_electronGsfTracks_*_*',
                  'keep *_muons_*_*',
                  'keep *_globalMuons_*_*',
                  'keep *_standAloneMuons_*_*',
                  'keep recoPFCandidates_particleFlow_*_*',
# genParticles & genJets
                  'keep *_genParticles_*_*',
                  'keep recoGenJets_ak5GenJets_*_*',
                  'keep recoGenJets_kt6GenJets_*_*',
# gen Info
                  'keep PileupSummaryInfos_*_*_*',
                  'keep GenEventInfoProduct_*_*_*',
                  'keep GenRunInfoProduct_*_*_*',
                  'keep LHEEventProduct_*_*_*',
                  'keep *_genEventScale_*_*',
###### MET products
                  'keep *_patMETs_*_*',
                  'keep *_patType1CorrectedPFMet_*_*',
### for HLT selection
                  'keep edmTriggerResults_TriggerResults_*_HLT',
#                  'keep *_hltTriggerSummaryAOD_*_*',
                  ),
)

process.out.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring("filterPath",
                                       )
)
process.out.dropMetaData = cms.untracked.string("DROPPED")

# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, sequence = 'p', hltProcess = '*' )

process.outPath = cms.EndPath(process.out)

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)

# process all the events
process.maxEvents.input = 200
process.options.wantSummary = True

