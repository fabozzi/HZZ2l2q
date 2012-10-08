## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# turn on when running on MC
runOnMC = True

# turn on when running on powheg signal MC 
isPowhegSignal = True

#add the L2L3Residual corrections only for data
if runOnMC:#MC
    jetCorrections=['L1FastJet','L2Relative','L3Absolute']
else:#Data
    jetCorrections=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

############ general options ####################
process.options.wantSummary = True
process.maxEvents.input = 400
process.MessageLogger.cerr.FwkReport.reportEvery = 100
########### global tag ############################
#from CMGTools.Common.Tools.getGlobalTag import getGlobalTag
#process.GlobalTag.globaltag = cms.string(getGlobalTag(runOnMC))
process.GlobalTag.globaltag = 'START44_V9C::All'
##################################################

############ PRINTOUT ###################

sep_line = "-" * 50
print sep_line
print 'running the following PFBRECO+PAT sequences:'
print '\tAK5'
print 'run on MC        : ', runOnMC
print sep_line
print 'Global tag       : ', process.GlobalTag.globaltag
print sep_line

######################################################

### INPUT COLLECTIONS ##########

process.source.fileNames = [
    'file:/data3/scratch/cms/mc/Fall11_44x/GluGluToHToZZTo2L2Q_M-400_PU_S6_START44_V9B/F8647FF8-C82E-E111-A2B1-00215E21D62A.root'
]

### DEFINITION OF THE PFBRECO+PAT SEQUENCES ##########
# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *

from PhysicsTools.PatAlgos.tools.pfTools import *

# ---------------- rho calculation for JEC ----------------------

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets

process.kt6PFJets = kt4PFJets.clone(
    rParam = cms.double(0.6),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
)

#compute rho correction for lepton isolation
process.kt6PFJetsForIso = process.kt6PFJets.clone(
    Rho_EtaMax = cms.double(2.5),
    Ghost_EtaMax = cms.double(2.5) )

# ---------------- Sequence AK5 ----------------------

postfixAK5 = "AK5"
jetAlgoAK5="AK5"

# Configure PAT to use PFBRECO instead of AOD sources
# this function will modify the PAT sequences. 
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgoAK5, runOnMC=runOnMC, postfix=postfixAK5,
          jetCorrections=('AK5PF', jetCorrections))

### DO NOT APPLY PFnoPU ###
getattr(process,"pfNoPileUp"+postfixAK5).enable = False
### setup additional projections ###
getattr(process,"pfNoMuon"+postfixAK5).enable = False 
getattr(process,"pfNoElectron"+postfixAK5).enable = False 
getattr(process,"pfNoTau"+postfixAK5).enable = False 
getattr(process,"pfNoJet"+postfixAK5).enable = True

# removing default cuts on muons 	 
getattr(process,"pfMuonsFromVertexAK5").dzCut = 99 	 
getattr(process,"pfMuonsFromVertexAK5").d0Cut = 99 	 
getattr(process,"pfSelectedMuonsAK5").cut="pt()>3" 	 
getattr(process,"pfIsolatedMuons"+postfixAK5).isolationCut = 999999 	 

# removing default cuts on electrons 	 
getattr(process,"pfElectronsFromVertexAK5").dzCut = 99 	 
getattr(process,"pfElectronsFromVertexAK5").d0Cut = 99 	 
getattr(process,"pfSelectedElectronsAK5").cut="pt()>5" 	 
getattr(process,"pfIsolatedElectrons"+postfixAK5).isolationCut = 999999 	 

# remove pfTau and pfPhoton from the sequence
process.PFBRECOAK5.remove( process.pfTauSequenceAK5 )
process.PFBRECOAK5.remove( process.pfNoTauAK5 )
process.PFBRECOAK5.remove( process.pfPhotonSequenceAK5 )

# make sure about patJets input
switchToPFJets(process, input=cms.InputTag('pfJetsAK5'), algo=jetAlgoAK5, postfix = postfixAK5, jetCorrections=('AK5PF', jetCorrections))

### we use "classic" muons and electrons (see below)
removeSpecificPATObjects(process, ['Taus'], postfix = "AK5")
removeSpecificPATObjects(process, ['Electrons'], postfix = "AK5")
removeSpecificPATObjects(process, ['Muons'], postfix = "AK5")
removeSpecificPATObjects(process, ['Photons'], postfix = "AK5")

############### remove useless modules #############
def removeUseless( modName ):
    getattr(process,"patDefaultSequence"+postfixAK5).remove(
        getattr(process, modName+postfixAK5)
        )

#removeUseless( "produceCaloMETCorrections" )
#removeUseless( "pfCandsNotInJet" )
#removeUseless( "pfJetMETcorr" )
#removeUseless( "pfCandMETcorr" )
#removeUseless( "pfchsMETcorr" )
#removeUseless( "pfType1CorrectedMet" )
#removeUseless( "pfType1p2CorrectedMet" )
removeUseless( "electronMatch" )
removeUseless( "muonMatch" )
removeUseless( "patPFTauIsolation" )
removeUseless( "tauMatch" )
removeUseless( "tauGenJets" )
removeUseless( "tauGenJetsSelectorAllHadrons" )
removeUseless( "tauGenJetMatch" )
removeUseless( "patShrinkingConePFTauDiscrimination" )
#########################################################

# curing a weird bug in PAT..
from CMGTools.Common.PAT.removePhotonMatching import removePhotonMatching
removePhotonMatching( process, postfixAK5 )

########## insert the PFMET significance calculation #############

process.load("CMGTools.Common.PAT.MetSignificance_cff")

setattr(process,"PFMETSignificance"+postfixAK5, process.pfMetSignificance.clone())
getattr(process,"patDefaultSequence"+postfixAK5).insert(getattr(process,"patDefaultSequence"+postfixAK5).index(getattr(process,"patMETs"+postfixAK5)),getattr(process,"PFMETSignificance"+postfixAK5))

####################################################################

########## add specific configuration for pat Jets ##############

getattr(process,"patJets"+postfixAK5).addTagInfos = True
getattr(process,"patJets"+postfixAK5).tagInfoSources  = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODAK5"),
    cms.InputTag("impactParameterTagInfosAODAK5")
    )
### non default embedding of AOD items for default patJets
getattr(process,"patJets"+postfixAK5).embedCaloTowers = False
getattr(process,"patJets"+postfixAK5).embedPFCandidates = True

### disable MC matching (will be done at analysis level)
getattr(process,"patJets"+postfixAK5).addGenPartonMatch = False
getattr(process,"patJets"+postfixAK5).addGenJetMatch = False

##############################################################
#add user variables to PAT-jets 
process.customPFJetsNoPUSub = cms.EDProducer(
    'PFJetUserData',
    JetInputCollection=cms.untracked.InputTag("selectedPatJetsAK5"),
    is2012Data=cms.untracked.bool(False),
    qgMap=cms.untracked.InputTag("qglAK5PF"),
    Verbosity=cms.untracked.bool(False)
    )

process.customPFJetsNoPUSubCentral = cms.EDFilter(
    "PATJetSelector",
    src = cms.InputTag("customPFJetsNoPUSub"),
    cut = cms.string("abs(eta) < 2.4")
    )

############## "Classic" PAT Muons and Electrons ########################
# (made from all reco muons, and all gsf electrons, respectively)
process.patMuons.embedTcMETMuonCorrs = False
process.patMuons.embedCaloMETMuonCorrs = False
process.patMuons.embedTrack = True

process.patElectrons.embedTrack = True
process.patElectrons.pfElectronSource = 'particleFlow'

# use PFIsolation
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons', 'PFIso')
process.muIsoSequence = setupPFMuonIso(process, 'muons', 'PFIso')
adaptPFIsoMuons( process, applyPostfix(process,"patMuons",""), 'PFIso')
adaptPFIsoElectrons( process, applyPostfix(process,"patElectrons",""), 'PFIso')

# setup recommended 0.3 cone for electron PF isolation
#process.pfIsolatedElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFIso"))
#process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFIso")
#process.pfIsolatedElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"), cms.InputTag("elPFIsoValueGamma03PFIdPFIso"))
process.patElectrons.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
    )
process.patElectrons.isolationValuesNoPFId = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
    )


process.stdMuonSeq = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.muIsoSequence +
    process.makePatMuons +
    process.selectedPatMuons
    )
process.stdElectronSeq = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.eleIsoSequence +
    process.makePatElectrons +
    process.selectedPatElectrons
    )

process.stdMuonSeq.remove( process.muonMatch )
process.stdElectronSeq.remove( process.electronMatch )
process.patMuons.addGenMatch = False
process.patElectrons.addGenMatch = False
process.patMuons.embedGenMatch = False
process.patElectrons.embedGenMatch = False

# Modules for Electron ID
# MVA Electron ID
process.load("EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi")
process.mvaeIdSequence = cms.Sequence(
    process.mvaTrigV0 +
    process.mvaNonTrigV0
)
# ElectronID in the VBTF prescription
#import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
# Switch to the official Electron VBTF Selection for 2011 Data (relax H/E cut in the Endcap):
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/VbtfEleID2011
import HiggsAnalysis.Higgs2l2b.simpleCutBasedElectronIDSummer11_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )
        
process.vbtfeIdSequence = cms.Sequence(
        process.eidVBTFRel95 +
        process.eidVBTFRel80 +
        process.eidVBTFCom95 +
        process.eidVBTFCom80
)

process.eidSequence = cms.Sequence(
    process.mvaeIdSequence +
    process.vbtfeIdSequence
)

process.patElectrons.electronIDSources = cms.PSet(
        eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
        eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
        eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
        eidVBTFCom80 = cms.InputTag("eidVBTFCom80"),
        #MVA 
        mvaTrigV0 = cms.InputTag("mvaTrigV0"),
        mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
)

process.stdLeptonSequence = cms.Sequence(
    process.stdMuonSeq +
    process.eidSequence +
    process.stdElectronSeq 
    )

# Classic Electrons with UserData

process.userDataSelectedElectrons = cms.EDProducer(
    "Higgs2l2bElectronUserData",
    src = cms.InputTag("selectedPatElectrons"),
    rho = cms.InputTag("kt6PFJetsForIso:rho"),
    primaryVertices=cms.InputTag("offlinePrimaryVertices")
)

# ID-selected electrons (only ID and conversion, no isolation)
process.selectedIDElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("userDataSelectedElectrons"),
    cut = cms.string("(electronID('eidVBTFCom95') == 7) ||"               +
                     " (electronID('eidVBTFCom95') == 5) "
                     )
)

# Isolated electrons: standard isolation
process.selectedIsoElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedIDElectrons"),
    cut = cms.string("electronID('eidVBTFCom95') == 7")
)

# Classic Muons with UserData
process.userDataSelectedMuons = cms.EDProducer(
    "Higgs2l2bMuonUserData",
    src = cms.InputTag("selectedPatMuons"),
    rho = cms.InputTag("kt6PFJetsForIso:rho"),
    primaryVertices=cms.InputTag("offlinePrimaryVertices")
)

# ID-selected muons
process.selectedIDMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("userDataSelectedMuons"),
    cut = cms.string(
            "isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
            "globalTrack().hitPattern().numberOfValidTrackerHits > 10 && "  +
            "globalTrack().hitPattern().numberOfValidPixelHits > 0 && "     +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "      +
            "abs(dB) < 0.2 && numberOfMatches > 1" )
)

# Isolated muons: standard isolation
process.selectedIsoMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedIDMuons"),
    cut = cms.string("trackIso + caloIso < 0.15 * pt")
)

process.userDataStandardLeptonSequence = cms.Sequence(
    process.userDataSelectedMuons *
    process.userDataSelectedElectrons *
    process.selectedIDMuons *
    process.selectedIDElectrons *
    process.selectedIsoMuons *
    process.selectedIsoElectrons 
    )


# Jet cleaning for patJets
process.cleanPatJetsNoPUIsoLept = cms.EDProducer(
    "PATJetCleaner",
    src = cms.InputTag("customPFJetsNoPUSubCentral"),
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

# ---------------- Common stuff ---------------

### PATH DEFINITION #############################################

# counters that can be used at analysis level to know the processed events
process.prePathCounter = cms.EDProducer("EventCountProducer")
process.postPathCounter = cms.EDProducer("EventCountProducer")

# trigger information (no selection)

process.p = cms.Path( process.prePathCounter )

#process.p += process.kt6PFJets

# PFBRECO+PAT ---

# used by metsignificance module
process.p += process.pfNoJet

process.p += getattr(process,"patPF2PATSequence"+postfixAK5)

#process.p += process.kt6PFJetsCHS
process.p += process.kt6PFJetsForIso
#process.p += process.kt6PFJetsCHSForIso

process.p += process.customPFJetsNoPUSub
#process.p += process.puJetIdSequenceAK5
process.p += process.customPFJetsNoPUSubCentral

process.p += process.stdLeptonSequence

process.p += process.userDataStandardLeptonSequence
process.p += process.cleanPatJetsNoPUIsoLept

# Select leptons
process.selectedPatMuons.cut = (
    "pt > 10 && abs(eta) < 2.4"
        )

process.selectedPatElectrons.cut = (
    "pt > 10.0 && abs(eta) < 2.5"
    )
# Select jets
process.selectedPatJetsAK5.cut = cms.string('pt > 25.0')

################# COMBINATORIAL ANALYSIS ###########################

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
                             checkOverlap = cms.bool(False),
                             cut = cms.string(''),
                             decay = cms.string("cleanPatJetsNoPUIsoLept cleanPatJetsNoPUIsoLept")
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

process.combinatorialSequence = cms.Sequence(
    process.zee +
    process.zmm +
    process.zem +
    process.zjj +
    process.hzzeejjBaseColl +
    process.hzzmmjjBaseColl +
    process.hzzemjjBaseColl +
    process.hzzeejj +
    process.hzzmmjj +
    process.hzzemjj 
)

process.p += process.combinatorialSequence

process.p += getattr(process,"postPathCounter") 


# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
                             src = cms.VInputTag("zee", "zmm", "zem")
)

process.zllFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("zll"),
                                 minNumber = cms.uint32(1),
)

process.jetFilterNoPUSub = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("customPFJetsNoPUSubCentral"),
                                 minNumber = cms.uint32(2),
)

process.filterPath1= cms.Path(
    process.zll *
    process.zllFilter *
    process.jetFilterNoPUSub
)

# event cleaning (in tagging mode, no event rejected)
process.load('CMGTools.Common.PAT.addFilterPaths_cff')

process.fullPath = cms.Schedule(
    process.p,
    process.filterPath1,
#    process.filterPath2,
    process.EcalDeadCellTriggerPrimitiveFilterPath,
    process.hcalLaserEventFilterPath,
    process.trackingFailureFilterPath,
    process.CSCTightHaloFilterPath,
    process.HBHENoiseFilterPath,
#    process.eeBadScFilterPath,
    process.primaryVertexFilterPath,
    process.noscrapingFilterPath
    )

#this is needed only for Madgraph MC:
if runOnMC:
    process.fullPath.append(process.totalKinematicsFilterPath)
else:
    del process.totalKinematicsFilterPath

### OUTPUT DEFINITION #############################################

# PFBRECO+PAT ---

# Add PFBRECO output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patTriggerEventContent, patTriggerStandAloneEventContent


process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('h2l2qSkimData.root'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring(
        'filterPath1')#,
#        'filterPath2')
      ),
    outputCommands =  cms.untracked.vstring(
      'drop *_*_*_*',
      ),
    )

process.out.dropMetaData = cms.untracked.string("DROPPED")

process.out.outputCommands.extend([
    'keep *_userDataSelectedElectrons_*_PAT',
    'keep *_userDataSelectedMuons_*_PAT',
#    'keep *_customPFJets_*_PAT',
    'keep *_selectedPatJetsAK5_pfCandidates_PAT',
    'keep *_customPFJetsNoPUSub_*_PAT',
    'keep *_cleanPatJetsNoPUIsoLept_*_PAT',
    # rho variables
    'keep *_*_rho_PAT',
    'keep *_kt6PFJetsCentralNeutral_rho_*',
    'keep *_kt6PFJets*_rho_*',
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
#    'keep *_secondaryVertexTagInfos*_*_*',
#    'keep *_impactParameterTagInfos*_*_*',
#    'keep *_*_*tagInfo*_*',
    # additional collections from AOD   
#    'keep *_generalTracks_*_*',
#    'keep *_electronGsfTracks_*_*',
#    'keep *_muons_*_*',
#    'keep *_globalMuons_*_*',
#    'keep *_standAloneMuons_*_*',
    'keep recoPFCandidates_particleFlow_*_*',
    # genParticles & genJets
    'keep *_genParticles_*_*',
    'keep recoGenJets_ak5GenJets_*_*',
#    'keep recoGenJets_kt6GenJets_*_*',
    # gen Info
    'keep PileupSummaryInfos_*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    'keep LHEEventProduct_*_*_*',
    'keep *_genEventScale_*_*',
    ###### MET products
    'keep *_patMETs*_*_*',
#    'keep *_patType1CorrectedPFMet_*_*', # NOT included for the moment
    ### for HLT selection
    'keep edmTriggerResults_TriggerResults_*_HLT'])

# additional products for event cleaning
process.out.outputCommands.extend([
    'keep *_TriggerResults_*_PAT',
    ])

process.out.outputCommands.extend(['keep edmMergeableCounter_*_*_*'])

### Ntuplization ###
process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2qedmNtuples_52_cff")

process.genInfoNtuple = cms.EDProducer(
    "GenInfoNtupleDump",
    src = cms.InputTag("genParticles")
    )

# select lepton daughters from the Higgs
process.genSelectorZDaughterLeptons = cms.EDFilter(
    "GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string(' ((abs(pdgId)==11 ) || (abs(pdgId)==13 )) && abs(mother.pdgId)==23 ')
    )

# select quark daughters from the Higgs
process.genSelectorZDaughterQuarks = cms.EDFilter(
    "GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string(' ((abs(pdgId)>=1 ) && (abs(pdgId)<=6 )) && abs(mother.pdgId)==23 ')
    )

# perform MC matching
process.hCandMatch = cms.EDProducer(
    "HiggsMatcher",
    elhiggs = cms.InputTag("hzzeejj:h"),
    muhiggs = cms.InputTag("hzzmmjj:h"),
    genLept = cms.InputTag("genSelectorZDaughterLeptons"),
    genQuarks = cms.InputTag("genSelectorZDaughterQuarks")
    )

process.PUInfoNtuple = cms.EDProducer(
    "GenPUNtupleDump",
    isData = cms.bool(False)
)

# Event rho dumper
process.rhoDumper = cms.EDProducer("EventRhoDumper",
                                    rho = cms.InputTag("kt6PFJets:rho"),
                                    restrictedRho = cms.InputTag("kt6PFJetsForIso:rho")
                                    )


# Met variables producer
process.metInfoProducer = cms.EDProducer("MetVariablesProducer",
                                    metTag = cms.InputTag("patMETsAK5"),
                                    t1CorrMetTag = cms.InputTag("patType1CorrectedPFMetAK5")
                                    )

process.analysisPath = cms.Sequence(
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.rhoDumper+
    process.metInfoProducer+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple+
    process.Higgsemu2bEdmNtuple+
    process.jetinfos
)

if runOnMC and isPowhegSignal:
    process.analysisPath += process.genInfoNtuple
    process.analysisPath += process.genSelectorZDaughterLeptons
    process.analysisPath += process.genSelectorZDaughterQuarks
    process.analysisPath += process.hCandMatch

process.p += process.analysisPath

process.edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('h2l2q_ntuple.root'),
    outputCommands = cms.untracked.vstring(
      "drop *",
      "keep *_hCandMatch_*_*",
      "keep *_genInfoNtuple_*_*",
      "keep *_eventVtxInfoNtuple_*_*",
      "keep *_PUInfoNtuple_*_*",
      "keep *_rhoDumper_*_*",
      "keep *_metInfoProducer_*_*",
      'keep *_kt6PFJetsCentralNeutral_rho_*',
      "keep *_Higgs2e2bEdmNtuple_*_*",
      "keep *_Higgs2mu2bEdmNtuple_*_*",
      "keep *_Higgsemu2bEdmNtuple_*_*",
      "keep *_jetinfos_*_*"
      ),
    dropMetaData = cms.untracked.string('ALL'),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring(
        'filterPath1')#,
#        'filterPath2')
      ),
    )

process.edmNtuplesOut.outputCommands.extend([
    'keep edmTriggerResults_TriggerResults_*_HLT',
    'keep *_TriggerResults_*_PAT',
    ])


process.endPath = cms.EndPath(process.edmNtuplesOut)
 

