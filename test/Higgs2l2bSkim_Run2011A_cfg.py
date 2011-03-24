# This is an example PAT configuration showing the usage of PAT on minbias data

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

## global tag for data
process.GlobalTag.globaltag = 'GR_R_311_V2::All'


# Triggers for the /Jet PD are from:
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HLTrigger/Configuration/python/HLT_FULL_data_cff.py?revision=1.63&view=markup
# The list is here, let's take all of them :
## Jet = cms.vstring( 'HLT_DiJetAve100U_v4',
##   'HLT_DiJetAve140U_v4',
##   'HLT_DiJetAve15U_v4',
##   'HLT_DiJetAve180U_v4',
##   'HLT_DiJetAve300U_v4',
##   'HLT_DiJetAve30U_v4',
##   'HLT_DiJetAve50U_v4',
##   'HLT_DiJetAve70U_v4',
##   'HLT_Jet110_v1',
##   'HLT_Jet150_v1',
##   'HLT_Jet190_v1',
##   'HLT_Jet240_v1',
##   'HLT_Jet30_v1',
##   'HLT_Jet370_NoJetID_v1',
##   'HLT_Jet370_v1',
##   'HLT_Jet60_v1',
##   'HLT_Jet80_v1' ),
mytrigs = ['*']

# Jet energy corrections to use: 
inputJetCorrLabel = ('AK5PF', ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

# add pf met
from PhysicsTools.PatAlgos.tools.metTools import *
removeMCMatching(process, ['All'])
addPfMET(process, 'PF')

# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
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


# Taus are currently broken in 4.1.x
removeSpecificPATObjects( process, ['Taus'] )
process.patDefaultSequence.remove( process.patTaus )

# require physics declared
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )
# HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if mytrigs is not None :
    process.hltSelection = hltHighLevel.clone(TriggerResultsTag = 'TriggerResults::HLT', HLTPaths = mytrigs)
    process.hltSelection.throw = False


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )


# Select jets
process.selectedPatJets.cut = cms.string('pt > 25')

# Add the files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()



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
                           src = cms.InputTag("patJets"),
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

process.hzzemjjBaseColl = cms.EDProducer("CandViewCombiner",
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

process.hzzemjj = cms.EDProducer("Higgs2l2bUserDataNoMC",
                                     higgs = cms.InputTag("hzzemjjBaseColl"),
                                     #gensTag = cms.InputTag("genParticles"),
                                     metTag = cms.InputTag("patMETs")
                                     )


readFiles.extend( [
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/547/DAFCA3B7-B850-E011-9ADF-0030487C778E.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/502/7297DBD1-6B50-E011-9FA5-0030487CD162.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/500/9E1BCB0E-6C50-E011-BC4D-0030487CD7EE.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/499/3850D990-9450-E011-BF72-0030487CD716.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/498/AC8FE79D-8450-E011-A58B-0030487CAEAC.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/495/88BEEA68-3850-E011-9F13-0030487CD17C.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/494/1C35E662-3350-E011-BAA0-003048F1183E.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/488/66D0A0F1-1850-E011-97C0-003048F1BF68.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/486/E2FB2282-1850-E011-BA61-001D09F2527B.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/484/6E1A5CBE-1650-E011-B823-003048D2C108.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/472/CA224AEB-1450-E011-84A1-001D09F2527B.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/471/5C765497-1450-E011-AD91-001D09F248F8.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/469/204A3C03-1A50-E011-B8D7-0030487C8CB6.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/467/6ABD98CE-1B50-E011-9D5F-003048F11114.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/466/F296557B-1E50-E011-8C7F-0030487CD6B4.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/463/72D3943E-2350-E011-B098-001617C3B654.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/462/E650148E-2650-E011-96BD-001617C3B654.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/456/863BF1F0-4C50-E011-9ABC-001617C3B5E4.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/455/46A4141B-4A50-E011-8ED4-001D09F24934.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/454/C2041925-DC50-E011-B0C2-0030487C7828.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/450/2EE61936-5050-E011-9731-0030487C7392.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/449/5C5D099B-4C50-E011-9779-000423D98B6C.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/447/B8CFC010-5650-E011-A6CA-0030487CD162.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/446/6810F3AF-5C50-E011-9EA3-003048F11114.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/445/5E83EDE8-5650-E011-B45D-001D09F23A3E.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/444/0E40C747-AC50-E011-8A91-0030487CD718.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/443/EAA93DE5-5750-E011-A127-001617C3B76E.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/442/70B27507-5B50-E011-99A9-0030487CD162.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/439/6658D888-B34F-E011-BE9D-001D09F24D67.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/433/12B7CF7F-B54F-E011-94B6-001D09F23C73.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/432/96E7E8AA-1150-E011-878B-001D09F28F1B.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/431/44A96589-2D50-E011-861C-003048F118D4.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/428/00E2E174-B34F-E011-9858-003048F110BE.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/427/52C35785-814F-E011-A17C-003048F117EC.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/425/F4681D39-844F-E011-B3E5-003048F024F6.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/423/2A0114B3-7C4F-E011-B719-001D09F24303.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/421/0A0F17E0-774F-E011-93C1-000423D94E70.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/413/28A09525-EB4F-E011-9AE5-003048F0258C.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/410/8094C4AE-5F4F-E011-B390-003048F118D2.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/406/4E1472A9-5F4F-E011-B062-0030487C7828.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/405/3AB185A4-D64F-E011-8394-0030487C2B86.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/404/38F2FFDA-374F-E011-AD2A-003048F024F6.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/403/64B3598B-364F-E011-8683-001D09F29849.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/386/5059F123-0B4F-E011-8DDB-0030487C2B86.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/384/463D8003-094F-E011-A4B1-0030487CD17C.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/383/7EDAB706-0D4F-E011-BBBF-0030487C6A66.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/379/14D52E14-0C4F-E011-B5AF-003048F11C58.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/377/5CBF6B69-0B4F-E011-855D-0030487CD906.root',
'/store/data/Run2011A/Jet/AOD/PromptReco-v1/000/160/329/6631A77A-2F4E-E011-AE5A-003048F118C2.root',
 ] )

process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

# let it run

process.p = cms.Path(
    process.hltSelection*
    process.scrapingVeto*
    process.primaryVertexFilter*
    process.HBHENoiseFilter*
    process.eidSequence *
    process.patDefaultSequence*
    process.cleanPatJets *
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
#    process.goodPatJets
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
    process.hltSelection *
    process.scrapingVeto *
    process.primaryVertexFilter *
    process.HBHENoiseFilter *
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
                  'keep *_selectedPatElectrons_*_PAT',
                  'keep *_selectedPatMuons_*_PAT',
                  'keep *_cleanPatJets_*_PAT',
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
                  'keep *_offlinePrimaryVertices_*_*'
                ),
)

process.out.SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring("filterPath",
                                       )
)

process.outPath = cms.EndPath(process.out)

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = 1000
process.options.wantSummary = True

