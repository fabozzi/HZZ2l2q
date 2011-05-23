
import FWCore.ParameterSet.Config as cms
import copy

### HIGGS variables


higgs =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("h2l2b"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("H"),
    eventInfo=cms.untracked.bool(True),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Mass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta")
    ),
    cms.PSet(
    tag = cms.untracked.string("Y"),
    quantity = cms.untracked.string("y")
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi")
    ),
    cms.PSet(
    tag = cms.untracked.string("Met"),
    quantity = cms.untracked.string("userFloat('met')")
    ),
    cms.PSet(
    tag = cms.untracked.string("Met2"),
    quantity = cms.untracked.string("userFloat('met2')")
    ),
    cms.PSet(
    tag = cms.untracked.string("MetSig"),
    quantity = cms.untracked.string("userFloat('metSig')")
    ),
    cms.PSet(
    tag = cms.untracked.string("MetPhi"),
    quantity = cms.untracked.string("userFloat('metPhi')")
    ),
    cms.PSet(
    tag = cms.untracked.string("AzimuthalAngle"),
    quantity = cms.untracked.string("userFloat('azimuthalAngle')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jminbmatch"),
    quantity = cms.untracked.string("userFloat('jminbmatch')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jmincmatch"),
    quantity = cms.untracked.string("userFloat('jmincmatch')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jmaxbmatch"),
    quantity = cms.untracked.string("userFloat('jmaxbmatch')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jmaxcmatch"),
    quantity = cms.untracked.string("userFloat('jmaxcmatch')")
    ),
    cms.PSet(
    tag = cms.untracked.string("zzdPhi"),
    quantity = cms.untracked.string("userFloat('zzdPhi')")
    ),
    cms.PSet(
    tag = cms.untracked.string("zzdEta"),
    quantity = cms.untracked.string("userFloat('zzdEta')")
    ),
    cms.PSet(
    tag = cms.untracked.string("zzdr"),
    quantity = cms.untracked.string("userFloat('zzdr')")
    ),
    cms.PSet(
    tag = cms.untracked.string("lldPhi"),
    quantity = cms.untracked.string("userFloat('lldPhi')")
    ),
    cms.PSet(
    tag = cms.untracked.string("lldEta"),
    quantity = cms.untracked.string("userFloat('lldEta')")
    ),
    cms.PSet(
    tag = cms.untracked.string("lldr"),
    quantity = cms.untracked.string("userFloat('lldr')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jjdPhi"),
    quantity = cms.untracked.string("userFloat('jjdPhi')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jjdEta"),
    quantity = cms.untracked.string("userFloat('jjdEta')")
    ),
    cms.PSet(
    tag = cms.untracked.string("jjdr"),
    quantity = cms.untracked.string("userFloat('jjdr')")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau1LooseID"),
    quantity = cms.untracked.string("userFloat('jet1LooseID')")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau2LooseID"),
    quantity = cms.untracked.string("userFloat('jet2LooseID')")
    ),
    cms.PSet(
    tag = cms.untracked.string("costhetaNT1"),
    quantity = cms.untracked.string("userFloat('costhetaNT1')")
    ),
    cms.PSet(
    tag = cms.untracked.string("costhetaNT2"),
    quantity = cms.untracked.string("userFloat('costhetaNT2')")
    ),
    cms.PSet(
    tag = cms.untracked.string("phistarNT1"),
    quantity = cms.untracked.string("userFloat('phistarNT1')")
    ),
    cms.PSet(
    tag = cms.untracked.string("phistarNT2"),
    quantity = cms.untracked.string("userFloat('phistarNT2')")
    ),
    cms.PSet(
    tag = cms.untracked.string("phistarNT12"),
    quantity = cms.untracked.string("userFloat('phistarNT12')")
    ),
    cms.PSet(
    tag = cms.untracked.string("phiNT"),
    quantity = cms.untracked.string("userFloat('phiNT')")
    ),
    cms.PSet(
    tag = cms.untracked.string("phiNT1"),
    quantity = cms.untracked.string("userFloat('phiNT1')")
    ),
    cms.PSet(
    tag = cms.untracked.string("phiNT2"),
    quantity = cms.untracked.string("userFloat('phiNT2')")
    ),
    cms.PSet(
    tag = cms.untracked.string("costhetastarNT"),
    quantity = cms.untracked.string("userFloat('costhetastarNT')")
    ),
    cms.PSet(
    tag = cms.untracked.string("HelyLD"),
    quantity = cms.untracked.string("userFloat('helyLD')")
    ),
    )
    )

### Zll variables

zll = (

    cms.PSet(
    tag = cms.untracked.string("zllMass"),
    quantity = cms.untracked.string("daughter(0).mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("zllPt"),
    quantity = cms.untracked.string("daughter(0).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("zllEta"),
    quantity = cms.untracked.string("daughter(0).eta")
    ),
    cms.PSet(
    tag = cms.untracked.string("zllY"),
    quantity = cms.untracked.string("daughter(0).y")
    ),
    cms.PSet(
    tag = cms.untracked.string("zllPhi"),
    quantity = cms.untracked.string("daughter(0).phi")
    ),
    ## pt
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Pt"),
    quantity = cms.untracked.string("daughter(0).daughter(0).pt ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Pt"),
    quantity = cms.untracked.string("daughter(0).daughter(1).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Px"),
    quantity = cms.untracked.string("daughter(0).daughter(0).px ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Px"),
    quantity = cms.untracked.string("daughter(0).daughter(1).px")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Py"),
    quantity = cms.untracked.string("daughter(0).daughter(0).py ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Py"),
    quantity = cms.untracked.string("daughter(0).daughter(1).py")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Pz"),
    quantity = cms.untracked.string("daughter(0).daughter(0).pz ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Pz"),
    quantity = cms.untracked.string("daughter(0).daughter(1).pz")
    ),
## eta
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Eta"),
    quantity = cms.untracked.string("daughter(0).daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Eta"),
    quantity = cms.untracked.string("daughter(0).daughter(1).eta ")
    ),
    ## phi
    cms.PSet(
    tag = cms.untracked.string("LeptDau1Phi"),
    quantity = cms.untracked.string("daughter(0).daughter(0).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2Phi"),
    quantity = cms.untracked.string("daughter(0).daughter(1).phi ")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(0).daughter(0).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2GlobalMuonBit"),
    quantity = cms.untracked.string("daughter(0).daughter(1).isGlobalMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1StandAloneBit"),
    quantity = cms.untracked.string("daughter(0).daughter(0).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2StandAloneBit"),
    quantity = cms.untracked.string("daughter(0).daughter(1).isStandAloneMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(0).daughter(0).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2TrackerMuonBit"),
    quantity = cms.untracked.string("daughter(0).daughter(1).isTrackerMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofMuonHits"),
    quantity = cms.untracked.string("?daughter(0).daughter(0).isGlobalMuon?daughter(0).daughter(0).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofMuonHits"),
    quantity = cms.untracked.string("?daughter(0).daughter(1).isGlobalMuon?daughter(0).daughter(1).masterClone.globalTrack.hitPattern.numberOfValidMuonHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofStripHits"),
    quantity = cms.untracked.string("?daughter(0).daughter(0).isGlobalMuon?daughter(0).daughter(0).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofStripHits"),
    quantity = cms.untracked.string("?daughter(0).daughter(1).isGlobalMuon?daughter(0).daughter(1).masterClone.globalTrack.hitPattern.numberOfValidStripHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofPixelHits"),
    quantity = cms.untracked.string("?daughter(0).daughter(0).isGlobalMuon?daughter(0).daughter(0).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofPixelHits"),
    quantity = cms.untracked.string("?daughter(0).daughter(1).isGlobalMuon?daughter(0).daughter(1).masterClone.globalTrack.hitPattern.numberOfValidPixelHits: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NormChi2"),
    quantity = cms.untracked.string("?daughter(0).daughter(0).isGlobalMuon?daughter(0).daughter(0).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NormChi2"),
    quantity = cms.untracked.string("?daughter(0).daughter(1).isGlobalMuon?daughter(0).daughter(1).masterClone.globalTrack.normalizedChi2: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1NofChambers"),
    quantity = cms.untracked.string("?daughter(0).daughter(0).isGlobalMuon?daughter(0).daughter(0).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2NofChambers"),
    quantity = cms.untracked.string("?daughter(0).daughter(1).isGlobalMuon?daughter(0).daughter(1).masterClone.numberOfChambers: -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1dB"),
    quantity = cms.untracked.string("daughter(0).daughter(0).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2dB"),
    quantity = cms.untracked.string("daughter(0).daughter(1).masterClone.dB")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1TrkIso"),
    quantity = cms.untracked.string("daughter(0).daughter(0).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2TrkIso"),
    quantity = cms.untracked.string("daughter(0).daughter(1).masterClone.trackIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau1EcalIso"),
    quantity = cms.untracked.string("daughter(0).daughter(0).masterClone.ecalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2EcalIso"),
    quantity = cms.untracked.string("daughter(0).daughter(1).masterClone.ecalIso")
    ),
      cms.PSet(
    tag = cms.untracked.string("LeptDau1HcalIso"),
    quantity = cms.untracked.string("daughter(0).daughter(0).masterClone.hcalIso")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2HcalIso"),
    quantity = cms.untracked.string("daughter(0).daughter(1).masterClone.hcalIso")
    ),
     cms.PSet(
    tag = cms.untracked.string("LeptDau1CombRelIso"),
    quantity = cms.untracked.string("(daughter(0).daughter(0).masterClone.hcalIso + daughter(0).daughter(0).masterClone.ecalIso + daughter(0).daughter(0).masterClone.trackIso )/ daughter(0).daughter(0).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("LeptDau2CombRelIso"),
    quantity = cms.untracked.string("(daughter(0).daughter(1).masterClone.hcalIso + daughter(0).daughter(1).masterClone.ecalIso + daughter(0).daughter(1).masterClone.trackIso )/ daughter(0).daughter(1).pt")
    ),

    )

zee =(

    cms.PSet(
    tag = cms.untracked.string("EleDau1VBTF80CombID"),
    quantity = cms.untracked.string("daughter(0).daughter(0).masterClone.electronID(\"eidVBTFCom80\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("EleDau2VBTF80CombID"),
    quantity = cms.untracked.string("daughter(0).daughter(1).masterClone.electronID(\"eidVBTFCom80\")")
    )
    )

zem =(

    cms.PSet(
    tag = cms.untracked.string("EleDau1VBTF80CombID"),
    quantity = cms.untracked.string("daughter(0).daughter(0).masterClone.electronID(\"eidVBTFCom80\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("EleDau2VBTF80CombID"),
    quantity = cms.untracked.string("-1")
    )    
    )

###  zjj standard and bDiscriminator variables

zjj = (
    
    cms.PSet(
    tag = cms.untracked.string("zjjMass"),
    quantity = cms.untracked.string("daughter(1).mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("zjjPt"),
    quantity = cms.untracked.string("daughter(1).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("zjjEta"),
    quantity = cms.untracked.string("daughter(1).eta")
    ),
    cms.PSet(
    tag = cms.untracked.string("zjjY"),
    quantity = cms.untracked.string("daughter(1).y")
    ),
    cms.PSet(
    tag = cms.untracked.string("zjjPhi"),
    quantity = cms.untracked.string("daughter(1).phi")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau1Pt"),
    quantity = cms.untracked.string("daughter(1).daughter(0).pt ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau2Pt"),
    quantity = cms.untracked.string("daughter(1).daughter(1).pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau1Px"),
    quantity = cms.untracked.string("daughter(1).daughter(0).px ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau2Px"),
    quantity = cms.untracked.string("daughter(1).daughter(1).px")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau1Py"),
    quantity = cms.untracked.string("daughter(1).daughter(0).py ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau2Py"),
    quantity = cms.untracked.string("daughter(1).daughter(1).py")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau1Pz"),
    quantity = cms.untracked.string("daughter(1).daughter(0).pz ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau2Pz"),
    quantity = cms.untracked.string("daughter(1).daughter(1).pz")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau1Eta"),
    quantity = cms.untracked.string("daughter(1).daughter(0).eta ")
    ),
    cms.PSet(
    tag = cms.untracked.string("JetDau2Eta"),
    quantity = cms.untracked.string("daughter(1).daughter(1).eta ")
    ),cms.PSet(
    tag = cms.untracked.string("JetDau1Phi"),
    quantity = cms.untracked.string("daughter(1).daughter(0).phi ")
    ),
cms.PSet(
    tag = cms.untracked.string("JetDau2Phi"),
    quantity = cms.untracked.string("daughter(1).daughter(1).phi ")
    ),

    ### b tagging ###

    cms.PSet(
    tag = cms.untracked.string("Jet1CSV"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2CSV"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"combinedSecondaryVertexBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1CSVMVA"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2CSVMVA"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"combinedSecondaryVertexMVABJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1JProb"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2JProb"),
    quantity = cms.untracked.string(" daughter(1).daughter(1).masterClone.bDiscriminator(\"jetProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1JbProb"),
    quantity = cms.untracked.string(" daughter(1).daughter(0).masterClone.bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2JbProb"),
    quantity = cms.untracked.string(" daughter(1).daughter(1).masterClone.bDiscriminator(\"jetBProbabilityBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1SSVHE"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2SSVHE"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"simpleSecondaryVertexHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1SSVHP"),
    quantity = cms.untracked.string(" daughter(1).daughter(0).masterClone.bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2SSVHP"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"simpleSecondaryVertexHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1ElPt"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"softElectronByPtBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2ElPt"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"softElectronByPtBJetTags\") ")
    ),
        cms.PSet(
    tag = cms.untracked.string("Jet1ElIp"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"softElectronByIP3dBJetTags\")")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2ElIp"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"softElectronByIP3dBJetTags\") ")
    ),
        cms.PSet(
    tag = cms.untracked.string("Jet1Mu"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2Mu"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"softMuonBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1MuPt"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2MuPt"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"softMuonByPtBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1MuIp"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2MuIp"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"softMuonByIP3dBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1TKHE"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2TKHE"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"trackCountingHighEffBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet1TKHP"),
    quantity = cms.untracked.string("daughter(1).daughter(0).masterClone.bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),
    cms.PSet(
    tag = cms.untracked.string("Jet2TKHP"),
    quantity = cms.untracked.string("daughter(1).daughter(1).masterClone.bDiscriminator(\"trackCountingHighPurBJetTags\") ")
    ),


##     ## jET ID
    
##     cms.PSet(
##     tag = cms.untracked.string("JetMinID"),
##     quantity = cms.untracked.string("? (daughter(0).masterClone.pt < daughter(1).masterClone.pt) && ((daughter(0).masterClone.neutralEmEnergy/daughter(0).masterClone.energy ) < 1.00) && ((daughter(0).masterClone.chargedEmEnergy/daughter(0).masterClone.energy)<1.00) && ((daughter(0).masterClone.chargedHadronEnergy/daughter(0).masterClone.energy)>0) ? true : false")
##     ),
##     cms.PSet(
##     tag = cms.untracked.string("JetMaxID"),
##     quantity = cms.untracked.string("? (daughter(0).masterClone.pt < daughter(1).masterClone.pt) && ? daughter(1).masterClone.bDiscriminator(\"trackCountingHighPurBJetTags\") :daughter(0).masterClone.bDiscriminator(\"trackCountingHighPurBJetTags\") ")
##     ),


    )


## zll.variables += baseKinematics


Higgs2mu2bEdmNtuple = copy.deepcopy(higgs)
Higgs2mu2bEdmNtuple.variables += zll
Higgs2mu2bEdmNtuple.variables += zjj
Higgs2mu2bEdmNtuple.src = cms.InputTag("hzzmmjjLD:h")
Higgs2mu2bEdmNtuple.prefix = cms.untracked.string("muHiggs")


Higgs2e2bEdmNtuple = copy.deepcopy(Higgs2mu2bEdmNtuple)
Higgs2e2bEdmNtuple.variables += zee
Higgs2e2bEdmNtuple.src = cms.InputTag("hzzeejjLD:h")
Higgs2e2bEdmNtuple.prefix = cms.untracked.string("elHiggs")

Higgsemu2bEdmNtuple = copy.deepcopy(Higgs2mu2bEdmNtuple)
Higgsemu2bEdmNtuple.variables += zem
Higgsemu2bEdmNtuple.src = cms.InputTag("hzzemjjLD:h")
Higgsemu2bEdmNtuple.prefix = cms.untracked.string("elmuHiggs")

eventVtxInfoNtuple = cms.EDProducer(
    "EventVtxInfoNtupleDumper",
    primaryVertices=cms.InputTag("offlinePrimaryVertices")
)


edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('Higgs2l2bEdmNtuples.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_eventVtxInfoNtuple_*_*",
    'keep *_kt6PFJets_rho_PAT',
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
#    "keep *_Higgsemu2bEdmNtuple_*_*"
    
    )
    )
