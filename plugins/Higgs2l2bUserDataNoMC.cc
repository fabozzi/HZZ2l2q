#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"

// PF candidates
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Helicity angles 
#include "HiggsAnalysis/Higgs2l2b/interface/Helicity.h"
// KinFit
#include "HiggsAnalysis/Higgs2l2b/interface/JetKinFitter.h"
// Likelihood discriminant
#include "HiggsAnalysis/Higgs2l2b/interface/HelicityLikelihoodDiscriminant.h"

#include <Math/VectorUtil.h>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bUserDataNoMC : public edm::EDProducer {
public:
  Higgs2l2bUserDataNoMC( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void helicityAngles(const reco::Candidate *, const reco::Candidate *);
  void helicityAnglesRefit(const reco::Candidate *, TLorentzVector, TLorentzVector);
  int runJetKinFit(TLorentzVector &, TLorentzVector &, 
  		   const TLorentzVector &, TLorentzVector &, TLorentzVector &,
  		   float &, float &);
  void getLDVariables(float, float, float, float, float, float, 
		      float &, float &, float &);  
 
  InputTag higgsTag, metTag, PFCandTag, vtxTag;
  double deltaZCut_;
  PFJetIDSelectionFunctor jetIDLoose;
  pat::strbitset ret; 
  // angular variables with unrefitted jet quantities
  double costhetaNT1, costhetaNT2, costhetastarNT, phiNT, phiNT1;
  // angular variables with refitted jet quantities
  double costhetaNT1Refit, costhetaNT2Refit, costhetastarNTRefit, phiNTRefit, phiNT1Refit;
  double zNominalMass_;
  JetKinFitter kinFitter_;
};

Higgs2l2bUserDataNoMC::Higgs2l2bUserDataNoMC( const ParameterSet & cfg ):
  higgsTag( cfg.getParameter<InputTag>( "higgs" ) ),
  PFCandTag( cfg.getParameter<edm::InputTag>("PFCandidates")),
  vtxTag(cfg.getParameter<InputTag>("primaryVertices")),
  deltaZCut_( cfg.getParameter<double>("dzCut")),
  jetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA,
	      PFJetIDSelectionFunctor::LOOSE ),
  zNominalMass_( 91.1876 ),
  kinFitter_( zNominalMass_, 0.0 )
{
  ret = jetIDLoose.getBitTemplate();
  produces<vector<pat::CompositeCandidate> >("h").setBranchAlias( "h" );
}

void Higgs2l2bUserDataNoMC::produce( Event & evt, const EventSetup & ) {

  Handle<std::vector<reco::CompositeCandidate> > higgsH;
  evt.getByLabel(higgsTag, higgsH);

  // get PFCandidates
  Handle<PFCandidateCollection> pfCandidates;
  evt.getByLabel(PFCandTag, pfCandidates);

  // get Primary vtx
  //  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
  //  evt.getByLabel(vtxTag, primaryVertices);
  //  const reco::Vertex &pv = (*primaryVertices)[0];

  auto_ptr<vector<pat::CompositeCandidate> > higgsColl( new vector<pat::CompositeCandidate> () );

  float zzdPhi, zzdEta, zzdr, lldPhi, lldEta,lldr, jjdPhi, jjdEta, jjdr; 
  float neutralEmEnergy, chargedEmEnergy, chargedHadronEnergy, energy;
  float jminid, jmaxid;
  float j0LooseID, j1LooseID;
  //  float  lminpt, lmineta, lminphi, lmaxpt, lmaxeta, lmaxphi, jminpt, jmineta, jminphi, jmaxpt,  jmaxeta,  jmaxphi;
  float j1RefitPt, j2RefitPt;
  float j1RefitEta, j2RefitEta;
  float j1RefitPhi, j2RefitPhi;
  float j1RefitE, j2RefitE;
  float ZjjRefitMass;
  float HZZRefitMass;
  float KFchiSquare, KFchiSquareProb;
  TLorentzVector j1corr;
  TLorentzVector j2corr;
  TLorentzVector HZZKinFit4mom, ZLL4mom, Zjj4mom; //initialized to (0, 0, 0 ,0)
  float helyLD;
  float ldSig, ldBkg;
  float helyLDRefit;
  float ldSigRefit, ldBkgRefit;
  //  float trkMetX, trkMetY, trkMet, trkPlusNeuMet;
  //  float neutralContributionX, neutralContributionY;
  //  float trkCorrectedMetX, trkCorrectedMetY;

  for (unsigned int i = 0; i< higgsH->size();++i){
    const reco::CompositeCandidate & H = (*higgsH)[i];
    edm::Ref<std::vector<reco::CompositeCandidate> > hRef(higgsH, i);
    pat::CompositeCandidate h(H);

    const Candidate * zDauRefl0 = H.daughter(0)->daughter(0);
    const Candidate * zDauRefl1 = H.daughter(0)->daughter(1);
    const Candidate * zDauRefj0 = H.daughter(1)->daughter(0);
    const Candidate * zDauRefj1 = H.daughter(1)->daughter(1);  
    const pat::Jet & j0 = dynamic_cast<const pat::Jet &>(*(zDauRefj0->masterClone()));
    const pat::Jet & j1 = dynamic_cast<const pat::Jet &>(*(zDauRefj1->masterClone()));

    const reco::Candidate * Zll = h.daughter(0);
    const reco::Candidate * Zjj = h.daughter(1);
    
    // compute helicity angles with unrefitted jet quantities
    helicityAngles( Zll, Zjj );

    // dPhi, dEta, dr between H and Zs daughters
    zzdPhi = fabs( deltaPhi(Zll->phi(), Zjj->phi()) ) ;
    zzdEta = fabs( (Zll->eta() - Zjj->eta()) );
    zzdr = deltaR(Zll->eta(), Zll->phi(), Zjj->eta(), Zjj->phi() );
    
    lldPhi = fabs(deltaPhi(zDauRefl0->phi(),zDauRefl1->phi() ) ) ;
    lldEta = fabs(zDauRefl0->eta() - zDauRefl1->eta());
    lldr = deltaR(zDauRefl0->eta(), zDauRefl0->phi(), zDauRefl1->eta(), zDauRefl1->phi() );
      
    jjdPhi = fabs(deltaPhi(zDauRefj0->phi(),zDauRefj1->phi() ) ) ;
    jjdEta = fabs(zDauRefj0->eta() - zDauRefj1->eta());
    jjdr = deltaR(zDauRefj0->eta(), zDauRefj0->phi(), zDauRefj1->eta(), zDauRefj1->phi() );

    // store jetID for Z->jj daughters
    j0LooseID = (float) jetIDLoose( j0, ret );
    j1LooseID = (float) jetIDLoose( j1, ret );
   
    if(j0.pt() < j1.pt() ){      
      neutralEmEnergy = j0.neutralEmEnergy();
      chargedEmEnergy = j0.chargedEmEnergy() ;
      chargedHadronEnergy =j0.chargedHadronEnergy() ;
      energy = j0.energy() ;
      if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jminid = true;
      else jminid = false;
  
      neutralEmEnergy = j1.neutralEmEnergy();
      chargedEmEnergy = j1.chargedEmEnergy() ;
      chargedHadronEnergy =j1.chargedHadronEnergy() ;
      energy = j1.energy() ;
      if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jmaxid = true;
      else jmaxid = false;
    }
    else{
      neutralEmEnergy = j1.neutralEmEnergy();
      chargedEmEnergy = j1.chargedEmEnergy() ;
      chargedHadronEnergy = j1.chargedHadronEnergy() ;
      energy = j1.energy() ;
      if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jminid = true;
      else jminid = false;

      neutralEmEnergy = j0.neutralEmEnergy();
      chargedEmEnergy = j0.chargedEmEnergy() ;
      chargedHadronEnergy = j0.chargedHadronEnergy() ;
      energy = j0.energy() ;
      if (neutralEmEnergy/energy < 1.00 && chargedEmEnergy/energy < 1.00 && chargedHadronEnergy/energy > 0) jmaxid = true;
      else jmaxid = false;
    }

    // gen info    
    //    if (zDauRefl0->pt() < zDauRefl1->pt()) {
    //      lminpt  = zDauRefl0->pt();
    //      lmineta = zDauRefl0->eta();
    //      lminphi = zDauRefl0->phi();
    //      lmaxpt  = zDauRefl1->pt();
    //      lmaxeta = zDauRefl1->eta();
    //      lmaxphi = zDauRefl1->phi();
    //    }
    //    
    //    else {
    //      lminpt  = zDauRefl1->pt();
    //      lmineta = zDauRefl1->eta();
    //      lminphi = zDauRefl1->phi();
    //      lmaxpt  = zDauRefl0->pt();
    //      lmaxeta = zDauRefl0->eta();
    //      lmaxphi = zDauRefl0->phi();
    //    }
    
    //    if (zDauRefj0->pt() < zDauRefj1->pt()) {
    //      jminpt  = zDauRefj0->pt();
    //      jmineta = zDauRefj0->eta();
    //      jminphi = zDauRefj0->phi();
    //      jmaxpt  = zDauRefj1->pt();
    //      jmaxeta = zDauRefj1->eta();
    //      jmaxphi = zDauRefj1->phi();
    //    }
    //    
    //    else {
    //      jminpt  = zDauRefj1->pt();
    //      jmineta = zDauRefj1->eta();
    //      jminphi = zDauRefj1->phi();
    //      jmaxpt  = zDauRefj0->pt();
    //      jmaxeta = zDauRefj0->eta();
    //      jmaxphi = zDauRefj0->phi();
    //    }
    
    // prepare input for KinFit
    double j1en = j0.energy(); 
    double j1pt = j0.pt();
    double j1eta = j0.eta(); 
    double j1phi = j0.phi();
    double j2en = j1.energy(); 
    double j2pt = j1.pt(); 
    double j2eta = j1.eta(); 
    double j2phi = j1.phi();
    j1corr.SetPtEtaPhiE(j1pt,j1eta,j1phi,j1en);
    j2corr.SetPtEtaPhiE(j2pt,j2eta,j2phi,j2en);
    ZLL4mom.SetPtEtaPhiM(Zll->pt(), Zll->eta(), Zll->phi(), Zll->mass());

    // run KinFit
    int kinfitstatus = runJetKinFit(j1corr, j2corr, ZLL4mom, Zjj4mom, HZZKinFit4mom, KFchiSquare, KFchiSquareProb);

    if (kinfitstatus==0) {
      j1RefitPt = j1corr.Pt();
      j2RefitPt = j2corr.Pt();
      j1RefitEta = j1corr.Eta();
      j2RefitEta = j2corr.Eta();
      j1RefitPhi = j1corr.Phi(); 
      j2RefitPhi = j2corr.Phi();
      j1RefitE = j1corr.E(); 
      j2RefitE = j2corr.E();
      ZjjRefitMass = Zjj4mom.M();
      HZZRefitMass = HZZKinFit4mom.M();
    } else {
      //kinematic fit failed
      j1RefitPt = 0;
      j2RefitPt = 0;
      j1RefitEta = 0;
      j2RefitEta = 0;
      j1RefitPhi = 0; 
      j2RefitPhi = 0;
      j1RefitE = 0; 
      j2RefitE = 0;
      ZjjRefitMass = 0;
      HZZRefitMass = 0;
      KFchiSquare = -1. ; 
      KFchiSquareProb = -1.;
    }

    // compute helicity angles with refitted jet quantities
    helicityAnglesRefit( Zll, j1corr, j2corr);

    // Get LD variable
    getLDVariables( costhetaNT1, costhetaNT2, costhetastarNT,
		    phiNT, phiNT1, H.mass(), 
		    ldSig, ldBkg, helyLD );
    
    if (kinfitstatus==0)
      getLDVariables( costhetaNT1Refit, costhetaNT2Refit, costhetastarNTRefit,
		      phiNTRefit, phiNT1Refit, HZZRefitMass, 
		      ldSigRefit, ldBkgRefit, helyLDRefit );
    else{
      ldSigRefit = -100.0;
      ldBkgRefit = -100.0;
      helyLDRefit = -100.0;
    }


    h.addUserFloat("zzdPhi", zzdPhi);
    h.addUserFloat("zzdEta", zzdEta);
    h.addUserFloat("zzdr", zzdr);
    h.addUserFloat("lldPhi", lldPhi);
    h.addUserFloat("lldEta", lldEta);
    h.addUserFloat("lldr", lldr);
    h.addUserFloat("jjdPhi", jjdPhi);
    h.addUserFloat("jjdEta", jjdEta);
    h.addUserFloat("jjdr", jjdr);
    h.addUserFloat("jminid",jminid);
    h.addUserFloat("jmaxid",jmaxid);
    h.addUserFloat("jet1LooseID",j0LooseID);
    h.addUserFloat("jet2LooseID",j1LooseID);
    h.addUserFloat("costhetaNT1",costhetaNT1);
    h.addUserFloat("costhetaNT2",costhetaNT2);
    h.addUserFloat("phiNT",phiNT);
    h.addUserFloat("phiNT1",phiNT1);
    h.addUserFloat("costhetastarNT",costhetastarNT);
    h.addUserFloat("costhetaNT1Refit",costhetaNT1Refit);
    h.addUserFloat("costhetaNT2Refit",costhetaNT2Refit);
    h.addUserFloat("phiNTRefit",phiNTRefit);
    h.addUserFloat("phiNT1Refit",phiNT1Refit);
    h.addUserFloat("costhetastarNTRefit",costhetastarNTRefit);
    h.addUserFloat("j1RefitPt", j1RefitPt);
    h.addUserFloat("j2RefitPt", j2RefitPt);
    h.addUserFloat("j1RefitEta", j1RefitEta);
    h.addUserFloat("j2RefitEta", j2RefitEta);
    h.addUserFloat("j1RefitPhi", j1RefitPhi);
    h.addUserFloat("j2RefitPhi", j2RefitPhi);
    h.addUserFloat("j1RefitE", j1RefitE);
    h.addUserFloat("j2RefitE", j2RefitE);
    h.addUserFloat("ZjjRefitMass", ZjjRefitMass);
    h.addUserFloat("HZZRefitMass", HZZRefitMass);
    h.addUserFloat("KFchiSquare", KFchiSquare);
    h.addUserFloat("KFchiSquareProb", KFchiSquareProb);
    h.addUserFloat("helyLD", helyLD);
    h.addUserFloat("ldSig", ldSig);
    h.addUserFloat("ldBkg", ldBkg);
    h.addUserFloat("helyLDRefit", helyLDRefit);
    h.addUserFloat("ldSigRefit", ldSigRefit);
    h.addUserFloat("ldBkgRefit", ldBkgRefit);
    // new met variables
    //    h.addUserFloat("trkMetX", trkMetX);
    //    h.addUserFloat("trkMetY", trkMetY);
    //    h.addUserFloat("trkCorrectedMetX", trkCorrectedMetX);
    //    h.addUserFloat("trkCorrectedMetY", trkCorrectedMetY);
    //    h.addUserFloat("trkMet", trkMet);
    //    h.addUserFloat("trkPlusNeuMet", trkPlusNeuMet);

    higgsColl->push_back(h);
  }

  
  evt.put( higgsColl, "h");

}


void Higgs2l2bUserDataNoMC::getLDVariables( float costhetaNT1, float costhetaNT2, float costhetastarNT,
					float phiNT, float phiNT1, float Hmass, 
					float & ldSig, float & ldBkg, float & helyLD ){
  HelicityLikelihoodDiscriminant LD_;
  HelicityLikelihoodDiscriminant::HelicityAngles myha;
  ldSig = -100.0; ldBkg=-100.0; helyLD=-100.0;
  myha.helCosTheta1    = costhetaNT1;
  myha.helCosTheta2    = costhetaNT2;
  myha.helCosThetaStar = costhetastarNT;
  myha.helPhi          = phiNT;
  myha.helPhi1         = phiNT1;
  myha.mzz             = Hmass;
  LD_.setMeasurables(myha);
  ldSig = LD_.getSignalProbability();
  ldBkg = LD_.getBkgdProbability();
  helyLD = ldSig / (ldSig + ldBkg);
}

void Higgs2l2bUserDataNoMC::helicityAngles (const reco::Candidate *Zll, const reco::Candidate *Zjj) {

    // prepare for helicity angles computation
    costhetaNT1 = -8.80; costhetaNT2 = -8.80; costhetastarNT = -8.80;
    phiNT       = -8.80; phiNT1      = -8.80;
    TLorentzVector p4lept1(0.0,0.0,0.0,0.0);
    TLorentzVector p4lept2(0.0,0.0,0.0,0.0);
    TLorentzVector p4jet1(0.0,0.0,0.0,0.0);
    TLorentzVector p4jet2(0.0,0.0,0.0,0.0);
    //set as lepton #1 the negative one
    int lM, lP;
    if(Zll->daughter(0)->charge()<0.0) 
      lM = 0;
    else   
      lM = 1;
    lP = 1-lM;

    p4lept1.SetPxPyPzE(Zll->daughter(lM)->p4().x(),Zll->daughter(lM)->p4().y(),Zll->daughter(lM)->p4().z(),Zll->daughter(lM)->p4().e());
    p4lept2.SetPxPyPzE(Zll->daughter(lP)->p4().x(),Zll->daughter(lP)->p4().y(),Zll->daughter(lP)->p4().z(),Zll->daughter(lP)->p4().e());

    p4jet1.SetPxPyPzE(Zjj->daughter(0)->p4().x(),Zjj->daughter(0)->p4().y(),Zjj->daughter(0)->p4().z(),Zjj->daughter(0)->p4().e());
    p4jet2.SetPxPyPzE(Zjj->daughter(1)->p4().x(),Zjj->daughter(1)->p4().y(),Zjj->daughter(1)->p4().z(),Zjj->daughter(1)->p4().e());
    //compute helicity angles
   Helicity myAngles;
   myAngles.calculateAngles(p4lept1, p4lept2, p4jet1, p4jet2, costhetaNT1, costhetaNT2, costhetastarNT, phiNT, phiNT1);

}


void Higgs2l2bUserDataNoMC::helicityAnglesRefit(const reco::Candidate *Zll, TLorentzVector p4jet1, TLorentzVector p4jet2) 
{  
  // prepare for helicity angles computation
  costhetaNT1Refit = -8.80; costhetaNT2Refit = -8.80; costhetastarNTRefit = -8.80;
  phiNTRefit       = -8.80; phiNT1Refit      = -8.80;
  TLorentzVector p4lept1(0.0,0.0,0.0,0.0);
  TLorentzVector p4lept2(0.0,0.0,0.0,0.0);
  //set as lepton #1 the negative one
  int lM, lP;
  if(Zll->daughter(0)->charge()<0.0) 
    lM = 0;
  else   
    lM = 1;
  lP = 1-lM;
  
  p4lept1.SetPxPyPzE(Zll->daughter(lM)->p4().x(),Zll->daughter(lM)->p4().y(),Zll->daughter(lM)->p4().z(),Zll->daughter(lM)->p4().e());
  p4lept2.SetPxPyPzE(Zll->daughter(lP)->p4().x(),Zll->daughter(lP)->p4().y(),Zll->daughter(lP)->p4().z(),Zll->daughter(lP)->p4().e());
  
  //compute helicity angles
  Helicity myAngles;
  myAngles.calculateAngles(p4lept1, p4lept2, p4jet1, p4jet2, costhetaNT1Refit, costhetaNT2Refit, costhetastarNTRefit, phiNTRefit, phiNT1Refit);
}

int Higgs2l2bUserDataNoMC::runJetKinFit(TLorentzVector &j1,TLorentzVector &j2,
				    const TLorentzVector &ZLL, TLorentzVector & Zjj, 
				    TLorentzVector &XZZ, float & chiSquare, 
				    float & chiSquareProb) {
  
  int status=0;

  //pass the two four momenta and initialize the Kinfit object
  kinFitter_.setJet4Mom(j1,j2);

  //ask the kinfit object to correct the four momenta
  status=kinFitter_.Refit();
  if(status==0){
    j1=kinFitter_.getCorrJets().at(0);
    j2=kinFitter_.getCorrJets().at(1);
    chiSquare     = kinFitter_.chiSquare();
    chiSquareProb = kinFitter_.chiSquareProb();
  }

  //update also 4-mom of XZZ and of ZJJ
  Zjj = j1+j2;
  XZZ = ZLL+Zjj;

  return status;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bUserDataNoMC );


