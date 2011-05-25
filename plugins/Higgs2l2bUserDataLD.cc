#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
// Likelihood discriminant
#include "HiggsAnalysis/Higgs2l2b/interface/HelicityLikelihoodDiscriminant.h"
// KinFit
#include "HiggsAnalysis/Higgs2l2b/interface/JetKinFitter.h"

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bUserDataLD : public edm::EDProducer {
public:
  Higgs2l2bUserDataLD( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  int runJetKinFit(TLorentzVector &, TLorentzVector &, 
		   const TLorentzVector &, TLorentzVector &, TLorentzVector &,
		   float &, float &);

  InputTag higgsTag;
  HelicityLikelihoodDiscriminant LD_;
  double zNominalMass_;
  JetKinFitter kinFitter_;
};



Higgs2l2bUserDataLD::Higgs2l2bUserDataLD( const ParameterSet & cfg ):
  higgsTag( cfg.getParameter<InputTag>( "higgs" ) ),
  zNominalMass_( 91.1876 ),
  kinFitter_( zNominalMass_, 0.0 )
{
  produces<vector<pat::CompositeCandidate> >("h").setBranchAlias( "h" );
}

void Higgs2l2bUserDataLD::produce( Event & evt, const EventSetup & ) {

  Handle<std::vector<pat::CompositeCandidate> > higgsH;
  evt.getByLabel(higgsTag, higgsH);

  auto_ptr<vector<pat::CompositeCandidate> > higgsColl( new vector<pat::CompositeCandidate> () );

  float helyLD;
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

  for (unsigned int i = 0; i< higgsH->size();++i){
    const pat::CompositeCandidate & H = (*higgsH)[i];
    edm::Ref<std::vector<pat::CompositeCandidate> > hRef(higgsH, i);
    pat::CompositeCandidate h(H);

    const reco::Candidate * Zll = H.daughter(0);
    const pat::Jet *  j1 = dynamic_cast<const pat::Jet*>(H.daughter(1)->daughter(0)->masterClone().get());
    const pat::Jet *  j2 = dynamic_cast<const pat::Jet*>(H.daughter(1)->daughter(1)->masterClone().get());
    
    double j1en = j1->energy(); 
    double j1pt = j1->pt();
    double j1eta = j1->eta(); 
    double j1phi = j1->phi();
    double j2en = j2->energy(); 
    double j2pt = j2->pt(); 
    double j2eta = j2->eta(); 
    double j2phi = j2->phi();
    
    j1corr.SetPtEtaPhiE(j1pt,j1eta,j1phi,j1en);
    j2corr.SetPtEtaPhiE(j2pt,j2eta,j2phi,j2en);
    
    ZLL4mom.SetPtEtaPhiM(Zll->pt(), Zll->eta(), Zll->phi(), Zll->mass());
    
    int kinfitstatus = runJetKinFit(j1corr, j2corr, ZLL4mom, Zjj4mom, HZZKinFit4mom, KFchiSquare, KFchiSquareProb);
    
    if (kinfitstatus==0) {
      j1RefitPt = j1corr.Pt();
      j2RefitPt = j2corr.Pt();
      j1RefitEta = j1corr.Eta();
      j2RefitEta = j2corr.Eta();
      j1RefitPhi = j1corr.Phi(); 
      j2RefitPhi = j2corr.Phi();
      j1RefitE = j1corr.E(); 
      j2RefitE = j1corr.E();
      ZjjRefitMass = Zjj4mom.M();
      HZZRefitMass = HZZKinFit4mom.M();
      //      cout << "### Jet Kine After KinFit: " <<endl
      //	   << "J1 ("<<j1corr.Pt()<<" , "<<j1corr.Eta()<<" , "<<j1corr.Phi()<<" , "<<j1corr.E()<<" )" <<endl
      //	   << "J2 ("<<j2corr.Pt()<<" , "<<j2corr.Eta()<<" , "<<j2corr.Phi()<<" , "<<j2corr.E()<<" )" <<endl;
      //      cout << "### XZZ mass After KinFit: "<< HZZRefitMass <<endl;
      //      cout << "### Zjj mass Before KinFit: "<< H.daughter(1)->mass() <<endl;
      //      cout << "### Zjj mass After KinFit: "<< ZjjRefitMass <<endl;
      //      cout << "### X^2/ndf and probability: "<< KFchiSquare <<" , "<< KFchiSquareProb <<endl;
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

    HelicityLikelihoodDiscriminant::HelicityAngles myha;
    myha.helCosTheta1    = H.userFloat("costhetaNT1");
    //    cout << "helCosTheta1 = " << myha.helCosTheta1 << endl;
    myha.helCosTheta2    = H.userFloat("costhetaNT2");
    //    cout << "helCosTheta2 = " << myha.helCosTheta2 << endl;
    myha.helCosThetaStar = H.userFloat("costhetastarNT");
    //    cout << "helCosThetaStar = " << myha.helCosThetaStar << endl;
    myha.helPhi          = H.userFloat("phiNT");
    //    cout << "helPhi = " << myha.helPhi << endl;
    myha.helPhi1         = H.userFloat("phiNT1");
    //    cout << "helPhi1 = " << myha.helPhi1 << endl;
    myha.mzz             = H.mass();
    //    cout << "mzz = " << myha.mzz << endl;
    if (kinfitstatus==0)
      myha.mzz = HZZRefitMass;

    LD_.setMeasurables(myha);
    float signProb = LD_.getSignalProbability();
    float bkgdProb = LD_.getBkgdProbability();
    helyLD = signProb / (signProb + bkgdProb);

    //    cout << "signProb = " << signProb << endl;
    //    cout << "bkgdProb = " << bkgdProb << endl;
    //    cout << "helyLD = " << helyLD << endl;


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

    higgsColl->push_back(h);
  }

  
  evt.put( higgsColl, "h");

}


int Higgs2l2bUserDataLD::runJetKinFit(TLorentzVector &j1,TLorentzVector &j2,
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

DEFINE_FWK_MODULE( Higgs2l2bUserDataLD );

