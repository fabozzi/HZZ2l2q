#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"

// #include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
// #include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
// #include "DataFormats/PatCandidates/interface/Isolation.h"
// #include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
// #include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"

// #include "DataFormats/BeamSpot/interface/BeamSpot.h"
// #include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "HiggsAnalysis/Higgs2l2b/interface/Helicity.h"

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
  void helicityAngles(const pat::CompositeCandidate &, const reco::Candidate *, const reco::Candidate *);
  
  InputTag higgsTag, metTag;
  PFJetIDSelectionFunctor jetIDLoose;
  pat::strbitset ret; 
  double costhetaNT1, costhetaNT2, costhetastarNT, phiNT,
    phistarNT1, phistarNT2, phistarNT12, phiNT1, phiNT2;
};



Higgs2l2bUserDataNoMC::Higgs2l2bUserDataNoMC( const ParameterSet & cfg ):
  higgsTag( cfg.getParameter<InputTag>( "higgs" ) ),
  metTag( cfg.getParameter<edm::InputTag>("metTag")),
  jetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA,
	      PFJetIDSelectionFunctor::LOOSE )
{
  ret = jetIDLoose.getBitTemplate();

  produces<vector<pat::CompositeCandidate> >("h").setBranchAlias( "h" );

}

void Higgs2l2bUserDataNoMC::produce( Event & evt, const EventSetup & ) {

  Handle<std::vector<reco::CompositeCandidate> > higgsH;
  evt.getByLabel(higgsTag, higgsH);
  
  Handle<pat::METCollection> metH;
  evt.getByLabel(metTag, metH);
  pat::METCollection met_h = *metH;

  

  auto_ptr<vector<pat::CompositeCandidate> > higgsColl( new vector<pat::CompositeCandidate> () );

  float phi; 
  float met, metSig, metPhi;
  float met2, Px, Py; 
  float zzdPhi, zzdEta, zzdr, lldPhi, lldEta,lldr, jjdPhi, jjdEta,jjdr; 
  float neutralEmEnergy, chargedEmEnergy, chargedHadronEnergy, energy;
  float jminid, jmaxid;
  float j0LooseID, j1LooseID;
  float  lminpt, lmineta, lminphi, lmaxpt, lmaxeta, lmaxphi, jminpt, jmineta, jminphi, jmaxpt,  jmaxeta,  jmaxphi;
  

  for (unsigned int i = 0; i< higgsH->size();++i){
    const reco::CompositeCandidate & H = (*higgsH)[i];
    edm::Ref<std::vector<reco::CompositeCandidate> > hRef(higgsH, i);
    pat::CompositeCandidate h(H);

    //Boost in Higgs rest Frame

    Booster hFrameBoost( H.boostToCM() );
    const Candidate * zDauRefl0 = H.daughter(0)->daughter(0);
    Candidate * boostedL0_HFrame = zDauRefl0->clone();
    const Candidate * zDauRefl1 = H.daughter(0)->daughter(1);
    Candidate * boostedL1_HFrame = zDauRefl1->clone();
    const Candidate * zDauRefj0 = H.daughter(1)->daughter(0);
    Candidate * boostedJ0_HFrame = zDauRefj0->clone();
    const pat::Jet & j0 = dynamic_cast<const pat::Jet &>(*(zDauRefj0->masterClone()));
    const Candidate * zDauRefj1 = H.daughter(1)->daughter(1);  
    const pat::Jet & j1 = dynamic_cast<const pat::Jet &>(*(zDauRefj1->masterClone()));
    Candidate * boostedJ1_HFrame = zDauRefj1->clone();

    const reco::Candidate * Zll = h.daughter(0);
    const reco::Candidate * Zjj = h.daughter(1);
    
    helicityAngles(h, Zll, Zjj);
    
    // dPhi, dEta, dr between H and Zs daughters
    
    zzdPhi = fabs(deltaPhi(H.daughter(0)->phi(),H.daughter(1)->phi() ) ) ;
    zzdEta = fabs( (H.daughter(0)->eta() - H.daughter(1)->eta()) );
    zzdr = deltaR(H.daughter(0)->eta(), H.daughter(0)->phi(), H.daughter(1)->eta(), H.daughter(1)->phi() );
    
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
    
    if (zDauRefl0->pt() < zDauRefl1->pt()) {
      lminpt  = zDauRefl0->pt();
      lmineta = zDauRefl0->eta();
      lminphi = zDauRefl0->phi();
      lmaxpt  = zDauRefl1->pt();
      lmaxeta = zDauRefl1->eta();
      lmaxphi = zDauRefl1->phi();
    }
    
    else {
      lminpt  = zDauRefl1->pt();
      lmineta = zDauRefl1->eta();
      lminphi = zDauRefl1->phi();
      lmaxpt  = zDauRefl0->pt();
      lmaxeta = zDauRefl0->eta();
      lmaxphi = zDauRefl0->phi();
    }
    
    if (zDauRefj0->pt() < zDauRefj1->pt()) {
      jminpt  = zDauRefj0->pt();
      jmineta = zDauRefj0->eta();
      jminphi = zDauRefj0->phi();
      jmaxpt  = zDauRefj1->pt();
      jmaxeta = zDauRefj1->eta();
      jmaxphi = zDauRefj1->phi();
    }
    
    else {
      jminpt  = zDauRefj1->pt();
      jmineta = zDauRefj1->eta();
      jminphi = zDauRefj1->phi();
      jmaxpt  = zDauRefj0->pt();
      jmaxeta = zDauRefj0->eta();
      jmaxphi = zDauRefj0->phi();
    }
    
    
   
   // Phi in H rest frame
    
    hFrameBoost.set( *boostedL0_HFrame );
    hFrameBoost.set( *boostedL1_HFrame);
    hFrameBoost.set( *boostedJ0_HFrame );
    hFrameBoost.set( *boostedJ1_HFrame );


    phi =  ROOT::Math::VectorUtil::Angle( (boostedL0_HFrame->momentum()).Cross(boostedL1_HFrame->momentum()), (boostedJ0_HFrame->momentum()).Cross(boostedJ1_HFrame->momentum()) );
       
    if (phi>M_PI/2) phi = M_PI -phi;

    met = met_h.front().et();
    //    cout<<"met: "<<met<<endl;
    metSig = met_h.front().mEtSig();
    metPhi = met_h.front().phi();
    Px = zDauRefl0->px() + zDauRefl1->px() + zDauRefj0->px() + zDauRefj1->px();
    Py = zDauRefl0->py() + zDauRefl1->py() + zDauRefj0->py() + zDauRefj1->py();
    met2 = sqrt(pow(Px,2)+pow(Py,2));

    h.addUserFloat("azimuthalAngle", phi);
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
    h.addUserFloat("met",met);
    h.addUserFloat("metSig",metSig);
    h.addUserFloat("metPhi",metPhi);
    h.addUserFloat("met2",met2);
    h.addUserFloat("jet1LooseID",j0LooseID);
    h.addUserFloat("jet2LooseID",j1LooseID);
    h.addUserFloat("costhetaNT1",costhetaNT1);
    h.addUserFloat("costhetaNT2",costhetaNT2);
    h.addUserFloat("phistarNT1", phistarNT1);
    h.addUserFloat("phistarNT2", phistarNT2);
    h.addUserFloat("phistarNT12", phistarNT12);  
    h.addUserFloat("phiNT",phiNT);
    h.addUserFloat("phiNT1",phiNT1);
    h.addUserFloat("phiNT2",phiNT2);
    h.addUserFloat("costhetastarNT",costhetastarNT);

    higgsColl->push_back(h);
    
  }
  
  evt.put( higgsColl, "h");

}



void Higgs2l2bUserDataNoMC::helicityAngles (const pat::CompositeCandidate &X, const reco::Candidate *Z1, const reco::Candidate *Z2) {

   TLorentzVector myp4H(X.p4().x(),X.p4().y(),X.p4().z(),X.p4().e() );
   TLorentzVector myp4Z1(Z1->p4().x(),Z1->p4().y(),Z1->p4().z(),Z1->p4().e());
   TLorentzVector myp4Z2(Z2->p4().x(),Z2->p4().y(),Z2->p4().z(),Z2->p4().e());

   TLorentzVector myp4M11(0.0,0.0,0.0,0.0);
   TLorentzVector myp4M12(0.0,0.0,0.0,0.0);
   TLorentzVector myp4M21(0.0,0.0,0.0,0.0);
   TLorentzVector myp4M22(0.0,0.0,0.0,0.0);

   //output
   costhetaNT1 = -8.80; costhetaNT2 = -8.80; costhetastarNT = -8.80;
   phistarNT1  = -8.80; phistarNT2  = -8.80; phistarNT12    = -8.80;
   phiNT       = -8.80; phiNT1      = -8.80; phiNT2         = -8.80;

   // Zll is always Z1 (in our case)
   bool Zllis1=true;

   // pick the jet with pos azim angle in the Zjj rest frame
   reco::Candidate * Zjjb = Z2->clone();
   CenterOfMassBooster boostZjj(*Zjjb);
   boostZjj.set(*Zjjb);

   int posphiZjjdau=-1;
   for(size_t dau=0;dau<2;dau++){    //Zjjb->numberOfDaughters();dau++){ //careful! there 4 daughters!
     if(Zjjb->daughter(dau)->phi()>0) posphiZjjdau = dau;
   }
   // std::cout << " ### posphiZjjdau " << posphiZjjdau << " ###" << std::endl;

   if(Zllis1){
     //set as lepton #1 the negative one
     if(Z1->daughter(0)->charge()<0.0){
        myp4M11.SetPxPyPzE(Z1->daughter(0)->p4().x(),Z1->daughter(0)->p4().y(),Z1->daughter(0)->p4().z(),Z1->daughter(0)->p4().e());
        myp4M12.SetPxPyPzE(Z1->daughter(1)->p4().x(),Z1->daughter(1)->p4().y(),Z1->daughter(1)->p4().z(),Z1->daughter(1)->p4().e());
     }
     else{
        myp4M11.SetPxPyPzE(Z1->daughter(1)->p4().x(),Z1->daughter(1)->p4().y(),Z1->daughter(1)->p4().z(),Z1->daughter(1)->p4().e());
        myp4M12.SetPxPyPzE(Z1->daughter(0)->p4().x(),Z1->daughter(0)->p4().y(),Z1->daughter(0)->p4().z(),Z1->daughter(0)->p4().e());
     }
     //set as jet #1 the pos phi one
     if( posphiZjjdau==0){
        myp4M21.SetPxPyPzE(Z2->daughter(0)->p4().x(),Z2->daughter(0)->p4().y(),Z2->daughter(0)->p4().z(),Z2->daughter(0)->p4().e());
        myp4M22.SetPxPyPzE(Z2->daughter(1)->p4().x(),Z2->daughter(1)->p4().y(),Z2->daughter(1)->p4().z(),Z2->daughter(1)->p4().e());
     }
     else{
        myp4M21.SetPxPyPzE(Z2->daughter(1)->p4().x(),Z2->daughter(1)->p4().y(),Z2->daughter(1)->p4().z(),Z2->daughter(1)->p4().e());
        myp4M22.SetPxPyPzE(Z2->daughter(0)->p4().x(),Z2->daughter(0)->p4().y(),Z2->daughter(0)->p4().z(),Z2->daughter(0)->p4().e());
     }
   }//end if Z1 == Z->ll
   else{
     //set as lepton #1 the negative one
     if(Z2->daughter(0)->charge()<0.0){
        myp4M21.SetPxPyPzE(Z2->daughter(0)->p4().x(),Z2->daughter(0)->p4().y(),Z2->daughter(0)->p4().z(),Z2->daughter(0)->p4().e());
        myp4M22.SetPxPyPzE(Z2->daughter(1)->p4().x(),Z2->daughter(1)->p4().y(),Z2->daughter(1)->p4().z(),Z2->daughter(1)->p4().e());
     }
     else{
        myp4M21.SetPxPyPzE(Z2->daughter(1)->p4().x(),Z2->daughter(1)->p4().y(),Z2->daughter(1)->p4().z(),Z2->daughter(1)->p4().e());
        myp4M22.SetPxPyPzE(Z2->daughter(0)->p4().x(),Z2->daughter(0)->p4().y(),Z2->daughter(0)->p4().z(),Z2->daughter(0)->p4().e());
     }
     //set as jet #1 the pos phi one
     if( posphiZjjdau==0 ){
        myp4M11.SetPxPyPzE(Z1->daughter(0)->p4().x(),Z1->daughter(0)->p4().y(),Z1->daughter(0)->p4().z(),Z1->daughter(0)->p4().e());
        myp4M12.SetPxPyPzE(Z1->daughter(1)->p4().x(),Z1->daughter(1)->p4().y(),Z1->daughter(1)->p4().z(),Z1->daughter(1)->p4().e());
     }
     else{
        myp4M11.SetPxPyPzE(Z1->daughter(1)->p4().x(),Z1->daughter(1)->p4().y(),Z1->daughter(1)->p4().z(),Z1->daughter(1)->p4().e());
        myp4M12.SetPxPyPzE(Z1->daughter(0)->p4().x(),Z1->daughter(0)->p4().y(),Z1->daughter(0)->p4().z(),Z1->daughter(0)->p4().e());
     }
   }

   bool swapped=false;

   Helicity dataAngles;

   dataAngles.calculateAngles(myp4H, myp4Z1, myp4M11, myp4M12, myp4Z2, myp4M21, myp4M22,
                            costhetaNT1, costhetaNT2, phiNT, costhetastarNT,
                            phistarNT1, phistarNT2, phistarNT12, phiNT1, phiNT2,swapped);

// std::cout << "Helicity angles: "<<
//              costhetaNT1 <<" , "<< costhetaNT2 <<" , "<< costhetastarNT <<" , "<<
//              phistarNT1  <<" , "<< phistarNT2  <<" , "<< phistarNT12    <<" , "<<
//              phiNT       <<" , "<< phiNT1      <<" , "<< phiNT2         << std::endl;

}





#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bUserDataNoMC );

