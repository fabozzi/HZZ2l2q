#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TVector3.h>
#include <algorithm>
#include <Math/VectorUtil.h>



class Higgs2l2bAnalysis : public edm::EDAnalyzer {
    public:
        explicit Higgs2l2bAnalysis(const edm::ParameterSet&);
        ~Higgs2l2bAnalysis();
    
    
    private:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob(); 

        edm::InputTag higgsTag;
  TH1D * l_theta, *l_cosTheta,* j_theta, *j_cosTheta, *h_phi, *h_cosPhi;

        
};

Higgs2l2bAnalysis::Higgs2l2bAnalysis(const edm::ParameterSet& iConfig):
    higgsTag(iConfig.getParameter<edm::InputTag>("higgsTag"))
{
  edm::Service<TFileService> fs;
  l_theta = fs->make<TH1D>( "Lepton_Theta", "Lepton_Theta", 50,  0, M_PI );
  l_cosTheta = fs->make<TH1D>( "Lepton_CosTheta", "Lepton_CosTheta", 50,  0., 1. );
  j_theta = fs->make<TH1D>( "Jet_Theta", "Jet_Theta", 50,  0, M_PI );
  j_cosTheta = fs->make<TH1D>( "Jet_CosTheta", "Jet_CosTheta", 50,  0., 1. );
  h_phi = fs->make<TH1D>( "AzimuthalAngle", "AzimuthalAngle", 50,  0, M_PI/2 );
  h_cosPhi = fs->make<TH1D>( "Cos2AzimuthalAngle", "Cos2AzimuthalAngle", 50,  -1, 1. );
 
}


Higgs2l2bAnalysis::~Higgs2l2bAnalysis() {
}


void Higgs2l2bAnalysis::beginJob() {

}


void Higgs2l2bAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    // Get hzzjjlls collections

    Handle<CandidateView> higgsH;
    iEvent.getByLabel(higgsTag, higgsH);
    double l0_theta, l1_theta, j0_theta, j1_theta, phi;


  for (size_t i = 0; i < higgsH->size(); ++i) {
    
    const Candidate &  h = (*higgsH)[i];
    const Candidate* zllptr = h.daughter(0);    
    const Candidate* zjjptr = h.daughter(1);
    
    //create booster
    ////CenterOfMassBooster boost_h(h);
    CenterOfMassBooster boost_ll(*zllptr);
    CenterOfMassBooster boost_jj(*zjjptr);

	
    // Boost to Higgs rest frame
    
    ////Candidate * hBoost = (h.hasMasterClone() ? &(*(h.masterClone())) : &h)->clone();
    Candidate * zllBoost = (zllptr->hasMasterClone() ? &(*(zllptr->masterClone())) : zllptr)->clone();
    Candidate * zjjBoost = (zjjptr->hasMasterClone() ? &(*(zjjptr->masterClone())) : zjjptr)->clone();

    ////boost_ll.set(*hBoost);

    ////const Candidate * z0 = const_cast<const Candidate *>(hBoost)->daughter(0);
    ////const Candidate * z1 = const_cast<const Candidate *>(hBoost)->daughter(1);

    /////phi = ROOT::Math::VectorUtil::Angle( ( ((z0->daughter(0))->momentum()).Cross((z0->daughter(1))->momentum()) ),( ((z1->daughter(0))->momentum()).Cross((z1->daughter(1))->momentum()) ) );

    ////if (phi>M_PI/2) phi = M_PI -phi;
    ///h_phi->Fill(phi);
    ////h_cosPhi->Fill(cos(2*phi));
    


    Booster hFrameBoost( h.boostToCM() );
    const Candidate * zDauRefl0 = h.daughter(0)->daughter(0);
    Candidate * boostedL0_HFrame = zDauRefl0->clone();
    hFrameBoost.set( *boostedL0_HFrame );
    const Candidate * zDauRefl1 = h.daughter(0)->daughter(1);
    Candidate * boostedL1_HFrame = zDauRefl1->clone();
    hFrameBoost.set( *boostedL1_HFrame);
    const Candidate * zDauRefj0 = h.daughter(1)->daughter(0);
    Candidate * boostedJ0_HFrame = zDauRefj0->clone();
    hFrameBoost.set( *boostedJ0_HFrame );
    const Candidate * zDauRefj1 = h.daughter(1)->daughter(1);  
    Candidate * boostedJ1_HFrame = zDauRefj1->clone();
    hFrameBoost.set( *boostedJ1_HFrame );

    phi =  ROOT::Math::VectorUtil::Angle( (boostedL0_HFrame->momentum()).Cross(boostedL1_HFrame->momentum()), (boostedJ0_HFrame->momentum()).Cross(boostedJ1_HFrame->momentum()) );


    if (phi>M_PI/2) phi = M_PI -phi;
    h_phi->Fill(phi);
    h_cosPhi->Fill(cos(2*phi));

    // Boost to Zs rest frame
    boost_ll.set(*zllBoost);
    boost_jj.set(*zjjBoost);
    
    // Leptons

    const Candidate * l0 = const_cast<const Candidate *>(zllBoost)->daughter(0);
    const Candidate * l1 = const_cast<const Candidate *>(zllBoost)->daughter(1);

    //std::cout<<l0->hasMasterClone()<<", "<<l1->hasMasterClone()<<endl;
    //std::cout<<"Before Boost: l0 eta "<<(zllptr)->daughter(0)->eta()<<"l1 eta"<<(zllptr)->daughter(1)->eta()<<endl;

    l0_theta = ROOT::Math::VectorUtil::Angle(l0->p4(), zllptr->p4());
    l1_theta = ROOT::Math::VectorUtil::Angle(l1->p4(), zllptr->p4());

    l_theta->Fill(l0_theta);
    l_theta->Fill(l1_theta);
    l_cosTheta->Fill(abs(cos(l0_theta)));

    // Jets

    const Candidate * j0 = const_cast<const Candidate *>(zjjBoost)->daughter(0);
    const Candidate * j1 = const_cast<const Candidate *>(zjjBoost)->daughter(1);

    j0_theta = ROOT::Math::VectorUtil::Angle(j0->p4(), zjjptr->p4());
    j1_theta = ROOT::Math::VectorUtil::Angle(j1->p4(), zjjptr->p4());

    j_theta->Fill(j0_theta);
    j_theta->Fill(j1_theta);
    j_cosTheta->Fill(abs(cos(j0_theta)));

      
  }
}

void Higgs2l2bAnalysis::endJob() {
}

DEFINE_FWK_MODULE(Higgs2l2bAnalysis);
