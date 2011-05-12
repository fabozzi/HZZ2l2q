#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include "DataFormats/Math/interface/deltaR.h"



class H2l2bFilter : public edm::EDFilter {
public:
  explicit H2l2bFilter(const edm::ParameterSet&);
  ~H2l2bFilter();
  
private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  edm::InputTag hzzeejj_;
  edm::InputTag hzzmmjj_;
  double zLepPtCut_;
  double zJetPtCut_;
  double zMuRelIsoCut_;
  double zLepMassCut_;
  double zJetMassCut_;
  double zllPtCut_;
  double metCut_;
  double jjdrCut_;
  double zMudBCut_;
  
};


H2l2bFilter::H2l2bFilter(const edm::ParameterSet& iConfig):
  hzzeejj_(iConfig.getParameter<edm::InputTag>("hzzeejjTag")),
  hzzmmjj_(iConfig.getParameter<edm::InputTag>("hzzmmjjTag")),
  zLepPtCut_(iConfig.getParameter<double>("zLepPtCut")),
  zJetPtCut_(iConfig.getParameter<double>("zJetPtCut")),
  zMuRelIsoCut_(iConfig.getParameter<double>("zLepRelIsoCut")),
  zLepMassCut_(iConfig.getParameter<double>("zLepMassCut")),
  zJetMassCut_(iConfig.getParameter<double>("zJetMassCut")),
  zllPtCut_(iConfig.getParameter<double>("zllPtCut")),
  metCut_(iConfig.getParameter<double>("metCut")),
  jjdrCut_(iConfig.getParameter<double>("jjdrCut")),
  zMudBCut_(iConfig.getParameter<double>("zMudBCut"))
{
}


H2l2bFilter::~H2l2bFilter() {
}

bool H2l2bFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    float zNominalMass = 91.19;

     // hzzmmjj candidates
    Handle<std::vector<pat::CompositeCandidate> > hzzmmjj;
    iEvent.getByLabel(hzzmmjj_, hzzmmjj);
    for(unsigned i=0; i<hzzmmjj->size(); ++i){
      const pat::CompositeCandidate & h = (*hzzmmjj)[i];
      const Candidate * zll = h.daughter(0);    
      const Candidate * zjj = h.daughter(1);
      const Candidate * zDauRefl0 = h.daughter(0)->daughter(0);
      const Candidate * zDauRefl1 = h.daughter(0)->daughter(1);
      const Candidate * zDauRefj0 = h.daughter(1)->daughter(0);
      const Candidate * zDauRefj1 = h.daughter(1)->daughter(1);  
      const pat::Muon * lept0mu = dynamic_cast<const pat::Muon *>(zDauRefl0->masterClone().get());
      const pat::Muon * lept1mu = dynamic_cast<const pat::Muon *>(zDauRefl1->masterClone().get());
      //const pat::Jet & j0 = dynamic_cast<const pat::Jet &>(*zDauRefj0->masterClone());
      //const pat::Jet & j1 = dynamic_cast<const pat::Jet &>(*zDauRefj1->masterClone());
      float met = h.userFloat("met");
      float jjdr = deltaR(zDauRefj0->eta(), zDauRefj0->phi(), zDauRefj1->eta(), zDauRefj1->phi() );
      bool dB = fabs(lept0mu->dB()) < zMudBCut_ && fabs(lept1mu->dB()) < zMudBCut_;
      float RelIso0 = (lept0mu->ecalIso()+lept0mu->hcalIso()+lept0mu->trackIso())/lept0mu->pt();
      float RelIso1 = (lept1mu->ecalIso()+lept1mu->hcalIso()+lept1mu->trackIso())/lept1mu->pt();
      bool Iso = RelIso0<zMuRelIsoCut_ && RelIso1<zMuRelIsoCut_;

      if( zDauRefl0->pt()>zLepPtCut_ && zDauRefl1->pt()>zLepPtCut_ && zDauRefj0->pt()>zJetPtCut_ && zDauRefj1->pt()>zJetPtCut_){
	if (dB && Iso){
	  if(fabs(zll->mass() - zNominalMass)< zLepMassCut_ && fabs(zjj->mass() - zNominalMass)< zJetMassCut_){
	    if( zll->pt() > zllPtCut_ ) {
	      if( jjdr < jjdrCut_ ){
		if( met < metCut_) {
		  return true;
		}
	      }
	    }
	  }
	}
      }
    }

    // hzzeejj candidates
    Handle<std::vector<pat::CompositeCandidate> > hzzeejj;
    iEvent.getByLabel(hzzeejj_, hzzeejj);
    for(unsigned i=0; i<hzzeejj->size(); ++i){
      const pat::CompositeCandidate & h = (*hzzeejj)[i];
      const Candidate * zll = h.daughter(0);    
      const Candidate * zjj = h.daughter(1);
      const Candidate * zDauRefl0 = h.daughter(0)->daughter(0);
      const Candidate * zDauRefl1 = h.daughter(0)->daughter(1);
      const Candidate * zDauRefj0 = h.daughter(1)->daughter(0);
      const Candidate * zDauRefj1 = h.daughter(1)->daughter(1);  
      const pat::Electron * lept0el = dynamic_cast<const pat::Electron *>(zDauRefl0->masterClone().get());
      const pat::Electron * lept1el = dynamic_cast<const pat::Electron *>(zDauRefl1->masterClone().get());
      //const pat::Jet & j0 = dynamic_cast<const pat::Jet &>(*zDauRefj0->masterClone());
      //const pat::Jet & j1 = dynamic_cast<const pat::Jet &>(*zDauRefj1->masterClone());
      float met = h.userFloat("met");
      float jjdr = deltaR(zDauRefj0->eta(), zDauRefj0->phi(), zDauRefj1->eta(), zDauRefj1->phi() );     
      bool VBTF80CombID = lept0el->electronID("eidVBTFCom80")==7 || lept1el->electronID("eidVBTFCom80")==7;
 
      if( zDauRefl0->pt()>zLepPtCut_ && zDauRefl1->pt()>zLepPtCut_ && zDauRefj0->pt()>zJetPtCut_ && zDauRefj1->pt()>zJetPtCut_){
	if (VBTF80CombID){
	  if(fabs(zll->mass() - zNominalMass)< zLepMassCut_ && fabs(zjj->mass() - zNominalMass)< zJetMassCut_){
	    if( zll->pt() > zllPtCut_ ) {
	      if( jjdr < jjdrCut_ ){ 
		if( met < metCut_) {
		  return true; 
		}
	      }
	    }
	  }
	}
      }
    }

    return false; 
}

void H2l2bFilter::beginJob() {
}

void H2l2bFilter::endJob() {
}

DEFINE_FWK_MODULE(H2l2bFilter);

