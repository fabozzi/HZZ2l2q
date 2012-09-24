#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <Math/VectorUtil.h>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

class HiggsMatcher : public edm::EDProducer {
public:
  HiggsMatcher( const edm::ParameterSet & );   
  // match type (0-> matching jets; 1-> matching leptons with a certain absId)
  bool verifyMatch( int matchType, const Candidate & H, 
		    Handle<GenParticleCollection> genColl, float drCut, int absId = 11);

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag elHiggsTag_, muHiggsTag_;
  InputTag genLeptonsTag_, genQuarksTag_;
};

HiggsMatcher::HiggsMatcher( const ParameterSet & cfg ):
  elHiggsTag_( cfg.getParameter<InputTag>( "elhiggs" ) ),
  muHiggsTag_( cfg.getParameter<InputTag>( "muhiggs" ) ),
  genLeptonsTag_( cfg.getParameter<InputTag>("genLept")),
  genQuarksTag_( cfg.getParameter<InputTag>("genQuarks"))
{
  produces< vector<unsigned int> >("elHiggsMatch").setBranchAlias( "elHiggsMatch" );
  produces< vector<unsigned int> >("muHiggsMatch").setBranchAlias( "muHiggsMatch" );
}

void HiggsMatcher::produce( Event & evt, const EventSetup & ) {

  Handle<CandidateView> elH;
  evt.getByLabel(elHiggsTag_, elH);

  Handle<CandidateView>  muH;
  evt.getByLabel(muHiggsTag_, muH);

  Handle<GenParticleCollection> genL;
  evt.getByLabel(genLeptonsTag_, genL);

  Handle<GenParticleCollection> genQ;
  evt.getByLabel(genQuarksTag_, genQ);


  auto_ptr< vector<unsigned int> > elHMatch( new vector<unsigned int> () );
  auto_ptr< vector<unsigned int> > muHMatch( new vector<unsigned int> () );

  unsigned int elHSize = elH->size();
  unsigned int muHSize = muH->size();

  // looping on ElHiggs candidates
  for (unsigned int i = 0; i < elHSize;++i){
    cout << "elH candidate = " << i << endl;

    const Candidate & H = (*elH)[i];
    // matching electrons
    bool isMatchedEl = verifyMatch(1, H, genL, 0.3, 11) ;
    cout << "matchedEL??? = " << isMatchedEl << endl;
    // matching jets
    bool isMatchedJet = verifyMatch(0, H, genQ, 0.5) ;
    cout << "ismatchedJet??? = " << isMatchedJet << endl;

    /*

    for(unsigned int lind = 0; lind<2; ++lind) {
      const Candidate * zDauReflep = H.daughter(0)->daughter(lind);    
      for (unsigned int j = 0; j < genL->size(); ++j){
	const GenParticle & lept = (*genL)[j];
	cout << "GEN LEPT " << j << endl; 
	// check if H->ee
	if(fabs(lept.pdgId()) != 11) 
	  break;
	float dr = deltaR(lept.eta(), lept.phi(), zDauReflep->eta(), zDauReflep->phi() );
	if( dr < 0.3 ) {
	  isMatchedLep = true;
	  break;
	} else {
	  isMatchedLep = false;
	} 
      }

      if( !isMatchedLep )
	break;
    }

    // matching jets
    bool isMatchedJet = true;
    for(unsigned int jind = 0; jind<2; ++jind) {
      const Candidate * zDauRefj = H.daughter(1)->daughter(jind);
      for (unsigned int j = 0; j < genQ->size(); ++j){
	const GenParticle & qdau = (*genQ)[j];
	cout << "GEN QUARK " << j << endl; 
	float dr = deltaR(qdau.eta(), qdau.phi(), zDauRefj->eta(), zDauRefj->phi() );
	if( dr < 0.3 ) {
	  isMatchedJet = true;
	  break;
	} else {
	  isMatchedJet = false;
	} 
      }

      if(isMatchedJet = false)
	break;
    }

    */


    if(isMatchedEl && isMatchedJet)
      (*elHMatch).push_back(1);
    else
      (*elHMatch).push_back(0);

  }

  // looping on MuHiggs candidates
  for (unsigned int i = 0; i < muHSize;++i){

    cout << "muH candidate = " << i << endl;

    const Candidate & H = (*muH)[i];
    // matching muons
    bool isMatchedMu = verifyMatch(1, H, genL, 0.3, 13) ;
    cout << "matchedMu??? = " << isMatchedMu << endl;
    // matching jets
    bool isMatchedJet = verifyMatch(0, H, genQ, 0.5) ;
    cout << "ismatchedJet??? = " << isMatchedJet << endl;

    if(isMatchedMu && isMatchedJet)
      (*muHMatch).push_back(1);
    else
      (*muHMatch).push_back(0);
    
  }


  evt.put( elHMatch, "elHiggsMatch");
  evt.put( muHMatch, "muHiggsMatch");

}

// match type (0->jet, 1->lept)
bool HiggsMatcher::verifyMatch( int matchType, const Candidate & H, 
				Handle<GenParticleCollection> genColl, float drCut, int absId ) {
  bool isMatchedCand(true);

  for(unsigned int lind = 0; lind<2; ++lind) {
    const Candidate * zDauRef;
    if(matchType == 1)
      zDauRef = H.daughter(0)->daughter(lind);
    else
      zDauRef = H.daughter(1)->daughter(lind);
    
    for (unsigned int j = 0; j < genColl->size(); ++j){
      const GenParticle & genLeg = (*genColl)[j];
      if(matchType == 1){
	// for leptons check if electron or muon
	if(fabs(genLeg.pdgId()) != absId) 
	  break;
      } 
      float dr = deltaR(genLeg.eta(), genLeg.phi(), zDauRef->eta(), zDauRef->phi() );
      cout << "dr  = " << dr << endl;
      if( dr < drCut ) {
	cout << "MATCHED " << endl;
	isMatchedCand = true;
	break;
      } else {
	isMatchedCand = false;
	cout << "UNMATCHED " << endl;
      } 
    }
    
    if( !isMatchedCand )
      break;
  }
  
  return isMatchedCand;

}




#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HiggsMatcher );


