#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class GenInfoNtupleDump : public edm::EDProducer {
public:
  GenInfoNtupleDump( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag genpTag_;
};

GenInfoNtupleDump::GenInfoNtupleDump( const ParameterSet & cfg ) : 
  genpTag_(cfg.getParameter<edm::InputTag>("src"))
{
  produces<float>( "genHiggsMass" ).setBranchAlias( "genHiggsMass" );
}

void GenInfoNtupleDump::produce( Event & evt, const EventSetup & ) {
  
  edm::Handle<GenParticleCollection> genParticles;
  evt.getByLabel(genpTag_, genParticles);

  auto_ptr<float> hmass_( new float );
  *hmass_ = -1.0;

  unsigned int gensize = genParticles->size();
  for (unsigned int i = 0; i<gensize; ++i) {
    if ((*genParticles)[i].pdgId() == 25){ // found the Higgs
      *hmass_ = (*genParticles)[i].mass();
      break;
    }
  }

  evt.put( hmass_, "genHiggsMass" );
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( GenInfoNtupleDump );

