#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace lhef;

class LHENUPDump : public edm::EDProducer {
public:
  LHENUPDump( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag LHEsrc_;
};

LHENUPDump::LHENUPDump( const ParameterSet & cfg ) : 
  LHEsrc_("source")
{
  produces<int>( "lheNup" ).setBranchAlias( "lheNup" );
}

void LHENUPDump::produce( Event & evt, const EventSetup & ) {
  
  Handle<LHEEventProduct> LHEevt;
  evt.getByLabel( LHEsrc_, LHEevt );

  auto_ptr<int> lheNup( new int );
  *lheNup = -1;
  
  const lhef::HEPEUP hepeup_ = LHEevt->hepeup();
  *lheNup = hepeup_.NUP; 

  //  cout << "LHE NUP = " << *lheNup << endl;

  evt.put( lheNup, "lheNup" );

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( LHENUPDump );

