#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

using namespace edm;
using namespace std;

class EventRhoDumper : public edm::EDProducer {
public:
  EventRhoDumper( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag rhoAll_;
  edm::InputTag rhoRestricted_;

};

EventRhoDumper::EventRhoDumper( const ParameterSet & cfg ) : 
  rhoAll_(cfg.getParameter<InputTag>("rho")),
  rhoRestricted_(cfg.getParameter<InputTag>("restrictedRho")) {
  produces<float>( "rhoAllEta" ).setBranchAlias( "rhoAllEta" );
  produces<float>( "rhoRestrictedEta" ).setBranchAlias( "rhoRestrictedEta" );
}

void EventRhoDumper::produce( Event & evt, const EventSetup & ) {
  edm::Handle<double> rhoAll;
  evt.getByLabel(rhoAll_, rhoAll);
  edm::Handle<double> rhoRestricted;
  evt.getByLabel(rhoRestricted_, rhoRestricted);


  auto_ptr<float> rhoVariable( new float );
  auto_ptr<float> rhoVariableRestricted( new float );

  *rhoVariable = float(*rhoAll);
  *rhoVariableRestricted = float(*rhoRestricted);

  evt.put( rhoVariable, "rhoAllEta" );
  evt.put( rhoVariableRestricted, "rhoRestrictedEta" );

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( EventRhoDumper );

