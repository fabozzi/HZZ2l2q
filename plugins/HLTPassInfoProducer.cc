#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"

#include <string>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class HLTPassInfoProducer : public edm::EDProducer {
public:
  HLTPassInfoProducer( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag trigEvtTag_;
  std::vector<std::string> trigNamesMu_, trigNamesEl_,;
  std::vector<std::string> trigNamesDoubleMu_, trigNamesDoubleEl_,;
  bool verifyHLTPass(std::vector<std::string>, pat::TriggerPathRefVector);

};

HLTPassInfoProducer::HLTPassInfoProducer( const ParameterSet & cfg ) : 
  trigEvtTag_(cfg.getParameter<InputTag>("triggerEvent")),
  trigNamesMu_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleMu")),
  trigNamesEl_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleEl")), 
  trigNamesDoubleMu_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleMu")),
  trigNamesDoubleEl_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleEl"))
{
  produces<bool>( "passSingleMuTrig" ).setBranchAlias( "passSingleMuTrig" );
  produces<bool>( "passDoubleMuTrig" ).setBranchAlias( "passDoubleMuTrig" );
  produces<bool>( "passSingleElTrig" ).setBranchAlias( "passSingleElTrig" );
  produces<bool>( "passDoubleElTrig" ).setBranchAlias( "passDoubleElTrig" );
}



void HLTPassInfoProducer::produce( Event & evt, const EventSetup & ) {
  
  Handle<pat::TriggerEvent> trigEvt;
  evt.getByLabel(trigEvtTag_, trigEvt);
  pat::TriggerPathRefVector passedPaths = trigEvt->acceptedPaths();

  auto_ptr<bool> trigOKMu_( new bool );
  auto_ptr<bool> trigOKDoubleMu_( new bool );
  auto_ptr<bool> trigOKEl_( new bool );
  auto_ptr<bool> trigOKDoubleEl_( new bool );

  
  *trigOKMu_ = false;
  *trigOKEl_ = false;
  *trigOKDoubleMu_ = false;
  *trigOKDoubleEl_ = false;
  
  *trigOKMu_ = verifyHLTPass(trigNamesMu_, passedPaths);
  *trigOKDoubleMu_ = verifyHLTPass(trigNamesDoubleMu_, passedPaths);
  *trigOKEl_ = verifyHLTPass(trigNamesEl_, passedPaths);
  *trigOKDoubleEl_ = verifyHLTPass(trigNamesDoubleEl_, passedPaths);
  
  cout << "Single Mu FLAG VALUE = " << *trigOKMu_ << endl;
  cout << "Double Mu FLAG VALUE = " << *trigOKDoubleMu_ << endl;
  cout << "Single El FLAG VALUE = " << *trigOKEl_ << endl;
  cout << "Double El FLAG VALUE = " << *trigOKDoubleEl_ << endl;

  evt.put( trigOKMu_, "passSingleMuTrig" );
  evt.put( trigOKDoubleMu_, "passDoubleMuTrig" );
  evt.put( trigOKEl_, "passSingleElTrig" );
  evt.put( trigOKDoubleEl_, "passDoubleElTrig" );


}


bool HLTPassInfoProducer::verifyHLTPass(std::vector<std::string> trigNames, pat::TriggerPathRefVector passedPaths) {

  // simply pass if no trigger path is specified in input
  if(trigNames.size()==0) {
    cout << "ALL PATHs" << endl;
    return true;
  }
  for(vector<string>::iterator nameIt = trigNames.begin(); nameIt != trigNames.end();
      ++nameIt) {
    //    if(passedTrig == true)
    //      break;
    string requestedPath = *nameIt;
    cout << "Examining " << requestedPath << endl;

    for(pat::TriggerPathRefVector::const_iterator pathIt = passedPaths.begin(); pathIt!=passedPaths.end(); 
	++pathIt) {
      string passedPathName = (*pathIt)->name();
      if(requestedPath == passedPathName){
	cout << "SETTING flag to TRUE" << endl;
	return true;
      }
      
    }

  }

  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( HLTPassInfoProducer );

