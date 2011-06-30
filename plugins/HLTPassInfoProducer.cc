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
  std::vector<int> runLimits_;
  std::vector<std::string> trigNamesMu_MC_, trigNamesEl_MC_;
  std::vector<std::string> trigNamesDoubleMu_MC_, trigNamesDoubleEl_MC_;
  std::vector<std::string> trigNamesMu_5e32_, trigNamesEl_5e32_;
  std::vector<std::string> trigNamesDoubleMu_5e32_, trigNamesDoubleEl_5e32_;
  std::vector<std::string> trigNamesMu_1e33_, trigNamesEl_1e33_;
  std::vector<std::string> trigNamesDoubleMu_1e33_, trigNamesDoubleEl_1e33_;
  std::vector<std::string> trigNamesMu_1p4e33_, trigNamesEl_1p4e33_;
  std::vector<std::string> trigNamesDoubleMu_1p4e33_, trigNamesDoubleEl_1p4e33_;
  bool verifyHLTPass(std::vector<std::string>, pat::TriggerPathRefVector);
};

HLTPassInfoProducer::HLTPassInfoProducer( const ParameterSet & cfg ) : 
  trigEvtTag_(cfg.getParameter<InputTag>("triggerEvent")),
  runLimits_(cfg.getParameter< std::vector<int> >("runLimits")),
  trigNamesMu_MC_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleMu_MC")),
  trigNamesEl_MC_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleEl_MC")), 
  trigNamesDoubleMu_MC_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleMu_MC")),
  trigNamesDoubleEl_MC_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleEl_MC")),
  trigNamesMu_5e32_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleMu_5e32")),
  trigNamesEl_5e32_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleEl_5e32")), 
  trigNamesDoubleMu_5e32_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleMu_5e32")),
  trigNamesDoubleEl_5e32_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleEl_5e32")),
  trigNamesMu_1e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleMu_1e33")),
  trigNamesEl_1e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleEl_1e33")), 
  trigNamesDoubleMu_1e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleMu_1e33")),
  trigNamesDoubleEl_1e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleEl_1e33")),
  trigNamesMu_1p4e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleMu_1p4e33")),
  trigNamesEl_1p4e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesSingleEl_1p4e33")), 
  trigNamesDoubleMu_1p4e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleMu_1p4e33")),
  trigNamesDoubleEl_1p4e33_(cfg.getParameter< std::vector<std::string> >("triggerNamesDoubleEl_1p4e33"))
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

  int runNumber = evt.run();

  cout << "RUN = " << runNumber << endl;

  enum{ start5E32, start1E33, start1P4E33 };
    
  std::vector<std::string> trigNamesMu_, trigNamesEl_;
  std::vector<std::string> trigNamesDoubleMu_, trigNamesDoubleEl_;

  auto_ptr<bool> trigOKMu_( new bool );
  auto_ptr<bool> trigOKDoubleMu_( new bool );
  auto_ptr<bool> trigOKEl_( new bool );
  auto_ptr<bool> trigOKDoubleEl_( new bool );

  
  *trigOKMu_ = false;
  *trigOKEl_ = false;
  *trigOKDoubleMu_ = false;
  *trigOKDoubleEl_ = false;


  if(runLimits_.size() == 0 ) {
    // is MC

    trigNamesMu_ = trigNamesMu_MC_;
    trigNamesEl_ = trigNamesEl_MC_;
    trigNamesDoubleMu_ = trigNamesDoubleMu_MC_;
    trigNamesDoubleEl_ = trigNamesDoubleEl_MC_;

  } else {
    // is Data

    if( (runNumber>=runLimits_[start5E32]) && ( runNumber<runLimits_[start1E33]) ) {
      // use 5e32 paths
      trigNamesMu_ = trigNamesMu_5e32_;
      trigNamesEl_ = trigNamesEl_5e32_;
      trigNamesDoubleMu_ = trigNamesDoubleMu_5e32_;
      trigNamesDoubleEl_ = trigNamesDoubleEl_5e32_;
    }
    if( (runNumber>=runLimits_[start1E33]) && ( runNumber<runLimits_[start1P4E33]) ) {
      // use 1e33 paths
      trigNamesMu_ = trigNamesMu_1e33_;
      trigNamesEl_ = trigNamesEl_1e33_;
      trigNamesDoubleMu_ = trigNamesDoubleMu_1e33_;
      trigNamesDoubleEl_ = trigNamesDoubleEl_1e33_;
    }
    if( (runNumber>=runLimits_[start1P4E33]) )  {
      // use 1.4e33 paths
      trigNamesMu_ = trigNamesMu_1p4e33_;
      trigNamesEl_ = trigNamesEl_1p4e33_;
      trigNamesDoubleMu_ = trigNamesDoubleMu_1p4e33_;
      trigNamesDoubleEl_ = trigNamesDoubleEl_1p4e33_;
    }

  }
    
  
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

  // if no trigger path is specified in input -> pass OK
  if(trigNames.size()==0) 
    return true;

  for(vector<string>::iterator nameIt = trigNames.begin(); nameIt != trigNames.end();
      ++nameIt) {
    string requestedPath = *nameIt;
    cout << "Examining " << requestedPath << endl;

    for(pat::TriggerPathRefVector::const_iterator pathIt = passedPaths.begin(); pathIt!=passedPaths.end(); 
	++pathIt) {
      string passedPathName = (*pathIt)->name();      
      int vpos = passedPathName.find_last_of("_");
      if(vpos < 0 )
	continue;
      passedPathName.erase(vpos);

      //      cout << "PASSED PATH  = " << passedPathName << endl;

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

