#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

using namespace edm;
using namespace std;
//using namespace reco;

class FilterByRun : public edm::EDFilter {
public:
  explicit FilterByRun( const edm::ParameterSet & );
  ~FilterByRun();

private:

  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  int runLimit_;
  string selMode_;
  // can be "g", "ge", "l" "le"
};

FilterByRun::FilterByRun( const ParameterSet & cfg ) : 
  runLimit_(cfg.getParameter<int>("run")),
  selMode_(cfg.getParameter<string>("selMode"))
{
}

FilterByRun::~FilterByRun() {
}

bool FilterByRun::filter( Event & evt, const EventSetup & ) {
  
  bool eventOK(false);
  int runNumber = evt.run();

  if(selMode_ == "g")
    eventOK = (runNumber > runLimit_) ;
  else if(selMode_ == "ge")
    eventOK = (runNumber >= runLimit_) ;
  else if(selMode_ == "l")
    eventOK = (runNumber < runLimit_) ;
  else if(selMode_ == "le")
    eventOK = (runNumber <= runLimit_) ;
  else
    cout << "Unknown selection mode! Event not selected" << endl;

  return eventOK;


}


void FilterByRun::beginJob() {
}

void FilterByRun::endJob() {
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FilterByRun);




