#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

using namespace edm;
using namespace std;
//using namespace reco;

class DYExclusiveFilter : public edm::EDFilter {
public:
  explicit DYExclusiveFilter( const edm::ParameterSet & );
  ~DYExclusiveFilter();

private:

  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  InputTag src_;
};

DYExclusiveFilter::DYExclusiveFilter( const ParameterSet & cfg ) :
  src_(cfg.getParameter<InputTag>("src"))
{
}

DYExclusiveFilter::~DYExclusiveFilter() {
}

bool DYExclusiveFilter::filter( Event & evt, const EventSetup & ) {
  
  bool eventOK(false);
  edm::Handle<int> nup;
  evt.getByLabel(src_, nup);

  if(*nup == 5)
    eventOK = true;
      
  return eventOK;

}


void DYExclusiveFilter::beginJob() {
}

void DYExclusiveFilter::endJob() {
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DYExclusiveFilter);




