#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include <vector>

using namespace edm;
using namespace std;


class GenPUNtupleDump : public edm::EDProducer {
public:
  GenPUNtupleDump( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag PileupSrc_;
};

GenPUNtupleDump::GenPUNtupleDump( const ParameterSet & cfg ) : 
  PileupSrc_("addPileupInfo") {
  produces<int>( "nGenInt" ).setBranchAlias( "nGenInt" );
}



void GenPUNtupleDump::produce( Event & evt, const EventSetup & ) {
  
  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  evt.getByLabel(PileupSrc_, PupInfo);

  auto_ptr<int> nGenInt( new int );
  *nGenInt = -1;
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    
    int BX = PVI->getBunchCrossing();
    
    if(BX == 0) { 
      *nGenInt = PVI->getPU_NumInteractions();
      continue;
    }
    
  }
  
  evt.put( nGenInt, "nGenInt" );

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( GenPUNtupleDump );

