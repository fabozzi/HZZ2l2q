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
  bool isData_;
};

GenPUNtupleDump::GenPUNtupleDump( const ParameterSet & cfg ) : 
  PileupSrc_("addPileupInfo"),
  isData_(cfg.getParameter<bool>("isData"))
{
  produces<float>( "nGenInt" ).setBranchAlias( "nGenInt" );
}



void GenPUNtupleDump::produce( Event & evt, const EventSetup & ) {
  
  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  evt.getByLabel(PileupSrc_, PupInfo);

  auto_ptr<float> nGenInt( new float );
  *nGenInt = -1.;
  
  if(!isData_) {
    std::vector<PileupSummaryInfo>::const_iterator PVI;

    float sum_nvtx = 0;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      
      int npv = PVI->getPU_NumInteractions();
      sum_nvtx += float(npv);

      //      int BX = PVI->getBunchCrossing();
      
      //      if(BX == 0) { 
      //	*nGenInt = PVI->getPU_NumInteractions();
      //	continue;
      //      }
      
    }

    *nGenInt = sum_nvtx/3.;

  }
  
  evt.put( nGenInt, "nGenInt" );
  
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( GenPUNtupleDump );

