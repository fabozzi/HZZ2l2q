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
  produces<int>( "nGenIntBXm1" ).setBranchAlias( "nGenIntBXm1" );
  produces<int>( "nGenIntBX0" ).setBranchAlias( "nGenIntBX0" );
  produces<int>( "nGenIntBXp1" ).setBranchAlias( "nGenIntBXp1" );
  // true number of interactions (available in Fall11 MC)
  produces<float>( "trueNInt" ).setBranchAlias( "trueNInt" );
}

void GenPUNtupleDump::produce( Event & evt, const EventSetup & ) {
  
  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  evt.getByLabel(PileupSrc_, PupInfo);

  auto_ptr<float> nGenInt( new float );
  *nGenInt = -1.;
  
  auto_ptr<int> nGenIntBXm1( new int );
  *nGenIntBXm1 = -1;
  
  auto_ptr<int> nGenIntBX0( new int );
  *nGenIntBX0 = -1;
  
  auto_ptr<int> nGenIntBXp1( new int );
  *nGenIntBXp1 = -1;
  
  auto_ptr<float> trueNInt( new float );
  *trueNInt = -1.;

  if(!isData_) {

    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();

      if(BX == -1) { 
	*nGenIntBXm1 = PVI->getPU_NumInteractions();
      }
      if(BX == 0) { 
	*nGenIntBX0 = PVI->getPU_NumInteractions();
	*trueNInt = PVI->getTrueNumInteractions();
      }
      if(BX == 1) { 
	*nGenIntBXp1 = PVI->getPU_NumInteractions();
      }

    }

    *nGenInt = float(*nGenIntBXm1 + *nGenIntBX0 + *nGenIntBXp1)/3.;

  }
  
  evt.put( nGenInt, "nGenInt" );
  evt.put( nGenIntBXm1, "nGenIntBXm1" );
  evt.put( nGenIntBX0, "nGenIntBX0" );
  evt.put( nGenIntBXp1, "nGenIntBXp1" );
  evt.put( trueNInt, "trueNInt" );
  
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( GenPUNtupleDump );

