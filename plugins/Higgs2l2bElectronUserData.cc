#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include <vector>
#include <TMath.h>

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bElectronUserData : public edm::EDProducer {
public:
  Higgs2l2bElectronUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag src_, rho_;
  const float R03;
};

Higgs2l2bElectronUserData::Higgs2l2bElectronUserData( const ParameterSet & cfg ):
  src_( cfg.getParameter<InputTag>("src") ),
  rho_( cfg.getParameter<edm::InputTag>("rho")),
  R03(0.3)
{
  produces<std::vector<pat::Electron> >();
}

void Higgs2l2bElectronUserData::produce( Event & evt, const EventSetup & ) {

  Handle<vector<pat::Electron>  > electrons;
  evt.getByLabel(src_, electrons);

  Handle<double> rhoHandle;
  evt.getByLabel(rho_,rhoHandle);

  double rho = *rhoHandle; 
  float PUEnergyInCone = (TMath::Pi()) * R03 * R03 * rho;  

  auto_ptr<vector<pat::Electron> > electronColl( new vector<pat::Electron> (*electrons) );
  for (unsigned int i = 0; i< electronColl->size();++i){
    pat::Electron & el = (*electronColl)[i];
    float absCombIsoPUCorrected = -1.0;
    if( el.isEB() )
      // pedestal subtraction for barrel 
      absCombIsoPUCorrected = el.dr03TkSumPt() + max(0., el.dr03EcalRecHitSumEt() - 1.) + el.dr03HcalTowerSumEt() - PUEnergyInCone;
    else
      absCombIsoPUCorrected = el.dr03TkSumPt() + el.dr03EcalRecHitSumEt() + el.dr03HcalTowerSumEt() - PUEnergyInCone;

    el.addUserFloat("absCombIsoPUCorrected", absCombIsoPUCorrected);
  }

  evt.put(electronColl);

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bElectronUserData );


