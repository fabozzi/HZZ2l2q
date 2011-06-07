#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


class MetVariablesProducer : public edm::EDProducer {
public:
  MetVariablesProducer( const edm::ParameterSet & );
   
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag metTag_;

};

MetVariablesProducer::MetVariablesProducer( const ParameterSet & cfg ) : 
  metTag_(cfg.getParameter<InputTag>("metTag")) {
  produces<float>( "met" ).setBranchAlias( "met" );
  produces<float>( "metSumEt" ).setBranchAlias( "metSumEt" );
  produces<float>( "metSig" ).setBranchAlias( "metSig" );
  produces<float>( "metSignificance" ).setBranchAlias( "metSignificance" );
  produces<float>( "metPhi" ).setBranchAlias( "metPhi" );
}



void MetVariablesProducer::produce( Event & evt, const EventSetup & ) {
  
  Handle<pat::METCollection> metHandle;
  evt.getByLabel(metTag_, metHandle);
  pat::METCollection met_h = *metHandle;

  auto_ptr<float> met_( new float );
  auto_ptr<float> metSumEt_( new float );
  auto_ptr<float> metSig_( new float );
  auto_ptr<float> metSignificance_( new float );
  auto_ptr<float> metPhi_( new float );

  *met_ = -100.;
  *metSumEt_ = -100.;
  *metSig_ = -100.;
  *metSignificance_ = -100.;
  *metPhi_ = -100.;

  *met_ = met_h.front().et();
  *metSumEt_ = met_h.front().sumEt();
  // rough met significance: met/sqrt(sumEt)
  *metSig_ = met_h.front().mEtSig();
  // met significance
  TMatrixD metmat = met_h.front().getSignificanceMatrix();
  //  metmat.Print();
  if( (metmat < 1.0e10) && (metmat > -1.0e10) ) {
    //    float determ_mat = metmat.Determinant();
    //    cout << "determinant = " << determ_mat << endl;
    *metSignificance_ = met_h.front().significance();
  } 
 
  //  cout << "metSignificance = " << *metSignificance_ << endl;

  *metPhi_ = met_h.front().phi();

  evt.put( met_, "met" );
  evt.put( metSumEt_, "metSumEt" );
  evt.put( metSig_, "metSig" );
  evt.put( metSignificance_, "metSignificance" );
  evt.put( metPhi_, "metPhi" );

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( MetVariablesProducer );

