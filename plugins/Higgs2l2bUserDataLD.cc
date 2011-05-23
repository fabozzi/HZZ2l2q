#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
// Likelihood discriminant
#include "HiggsAnalysis/Higgs2l2b/interface/HelicityLikelihoodDiscriminant.h"

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bUserDataLD : public edm::EDProducer {
public:
  Higgs2l2bUserDataLD( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  
  InputTag higgsTag;
  HelicityLikelihoodDiscriminant LD_;

};



Higgs2l2bUserDataLD::Higgs2l2bUserDataLD( const ParameterSet & cfg ):
  higgsTag( cfg.getParameter<InputTag>( "higgs" ) )
{
  produces<vector<pat::CompositeCandidate> >("h").setBranchAlias( "h" );
}

void Higgs2l2bUserDataLD::produce( Event & evt, const EventSetup & ) {

  Handle<std::vector<pat::CompositeCandidate> > higgsH;
  evt.getByLabel(higgsTag, higgsH);

  auto_ptr<vector<pat::CompositeCandidate> > higgsColl( new vector<pat::CompositeCandidate> () );

  float helyLD;

  for (unsigned int i = 0; i< higgsH->size();++i){
    const pat::CompositeCandidate & H = (*higgsH)[i];
    edm::Ref<std::vector<pat::CompositeCandidate> > hRef(higgsH, i);
    pat::CompositeCandidate h(H);


    HelicityLikelihoodDiscriminant::HelicityAngles myha;
    myha.helCosTheta1    = H.userFloat("costhetaNT1");
    //    cout << "helCosTheta1 = " << myha.helCosTheta1 << endl;
    myha.helCosTheta2    = H.userFloat("costhetaNT2");
    //    cout << "helCosTheta2 = " << myha.helCosTheta2 << endl;
    myha.helCosThetaStar = H.userFloat("costhetastarNT");
    //    cout << "helCosThetaStar = " << myha.helCosThetaStar << endl;
    myha.helPhi          = H.userFloat("phiNT");
    //    cout << "helPhi = " << myha.helPhi << endl;
    myha.helPhi1         = H.userFloat("phiNT1");
    //    cout << "helPhi1 = " << myha.helPhi1 << endl;
    myha.mzz             = H.mass();
    //    cout << "mzz = " << myha.mzz << endl;

    LD_.setMeasurables(myha);
    float signProb = LD_.getSignalProbability();
    float bkgdProb = LD_.getBkgdProbability();
    helyLD = signProb / (signProb + bkgdProb);

    //    cout << "signProb = " << signProb << endl;
    //    cout << "bkgdProb = " << bkgdProb << endl;
    //    cout << "helyLD = " << helyLD << endl;


    h.addUserFloat("helyLD", helyLD);

    higgsColl->push_back(h);
  }

  
  evt.put( higgsColl, "h");

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bUserDataLD );

