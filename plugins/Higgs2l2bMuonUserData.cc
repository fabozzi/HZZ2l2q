#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <vector>
#include <TMath.h>

using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bMuonUserData : public edm::EDProducer {
public:
  Higgs2l2bMuonUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag src_, rho_, primaryVertices_;
  const float R03;
};

Higgs2l2bMuonUserData::Higgs2l2bMuonUserData( const ParameterSet & cfg ):
  src_( cfg.getParameter<InputTag>("src") ),
  rho_( cfg.getParameter<edm::InputTag>("rho")),
  primaryVertices_(cfg.getParameter<InputTag>("primaryVertices")),
  R03(0.3)
{
  produces<std::vector<pat::Muon> >();
}

void Higgs2l2bMuonUserData::produce( Event & evt, const EventSetup & ) {

  Handle<vector<pat::Muon>  > muons;
  evt.getByLabel(src_,muons);

  Handle<double> rhoHandle;
  evt.getByLabel(rho_,rhoHandle);

  double rho = *rhoHandle; 
  float PUEnergyInCone = (TMath::Pi()) * R03 * R03 * rho;  

  Handle<reco::VertexCollection> primaryVertices;  // Collection of primary Vertices
  evt.getByLabel(primaryVertices_, primaryVertices);
  const reco::Vertex &pv = (*primaryVertices)[0];

  auto_ptr<vector<pat::Muon> > muonColl( new vector<pat::Muon> (*muons) );
  for (unsigned int i = 0; i< muonColl->size();++i){
    pat::Muon & m = (*muonColl)[i];
    float absCombIsoPUCorrected = m.trackIso() + m.caloIso() - PUEnergyInCone;
    m.addUserFloat("absCombIsoPUCorrected", absCombIsoPUCorrected);
    float dzVtx = m.innerTrack()->dz(pv.position());
    m.addUserFloat("dzVtx", dzVtx);

  }

  evt.put( muonColl);

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bMuonUserData );


