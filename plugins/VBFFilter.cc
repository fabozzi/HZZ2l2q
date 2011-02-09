#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

class VBFFilter : public edm::EDFilter {
    public:
        explicit VBFFilter(const edm::ParameterSet&);
        ~VBFFilter();

    private:
        virtual void beginJob();
        virtual bool filter(edm::Event&, const edm::EventSetup&);
        virtual void endJob();

        edm::InputTag src;

        int higgs;
        int gluon;
};


VBFFilter::VBFFilter(const edm::ParameterSet& iConfig):
  src(iConfig.getParameter<edm::InputTag>("src")),
  higgs(25),
  gluon(21)
{
}


VBFFilter::~VBFFilter() {
}

bool VBFFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<GenParticleCollection> gensH;
    iEvent.getByLabel(src, gensH);
    GenParticleCollection gens = *gensH;

    for (size_t i = 0; i < gens.size(); i++) {
        if (gens[i].pdgId() == higgs && gens[i].status() == 3 && gens[i].mother()->pdgId() != gluon) return true;
    }

    return false;
}

void VBFFilter::beginJob() {
}

void VBFFilter::endJob() {
}

DEFINE_FWK_MODULE(VBFFilter);

