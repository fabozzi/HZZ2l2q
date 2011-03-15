#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "HiggsAnalysis/Higgs2l2b/interface/AlpgenBBSamplesFlag.hh"

//-------------------------------------------------------
//-------------------------------------------------------


AlpgenBBSamplesFlag::AlpgenBBSamplesFlag(const edm::ParameterSet& iConfig)
{
  produces<unsigned int>( "flavorflag" ).setBranchAlias( "flavorflag" );
  
  _singleCPtCut=20;  // Let's try with 30 GeV. 
  _deltaRPairs=0.5;  // Default is 0.4. Recommened alternatices: 0.3 and 0.5
  _pairsPtCut=-1;    // By default we use all the pairs (15 is a clean alternative).
}

AlpgenBBSamplesFlag::~AlpgenBBSamplesFlag()
{
}

//-------------------------------------------------------

void AlpgenBBSamplesFlag::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  bool eventIsRejected =false;  // Event is accepted by default                                           
  
  int ncharm = 0;
  int nbottom = 0;

  Handle<GenParticleCollection> genpartH;
  iEvent.getByLabel("genParticles",genpartH);
  GenParticleCollection genpart = *genpartH;
  
  ///< To identify the type of events (accepted or not) 
  ///<and to store it during the skimming (1,....,15)
  std::auto_ptr<unsigned int> flavorflag ( new unsigned int );
  *flavorflag=0;

  for (size_t i=0; i<genpart.size() ;++i) {
    if (genpart[i].pdgId()==4 || genpart[i].pdgId()==-4)  ++ncharm;
    if (genpart[i].pdgId()==5 || genpart[i].pdgId()==-5)  ++nbottom;
  }
  
 
  // Checks:                                                                          
  if (nbottom==1) {
    cout<<"ERROR: Single-bottom production! "<<endl;
  }
  
  // Checking if we have a single-C                                                     
  
  if (ncharm==1) {
    eventIsRejected = processSingleC(genpart);  // It should be true if we reject the event              
    
    if (!eventIsRejected){
      *flavorflag=3;
    }
  }
  
  // Or more than one C, leading to CC pairs                                          
  
  else if(ncharm>1) {
    eventIsRejected = processCCPair(genpart);  // It should be true if we reject the event              

    if (eventIsRejected){
    }else{
      *flavorflag=2;
    }
  }
  
  // Checking if there are BB pairs in the event (if event is to be accepted)         
  
  if (!eventIsRejected && nbottom>1) {
    eventIsRejected = processBBPair(genpart);  // It should be true if we reject the event  

    if (eventIsRejected){ 
      *flavorflag=9;
    }else{
      *flavorflag=1;
    }
  }

  iEvent.put( flavorflag, "flavorflag" );

  //return !eventIsRejected; 
 
}

//-------------------------------------------------------
//-------------------------------------------------------
bool AlpgenBBSamplesFlag::processCCPair (GenParticleCollection genpart)
  // Make the processing for events with CC pairs, returning                          
  // true if the event has to be rejected.                                             
  // It also fills the corresponding plots for that kind of event.                    
{
  using namespace edm;
  using namespace reco;
  using namespace std;


  GenParticleCollection charm;
  GenParticleCollection  anticharm;

  float maxpt=-999;  // Charm with the highest pt:

  for (size_t i=0; i<genpart.size() ;++i) {
    //const PurdHepgPart *x = _genReco->getParticleAtPartonLevel(i);

    if (genpart[i].pdgId()==4) {
      charm.push_back(genpart[i]);
      if (genpart[i].pt()>maxpt) maxpt=genpart[i].pt();
    }
    else if (genpart[i].pdgId()==-4) {
      anticharm.push_back(genpart[i]);
      if (genpart[i].pt()>maxpt) maxpt=genpart[i].pt();
    }
  }

  if (charm.size()==0 || anticharm.size()==0) {
    cout<<"ERROR: CC Pair event and charm or anicharm not found!!! "
        <<charm.size()<<" "<<anticharm.size()<<endl;
    return true;  // considered problematic
  }

  // We look for the maximum of the minimum for all possible valid pairs:             

  float maxdr=-0.2;
  int npairs=0;

  for (size_t ic = 0; ic < charm.size(); ++ic) {

    float mindr=9999;

    for (size_t iac = 0 ; iac < anticharm.size(); ++iac) {
      if (charm[ic].pt()>=getPtCutforPairs() ||
          anticharm[iac].pt()>=getPtCutforPairs() ||
          charm[ic].pt()+anticharm[iac].pt()>=getPtCutforPairs()) {

	float dr = deltaR(charm[ic].eta(), charm[ic].phi(), anticharm[iac].eta(), anticharm[iac].phi());

	++npairs;

        if (dr<mindr) mindr=dr;
      }
    }

    if (mindr>9998) cout<<"ERROR: Problems with the charm-anticharm distance:  "
			<<mindr<<endl;
    else if (mindr>maxdr) maxdr=mindr;
  }

  // If there are no valid pairs, we assume it is like single-charm
  // and use the one with the highest pt.

  if (maxdr<-0.1) {
    return cuttingSingleC(maxpt);
  }

 // Applying the cut:

  return cuttingCCPair(maxdr);
}
//----------------------------------------------------------------------
bool AlpgenBBSamplesFlag::processBBPair (GenParticleCollection genpart)
  // Make the processing for events with BB pairs, returning
  // true if the event has to be rejected.
  // It also fills the corresponding plots for that kind of event.
{
  using namespace edm;
  using namespace reco;
  using namespace std;


  GenParticleCollection bottom;
  GenParticleCollection antibottom;

  for (size_t i=0; i<genpart.size() ;++i) {
    if (genpart[i].pdgId()==5) bottom.push_back(genpart[i]);
    else if (genpart[i].pdgId()==-5) antibottom.push_back(genpart[i]);
  }

  if (bottom.size()==0 || antibottom.size()==0) {
    cout<<"ERROR: BB Pair event but charm or anicharm not found !!!!  "
        <<bottom.size()<<" "<<antibottom.size()<<endl;
    return true;  // considered problematic                                           
  }
  // We look for the maximum of the minimum for all possible valid pairs:             

  double maxdr=0;
  int npairs=0;

  for (size_t ib = 0; ib < bottom.size();++ib) {

    float mindr=9999;

    for (size_t iab =0; iab < antibottom.size();++iab) {
      if (bottom[ib].pt()>=getPtCutforPairs() ||
          antibottom[iab].pt()>=getPtCutforPairs() ||
          bottom[ib].pt()+antibottom[iab].pt()>=getPtCutforPairs()) {

        float dr = deltaR(bottom[ib].eta(), bottom[ib].phi(), antibottom[iab].eta(), antibottom[iab].phi());

	++npairs;
	
        if (dr<mindr) mindr=dr;
      }
    }

    if (mindr>9998) cout<<"ERROR: Problems with the bottom-antibottom distance:  "
		       <<mindr<<endl;
    else if (mindr>maxdr) maxdr=mindr;
  }

  // Applying the cut:

  return cuttingBBPair(maxdr);
}

//----------------------------------------------------------------------
bool AlpgenBBSamplesFlag::processSingleC (GenParticleCollection genpart)
  // Make the processing for events with a single C, returning
  // true if the event is to be rejected.
  // It also fills the corresponding plots for that kind of event.
{
  
  double pt=-1;

  using namespace edm;
  using namespace reco;
  using namespace std;

  for (size_t i=0; i<genpart.size() ;++i) {
    if (genpart[i].pdgId()==4 || genpart[i].pdgId()==-4) {pt=genpart[i].pt(); break;}
  }

  if (pt<0) {
    cout<<"ERROR: Single-charm event but charm quark was not found !!!!"
	<<endl;
    return true;  // considered problematic
  }

  // Applying the cut:

  return cuttingSingleC(pt);
}
//-------------------------------------------------------
//-------------------------------------------------------

bool AlpgenBBSamplesFlag::cuttingSingleC (float)
  // Makes the selection for single-C events, using the pt of the charm quark.
  // Returns true if the event is rejected, that is always because events with
  // no valid charm pairs are not considered in W+CC but in W+C (or W+Np).
{
  return false; // Events with C are accepted
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
bool AlpgenBBSamplesFlag::cuttingCCPair (float dr) 
  // Makes the selection for CC-pair events, using the maximum of the minimum
  // Dr for the valid pairs.
  // Returns true if the event is rejected.
{
  return false; // Event with CC are accepted
}
//-----------------------------------------------------------------------
bool AlpgenBBSamplesFlag::cuttingBBPair (float dr) 
  // Makes the selection for BB-pair events, using the maximum of the minimum
  // Dr for the valid pairs.
  // Returns true if the event is rejected.
{
  if (dr>=getPairDeltaRCut()) return false; // Event is fine for W+BB samples
  return true; // Event must come from other samples.
}
//-------------------------------------------------------
//-------------------------------------------------------
void AlpgenBBSamplesFlag::beginJob() {
}

void AlpgenBBSamplesFlag::endJob() {
}

DEFINE_FWK_MODULE(AlpgenBBSamplesFlag);
