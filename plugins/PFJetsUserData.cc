#include <memory>
#include <Riostream.h>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
//#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
//#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/CaloTowers/interface/CaloTower.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

typedef std::vector< pat::Jet > PFJetCollectionAB;
typedef std::vector< reco::PFCandidate > PFCandCollectionAB;

class PFJetUserData : public edm::EDProducer{

public:
  PFJetUserData(const edm::ParameterSet&);
  ~PFJetUserData(){
    //
  }

private:
 virtual void beginJob();
  virtual void endJob();
  void produce( edm::Event &, const edm::EventSetup &);

  //data members
  edm::InputTag   jetLabel_;
  bool verbose_;
};

PFJetUserData::PFJetUserData(const edm::ParameterSet &pSet){
  jetLabel_ =pSet.getUntrackedParameter<edm::InputTag>("JetInputCollection");
  verbose_=pSet.getUntrackedParameter<bool>("Verbosity");
  // issue the produce<>
 produces< std::vector< pat::Jet > >();
 
}

void PFJetUserData::beginJob(){
  //
}

void PFJetUserData::endJob(){
}

void PFJetUserData::produce(edm::Event &iEvt,  const edm::EventSetup &iSetup){
  if(verbose_)std::cout<<"Processing run "<<iEvt.id().run()<<", event "<<iEvt.id().event()<<std::endl;

 //pick from the event the input jet collection 
  Handle< PFJetCollectionAB > jetColl;
  iEvt.getByLabel(jetLabel_, jetColl);

  /*
    Handle< std::vector<pat::Muon> > muColl;
    iEvt.getByLabel(muSrc_, muColl);
    std::auto_ptr< std::vector<pat::Muon> > muCollNew(new std::vector<pat::Muon>(*muColl))
    for(unsigned int i=0; i<muColl->size();++i ){
    //  const pat::Muon& myMu=(*imu);
    pat::Muon& myMu=(*muColl)[i];
    myMu.addUserFloat("dummyFloat", 9999.0);
    }
  */


  int nTotInJets=0;
  nTotInJets += jetColl->size();

  if(nTotInJets>0){
    if(!jetColl->at(0).isPFJet() ){
      throw cms::Exception("Bad Input") <<"ERROR in PFJetUserData::produce ! Jet collection is NOT made by PF Jets. "<<std::endl;
    }
  }

  if(verbose_)std::cout<<"In the EVent are present "<<nTotInJets<<" PF jets."<<std::endl;

   //create output collection with PF Candidates in the jets
  std::auto_ptr< std::vector< pat::Jet > > outputPFJets(new std::vector< pat::Jet >(*jetColl));
 
  for(PFJetCollectionAB::iterator ijet=outputPFJets->begin(); ijet!=outputPFJets->end();++ijet ){
    
    if(verbose_){
      std::cout<<"From PFJetUserData::produce: PF Jet substructure (q/g separation) -> "<<
	"pt= "<<ijet->pt()<<"  eta= "<<ijet->eta()<<"   phi= "<<ijet->phi()<<flush;
    }
    
	   
    
    std::vector<reco::PFCandidatePtr> pfCandidates = ijet->getPFConstituents();
    
    float sumPt_cands=0.;
    float sumPt2_cands=0.;
    float rms_cands=0.;
    
    //loop on jet constituents       
    for (vector<reco::PFCandidatePtr>::const_iterator jt = pfCandidates.begin(); jt != pfCandidates.end(); ++jt) {
      
      //PFCandidate::ParticleType id = (*jt)->particleId();
      // Convert particle momentum to normal TLorentzVector, wrong type :(
      math::XYZTLorentzVectorD const& p4t = (*jt)->p4();
      TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
      TLorentzVector jetp4;
      //	 jetp4.SetPtEtaPhiE(ijet->pt(), ijet->eta(), ijet->phi(), ijet->energy());
      jetp4.SetPtEtaPhiE(ijet->pt(), ijet->eta(), ijet->phi(), ijet->energy());
      
      sumPt_cands += p4.Pt();
      sumPt2_cands += (p4.Pt()*p4.Pt());
      //float deltaR = ijet->p4().DeltaR(p4);
      float deltaR = jetp4.DeltaR(p4);
      rms_cands += (p4.Pt()*p4.Pt()*deltaR*deltaR);
      
    }//end loop on PF cands
    
    float nChrgdMult=float(ijet->chargedHadronMultiplicity());
    float nNeutrMult=float(ijet->neutralHadronMultiplicity());
    float ptDJet = sqrt( sumPt2_cands )/sumPt_cands;
    float rmsCandJet = rms_cands/(sumPt_cands*sumPt_cands);
    
    if(verbose_){
      std::cout<<"  NchgdHadrMult="<<nChrgdMult<<"  NneutrHadrMult= "<<nNeutrMult<<"   ptD= "<<ptDJet<<"   rmsJET= "<<rmsCandJet<<std::endl;
    }
    
    ijet->addUserFloat("nChrgdHadrMult", float(ijet->chargedHadronMultiplicity()));
    ijet->addUserFloat("nNeutrHadrMult", float(ijet->neutralHadronMultiplicity()));
    ijet->addUserFloat("nPhotMult", float(ijet->photonMultiplicity()));
    ijet->addUserFloat("nElecMult", float(ijet->electronMultiplicity()));
    ijet->addUserFloat("nMuonMult", float(ijet->muonMultiplicity()));
    ijet->addUserFloat("nHFHadrMult", float(ijet->HFHadronMultiplicity ()));
    ijet->addUserFloat("nHFEMMult", float(ijet->HFEMMultiplicity()));
    ijet->addUserFloat("pfChrgdHadrEnergy", float(ijet->chargedHadronEnergy()));
    ijet->addUserFloat("pfNeutrHadrEnergy", float(ijet->neutralHadronEnergy()));
    ijet->addUserFloat("pfPhotEnergy", float(ijet->photonEnergy()));
    ijet->addUserFloat("pfElecEnergy", float(ijet->electronEnergy()));
    ijet->addUserFloat("pfHFHadrEnergy", float(ijet->HFHadronEnergy()));
    ijet->addUserFloat("pfHFEMEnergy", float(ijet->HFEMEnergy()));
    ijet->addUserFloat("ptDJet", ptDJet);
    ijet->addUserFloat("RMSJet", rmsCandJet);

    
  }//end loop on jets
  
  //  for(std::vector<pat::Muon>::const_iterator imu=muColl->begin(); imu!=muColl->end();++imu ){
  if(verbose_) std::cout<<"Ended run "<<iEvt.id().run()<<", event "<<iEvt.id().event()<<", processed "<< outputPFJets->size()<<" jets."<<std::endl;
  iEvt.put(outputPFJets);
  if(verbose_)std::cout<<"UserJetCollection added to edmEvent."<<std::endl;

}//end produce



// ========= MODULE DEF ==============
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFJetUserData);
