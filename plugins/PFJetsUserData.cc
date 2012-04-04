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

#include <memory>
#include <Riostream.h>
#include <string>
#include <vector>

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

  // Utilities for beta/beta*
  /// Computes the beta and betastar variables for the given PAT jet.
  void computeBeta (const pat::Jet &jjet,float *beta,float *betastar);
  /// Reads the vertices from the event record and process them to be used
  /// for the calculation of beta and betastar. Note that some selection is
  /// performed on what is a vertex.
  /// It also selects which is the main vertex of the analysis... taken to
  /// be the first that is valid and not fake (if it passes the cuts).
  void readVertices (const edm::Event &iEvent);

  //data members
  edm::InputTag   jetLabel_;

  //Z variable of the reconstructed vertices.
  std::vector<float> _verticesZ;  
  //Index of the main vertex (-1 if no main vertex).  
  int _mainVertex;  
  // Variables for statistics
  int _nValidVertices;   // Number of valid vertices.
  int _nSelectedVertices;  // Number of selected vertices.
  int _nEventsWithValidVtx;  // Number of events with a "valid" vertex.
  int _nEventsWithMainVtx;  // Number of events with a main vertex.

  bool verbose_;
};

PFJetUserData::PFJetUserData(const edm::ParameterSet &pSet){
  jetLabel_ =pSet.getUntrackedParameter<edm::InputTag>("JetInputCollection");
  verbose_=pSet.getUntrackedParameter<bool>("Verbosity");

  // for beta/beta* 
  _verticesZ.clear();
  _mainVertex=-1;
  _nValidVertices=0;
  _nSelectedVertices=0;
  _nEventsWithValidVtx=0;
  _nEventsWithMainVtx=0;

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


  // Initial setup for beta/beta* variables
  _verticesZ.clear();
  _mainVertex=-1;
  // read the vertices! 
  readVertices(iEvt);

  
  for(PFJetCollectionAB::iterator ijet=outputPFJets->begin(); ijet!=outputPFJets->end();++ijet ){


    // computing beta/beta* variables
    //    pat::Jet *jjet = &(*ijet);
    const pat::Jet & patjjet = *ijet;
    // We ask for the variables to store
    float beta=-1;
    float betastar=-1;
    //    computeBeta(*jjet,&beta,&betastar);
    computeBeta(patjjet,&beta,&betastar);
      


    
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
    
    // beta/beta* variables
    ijet->addUserFloat("puBeta",beta);
    ijet->addUserFloat("puBetaStar",betastar);
    //
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


//-----------------------------------------------------------------------
void PFJetUserData::computeBeta (const pat::Jet &jjet,float *beta,float *betastar)
// Computes the beta and betastar variables for the given PAT jet.
{
  // We set the values to 0 only if there will be information
  // associated to them... i.e. if there are vertices.
  if (_verticesZ.size()>0) {
    *betastar=0;
    if (_mainVertex!=-1) *beta=0;
  }

  float totalpt=0;  // Scalar sum of the pt of the charged PF constituent of the jet

  // We loop over the charged particles in the jet 

  for (std::vector<reco::PFCandidatePtr>::const_iterator xpart = jjet.getPFConstituents().begin();
         xpart!=jjet.getPFConstituents().end();++xpart) {

    reco::PFCandidate::ParticleType typ = (*xpart)->particleId();
    if (typ!=reco::PFCandidate::h
	&& typ!=reco::PFCandidate::e
	&& typ!=reco::PFCandidate::mu) continue;  // Neutral, ignored

    // We check the distance wrt to the vertices:

    float zpfo = (*xpart)->vz();

    totalpt += (*xpart)->pt();

    int minvtx=-1;
    float mindist=0.2;

    {int ivtx=0;
    for (std::vector<float>::const_iterator zvtx = _verticesZ.begin();
	 zvtx!=_verticesZ.end();++zvtx,++ivtx) {

      float d = fabs(zpfo-*zvtx);
      
      if (d<mindist) {
	minvtx=ivtx;
	mindist=d;
      }
    }}

    // Depending if the clostest vertex (if any) is the main or other
    // we fill beta or betastar
    if (minvtx==_mainVertex) {
      *beta += (*xpart)->pt();
    }
    else if (minvtx>=0) {
      *betastar += (*xpart)->pt();
    }
  }

  // We normalize the variables properly:

  if (totalpt>0) {
    if (*beta>0) *beta /= totalpt;
    if (*betastar>0) *betastar /= totalpt;
  }  

  //TEST  std::cout<<"JET: "<<jjet.pt()<<" "<<jjet.rapidity()<<" "<<jjet.phi()<<" "<<*beta<<" "<<*betastar<<std::endl;
}

//-----------------------------------------------------------------------
void PFJetUserData::readVertices (const edm::Event &iEvent)
// Reads the vertices from the event record and process them to be used
// for the calculation of beta and betastar. Note that some selection is
// performed on what is a vertex.
// It also selects which is the main vertex of the analysis... taken to
// be the first that is valid and not fake (if it passes the cuts).
{
  // Reading the vertices
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  
  // Processing the information

  int formain=0;

  for (reco::VertexCollection::const_iterator xvtx = recVtxs->begin();
       xvtx!=recVtxs->end();++xvtx) {

    if (!xvtx->isValid() || xvtx->isFake()) continue;  // Only valid and not fake vertices

    if (formain==0) ++_nEventsWithValidVtx;
    
    ++formain;
    /// We are selecting the first valid and not fake vertex as the
    /// main one... there is no main vertex if it fails the selection.

    // Cuts to select the vertices:

    if (xvtx->ndof()<4) continue;  // Number of dof

    if (fabs(xvtx->z())>24) continue;  // Z of the vertex

    {float rho = sqrt(xvtx->x()*xvtx->x()+xvtx->y()*xvtx->y());
      if (rho>2) continue;  // Rho of the vertex
    }
    
    // Here we have a good vertex:

    if (formain==1) {  // It is the first that is valid... main vertex.
      _mainVertex=_verticesZ.size();
      ++_nEventsWithMainVtx;
    }
    _verticesZ.push_back(xvtx->z());
  }

  // Counting the valid and not fake vertices.
  _nValidVertices+=formain;

  // Counting the selected vertices.
  _nSelectedVertices += _verticesZ.size();

}


// ========= MODULE DEF ==============
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFJetUserData);