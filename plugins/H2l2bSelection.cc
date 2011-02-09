#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "PhysicsTools/CandUtils/interface/CenterOfMassBooster.h"
#include "PhysicsTools/CandUtils/interface/Booster.h"
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TVector3.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <Math/VectorUtil.h>



class H2l2bSelection : public edm::EDAnalyzer {
public:
  explicit H2l2bSelection(const edm::ParameterSet&);
  ~H2l2bSelection();
  
  
private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob(); 
  float massSel, btagSel, zllptSel, metSel, jjdrSel, hmassSel;
  float skimEvents;
  bool massSelected, btagSelected, zllptSelected, metSelected, jjdrSelected, hmassSelected;
  float met;
  float zNominalMass;
  edm::InputTag higgsTag_;
  double zLepMassCut_;
  double zJetMassCut_;
  double combSecVtxCut_;
  double bTagCSVCut_;
  double bTagJProbCut_;
  double zllPtCut_;
  double metCut_;
  double jjdrCut_;
  double hMassMinCut_;
  double hMassMaxCut_;
  double lumiNormalization_;
  std::string channel,output_name_;
  ofstream outFile;

        
};

H2l2bSelection::H2l2bSelection(const edm::ParameterSet& iConfig):
  higgsTag_(iConfig.getParameter<edm::InputTag>("higgsTag")),
  zLepMassCut_(iConfig.getParameter<double>("zLepMassCut")),
  zJetMassCut_(iConfig.getParameter<double>("zJetMassCut")),
  bTagCSVCut_(iConfig.getParameter<double>("bTagCSVCut")),
  bTagJProbCut_(iConfig.getParameter<double>("bTagJProbCut")),
  zllPtCut_(iConfig.getParameter<double>("zllPtCut")),
  metCut_(iConfig.getParameter<double>("metCut")),
  jjdrCut_(iConfig.getParameter<double>("jjdrCut")),
  hMassMinCut_(iConfig.getParameter<double>("hMassMinCut")),
  hMassMaxCut_(iConfig.getParameter<double>("hMassMaxCut")),
  lumiNormalization_(iConfig.getParameter<double>("lumiNormalization")),
  output_name_(iConfig.getParameter<std::string>("output_name"))

{
  edm::Service<TFileService> fs;
  massSel=0;
  massSelected=false;
  btagSel=0;
  btagSelected=false;
  zllptSel=0;
  zllptSelected=false;
  metSel=0;
  metSelected=false;
  jjdrSel=0;
  jjdrSelected=false;
  hmassSel=0;
  hmassSelected=false;
  channel = "";
  //std::cout<<"normalization "<<lumiNormalization_<<std::endl;
  
}


H2l2bSelection::~H2l2bSelection() {

}


void H2l2bSelection::beginJob() {
  zNominalMass = 91;
  skimEvents = 0;
  
}


void H2l2bSelection::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    // Get hzzjjlls collections

    std::cout<<"new event "<<std::endl;
    Handle<CandidateView> higgsH;
    iEvent.getByLabel(higgsTag_, higgsH);

 //    Handle<pat::METCollection> metH;
//     iEvent.getByLabel(metTag, metH);
//     pat::METCollection met_h = *metH;
//     met = met_h.front().et();

    Handle<float> met_elH;
    //iEvent.getByLabel(metTag,metH);
    iEvent.getByLabel("hzzeejj","met",met_elH);

    Handle<float> met_muH;
    iEvent.getByLabel("hzzmmjj","met",met_muH);
 

    skimEvents++;

 
    massSelected = false;
    btagSelected = false;
    zllptSelected = false;
    metSelected = false;
    jjdrSelected = false;
    hmassSelected = false;
    float jjdr;

    std::cout<<" inizialization "<<massSelected<<std::endl;    
    
    for (size_t i = 0; i < higgsH->size(); ++i) {
    
      const Candidate &  h = (*higgsH)[i];
      const float &  metEl = (*met_elH);
      const float &  metMu = (*met_muH);
      const Candidate * zll = h.daughter(0);    
      const Candidate * zjj = h.daughter(1);
      const Candidate * zDauRefl0 = h.daughter(0)->daughter(0);
      const Candidate * zDauRefl1 = h.daughter(0)->daughter(1);
      const Candidate * zDauRefj0 = h.daughter(1)->daughter(0);
      const Candidate * zDauRefj1 = h.daughter(1)->daughter(1);  




      const pat::Electron * lept0el = dynamic_cast<const pat::Electron *>(zDauRefl0->masterClone().get());
      const pat::Electron * lept1el = dynamic_cast<const pat::Electron *>(zDauRefl1->masterClone().get());
      const pat::Muon * lept0mu;
      const pat::Muon * lept1mu;

      float dB0, dB1;
      bool VBTF80CombID;
      bool dB;
      if(lept0el != 0){
	if(channel == ""){
	  channel = "Electron";
	  outFile.open((output_name_+channel+"Selection.txt").c_str());	  
	}	
	VBTF80CombID = lept0el->electronID("eidVBTFCom80")==7 || lept1el->electronID("eidVBTFCom80")==7;
	dB = true; 
      }else{
	if(channel == ""){
	  channel = "Muon";
	  outFile.open((output_name_+channel+"Selection.txt").c_str());	  
	}	
	lept0mu = dynamic_cast<const pat::Muon *>(zDauRefl0->masterClone().get());
	lept1mu = dynamic_cast<const pat::Muon *>(zDauRefl1->masterClone().get());      
	dB0 = TMath::Abs(lept0mu->dB());
	dB1 = TMath::Abs(lept1mu->dB());
	VBTF80CombID = true;
	dB = dB0< 0.02 && dB1<0.02;
      }

	

      const pat::Jet & j0 = dynamic_cast<const pat::Jet &>(*zDauRefj0->masterClone());
      const pat::Jet & j1 = dynamic_cast<const pat::Jet &>(*zDauRefj1->masterClone());

      jjdr = deltaR(zDauRefj0->eta(), zDauRefj0->phi(), zDauRefj1->eta(), zDauRefj1->phi() );

      /* base selection and progressive cuts -> put selection on VBTFCombID for electrons */ 
      
      if( zDauRefl0->pt()>20 && zDauRefl1->pt()>20 && zDauRefj0->pt()>30 && zDauRefj1->pt()>30 	&& dB  && VBTF80CombID){             

	if(TMath::Abs(zll->mass() - zNominalMass)< zLepMassCut_ && TMath::Abs(zjj->mass() - zNominalMass)< zJetMassCut_){
	  massSelected=true ;
	  //std::cout<<"zll mass: "<< zll->mass()<<" zjj mass: "<< zjj->mass()<<std::endl;
	  
	  if(  ((j0.bDiscriminator("combinedSecondaryVertexMVABJetTags")>bTagCSVCut_ && j1.bDiscriminator("jetBProbabilityBJetTags")>bTagJProbCut_)||(j0.bDiscriminator("jetBProbabilityBJetTags")>bTagJProbCut_ && j1.bDiscriminator("combinedSecondaryVertexMVABJetTags")>bTagCSVCut_) ) ) {
	    btagSelected=true ;
	    
	    if( zll->pt()>zllPtCut_ ) {
	      zllptSelected=true ;
	      
	      if( (channel == "Electron" && metEl < metCut_) || (channel == "Muon" && metMu < metCut_) ) {
		metSelected=true ;
		
		if( jjdr < jjdrCut_ ){ 
		  jjdrSelected=true ;
		  
		  if( h.mass()<hMassMaxCut_ && h.mass()>hMassMinCut_ ){
		    hmassSelected=true ;
		  }	       
		}
	      }
	    }
	  }
	}
      }
      
    }




    if(massSelected == true) ++massSel;
    //std::cout<<" event has been selected "<<massSelected*lumiNormalization_<<std::endl;

    if(btagSelected == true) ++btagSel;
    //std::cout<<" event has been selected after btag selection "<<btagSelected*lumiNormalization_<<std::endl;

    if(zllptSelected == true) ++zllptSel;
    //std::cout<<" event has been selected after zll ptselection "<<btagSelected*lumiNormalization_<<std::endl;

    if(metSelected == true) ++metSel;
    //std::cout<<" event has been selected after met selection "<<metSelected*lumiNormalization_<<std::endl;

    if(jjdrSelected == true) ++jjdrSel;
    //std::cout<<" event has been selected after jjdr selection "<<jjdrSelected*lumiNormalization_<<std::endl;

    if(hmassSelected == true) ++hmassSel;
    //std::cout<<" event has been selected after higgs mass selection "<<hmassSelected*lumiNormalization_<<std::endl;

}

void H2l2bSelection::endJob() {
  // insert counter and write in a txt file
  std::cout<<std::endl;
  std::cout<<" Selection yields for H2b2l in the "+ channel +" channel: "<< std::endl;
  std::cout<<" events after skimming: "<<skimEvents*lumiNormalization_<<std::endl;
  std::cout<<" events after mass selection: "<<massSel*lumiNormalization_<<std::endl;
  std::cout<<" events after btag selection: "<<btagSel*lumiNormalization_<<std::endl;
  std::cout<<" events after zllpt selection: "<<zllptSel*lumiNormalization_<<std::endl;
  std::cout<<" events after met selection: "<<metSel*lumiNormalization_<<std::endl;
  std::cout<<" events after jjdr selection: "<<jjdrSel*lumiNormalization_<<std::endl;
  std::cout<<" events after hmass selection: "<<hmassSel*lumiNormalization_<<std::endl;

  outFile<<" Selection yields for H2b2l in the "+channel+" channel: "<< std::endl;
  outFile<<" events after skimming: "<<skimEvents*lumiNormalization_<<std::endl;
  outFile<<" events after mass selection: "<<massSel*lumiNormalization_<<std::endl;
  outFile<<" events after btag selection: "<<btagSel*lumiNormalization_<<std::endl;
  outFile<<" events after zllpt selection: "<<zllptSel*lumiNormalization_<<std::endl;
  outFile<<" events after met selection: "<<metSel*lumiNormalization_<<std::endl;
  outFile<<" events after jjdr selection: "<<jjdrSel*lumiNormalization_<<std::endl;
  outFile<<" events after hmass selection: "<<hmassSel*lumiNormalization_<<std::endl;

  outFile.close();

}

DEFINE_FWK_MODULE(H2l2bSelection);
