/// @file
/// File containing the definition of the methods associated to the class.
///

#include "../interface/LeptonTypeCounting.h"

// CMS-based include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//OLD #include "FWCore/ServiceRegistry/interface/Service.h"

//OLD #include "CommonTools/UtilAlgos/interface/TFileService.h"

//OLD #include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT classes

//OLD #include "TH1.h"



// C++ classes

#include <memory>

#include <vector>
using std::vector;

#include <iostream>
using std::endl;
using std::cout;
using std::cerr;
using std::setw;

#include <string>
using std::string;

// Defines:

#define outbool(x) (x?"true":"false")

//-----------------------------------------------------------------------
LeptonTypeCounting::LeptonTypeCounting (const edm::ParameterSet& iConfig) :
  edm::EDFilter(),
  
  requireElectron_(iConfig.getUntrackedParameter<bool>("requireElectron",false)),
  requireMuon_(iConfig.getUntrackedParameter<bool>("requireMuon",false)),
  requireTau_(iConfig.getUntrackedParameter<bool>("requireTau",false))

// Constructor of the class
{
  _nEvents=0;

  _nElectronEvents=0;
  _nMuonEvents=0;
  _nTauEvents=0;
}

//-----------------------------------------------------------------------
LeptonTypeCounting::~LeptonTypeCounting (void)
// Destructor of the class
{


}

//-----------------------------------------------------------------------
void LeptonTypeCounting::beginJob (void)
// Method runs for the EDAnalyzer at the beginning of the job.
{
  cout<<"Start-up of LeptonTypeCounting: "<<endl;
}

//-----------------------------------------------------------------------
void LeptonTypeCounting::endJob (void)
// Method run for the EDAnalyzer at the end of the job.
{
  cout<<"--------------------------------------------------------------------------"<<endl;
  cout<<"Report from LeptonTypeCounting: "<<endl;
  cout<<"   - Number of processed events: "<<_nEvents<<endl<<endl;
  cout<<"   - Events with electrons (Hard Process): "<<_nElectronEvents<<" [Required: "<<outbool(requireElectron_)<<"]"<<endl;
  cout<<"   - Events with muons (Hard Process): "<<_nMuonEvents<<" [Required: "<<outbool(requireMuon_)<<"]"<<endl;
  cout<<"   - Events with taus (Hard Process): "<<_nTauEvents<<" [Required: "<<outbool(requireTau_)<<"]"<<endl;
  cout<<"\n   - Total with matches: "<<_nElectronEvents+_nMuonEvents+_nTauEvents<<endl;
  cout<<"--------------------------------------------------------------------------"<<endl;
}

//-----------------------------------------------------------------------
bool LeptonTypeCounting::filter (edm::Event& iEvent, const edm::EventSetup& iSetup)
// Method run for each event in the analysis.
{
  ++_nEvents;

  bool hasElectron=false;
  bool hasMuon=false;
  bool hasTau=false;

  // We process the information

  edm::Handle<reco::GenParticleCollection> genpartH;
  iEvent.getByLabel("genParticles",genpartH);

  for (reco::GenParticleCollection::const_iterator xpart = genpartH->begin();
       xpart!=genpartH->end();++xpart) {

    // We just scan the hard-process of the table. As soon as status() is not 3, we abort.

    if (xpart->status()!=3) break;

    int hepid = abs(xpart->pdgId());

    // checking the lepton types:

    if (!hasElectron && hepid==11) {hasElectron=true; ++_nElectronEvents;}
    else if (!hasMuon && hepid==13) {hasMuon=true; ++_nMuonEvents;}
    else if (!hasTau && hepid==15) {hasTau=true; ++_nTauEvents;}
  }

  // Filtering events:

  if (requireElectron_ && !hasElectron) return false;
  if (requireMuon_ && !hasMuon) return false;
  if (requireTau_ && !hasTau) return false;

  return true;  // Event is accepted!
}

//-----------------------------------------------------------------------
//OLD //-----------------------------------------------------------------------
//OLD float LeptonTypeCounting::getDistanceForPairs (const reco::GenParticleRefVector &quark,
//OLD 						 const reco::GenParticleRefVector &antiquark) const
//OLD // Gets the maximum distance for quark/antiquark pairs (matched with the minimum
//OLD // distance.
//OLD {
//OLD   float drmax=-1; 
//OLD 
//OLD   // We loop over the quark list and get the minimum distance to an antiquark:
//OLD 
//OLD   //  cout<<"HOLA "<<quark.size()<<" "<<antiquark.size()<<endl;
//OLD 
//OLD   for (reco::GenParticleRefVector::const_iterator xq = quark.begin();
//OLD        xq!=quark.end();++xq) {
//OLD     float drmin=-1;
//OLD 
//OLD     for (reco::GenParticleRefVector::const_iterator xaq = antiquark.begin();
//OLD 	 xaq!=antiquark.end();++xaq) {
//OLD 
//OLD       // One of the two quarks must pass the cut on pt
//OLD 
//OLD       if ((*xq)->pt()<ptPairCut_ && (*xaq)->pt()<ptPairCut_) continue;
//OLD 
//OLD //      float dr = (xq->rapidity()-xaq->rapidity());
//OLD //
//OLD //      float dphi = fabs(xq->phi()-xaq->phi());
//OLD //      if (dphi>M_PI) dphi = 2*M_PI-dphi;
//OLD //
//OLD //      dr = dr*dr + dphi*dphi;  // We use the squared
//OLD       float dr = deltaR2 ((*xq)->eta(),(*xq)->phi(),
//OLD 			  (*xaq)->eta(),(*xaq)->phi());  // Using the squared
//OLD       
//OLD       //      cout<<"       "<<sqrt(dr)<<" "<<(*xq)->pt()<<" "<<(*xq)->eta()<<" "<<(*xq)->phi()<<" "<<(*xaq)->pt()<<" "<<(*xaq)->eta()<<" "<<(*xaq)->phi()<<endl;
//OLD       
//OLD       if (dr<drmin || drmin<0) drmin=dr;
//OLD     }
//OLD 
//OLD     if (drmin>drmax) drmax=drmin;
//OLD   }
//OLD 
//OLD   if (drmax<0) return drmax;
//OLD   return sqrt(drmax);   // We return the actual DR value.
//OLD }
//OLD 
//OLD //-----------------------------------------------------------------------
//OLD //-----------------------------------------------------------------------
//OLD void LeptonTypeCounting::bookHistograms (void)
//OLD // Booking of the histograms to be used in the class when doing analysis.
//OLD {
//OLD   edm::Service<TFileService> fs;
//OLD 
//OLD   //OLD  TFileDirectory mainanal = fs->mkdir(string("LeptonTypeCounting"));
//OLD 
//OLD   // Control histograms:
//OLD 
//OLD   _control.ptmax_c = fs->make<TH1D>("ptmax_charm","Pt of the charm (anti)quark with highest pt",120,0.0,600.0);
//OLD   _control.ymin_c = fs->make<TH1D>("ymin_charm","Rapidity (absolute) of the charm (anti)quark with smallest abs(y)",60,0.0,6.0);
//OLD 
//OLD   _control.ptmax_b = fs->make<TH1D>("ptmax_bottom","Pt of the bottom (anti)quark with highest pt",120,0.0,600.0);
//OLD   _control.ymin_b = fs->make<TH1D>("ymin_bottom","Rapidity (absolute) of the bottom (anti)quark with smallest abs(y)",60,0.0,6.0);
//OLD 
//OLD   // Variables for the matching:
//OLD 
//OLD   _control.pt_singlec = fs->make<TH1D>("pt_singlec","Pt of the leading charm in single-c production events",120,0.0,300.0);
//OLD   _control.dr_cc = fs->make<TH1D>("dr_cc","Dr of the hardest charm-anticharm pair in events with CC",120,0.0,12.0);
//OLD   _control.dr_bb = fs->make<TH1D>("dr_bb","Dr of the hardest bottom-antibottom pair in events with BB",120,0.0,12.0);
//OLD }

//-----------------------------------------------------------------------
//define this as a plug-in
DEFINE_FWK_MODULE(LeptonTypeCounting);
//=======================================================================
