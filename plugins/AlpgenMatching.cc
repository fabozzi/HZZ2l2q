/// @file
/// File containing the definition of the methods associated to the class.
///

#include "../interface/AlpgenMatching.h"

// CMS-based include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

// ROOT classes

#include "TH1.h"



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
AlpgenMatching::AlpgenMatching (const edm::ParameterSet& iConfig) :
  edm::EDProducer(),

  doPlots_(iConfig.getUntrackedParameter<bool>("doPlots",false)),
  plotWeight_((float) iConfig.getUntrackedParameter<double>("plotWeight",1)),

  singleCPtCut_((float) iConfig.getUntrackedParameter<double>("singleCPtCut",20)),
  deltaRCut_((float) iConfig.getUntrackedParameter<double>("deltaRCut",1.0)),
  ptPairCut_((float) iConfig.getUntrackedParameter<double>("ptPairCut",20))
// Constructor of the class
{
  _nEvents=0;
  _nPassingEvents=0;

  for (int i=0;i<4;++i) {
    _nEventsLF[i]=0;
    _nEventsSC[i]=0;
    _nEventsCC[i]=0;
    _nEventsBB[i]=0;
  }

  // Variables/branches to be produced:
  produces<unsigned int>("alpgenflavorflag").setBranchAlias("alpgenflavorflag");
  produces<float>("alpgenptsinglec").setBranchAlias("alpgenptsinglec");
  produces<float>("alpgendrcc").setBranchAlias("alpgendrcc");
  produces<float>("alpgendrbb").setBranchAlias("alpgendrbb");
}

//-----------------------------------------------------------------------
AlpgenMatching::~AlpgenMatching (void)
// Destructor of the class
{


}

//-----------------------------------------------------------------------
void AlpgenMatching::beginJob (void)
// Method runs for the EDAnalyzer at the beginning of the job.
{
  cout<<"Start-up of AlpgenMatching: "<<endl;

  if (doPlots_) bookHistograms();
}

//-----------------------------------------------------------------------
void AlpgenMatching::endJob (void)
// Method run for the EDAnalyzer at the end of the job.
{
  cout<<"--------------------------------------------------------------------------"<<endl;
  cout<<"Report from AlpgenMatching: "<<endl;
  cout<<"   - Producing plots: "<<outbool(doPlots_)<<" "<<plotWeight_<<endl;
  cout<<"   - Number of processed events: "<<_nEvents<<endl;
  cout<<"   - Number of events with good matching: "<<_nPassingEvents<<endl;
  cout<<endl;
  cout<<"                  Events-LF  Events-SC   Events-CC  Events-BB"<<endl;
  cout<<"- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"<<endl;
  cout<<" Matched as LF:  "<<setw(7)<<_nEventsLF[0]<<"    "<<setw(7)<<_nEventsSC[0]<<"     "<<setw(7)<<_nEventsCC[0]<<"     "<<setw(7)<<_nEventsBB[0]<<endl;
  cout<<" Matched as sC:  "<<setw(7)<<_nEventsLF[1]<<"    "<<setw(7)<<_nEventsSC[1]<<"     "<<setw(7)<<_nEventsCC[1]<<"     "<<setw(7)<<_nEventsBB[1]<<endl;
  cout<<" Matched as CC:  "<<setw(7)<<_nEventsLF[2]<<"    "<<setw(7)<<_nEventsSC[2]<<"     "<<setw(7)<<_nEventsCC[2]<<"     "<<setw(7)<<_nEventsBB[2]<<endl;
  cout<<" Matched as BB:  "<<setw(7)<<_nEventsLF[3]<<"    "<<setw(7)<<_nEventsSC[3]<<"     "<<setw(7)<<_nEventsCC[3]<<"     "<<setw(7)<<_nEventsBB[3]<<endl;
  cout<<"--------------------------------------------------------------------------"<<endl;
}

//-----------------------------------------------------------------------
void AlpgenMatching::produce (edm::Event& iEvent, const edm::EventSetup& iSetup)
// Method run for each event in the analysis.
{
  ++_nEvents;

  reco::GenParticleRefVector charm;
  reco::GenParticleRefVector anticharm;
  reco::GenParticleRefVector bottom;
  reco::GenParticleRefVector antibottom;

  float ptmax_b=-1;
  float ymin_c=-1;
  float ymin_b=-1;

  // Variables for storing in the Event
  std::auto_ptr<unsigned int> flavorflag(new unsigned int);
  *flavorflag=0;

  std::auto_ptr<float> pt_singlec(new float);
  *pt_singlec = -1;
  std::auto_ptr<float> dr_cc(new float);
  *dr_cc = -1;
  std::auto_ptr<float> dr_bb(new float);
  *dr_bb = -1;

  // For production in ALPGEN (status code==3)
  
  int ncharms=0;
  int nbottoms=0;

  float ptmin_prodc=-1;
  float ymax_prodc=-1;

  float ptmin_prodb=-1;
  float ymax_prodb=-1;

  // We process the information

  //cout<<"-------------------------------------------------------"<<endl;

  edm::Handle<reco::GenParticleCollection> genpartH;
  iEvent.getByLabel("genParticles",genpartH);

  //  GenParticleCollection genpart = *genpartH;
  for (reco::GenParticleCollection::const_iterator xpart = genpartH->begin();
       xpart!=genpartH->end();++xpart) {

    // For charm and bottom with status==3 we do some counting to try to
    // discover which class of event it is.

    if (xpart->status()==3) {
      int hepid = abs(xpart->pdgId());

      if (hepid==4) {
	if (ptmin_prodc<-0.5 || xpart->pt()<ptmin_prodc) ptmin_prodc=xpart->pt(); 
	if (ymax_prodc<-0.5 || fabs(xpart->rapidity())>ymax_prodc) ymax_prodc=fabs(xpart->rapidity());
	//cout<<"CHARM QUARK: "<<hepid<<" "<<xpart->mass()<<" "<<xpart->pt()<<" "<<xpart->rapidity()<<endl;
	if (xpart->mass()>1.0) ++ncharms;  // Only massive quarks
      }
      else if (hepid==5) {
	if (ptmin_prodb<-0.5 || xpart->pt()<ptmin_prodb) ptmin_prodb=xpart->pt();
        if (ymax_prodb<-0.5 || fabs(xpart->rapidity())>ymax_prodb) ymax_prodb=fabs(xpart->rapidity());

	//cout<<"BOTTOM QUARK: "<<xpart->pdgId()<<" "<<xpart->mass()<<endl;
	if (xpart->mass()>3.0) ++nbottoms;
      }
    }

    // We need to identify the charm and the bottom at parton level. The simplest
    // is to look for the strings/clusters and scan their mothers.

    if (xpart->pdgId()==92 || xpart->pdgId()==91) {
      reco::GenParticleRefVector mothers = xpart->motherRefVector();
      for (reco::GenParticleRefVector::const_iterator xmot = mothers.begin();
	   xmot!=mothers.end();++xmot) {

	// We store the found b's and c's:
	if ((*xmot)->pdgId()==4) {
	  //cout<<"charm: "<<(*xmot)->pt()<<" "<<(*xmot)->rapidity()<<endl;
	  charm.push_back(*xmot);
	  if ((*xmot)->pt()>(*pt_singlec)) (*pt_singlec)=(*xmot)->pt();
	  if (ymin_c<0 || fabs((*xmot)->rapidity())<ymin_c) ymin_c=fabs((*xmot)->rapidity());
	}
	else if ((*xmot)->pdgId()==-4) {
	  //cout<<"anticharm: "<<(*xmot)->pt()<<" "<<(*xmot)->rapidity()<<endl;
	  anticharm.push_back(*xmot);
          if ((*xmot)->pt()>(*pt_singlec)) (*pt_singlec)=(*xmot)->pt();
          if (ymin_c<0 || fabs((*xmot)->rapidity())<ymin_c) ymin_c=fabs((*xmot)->rapidity());
	}
	else if ((*xmot)->pdgId()==5) {
	  bottom.push_back(*xmot);
          if ((*xmot)->pt()>ptmax_b) ptmax_b=(*xmot)->pt();
          if (ymin_b<0 || fabs((*xmot)->rapidity())<ymin_b) ymin_b=fabs((*xmot)->rapidity());
	}
        else if ((*xmot)->pdgId()==-5) {
	  antibottom.push_back(*xmot);
	  if ((*xmot)->pt()>ptmax_b) ptmax_b=(*xmot)->pt();
          if (ymin_b<0 || fabs((*xmot)->rapidity())<ymin_b) ymin_b=fabs((*xmot)->rapidity());
	}
      }
    }
  }

  // check to identify HF samples
  if (nbottoms>2 || ncharms>2) {
    cerr<<"HEPTF-ERROR: Problems identifying the type of sample... matching may be wrong!"<<endl;
  }

  // Filling control histograms:

  if (doPlots_) {
    if ( (*pt_singlec)>0) {
      _control.ptmax_c->Fill((*pt_singlec),plotWeight_);
      _control.ymin_c->Fill(ymin_c,plotWeight_);
    }

    if (ptmax_b>0) {
      _control.ptmax_b->Fill(ptmax_b,plotWeight_);
      _control.ymin_b->Fill(ymin_b,plotWeight_);
    }

    if (ptmin_prodc>-0.5) {
      _control.ptmin_prodc->Fill(ptmin_prodc,plotWeight_);
      _control.ymax_prodc->Fill(ymax_prodc,plotWeight_);
    }

    if (ptmin_prodb>-0.5) {
      _control.ptmin_prodb->Fill(ptmin_prodb,plotWeight_);
      _control.ymax_prodb->Fill(ymax_prodb,plotWeight_);
    }
  }

  // Some checks:

  if (abs(charm.size()-anticharm.size())>1 || bottom.size()!=antibottom.size()) {
    cerr<<"HEPTF-ERROR: Problems with the number of charms and bottoms: "
	<<charm.size()<<" "<<anticharm.size()<<" "<<bottom.size()<<" "
	<<antibottom.size()<<endl;
  }

  // We process the variables we use in the matching:

  if (charm.size()!=anticharm.size()) {  // single-c production (with W samples)
    // These may come from the single-C or the LF depending on the pt of the
    // highest-pt charm, independently of the rest.

    // The variable to use is (*pt_singlec):

    if (doPlots_) _control.pt_singlec->Fill((*pt_singlec),plotWeight_);
    
    if ((*pt_singlec)>singleCPtCut_) *flavorflag |=0x1;  // Single-C must come from W+C samples.
  }
  else {
    (*pt_singlec)=-1;  // No single-c 

    // For the usual pair production the argument is for each quark to look
    // for the closest antiquark (in y-phi), and use the maximum of the distance.
    // For both, charm and bottom:

    *dr_cc = getDistanceForPairs(charm,anticharm);
    if ((*dr_cc)>0) {
      if (doPlots_) _control.dr_cc->Fill(*dr_cc,plotWeight_);
      if ((*dr_cc)>deltaRCut_) *flavorflag |=0x2;   // CC must come from W/Z+CC samples.
      //if ((*dr_cc)>6) cout<<"MIERDA "<<endl;
    }

    *dr_bb = getDistanceForPairs(bottom,antibottom);
    if ((*dr_bb)>0) {
      if (doPlots_) _control.dr_bb->Fill(*dr_bb,plotWeight_);
      if ((*dr_bb)>deltaRCut_) {
	*flavorflag |=0x4;   // BB must come from W/Z+BB samples.
      }
    }
  }

  // Usign the values, we decide the "code".
  unsigned int icode=0; // LF
  if ( (*flavorflag)&0x04 ) icode=3;  // Hard-BB
  else if ( (*flavorflag)&0x02 ) icode=2;  // Hard-CC
  else if ( (*flavorflag)&0x01 ) icode=1;  // single-C

  // We store the code which is a bit pattern:
  //    bit 0.... single-c with high pt.
  //    bit 1.... CC pair above Dyphi and pt threshold.
  //    bit 2.... BB pair above Dyphi and pt threshold.
  //
  // With this code: BB samples take events with code>=4 
  //                 CC samples take events with code>=2 and code<4
  //                 Single-c samples take events with code==1  (>=1 && <2)
  //                 Light-flavour samples take events with code==0 (>=0 && <1)
  // 
  // To simplify usage... bit 3 is set to 1 if for the sample we need to take the event).

  if (nbottoms==2) {  // It is as BB sample
    if ( (*flavorflag)>=4) *flavorflag |=0x8;
    ++_nEventsBB[icode];
  }
  else if (ncharms==2) {  // It is as CC sample
    if ( (*flavorflag)>=2 && (*flavorflag)<4 ) *flavorflag |=0x8;
    ++_nEventsCC[icode];
  }
  else if (ncharms==1) {  // It is a single-C sample
    if ( (*flavorflag)==1) *flavorflag |=0x8;
    ++_nEventsSC[icode];
  }
  else {  // It is a light-flavour sample
    if ( (*flavorflag)==0) *flavorflag |=0x8;
    ++_nEventsLF[icode];
  }

  if ((*flavorflag)>=8) ++_nPassingEvents;  // Matching accepted the event.
  //  else cout<<"MIERDA "<<endl;


  // Storing the final code:

  iEvent.put(flavorflag,"alpgenflavorflag");
  
  // We also store the variables we use in the matching
  iEvent.put(pt_singlec,"alpgenptsinglec");
  iEvent.put(dr_cc,"alpgendrcc");
  iEvent.put(dr_bb,"alpgendrbb");
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
float AlpgenMatching::getDistanceForPairs (const reco::GenParticleRefVector &quark,
						 const reco::GenParticleRefVector &antiquark) const
// Gets the maximum distance for quark/antiquark pairs (matched with the minimum
// distance.
{
  float drmax=-1; 

  // We loop over the quark list and get the minimum distance to an antiquark:

  //  cout<<"HOLA "<<quark.size()<<" "<<antiquark.size()<<endl;

  for (reco::GenParticleRefVector::const_iterator xq = quark.begin();
       xq!=quark.end();++xq) {
    float drmin=-1;

    bool valid=true;

    for (reco::GenParticleRefVector::const_iterator xaq = antiquark.begin();
	 xaq!=antiquark.end();++xaq) {

      //BUG      // One of the two quarks must pass the cut on pt
      //BUG if ((*xq)->pt()<ptPairCut_ && (*xaq)->pt()<ptPairCut_) continue;

//      float dr = (xq->rapidity()-xaq->rapidity());
//
//      float dphi = fabs(xq->phi()-xaq->phi());
//      if (dphi>M_PI) dphi = 2*M_PI-dphi;
//
//      dr = dr*dr + dphi*dphi;  // We use the squared
      float dr = deltaR2 ((*xq)->eta(),(*xq)->phi(),
			  (*xaq)->eta(),(*xaq)->phi());  // Using the squared
      
      //      cout<<"       "<<sqrt(dr)<<" "<<(*xq)->pt()<<" "<<(*xq)->eta()<<" "<<(*xq)->phi()<<" "<<(*xaq)->pt()<<" "<<(*xaq)->eta()<<" "<<(*xaq)->phi()<<endl;
      
      // For good pairs we will have to check whether the minimum is a pair
      // including one high-pt guy.
      // Note: setting the && to || makes the TWO to pass the pt requirement.

      if (dr<drmin || drmin<0) {
	valid=true;
	//if ((*xq)->pt()<ptPairCut_ && (*xaq)->pt()<ptPairCut_) valid=false;
	if ((*xq)->pt()<ptPairCut_ || (*xaq)->pt()<ptPairCut_) valid=false;
	drmin=dr;
      }
    }
    
    // Once that for each quark we have the minimum... we consider it for the
    // maximum if one of the two passes the cut on pt
    if (valid && drmin>drmax) drmax=drmin;
  }

  if (drmax<0) return drmax;
  return sqrt(drmax);   // We return the actual DR value.
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void AlpgenMatching::bookHistograms (void)
// Booking of the histograms to be used in the class when doing analysis.
{
  edm::Service<TFileService> fs;

  // Control histograms:

  _control.ptmax_c = fs->make<TH1D>("ptmax_charm","Pt of the charm (anti)quark with highest pt",120,0.0,600.0);
  _control.ymin_c = fs->make<TH1D>("ymin_charm","Rapidity (absolute) of the charm (anti)quark with smallest abs(y)",60,0.0,6.0);

  _control.ptmax_b = fs->make<TH1D>("ptmax_bottom","Pt of the bottom (anti)quark with highest pt",120,0.0,600.0);
  _control.ymin_b = fs->make<TH1D>("ymin_bottom","Rapidity (absolute) of the bottom (anti)quark with smallest abs(y)",60,0.0,6.0);

  _control.ptmin_prodc = fs->make<TH1D>("ptmin_prodcharm","Pt of the Produced charm (anti)quark with lowest pt",120,0.0,240.0);
  _control.ymax_prodc = fs->make<TH1D>("ymax_prodcharm","Rapidity (absolute) of the Produced charm (anti)quark with highest abs(y)",60,0.0,6.0);

  _control.ptmin_prodb = fs->make<TH1D>("ptmin_prodbottom","Pt of the Produced bottom (anti)quark with lowest pt",120,0.0,240.0);
  _control.ymax_prodb = fs->make<TH1D>("ymax_prodbottom","Rapidity (absolute) of the Produced bottom (anti)quark with highest abs(y)",60,0.0,6.0);

  // Variables for the matching:

  _control.pt_singlec = fs->make<TH1D>("pt_singlec","Pt of the leading charm in single-c production events",120,0.0,300.0);
  _control.dr_cc = fs->make<TH1D>("dr_cc","Dr of the hardest charm-anticharm pair in events with CC",120,0.0,12.0);
  _control.dr_bb = fs->make<TH1D>("dr_bb","Dr of the hardest bottom-antibottom pair in events with BB",120,0.0,12.0);
}

//-----------------------------------------------------------------------
//define this as a plug-in
DEFINE_FWK_MODULE(AlpgenMatching);
//=======================================================================
