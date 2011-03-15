/// @file
/// This file contains the declaration of the abstract class to make the selection and check
/// plots of the ALPGEN samples for HF matching (and avoid double-counting).
/// ---------------------------------------
/// Originally written by Oscar Gonzalez (26/IV/2009)
/// Modified for CMS usage by Miguel Vidal (21/X/2010)

#ifndef __ALPGENCCSAMPLESFLAG_HH
#define __ALPGENCCSAMPLESFLAG_HH

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

#include <string>
using namespace std;
using namespace edm;
using namespace reco;

class AlpgenCCSamplesFlag : public edm::EDProducer {

  // Parameters

  float _singleCPtCut;  ///< Cut on the pt of the single charm contributions.

  float _deltaRPairs; ///< Cut on the Delta-R for the CC and BB pair contributions.

  float _pairsPtCut;  ///< Cut on the pt of the quarks in CC and BB pairs to be 
                      ///< considered as "a pair".

  // Internal functions
  
  /// EDM Filter functions
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  /// Make the processing for events with a single C, returning
  /// true if the event has to be rejected.
  /// It also fills the corresponding plots for that kind of event.
  bool processSingleC (GenParticleCollection gen);

  /// Makes the selection for single-C events, using the pt of the charm quark.
  /// Returns true if the event is rejected.
  bool cuttingSingleC (float pt);


  /// Make the processing for events with CC pairs, returning
  /// true if the event has to be rejected.
  /// It also fills the corresponding plots for that kind of event.
  bool processCCPair (GenParticleCollection gen);

  /// Makes the selection for CC-pair events, using the maximum of the minimum
  /// Dr for the valid pairs.
  /// Returns true if the event is rejected.
  bool cuttingCCPair (float dr);

  /// Make the processing for events with CC pairs, returning
  /// true if the event has to be rejected.
  /// It also fills the corresponding plots for that kind of event.
  bool processBBPair (GenParticleCollection gen);

  /// Makes the selection for BB-pair events, using the maximum of the minimum
  /// Dr for the valid pairs.
  /// Returns true if the event is rejected.
  bool cuttingBBPair (float dr);

protected:


public:

  /// Destructor.
  ~AlpgenCCSamplesFlag (void);

  /// Constructor.
  explicit AlpgenCCSamplesFlag (const edm::ParameterSet&);


  // /////////////////////////////////////////////
  // Accesors
  // /////////////////////////////////////////////

  /// Returns the value of the cut used for the Pt of the single-charm
  /// for the matching.
  float getSingleCPtCut (void) const {return _singleCPtCut;}

  /// Returns the value of the pt cut for quarks in a pair for the pair
  /// to be a good one used in the calculation of the maximum in Delta-R.
  float getPtCutforPairs (void) const {return _pairsPtCut;}

  /// Returns the value of the Delta-R cuts for performing the matching
  /// of the Heavy-flavour samples using CC and BB pairs in the event.
  float getPairDeltaRCut (void) const {return _deltaRPairs;}

};

#endif
//=======================================================================
