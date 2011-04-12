/// @file
/// This file contains the declaration of the AlpgenMatching, used to get
/// the basic information (and plots) that can be used to perform the matching
/// of the light-flavour and Heavy-flavour samples for ALPGEN.
/// <PRE>
/// Written by O. Gonzalez (5/IV/2011)
/// </PRE>

#ifndef __ALPGENMATCHING__H_
#define __ALPGENMATCHING__H_

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

// CMS-based classes

class Event;
class ParameterSet;
class EventSetup;

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Root classes

class TH1D;


// C++ classes

#include <string>


//-----------------------------------------------------------------------
/// This class is used to get the basic information (and plots) that can be used
/// to perform the matching of the light-flavour and Heavy-flavour samples
/// for ALPGEN.

class AlpgenMatching : public edm::EDProducer {

  // Internal variables

  // Configuration variables through Python

  bool doPlots_;  ///< To fill the histograms.

  float singleCPtCut_;  ///< Matching value for the single-C samples.
  float deltaRCut_;  ///< Matching value for Delta-R for the CC and BB pairs.

  float ptPairCut_;  ///< Cut of the high-pt member of a pair for the pair to pass the cut.

  // Variables for statistics and final report

  int _nEventsLF[4];  ///< Number of events with LF, single-C, hard-CC o hard-BB.
  int _nEventsSC[4];  ///< Number of events with LF, single-C, hard-CC o hard-BB.
  int _nEventsCC[4];  ///< Number of events with LF, single-C, hard-CC o hard-BB.
  int _nEventsBB[4];  ///< Number of events with LF, single-C, hard-CC o hard-BB.

  int _nEvents;  ///< Number of processed events.

  int _nPassingEvents;  ///< Number of events that are accepted in matching.


  /// Pointer to histograms (for having some packing we use a structure).
  struct {
    TH1D *ptmax_c;
    TH1D *ymin_c;
    TH1D *ptmax_b;
    TH1D *ymin_b;

    TH1D *pt_singlec;
    TH1D *dr_cc;
    TH1D *dr_bb;
  } _control;


  // Internal methods

  /// Method runs for the EDAnalyzer at the beginning of the job.
  virtual void beginJob (void);

  /// Method run for each event in the analysis.
  virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup);

  /// Method run for the EDAnalyzer at the end of the job.
  virtual void endJob (void);

  /// Booking of the histograms to be used in the class when doing analysis.
  void bookHistograms (void);

  /// Gets the maximum distance for quark/antiquark pairs (matched with the minimum
  /// distance.
  float getDistanceForPairs (const reco::GenParticleRefVector &quark, 
			     const reco::GenParticleRefVector &antiquark) const;


public:

  /// Destructor of the class.
  virtual ~AlpgenMatching (void);

  /// Constructor of the class.
  explicit AlpgenMatching (const edm::ParameterSet &iConfig);


  // /////////////////////////////////////////////
  // Setters
  // /////////////////////////////////////////////



  // ////////////////////////////////////////////
  // Orders
  // /////////////////////////////////////////////

  // /////////////////////////////////////////////
  // Accesors
  // /////////////////////////////////////////////




};
#endif

//====================================================================
