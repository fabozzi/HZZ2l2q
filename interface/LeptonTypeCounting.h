/// @file
/// This file contains the declaration of the LeptonTypeCounting class
/// which counts the events with electrons, muons and taus 
/// in the Hard process (status=3).
/// <PRE>
/// Written by O. Gonzalez (12/IV/2011)
/// </PRE>

#ifndef __LEPTONTYPECOUNTING__H_
#define __LEPTONTYPECOUNTING__H_

#include "FWCore/Framework/interface/EDFilter.h"

// CMS-based classes

#include "FWCore/Utilities/interface/InputTag.h"

class Event;
class ParameterSet;
class EventSetup;

//OLD #include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Root classes

class TH1D;


// C++ classes

#include <string>


//-----------------------------------------------------------------------
/// This class is used to get the counting of events with electrons, muons
/// and taus at the hard process.

class LeptonTypeCounting : public edm::EDFilter {

  // Internal variables

  // Configuration variables through Python

  bool requireElectron_;   ///< To enable the requirement of electrons for the event to pass.
  bool requireMuon_;       ///< To enable the requirement of muons for the event to pass.
  bool requireTau_;        ///< To enable the requirement of taus for the event to pass.

  // Variables for statistics and final report

  int _nEvents;  ///< Number of processed events.

  int _nElectronEvents;  ///< Number of events with electrons.
  int _nMuonEvents;  ///< Number of events with muons.
  int _nTauEvents;  ///< Number of events with taus.


  // Internal methods

  /// Method runs for the EDAnalyzer at the beginning of the job.
  virtual void beginJob (void);

  /// Method run for each event in the analysis.
  virtual bool filter (edm::Event &iEvent, const edm::EventSetup &iSetup);

  /// Method run for the EDAnalyzer at the end of the job.
  virtual void endJob (void);

  /// Booking of the histograms to be used in the class when doing analysis.
  void bookHistograms (void);


public:

  /// Destructor of the class.
  virtual ~LeptonTypeCounting (void);

  /// Constructor of the class.
  explicit LeptonTypeCounting (const edm::ParameterSet &iConfig);


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
