#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "FWCore/Framework/interface/ESHandle.h"



using namespace edm;
using namespace std;
using namespace reco;

class Higgs2l2bUserDataSF : public edm::EDProducer {
public:
  Higgs2l2bUserDataSF( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag higgsTag;
  float Btageff_SF_;
};



Higgs2l2bUserDataSF::Higgs2l2bUserDataSF( const ParameterSet & cfg ):
  higgsTag( cfg.getParameter<InputTag>( "higgs" ) ),
  Btageff_SF_(1.0)
{
  produces<vector<pat::CompositeCandidate> >("h").setBranchAlias( "h" );
}

void Higgs2l2bUserDataSF::produce( Event & evt, const EventSetup & iSetup ) {

  ///////////////////////////////
  ///   Begin DB setup
  ///////////////////////////////
  
  //// This is needed for the DB
  
  std::map<std::string,PerformanceResult::ResultType> measureMap;
  measureMap["BTAGBEFF"]=PerformanceResult::BTAGBEFF;
  measureMap["BTAGBERR"]=PerformanceResult::BTAGBERR;
  measureMap["BTAGCEFF"]=PerformanceResult::BTAGCEFF;
  measureMap["BTAGCERR"]=PerformanceResult::BTAGCERR;
  measureMap["BTAGLEFF"]=PerformanceResult::BTAGLEFF;
  measureMap["BTAGLERR"]=PerformanceResult::BTAGLERR;
  measureMap["BTAGNBEFF"]=PerformanceResult::BTAGNBEFF;
  measureMap["BTAGNBERR"]=PerformanceResult::BTAGNBERR;
  measureMap["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR;
  measureMap["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  measureMap["BTAGCEFFCORR"]=PerformanceResult::BTAGCEFFCORR;
  measureMap["BTAGCERRCORR"]=PerformanceResult::BTAGCERRCORR;
  measureMap["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  measureMap["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;
  measureMap["BTAGNBEFFCORR"]=PerformanceResult::BTAGNBEFFCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["MUEFF"]=PerformanceResult::MUEFF;
  measureMap["MUERR"]=PerformanceResult::MUERR;
  measureMap["MUFAKE"]=PerformanceResult::MUFAKE; 
  measureMap["MUEFAKE"]=PerformanceResult::MUEFAKE;
  
  edm::ESHandle<BtagPerformance> perfH;
  std::vector<std::string> measureName;
  std::vector<std::string> measureType;
  
  // Define which Btag and Mistag algorithm you want to use. These are not user defined and need to be exact
  measureName.push_back("MISTAGTCHEM");
  measureName.push_back("BTAGTCHEM");
  measureName.push_back("MISTAGTCHEL");
  measureName.push_back("BTAGTCHEL");
  measureName.push_back("MISTAGTCHEM");
  measureName.push_back("MISTAGTCHEL");
  
  // Tell DB you want the SF. These are not user defined and need to be exact
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFF");
  measureType.push_back("BTAGLEFF");
  
  // These are user defined maps that we will use to store the SF 
  std::map<std::string,float> ScaleFactors_j0;   //store the Btag and Mistag SF for jet0
  std::map<std::string,float> ScaleFactors_j1;   //store the Btag and Mistag SF for jet1
  std::map<std::string,float> ScaleFactorsEff_j0;   //store the Mistag eff for jet0
  std::map<std::string,float> ScaleFactorsEff_j1;   //store the Mistag eff for jet1
  
  
  
  ///////////////////////////////
  ///   End DB setup
  ///////////////////////////////
  

  Handle<std::vector<pat::CompositeCandidate> > higgsH;
  evt.getByLabel(higgsTag, higgsH);

  auto_ptr<vector<pat::CompositeCandidate> > higgsColl( new vector<pat::CompositeCandidate> () );

   for (unsigned int i = 0; i< higgsH->size();++i){
    const pat::CompositeCandidate & H = (*higgsH)[i];
    edm::Ref<std::vector<pat::CompositeCandidate> > hRef(higgsH, i);
    pat::CompositeCandidate h(H);

    //    const reco::Candidate * Zll = H.daughter(0);
    const pat::Jet *  j1 = dynamic_cast<const pat::Jet*>(H.daughter(1)->daughter(0)->masterClone().get());
    const pat::Jet *  j2 = dynamic_cast<const pat::Jet*>(H.daughter(1)->daughter(1)->masterClone().get());
    
    double j1et = j1->et();
    double j1eta = j1->eta(); 
    double j2et = j2->et(); 
    double j2eta = j2->eta(); 

    //    cout << "j1et = " << j1et << endl;
    //    cout << "j1eta = " << j1eta << endl;
    //    cout << "j2et = " << j2et << endl;
    //    cout << "j2eta = " << j2eta << endl;

    for( size_t iMeasure = 0; iMeasure < measureName.size(); iMeasure++ ) {
      //Setup our measurement
      iSetup.get<BTagPerformanceRecord>().get( measureName[ iMeasure ],perfH);
      const BtagPerformance & perf = *(perfH.product());
      BinningPointByMap measurePoint;
      measurePoint.reset();
      measurePoint.insert(BinningVariables::JetEt, j1et );    ///// pass in the et of the jet
      measurePoint.insert(BinningVariables::JetAbsEta, abs( j1eta ) ); ///// pass in the absolute eta of the jet
      
      // Extract the mistag eff value
      if ( measureType[ iMeasure ] == "BTAGLEFF") {
	std::string suffix = "eff"; // add a suffix eff so we can distingiush it from other values
	ScaleFactorsEff_j0[ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
      }
      else{ // Extract the mistag and btag SF
	// The factor Btageff_SF_ is used for Btagging systematics and should be set to 1.0 as default
	ScaleFactors_j0[ measureName[ iMeasure ] ] = Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);  
      }
      
      // Reset point and do for second jet

      measurePoint.reset();
      measurePoint.insert(BinningVariables::JetEt, j2et );
      measurePoint.insert(BinningVariables::JetAbsEta, abs( j2eta ));
      
      if ( measureType[ iMeasure ] == "BTAGLEFF") {
	std::string suffix = "eff";
	ScaleFactorsEff_j1[ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
      }
      else{
	ScaleFactors_j1[ measureName[ iMeasure ] ] = Btageff_SF_*perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);  
      }

    }

    //    cout << "SFLj1= " << ScaleFactors_j0["BTAGTCHEL"] <<  endl;
    //    cout << "SFMj1= " << ScaleFactors_j0["BTAGTCHEM"] <<  endl;
    //    cout << "SFLj2= " << ScaleFactors_j1["BTAGTCHEL"] <<  endl;
    //    cout << "SFMj2= " << ScaleFactors_j1["BTAGTCHEM"] <<  endl;
    //    cout << "SFmisLj1= " << ScaleFactors_j0["MISTAGTCHEL"] <<  endl;
    //    cout << "SFmisMj1= " << ScaleFactors_j0["MISTAGTCHEM"] <<  endl;
    //    cout << "SFmisLj2= " << ScaleFactors_j1["MISTAGTCHEL"] <<  endl;
    //    cout << "SFmisMj2= " << ScaleFactors_j1["MISTAGTCHEM"] <<  endl;

    h.addUserFloat("j1SFBtagTCHEL", ScaleFactors_j0["BTAGTCHEL"]);
    h.addUserFloat("j1SFBtagTCHEM", ScaleFactors_j0["BTAGTCHEM"]);
    h.addUserFloat("j1SFMistagTCHEL", ScaleFactors_j0["MISTAGTCHEL"]);
    h.addUserFloat("j1SFMistagTCHEM", ScaleFactors_j0["MISTAGTCHEM"]);
    h.addUserFloat("j1SFMistagTCHELeff", ScaleFactorsEff_j0["MISTAGTCHELeff"]);
    h.addUserFloat("j1SFMistagTCHEMeff", ScaleFactorsEff_j0["MISTAGTCHEMeff"]);

    h.addUserFloat("j2SFBtagTCHEL", ScaleFactors_j1["BTAGTCHEL"]);
    h.addUserFloat("j2SFBtagTCHEM", ScaleFactors_j1["BTAGTCHEM"]);
    h.addUserFloat("j2SFMistagTCHEL", ScaleFactors_j1["MISTAGTCHEL"]);
    h.addUserFloat("j2SFMistagTCHEM", ScaleFactors_j1["MISTAGTCHEM"]);
    h.addUserFloat("j2SFMistagTCHELeff", ScaleFactorsEff_j1["MISTAGTCHELeff"]);
    h.addUserFloat("j2SFMistagTCHEMeff", ScaleFactorsEff_j1["MISTAGTCHEMeff"]);


    higgsColl->push_back(h);
  }

  
  evt.put( higgsColl, "h");

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( Higgs2l2bUserDataSF );

