#include "HiggsAnalysis/Higgs2l2b/interface/DiJetKinFitter.h"

#include "HiggsAnalysis/Higgs2l2b/interface/TFitConstraintM.h"
#include "HiggsAnalysis/Higgs2l2b/interface/TFitParticleEtEtaPhi.h"
#include "HiggsAnalysis/Higgs2l2b/interface/TKinFitter.h"



DiJetKinFitter::DiJetKinFitter( const TString&  name, const TString& title, double mass ) : TKinFitter( name, title ) {

  name_ = name;
  title_ = title;
  mass_ = mass;

}



std::pair<TLorentzVector,TLorentzVector> DiJetKinFitter::fit(TLorentzVector jet1, TLorentzVector jet2 ) {

  this->reset();

  TMatrixD m_jet1(3,3);
  TMatrixD m_jet2(3,3);

  m_jet1(0,0) = 0.5*this->ErrEt (jet1.Et(), jet1.Eta()); // et
  m_jet1(1,1) = 0.5*this->ErrEta(jet1.Et(), jet1.Eta()); // eta
  m_jet1(2,2) = 0.5*this->ErrPhi(jet1.Et(), jet1.Eta()); // phi
  m_jet2(0,0) = 0.5*this->ErrEt (jet2.Et(), jet2.Eta()); // et
  m_jet2(1,1) = 0.5*this->ErrEta(jet2.Et(), jet2.Eta()); // eta
  m_jet2(2,2) = 0.5*this->ErrPhi(jet2.Et(), jet2.Eta()); // phi

  TFitParticleEtEtaPhi *fitJet1 = new TFitParticleEtEtaPhi( "Jet1", "Jet1", &jet1, &m_jet1 );
  TFitParticleEtEtaPhi *fitJet2 = new TFitParticleEtEtaPhi( "Jet2", "Jet2", &jet2, &m_jet2 );
  
  TFitConstraintM *mCons_jets = new TFitConstraintM( "ZMassConstraint_jets", "ZMass-Constraint", 0, 0 , mass_);
  mCons_jets->addParticles1( fitJet1, fitJet2 );

  
  this->addMeasParticle( fitJet1 );
  this->addMeasParticle( fitJet2 );
  this->addConstraint( mCons_jets );
    

  //Set convergence criteria
  this->setMaxNbIter( 30 );
  this->setMaxDeltaS( 1e-2 );
  this->setMaxF( 1e-1 );
  this->setVerbosity(0);

  //Perform the fit
  TKinFitter::fit();

  TLorentzVector jet1_kinfit(*fitJet1->getCurr4Vec());
  TLorentzVector jet2_kinfit(*fitJet2->getCurr4Vec());

  std::pair<TLorentzVector,TLorentzVector> jetpair;
  jetpair.first  = jet1_kinfit;
  jetpair.second = jet2_kinfit;

  delete fitJet1 ; delete fitJet2;delete mCons_jets ;

  return jetpair;

}



double DiJetKinFitter::ErrEt(double Et, double Eta) {

  double InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 5.6;
    b = 1.25;
    c = 0.033;
  }
  else{
    a = 4.8;
    b = 0.89;             
    c = 0.043;
  }   
  InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
  return InvPerr2;        
}



double DiJetKinFitter::ErrEta(double Et, double Eta) {
  double InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 1.215;
    b = 0.037;
    c = 7.941 * 0.0001;
  }
  else{
    a = 1.773;
    b = 0.034;
    c = 3.56 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}


double DiJetKinFitter::ErrPhi(double Et, double Eta) {
  double InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }
  else{
    a = 2.908;
    b = 0.021;
    c = 2.59 * 0.0001;
  }
  InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}


