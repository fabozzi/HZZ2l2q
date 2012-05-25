#! /bin/csh -f

cd ../../

cvs co -r cbern_03May12 -d AnalysisDataFormats/CMGTools UserCode/CMG/AnalysisDataFormats/CMGTools
                    
cvs co -r cbern_03May12 -d CMGTools/Common  UserCode/CMG/CMGTools/Common                               
cvs co -r V00-00-09 -d CMGTools/External  UserCode/CMG/CMGTools/External                              
cvs co -r cbern_03May12 -d CMGTools/Production UserCode/CMG/CMGTools/Production
                             
cvs co -r cbern_03May12  -d CMGTools/RootTools UserCode/CMG/CMGTools/RootTools                              
cvs co -r V00-03-07-01 CommonTools/ParticleFlow
                   
cvs co -r cbern_25Apr12 DataFormats/PatCandidates

cvs co -r V00-00-08 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools

cvs co -r Shervin20022012_2011Jan16ReRec0_and_Shervin05032012_Fall11MC_smearing_V01 -d EgammaCalibratedGsfElectrons/CalibratedElectronAlgos UserCode/EGamma/EgammaCalibratedGsfElectrons/CalibratedElectronAlgos

cvs co -r Shervin20022012_2011Jan16ReRec0_and_Shervin05032012_Fall11MC_smearing_V01 -d EgammaCalibratedGsfElectrons/CalibratedElectronProducers UserCode/EGamma/EgammaCalibratedGsfElectrons/CalibratedElectronProducers

cvs co -r V00-00-60  FWCore/GuiBrowsers                               

cvs co -r Colin_TaggingMode_June30 JetMETAnalysis/ecalDeadCellTools 

cvs co -r V08-09-02 PhysicsTools/PatAlgos                            

cvs co -r V03-09-18-03 PhysicsTools/PatUtils                            

cvs co -r V08-03-16 PhysicsTools/Utilities                           

cvs co -r lhx_14APR2012_v1 RecoMET/METFilters                               

cvs co -r V15-01-05 RecoParticleFlow/PFProducer                      

cvs co -r wreece_020512 RecoParticleFlow/PostProcessing UserCode/RecoParticleFlow/PostProcessing

cd EGamma/EGammaAnalysisTools/data

cat download.url | xargs wget

cd ../../../


