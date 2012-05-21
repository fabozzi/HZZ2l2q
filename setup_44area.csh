#! /bin/csh -f

cd ../../

cvs co -r cbern_metreg_10May12 -d AnalysisDataFormats/CMGTools UserCode/CMG/AnalysisDataFormats/CMGTools
                    
cvs co -r colin_metregclean_10May12 -d CMGTools/Common  UserCode/CMG/CMGTools/Common                               

cvs co -r V00-02-02 -d CMGTools/External  UserCode/CMG/CMGTools/External                              
cvs co -r colin_metregclean_10May12 -d CMGTools/H2TauTau UserCode/CMG/CMGTools/H2TauTau

cvs co -r cbern_09May12 -d CMGTools/Production UserCode/CMG/CMGTools/Production

cvs co -r cbern_init_06May12 -d CMGTools/RootTools UserCode/CMG/CMGTools/RootTools                              

cvs co -r V00-03-05-08 CommonTools/ParticleFlow

cvs co -r V00-04-01 CondFormats/EgammaObjects                        

cvs co -r lhx_12JAN2012_v1 DataFormats/METReco                              

cvs co -r V06-05-01 DataFormats/PatCandidates

cvs co -r V00-00-08 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools

cvs co -r cbern_09May12 -d EGamma/EGammaAnalysisToolsSiXie UserCode/sixie/EGamma/EGammaAnalysisTools

cvs co -r Shervin20022012_2011Jan16ReRec0_and_Shervin05032012_Fall11MC_smearing_V01 -d EgammaCalibratedGsfElectrons/CalibratedElectronAlgos UserCode/EGamma/EgammaCalibratedGsfElectrons/CalibratedElectronAlgos

cvs co -r Shervin20022012_2011Jan16ReRec0_and_Shervin05032012_Fall11MC_smearing_V01 -d EgammaCalibratedGsfElectrons/CalibratedElectronProducers UserCode/EGamma/EgammaCalibratedGsfElectrons/CalibratedElectronProducers

cvs co -r V00-00-63 FWCore/GuiBrowsers                               

cvs co -r V00-11-00-01 GeneratorInterface/GenFilters  

cvs co -r Colin_TaggingMode_June30 JetMETAnalysis/ecalDeadCellTools 

cvs co -r V00-00-09 -d Muon/MuonAnalysisTools UserCode/sixie/Muon/MuonAnalysisTools

cvs co -r V11-03-16 PhysicsTools/HepMCCandAlgos                      

cvs co -r V01-05-07 PhysicsTools/IsolationAlgos 

cvs co -r V08-07-44 PhysicsTools/PatAlgos                            

cvs co -r V03-09-18-03 PhysicsTools/PatUtils                            

cvs co -r V08-03-16 PhysicsTools/Utilities                           

cvs co -r V00-03-32 RecoEgamma/ElectronIdentification

cvs co -r V03-03-07 RecoLuminosity/LumiDB

cvs co -r V00-00-08 RecoMET/METAnalyzers 

cvs co -r lhx_14APR2012_v1 RecoMET/METFilters                               

cvs co -r wreece_020512 -d RecoParticleFlow/PostProcessing UserCode/RecoParticleFlow/PostProcessing

cvs co -r joseNov14 -d TauAnalysis/SVFitStandAlone UserCode/TauAnalysis/SVFitStandAlone

cvs up -r 1.53 PhysicsTools/PatAlgos/python/tools/tauTools.py 

cd EGamma/EGammaAnalysisTools/data

cat download.url | xargs wget

cd ../../../

