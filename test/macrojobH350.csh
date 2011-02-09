#!/bin/csh

cd /afs/cern.ch/user/d/decosa/scratch0/Higgs/CMSSW_3_8_6/src/HiggsAnalysis/Higgs2l2b/test
eval `scramv1 runtime -csh`

bkgh2l2b0100



rfcp "Z0Jet.root" /castor/cern.ch/user/d/decosa/Higgs/Z0jet/histos/Z0Jet.root
rfcp "Z1Jet.root" /castor/cern.ch/user/d/decosa/Higgs/Z1jet/histos/Z1Jet.root
rfcp "Z2Jet.root" /castor/cern.ch/user/d/decosa/Higgs/Z2jet/histos/Z2Jet.root
rfcp "Z3Jet.root" /castor/cern.ch/user/d/decosa/Higgs/Z3jet/histos/Z3Jet.root
rfcp "Z4Jet.root" /castor/cern.ch/user/d/decosa/Higgs/Z4jet/histos/Z4Jet.root
rfcp "Z5Jet.root" /castor/cern.ch/user/d/decosa/Higgs/Z5jet/histos/Z5Jet.root





