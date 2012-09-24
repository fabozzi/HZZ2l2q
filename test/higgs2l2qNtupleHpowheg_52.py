import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIM")

applyFilter = False

# turn on when running on powheg signal MC (-> to get generated Higgs mass)
isPowhegSignal = True
#hMassHyp = "700"
#comEn = "8TeV"
#fileWeight = "MMozer/powhegweight/data/mZZ_Higgs"+hMassHyp+"_"+comEn+"_W.txt_I.txt"

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#from MMozer.powhegweight.tongguang600 import *
#process.powWeightProducer = tongguangweights600.clone(
#    filename = cms.FileInPath(fileWeight) )

process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2qedmNtuples_52_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.MessageLogger.cerr.threshold = ''
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource")

process.source.fileNames=cms.untracked.vstring(
    'file:h2l2qSkimData.root'
)


#### Event cleaning 
process.badEventFilter = cms.EDFilter(
    "HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","PAT"),
    HLTPaths = cms.vstring('primaryVertexFilterPath',
#                           'CSCTightHaloFilterPath',
#                           'EcalDeadCellTriggerPrimitiveFilterPath',
                           'EcalDeadCellBoundaryEnergyFilterPath',
                           'noscrapingFilterPath',          
                           'hcalLaserEventFilterPath',
                           'HBHENoiseFilterPath'#,
#                           'totalKinematicsFilterPath' #only for Madgraph MC
                           ),
    eventSetupPathsKey = cms.string(''),
    # how to deal with multiple triggers: True (OR) accept if ANY is true, False
    # (AND) accept if ALL are true
    andOr = cms.bool(False),
    throw = cms.bool(True)  # throw exception on unknown path names
    )

process.genInfoNtuple = cms.EDProducer(
    "GenInfoNtupleDump",
    src = cms.InputTag("genParticles")  
    )

# select lepton daughters from the Higgs
process.genSelectorZDaughterLeptons = cms.EDFilter(
    "GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string(' ((abs(pdgId)==11 ) || (abs(pdgId)==13 )) && abs(mother.pdgId)==23 ')
) 

# select quark daughters from the Higgs
process.genSelectorZDaughterQuarks = cms.EDFilter(
    "GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string(' ((abs(pdgId)>=1 ) && (abs(pdgId)<=6 )) && abs(mother.pdgId)==23 ')
) 

# perform MC matching
process.hCandMatch = cms.EDProducer(
    "HiggsMatcher",
    elhiggs = cms.InputTag("hzzeejj:h"),
    muhiggs = cms.InputTag("hzzmmjj:h"),
    genLept = cms.InputTag("genSelectorZDaughterLeptons"),
    genQuarks = cms.InputTag("genSelectorZDaughterQuarks")
) 

process.PUInfoNtuple = cms.EDProducer(
    "GenPUNtupleDump",
    isData = cms.bool(False)
)

# Event rho dumper
process.rhoDumper = cms.EDProducer("EventRhoDumper",
                                    rho = cms.InputTag("kt6PFJets:rho"),
                                    restrictedRho = cms.InputTag("kt6PFJetsForIso:rho")
                                    )

# Met variables producer
process.metInfoProducer = cms.EDProducer("MetVariablesProducer",
                                    metTag = cms.InputTag("patMETsAK5"),
                                    t1CorrMetTag = cms.InputTag("patType1CorrectedPFMetAK5")
                                    )

process.edmNtuplesOut.fileName = cms.untracked.string('h2l2q_ntuple.root')
process.edmNtuplesOut.outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_hCandMatch_*_*",
    "keep *_genInfoNtuple_*_*",
    "keep *_eventVtxInfoNtuple_*_*",
    "keep *_PUInfoNtuple_*_*",
    "keep *_rhoDumper_*_*",
    "keep *_metInfoProducer_*_*",
# keep rho recommended for muoniso
    'keep *_kt6PFJetsCentralNeutral_rho_*',
    "keep *_Higgs2e2bEdmNtuple_*_*",
    "keep *_Higgs2mu2bEdmNtuple_*_*",
    "keep *_Higgsemu2bEdmNtuple_*_*",
    "keep *_jetinfos_*_*"
)

if applyFilter:
    process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('analysisPath')
        )

process.edmNtuplesOut.dropMetaData = cms.untracked.string('ALL')

process.edmNtuplesOut.outputCommands.extend([
    'keep edmTriggerResults_TriggerResults_*_HLT',
    'keep *_TriggerResults_*_PAT',
    ])

process.analysisPath = cms.Path()

if applyFilter:
    process.analysisPath += process.badEventFilter

if isPowhegSignal:
    process.analysisPath += process.genInfoNtuple
    process.analysisPath += process.genSelectorZDaughterLeptons
    process.analysisPath += process.genSelectorZDaughterQuarks
    process.analysisPath += process.hCandMatch


#process.analysisPath = cms.Path(
#    process.badEventFilter +
process.analysisSequence = cms.Sequence(
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.rhoDumper+
    process.metInfoProducer+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple+
    process.Higgsemu2bEdmNtuple+
    process.jetinfos
)

process.analysisPath += process.analysisSequence

process.endPath = cms.EndPath(process.edmNtuplesOut)

process.schedule = cms.Schedule(
    process.analysisPath,
    process.endPath
    )
