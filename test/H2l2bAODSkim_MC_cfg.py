import FWCore.ParameterSet.Config as cms

process = cms.Process("H2l2bAODSkim")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
    'file:/data3/scratch/cms/mc/Summer11/GluGluToHToZZTo2L2Q_M-300/F6D80D8E-3794-E011-9800-00215E93D944.root'
     )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
## global tag for MC
process.GlobalTag.globaltag = 'START42_V17::All'

# Muon filter 
process.goodMuons = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("muons"),
  cut = cms.string('pt > 10'),
  filter = cms.bool(True)
)

# dimuon filters
process.dimuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodMuons"),
    minNumber = cms.uint32(2)
)

# Electron filter
process.goodElectrons = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("gsfElectrons"),
  cut = cms.string('pt > 10'),
  filter = cms.bool(True)                                
)

# dielectron filters
process.dielectronFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodElectrons"),
    minNumber = cms.uint32(2)
)

# Jet filter
process.goodJets = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("ak5PFJets"),
  cut = cms.string('pt > 20'),
  filter = cms.bool(True)                                
)

# dijet filters
process.dijetFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("goodJets"),
    minNumber = cms.uint32(2)
)

process.dijetSequence = cms.Sequence(
    process.goodJets *
    process.dijetFilter
    )

# Skim paths
process.dimuonPath = cms.Path(
    process.dijetSequence *
    process.goodMuons *
    process.dimuonFilter 
    )

process.dielectronPath = cms.Path(
    process.dijetSequence *
    process.goodElectrons *
    process.dielectronFilter
    )

process.emuPath = cms.Path(
    process.dijetSequence *
    process.goodElectrons *
    process.goodMuons
    )


# Output module configuration
from Configuration.EventContent.EventContent_cff import *


H2l2bAODSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring()
)
H2l2bAODSkimEventContent.outputCommands.extend(AODSIMEventContent.outputCommands)

H2l2bAODSkimEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
           'dimuonPath',
           'dielectronPath',
           'emuPath')
    )
)

process.H2l2bAODSkimOutputModule = cms.OutputModule("PoolOutputModule",
    H2l2bAODSkimEventContent,
    H2l2bAODSkimEventSelection,
    fileName = cms.untracked.string('h2l2bAODSkim.root')
)

process.outpath = cms.EndPath(process.H2l2bAODSkimOutputModule)


