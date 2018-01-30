import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('PhiAnalyzer'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
