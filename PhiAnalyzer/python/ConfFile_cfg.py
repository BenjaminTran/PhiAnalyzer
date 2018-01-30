import FWCore.ParameterSet.Config as cms

process = cms.Process("PhiAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhiAnalyzer.PhiAnalyzer.PhiSelector_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v4/170803_222621/0000/pPb_HM_111.root'
    )
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('histo.root')
)

process.test = cms.Sequence(process.PhiSelector)

process.p = cms.Path(process.test)
