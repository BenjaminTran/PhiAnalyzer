import FWCore.ParameterSet.Config as cms

process = cms.Process("PhiAnalyzer")

# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhiAnalyzer.PhiAnalyzer.PhiSelector_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_FullSkim_v4/170803_222621/0000/pPb_HM_111.root'
    #'root://cmsxrootd.fnal.gov//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_Pbp_V0Cascade_v1/170301_202152/0000/pPb_HM_105.root'
    #),
    #secondaryFileNames = cms.untracked.vstring(
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/3A1D69BC-01B7-E611-AD9A-FA163E98E135.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/3A8F14C0-F9B6-E611-9FDF-02163E0133B5.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/5878662E-F6B6-E611-A11A-FA163E02A339.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/AED413F1-F6B6-E611-AEBE-02163E0125FC.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/C408E634-F6B6-E611-AB11-FA163E592268.root',
    #    '/store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/286/010/00000/DECFA2B3-F7B6-E611-B733-FA163EB68B63.root'
        )
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('Phi_dedx_v1.root')
)

process.test = cms.Sequence(process.PhiSelector)

process.p = cms.Path(process.test)