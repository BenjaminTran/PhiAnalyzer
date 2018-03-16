from PhiAnalyzer.PhiAnalyzer.phianalysis_cfg import process, cms

# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            'root://cmsxrootd.fnal.gov//store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/AODSIM/MB_80X_mcRun2_pA_v4-v2/120000/003FC5B1-2009-E711-8AD2-0025904C6414.root'
            ),
        secondaryFileNames = cms.untracked.vstring(
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/525A69B0-0609-E711-BE75-0025905D1D52.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/58345DA2-FF08-E711-9068-0025904C51FE.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/AEDB5C26-0109-E711-BAE2-0025905D1D52.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/CCC2E6D2-0609-E711-8D78-0025905C54BA.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/D2552255-0009-E711-8C4E-0025904CDDEC.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/D4FBBEA2-FF08-E711-A81F-0025905C54B8.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/D6A2BBA2-FF08-E711-A5AE-0025905C54C4.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/E65FB2E9-0D09-E711-80BE-0025905C5500.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/F4A3FDA5-FF08-E711-ADE8-0025905C3E36.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/186A8744-0109-E711-8503-0025905C53A4.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/2CAE0D25-0109-E711-A4BE-0025905C5488.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/38C069A2-FF08-E711-8BDC-0025905C2CD0.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/504F8A24-0109-E711-94AC-0025904C641C.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120000/B495F73B-5A08-E711-9144-0025905C53D0.root',
'/store/himc/pPb816Summer16DR/ReggeGribovPartonMC_EposLHC_pPb_4080_4080_DataBS/GEN-SIM-RAW/MB_80X_mcRun2_pA_v4-v2/120001/16BBA6D5-0609-E711-B217-0025905C5430.root'
            )
        )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('PhiGenMatch_v1.root')
        )

process.test = cms.Sequence(process.PhiGenMatch)

process.p = cms.Path(process.test)
