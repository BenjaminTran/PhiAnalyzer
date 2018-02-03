import os
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.workArea = 'crab_dir/Phi_DeDx_2016pPb'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = os.path.expandvars('$CMSSW_BASE/src/PhiAnalyzer/PhiAnalyzer/test/phianalysis_cfg.py')

config.section_("Data")
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
#config.Data.primaryDataset = 'MinBias_TuneZ2star_7TeV_pythia6'
config.Data.splitting = 'EventAwareLumiBased'
#config.Data.totalUnits = 660000
config.Data.unitsPerJob = 50000
config.Data.outLFNDirBase = '/store/group/phys_heavyions/btran/Phi'
config.Data.useParent = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'PhiDeDx_Mass'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_US_MIT']
#config.Site.whitelist = ['T2_US_Vanderbilt']