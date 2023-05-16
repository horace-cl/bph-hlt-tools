from CRABClient.UserUtilities import config
config = config()

tag   = '_7Mar23'
#inDS = '/ZeroBias/Run2023B-PromptReco-v1/MINIAOD'
inDS  = '/ParkingDoubleMuonLowMass1/Run2023B-PromptReco-v1/MINIAOD'
#inDS  = '/ParkingDoubleElectronLowMass/Run2023B-PromptReco-v1/MINIAOD'


config.General.requestName = inDS.split('/')[1]+'_'+inDS.split('/')[2]

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName   = '../test/PsikaonRootupler.py'

config.Data.inputDataset = inDS
#config.Data.inputDataset = '/BcToJpsiPi_TuneCP5_14TeV_pythia8_evtgen/Run3Summer21MiniAOD-Pilot_120X_mcRun3_2021_realistic_v5_ext1-v2/MINIAODSIM'
config.Data.splitting   = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.publication = True
# This string is used to construct the output dataset name
config.General.workArea      = 'BPHTriggerNtuples'
config.Data.outputDatasetTag = inDS.split('/')[1]+tag
# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = './jsonFolder/Cert_Collisions2022_355100_362760_Golden.json'

# Select input data based on run-ranges
#config.Data.runRange = '362616'

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'
#T2_MX_CINVESTAV