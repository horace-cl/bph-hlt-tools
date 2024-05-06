L1_seeds = ['L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2',
                                'L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6',
                                'L1_DoubleMu5_SQ_OS_dR_Max1p6',
                                'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6',
                                'L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2',
                                'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6',
                                'L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6',
                                'L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5',
                                'L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4',
                                'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4',
                                'L1_DoubleMu4p5_SQ_OS_dR_Max1p2',
                                'L1_DoubleMu4_SQ_OS_dR_Max1p2', 
                                'L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4',
                                
                                'L1_DoubleMu0er1p5_SQ_dR_Max1p4',
                                'L1_DoubleMu0er2p0_SQ_dEta_Max1p6',
                                'L1_DoubleMu0er2p0_SQ_dEta_Max1p5',

                                'L1_SingleMu10_SQ14_BMTF',
                                'L1_SingleMu11_SQ14_BMTF',
                                'L1_SingleMu0_BMTF',

                                ]

HLT_Paths=[ "HLT_DoubleMu4_3_LowMass_v",
            "HLT_DoubleMu4_3_LowMass_SS_v",
            "HLT_DoubleMu4_LowMass_Displaced_v",
            "HLT_DoubleMu4_MuMuTrk_Displaced_v",
            "HLT_DoubleMu4_3_Bs_v",
            "HLT_Dimuon10_y1p4_v", # This may not be in the menu!
            "HLT_DoubleMu4_3_Jpsi_v",
            "HLT_Mu8_v", 
            "HLT_IsoMu24_v",
            "HLT_IsoMu27_v",
            "HLT_Mu4_L1DoubleMu_v",
            "HLT_Mu0_L1DoubleMu_v",
            "HLT_Mu3_PFJet40_v",
            "HLT_Mu15_v",
            "HLT_DoubleMu2_Jpsi_LowPt_v",
            
            "HLT_Mu0_Barrel_v",
            "HLT_Mu0_Barrel_L1HP10_v",
            "HLT_Mu0_Barrel_L1HP11_v",
            "HLT_Mu9_Barrel_L1HP10_IP6_v",
            "HLT_Mu10_Barrel_L1HP11_IP6_v",
            "HLT_Dimuon10_Upsilon_y1p4_v",            
]


# HLT_Paths=[ "HLT_DoubleMu4_3_LowMass_v",
#             "HLT_DoubleMu4_3_LowMass_SS_v",
#             "HLT_DoubleMu4_LowMass_Displaced_v",
#             "HLT_DoubleMu4_MuMuTrk_Displaced_v",
#             "HLT_DoubleMu4_3_Bs_v",
#             "HLT_DoubleMu2_Jpsi_LowPt",

#             #"HLT_Dimuon10_y1p4_v", # This may not be in the menu!
#             #"HLT_DoubleMu4_3_Jpsi_v",
#             #"HLT_Mu8_v", 
#             #"HLT_IsoMu24_v",
#             #"HLT_IsoMu27_v",
#             #"HLT_Mu4_L1DoubleMu_v",
#             #"HLT_Mu0_L1DoubleMu_v",
#             #"HLT_Mu3_PFJet40_v",
#             #"HLT_Mu15_v",            
            
#             "HLT_Mu0_Barrel_v",
#             "HLT_Mu0_Barrel_L1HP10_v",
#             "HLT_Mu0_Barrel_L1HP11_v",
#             "HLT_Mu9_Barrel_L1HP10_IP6_v",
#             "HLT_Mu10_Barrel_L1HP11_IP6_v",
#             "HLT_Dimuon10_Upsilon_y1p4_v8",            
# ]


import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)

options.register('maxE', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Maximum number of events"
)

options.register('debug', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Debugging"
)

options.parseArguments()





process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4') # for 2022
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_PromptAnalysis_v2') # for 2022
#process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_HLT_v2')      #2023 B
#process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_HLT_v2')  #for CMSSW >= 13_0_3 
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v2')  #for CMSSW >= 13_0_3 


## Message Logger and Event range
process.MessageLogger.cerr.FwkReport.reportEvery = 10000 
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxE))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#data_file = '/store/data/Run2023C/ParkingDoubleMuonLowMass7/MINIAOD/PromptReco-v4/000/367/770/00000/0c62c10b-6629-4caf-9f6f-ed6c0bbc65b0.root'
#data_file = "/store/data/Run2023C/Muon0/MINIAOD/PromptReco-v4/000/367/770/00000/2c7a455b-304f-4ece-b66a-dabf30e9c7b1.root"
#data_file = "/store/data/Run2023D/Muon1/MINIAOD/PromptReco-v2/000/370/776/00000/a92b94e8-1455-4a5d-b4c9-323e43d486f9.root"
#data_file = "/store/data/Run2022F/Muon/MINIAOD/PromptReco-v1/000/360/389/00000/ad0997b9-ff20-4b2c-9c51-1d6ef49100f4.root"
#data_file = "/store/data/Run2022F/Muon/MINIAOD/PromptReco-v1/000/360/335/00000/db3a7d95-2b78-4e72-86e4-8436005406bf.root"
#data_file  = "/store/data/Run2022F/Muon/MINIAOD/PromptReco-v1/000/360/390/00000/be5c66b6-fea2-48f4-879f-846f6de0e511.root"
data_file  = "/store/data/Run2024B/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/379/058/00000/d57fe8ca-ccb0-4df9-a027-d6fa9788b51d.root" #2024B
data_file  = "/store/data/Run2024B/ParkingSingleMuon0/MINIAOD/PromptReco-v1/000/379/252/00000/75e2e99e-80b2-429a-9d60-1e7e4c9682fd.root"


if options.isMC:
    data_file = '/store/mc/Run3Summer22EEMiniAODv3/ButoJpsiK_Jpsito2Mu_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2520000/04f4baff-0d98-4f64-a0e2-0421d62db508.root'
    data_file = '/store/mc/Run3Summer22MiniAODv4/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v1/2550000/28491152-7072-4c6d-ac38-6e20c57ae346.root'
    data_file = '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/0009bf67-3c7a-4e81-b50a-e3914b3d2ffa.root'
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #/ParkingDoubleMuonLowMass*/Run2023*-PromptReco-v4/MINIAOD
        #'/store/data/Run2023C/ParkingDoubleMuonLowMass7/MINIAOD/PromptReco-v4/000/367/770/00000/0c62c10b-6629-4caf-9f6f-ed6c0bbc65b0.root', # Data
       #'/store/mc/Run3Summer22EEMiniAODv3/ButoJpsiK_Jpsito2Mu_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2520000/04f4baff-0d98-4f64-a0e2-0421d62db508.root', #MonteCarlo
       data_file

        
 )
)

# process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
#                                         triggerConditions = cms.vstring(#'HLT_Dimuon25_Jpsi_v*',
#                                                                         #'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
#                                                                         #'HLT_DoubleMu4_MuMuTrk_Displaced_v*',
#                                                                         #'HLT_DoubleMu4_JpsiTrk_Bc_v*',
#                                                                         #'HLT_DoubleMu4_LowMass_Displaced_v*',
#                                                                         'HLT_DoubleMu4_Jpsi_Displaced_v*'                                   
#                                                                        ),
#                                         hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
#                                         l1tResults = cms.InputTag( "" ),
#                                         throw = cms.bool(False)
#                                         )

#process.load("myAnalyzers.JPsiKsPAT.PsikaonRootupler_cfi")
#process.rootuple.dimuons = cms.InputTag('slimmedMuons')
process.rootuple = cms.EDAnalyzer('MuMu',
                          #stageL1Trigger = 2,
                          
                          muons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          #Trak_lowpt = cms.InputTag("lostTracks"),
                          #GenParticles = cms.InputTag("genParticles"),
                          GenParticles = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          
                          TriggerInput = cms.InputTag("slimmedPatTrigger"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          algInputTag = cms.InputTag("gtStage2Digis", "", "RECO"),
                          l1Muons = cms.InputTag("gmtStage2Digis", "Muon", "RECO"),
                          HLTPaths = cms.vstring(HLT_Paths),
                          L1Seeds = cms.vstring(L1_seeds),
                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(options.isMC),
                          OnlyGen = cms.bool(False),
                          mumuMassConstraint = cms.bool(False),       
                          mumuMasscut     = cms.vdouble(0,15),                                  
                          Trkmass           = cms.double(0.493677),
                          #Trkmass           = cms.double(0.13957018),
                          BarebMasscut      = cms.vdouble(4.2,6.8),
                          bMasscut          = cms.vdouble(5.0,6.0),
                          debug = cms.bool(options.debug)        
                          )

from Configuration.Eras.Modifier_stage2L1Trigger_cff import stage2L1Trigger
stage2L1Trigger.toModify(process.rootuple, stageL1Trigger = 2)


dataset_name = data_file.split('/')[4][:11]
file_name = f'Rootuple_DiMu-MiniAOD_{dataset_name}.root'
if options.isMC:
    file_name = f'Rootuple_MC_DiMu-MiniAOD_{dataset_name}.root'

    
process.TFileService = cms.Service("TFileService",
       #fileName = cms.string('Rootuple_Butomumu_2023-MiniAOD.root'),
       fileName = cms.string(file_name),
)

process.mySequence = cms.Sequence(
                                   #process.triggerSelection *
                                   process.rootuple
				   )

#process.p = cms.Path(process.mySequence)

#process.p = cms.Path(process.triggerSelection*process.rootuple)
process.p = cms.Path(process.rootuple)
#process.p = cms.Path(process.rootuple)


