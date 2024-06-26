#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
source /cvmfs/cms.cern.ch/common/crab-setup.sh
"""
import os
from optparse import OptionParser

import CRABClient
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException #python2
from http.client import HTTPException #python3
#from CRABClient.UserUtilities import config

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()
    
    # The submit command needs special treatment.
    if options.crabCmd == 'submit':
        
        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()
        
        #jobname = 'Bsmumuphi_miniAOD_LowMass1_2023Cv0'
        
        config.General.requestName = None
        config.General.workArea = 'TriggerEfficiencies_MC_Try1'
        config.General.transferOutputs = True
        config.General.transferLogs = False

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = '../test/mumuRootupler.py'
        config.JobType.pyCfgParams = ['isMC=True']
        #config.JobType.allowUndistributedCMSSW = True
        
        config.Data.inputDataset = None
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 7
        config.Data.totalUnits  = 500
        #config.Data.lumiMask = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt'
        
        config.Data.outLFNDirBase = '/store/user/hcrottel/'+str(config.General.workArea)
        config.Data.publication = False
        config.Data.outputDatasetTag = None
        config.Site.storageSite = 'T3_CH_CERNBOX'
        #config.Site.whitelist = ['T2_US*']
        #config.Data.ignoreLocality = True
        #config.Site.storageSite = None # Choose your site.
        
        #--------------------------------------------------------
        # Will submit one task for each of these input datasets.
        inputDatasets = [
            "/ButoJpsiK_Jpsito2Mu_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM",
            "/BuToKJPsiMuMu_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v2-v2/MINIAODSIM",
            "/K0sToMuMu_K0sFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1/MINIAODSIM",
            "/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5/MINIAODSIM",
            #'/ParkingDoubleMuonLowMass0/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass1/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass2/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass3/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass4/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass5/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass6/Run2023C-PromptReco-v4/MINIAOD',
            #'/ParkingDoubleMuonLowMass7/Run2023C-PromptReco-v4/MINIAOD',
        ]
        
        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            #config.General.requestName = inDS.split('/')[1]+'-'+inDS.split('/')[2]
            config.General.requestName = inDS.split('/')[1]
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
            # Submit.
            try:
                print ( "Submitting for input dataset %s" % (inDS) )
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print ( "Submission for input dataset %s failed: %s" % (inDS, hte.headers) )
            except ClientException as cle:
                print ( "Submission for input dataset %s failed: %s" % (inDS, cle) )
                
                # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print ( "-"*len(msg) )
            print ( msg )
            print ( "-"*len(msg) )
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print ( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers) )
            except ClientException as cle:
                print ( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle) )


if __name__ == '__main__':
    main()
