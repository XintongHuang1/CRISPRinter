#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
CRISPRinterLOCAL - A CRISPR homologous enhancement-based host/phage interaction prediction tool, LOCAL
Author: Author: xintong HUANG, yishu MA, sijia MA
Beijing Normal University - Hong Kong Baptist University United International College, Zhuhai, China
Email: 2030026064@mail.uic.edu.hk; 2030026105@mail.uic.edu.hk; 2030026107@mail.uic.edu.hk
'''

import sys
running_python3 = False
if sys.version_info > (3, 0):
    running_python3 = True

import os
import time
import warnings
warnings.filterwarnings("ignore")


print('''
  ____ ____  ___ ____  ____  ____  _       _
 / ___| _  \|_ _/ ___||  _ \|  _ \(_)_ __ | |_ ___ _ __
| |   | |_) || |\___ \| |_) | |_) | | '_ \| __/ _ \ '__|
| |___| _ <  | | ___) |  __/| _ < | | | | | ||  __/ |
 \____|_| \_\___|____/|_|   |_| \_\_|_| |_|\__\___|_| 

''')

print("CRISPRinter")
print("Author: xintong HUANG, yishu MA, sijia MA")
print("Beijing Normal University - Hong Kong Baptist University United International College, Zhuhai, China")
print("===============================================================")
print("Start loading resources...")
resourcesTime = time.time()
runningTime = time.time()

from CRISPRinterLOCAL import CRISPRinterParams
from CRISPRinterLOCAL import CRISPRinterReviser
from CRISPRinterLOCAL import CRISPRinterResources
from CRISPRinterLOCAL import OnlineResources

print("Loading resources finished. Time used: %s seconds" % (time.time() - resourcesTime))

def main():
    print("===============================================================")
    print("Initializing parameters...")

    arg_parser = CRISPRinterParams.get_arg_parser()
    opts = arg_parser.parse_args()

    InputGenome = opts.input
    OutDir = opts.outdir
    BlastMode = opts.blastmode
    BacteriaDatabase = opts.bacteriaDB
    PhageDatabase = opts.phageDB
    ProteinCOV = opts.pcovs
    ProteinIDENT = opts.pident
    UniqueMode = opts.unique
    RepeatIDENT = opts.rident
    MaxProteinNum = opts.MaxProteinNum
    CRISPRCasTyperDB = opts.cctyperDB

    ### 文件检查

    CRISPRinterParams.check_environment()
    CRISPRinterParams.check_outdir(OutDir)

    ### 模式检查和间隔修改

    print("Initializing parameters finished. Time used: %s seconds" % (time.time() - runningTime))
    print("===============================================================")

    ### 预测输入细菌基因组的所有潜在CRISPR-Cas系统
    CRISPRCasTyperDir = os.path.join(OutDir,"CRISPRCasTyperResult/")
    CRISPRinterResources.run_CRISPRCasTyper(InputGenome,CRISPRCasTyperDB,CRISPRCasTyperDir)

    CRISPRinterResources.check_CRISPRCas_status(CRISPRCasTyperDir)
    CRISPRinfoLabeledDict = CRISPRinterResources.get_CRISPRCas_info(CRISPRCasTyperDir)

    ### [crisprid]:[ContigID,SpacerFile,RepeatSeq,ProteinSeq,Label,CasFlag]
    ### Label: WithProt & NoProt
    ### CasFlag: Cas9/Cas3/Cas12/Cas13/None ...

    ReviserDir = os.path.join(OutDir,"Reviser/")
    CRISPRinterParams.check_directory(ReviserDir)

    for CRISPRunit, CRISPRunitInfo in CRISPRinfoLabeledDict.items():
        print(CRISPRunitInfo)

    for CRISPRunit, CRISPRunitInfo in CRISPRinfoLabeledDict.items():
        CRISPRLabel = CRISPRunitInfo[4]
        CasFlag = CRISPRunitInfo[5]
        # print("This is CasFlag: "+CasFlag)
        targetDir = os.path.join(ReviserDir,CRISPRunit)
        if CRISPRLabel == "WithProt":
            ### Start CRISPR homologous enhancement

            CRISPRinterParams.check_directory(targetDir)
            finalspacer = CRISPRinterReviser.main(CRISPRunitInfo[1],CRISPRunitInfo[2],CRISPRunitInfo[3],targetDir,UniqueMode,ProteinCOV,ProteinIDENT,RepeatIDENT,MaxProteinNum,BacteriaDatabase,CasFlag)
            if finalspacer == False:
                finalspacer = CRISPRunitInfo[1]
        else:
            ### No enhancement, use spacer file to predict phage
            finalspacer = CRISPRunitInfo[1]
        ### 修改finalspacer文件

        spaceridLenDict = CRISPRinterResources.modify_spacer_file(finalspacer)
        ### run spacer blast
        ### 设置SPACER BLAST 目录
        SpacerBlastDir = os.path.join(targetDir,"SpacerBlast/")
        CRISPRinterParams.check_directory(SpacerBlastDir)
        SpacerInputDir = os.path.join(SpacerBlastDir,"Input/")
        CRISPRinterParams.check_directory(SpacerInputDir)
        SpacerBlastOutputDir = os.path.join(SpacerBlastDir,"RawOutput/")
        CRISPRinterParams.check_directory(SpacerBlastOutputDir)

        ### 检查间隔长度,如果长度> 20,切至20
        SpacerNumbers = os.popen("grep -c '>' %s" % finalspacer).read().strip()
        if int(SpacerNumbers) > 20:
            CRISPRinterResources.split_spacer_file(finalspacer,SpacerInputDir,int(SpacerNumbers))
        else:
            os.system("cp %s %s" % (finalspacer,SpacerInputDir))
        if(PhageDatabase):
            for InputSpacerFile in os.listdir(SpacerInputDir):
                InputSpacerFile = os.path.join(SpacerInputDir, InputSpacerFile)
                OutputSpacerFile = os.path.join(SpacerBlastOutputDir, os.path.basename(InputSpacerFile))
                # print(OutputSpacerFile)
                CRISPRinterResources.local_spacer_blast(InputSpacerFile, OutputSpacerFile, PhageDatabase)
        ### 对重复数据库BLAST SPACER
        else:
            for InputSpacerFile in os.listdir(SpacerInputDir):
                InputSpacerFile = os.path.join(SpacerInputDir,InputSpacerFile)
                OutputSpacerFile = os.path.join(SpacerBlastOutputDir,os.path.basename(InputSpacerFile))
                OnlineResources.run_spacer_blast(InputSpacerFile,OutputSpacerFile,BlastMode)

        ### 合并BLAST输出， 如果有多个BLAST输出
        MergedSpacerBlastOutput = os.path.join(SpacerBlastDir,"MergedSpacerBlastOutput.txt")
        if len(os.listdir(SpacerBlastOutputDir)) > 1:
            CRISPRinterResources.merge_blast_output(SpacerBlastOutputDir,MergedSpacerBlastOutput)
        else:
            os.system("cp %s %s" % (os.path.join(SpacerBlastOutputDir,os.listdir(SpacerBlastOutputDir)[0]),MergedSpacerBlastOutput))

        ### 获得SIGNIFICANT SPACER BLAST输出与NON-SIGNIFICANT HIT SPACER的ID
        SignificantSpacerBlastOutput = os.path.join(SpacerBlastDir,"SignificantSpacerBlastOutputRaw.txt")
        InSignificantSpacerID = CRISPRinterResources.get_significant_spacer_blast_output(MergedSpacerBlastOutput,SignificantSpacerBlastOutput,spaceridLenDict)
        
    ### 设置最终SpacerHit文件
    FinalSpacerHitFile = os.path.join(OutDir, "FinalSpacerHit.txt")
    CRISPRinterResources.merge_final_spacer_hit(ReviserDir, FinalSpacerHitFile)
    MergeFile = os.path.join(OutDir, "merge.txt")
    HostCountFile = os.path.join(OutDir, "host_count.txt")
    Host2CountFile = os.path.join(OutDir, "host2_count.txt")
    OnlineResources.final_prediction(FinalSpacerHitFile, MergeFile, HostCountFile, Host2CountFile)
        
if __name__ == "__main__":
    main()
