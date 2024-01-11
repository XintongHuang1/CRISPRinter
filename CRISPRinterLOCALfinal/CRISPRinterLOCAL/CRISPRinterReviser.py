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
from Bio import SeqIO

from CRISPRinterLOCAL import CRISPRinterParams
from CRISPRinterLOCAL import OnlineResources
from CRISPRinterLOCAL import CRISPRinterResources

def main(spacerFile,repeatseq,proteinSeq,outDir,inunique,incov,inident,inrident,inmaxprotlen,blastdb,effectorflag):

    ### 读取蛋白质长度并检查修改状态
    proteinLength = len(proteinSeq)
    ProteinBlastDir = os.path.join(outDir,'ProteinBlast/')
    ProteinBlastResult = os.path.join(ProteinBlastDir,'ProteinBlastResult.txt')
    CRISPRinterParams.check_directory(ProteinBlastDir)

    ### 运行蛋白胚预测相关基因组
    ### 建立临时蛋白质fasta文件
    tempProteinFile = os.path.join(ProteinBlastDir,'tempProteinSeq.fasta')
    CRISPRinterResources.build_temp_protein_fasta(proteinSeq,tempProteinFile)
    CRISPRinterResources.run_local_blastp(tempProteinFile,ProteinBlastResult,blastdb,effectorflag)

    if CRISPRinterParams.check_null_files(ProteinBlastResult):
        print("No related protein found. Spacer revise failed.")
        return False
        
    ### 获得significant hits，返回假定的来源的tax name
    SignificantHitResult = os.path.join(ProteinBlastDir,"SignificantHits.txt")
    LevelBox = CRISPRinterResources.get_local_significant_blastp_hits(ProteinBlastResult,SignificantHitResult,incov,inident)
    print(LevelBox)
    if CRISPRinterParams.check_null_files(SignificantHitResult):
        print("No significant hits found. Spacer revise failed. Use original spacer file")
        return False

    ### 根据MinCED结果选择正确的DRs
    CasFullInfoPickle = os.path.join(blastdb,"CasInfoPICKLE/",effectorflag+".pickle")
    spacerPickle = os.path.join(blastdb,"SpacerCollection.pkl")
    SelectedSpacerFilesDict = CRISPRinterResources.select_correct_dr(spacerPickle,repeatseq,inrident,CasFullInfoPickle,LevelBox)

    ### 检查选定的Spacer文件字典
    if len(SelectedSpacerFilesDict) == 0:
        print("No correct DR found. Spacer revise failed. Use original spacer file")
        return False
    
    ### 修改修改Spacer文件
    RevisedSpacerFile = os.path.join(outDir,'RevisedSpacer.txt')
    #CRISPRinterResources.revise_spacer_file(spacerFile,MinCEDDir,SelectedSpacerFilesDict,RevisedSpacerFile,inunique)
    CRISPRinterResources.revise_spacer_file(spacerFile,SelectedSpacerFilesDict,inunique,RevisedSpacerFile)

    ### Return revised spacer file
    return RevisedSpacerFile

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
