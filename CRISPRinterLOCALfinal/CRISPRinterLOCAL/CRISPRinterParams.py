#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
CRISPRinterLOCAL - A CRISPR homologous enhancement-based host/phage interaction prediction tool, LOCAL
Author: Author: xintong HUANG, yishu MA, sijia MA
Beijing Normal University - Hong Kong Baptist University United International College, Zhuhai, China
Email: 2030026064@mail.uic.edu.hk; 2030026105@mail.uic.edu.hk; 2030026107@mail.uic.edu.hk
'''
import argparse
import sys
import os
from Bio import SeqIO

def get_arg_parser():

    ap = argparse.ArgumentParser(description='PAMPHLET - PAM Prediction HomoLogous Enhancement Toolkit')

    ap.add_argument("-i","--input",help="Input bacteria genome, required",required=True)
    ap.add_argument("-o","--outdir",help="Output directory, required.",required=True)
    ap.add_argument("-b","--blastmode",help="Spacer blastn mode. Common mode means use default blastn parameters and strict mode means use specific parameters, which could get more putative protospacers but also could cause false positives. Default is [common].",choices=['relax','common'],default='common')
    ap.add_argument("-u","--unique",help="Unique mode. If set, only revise unique spacer. Default is [False].",action="store_true",default=False)
    ap.add_argument("-bDB",'--bacteriaDB',help="Bacteria database for CRISPRinterLOCAL, you can use BACinterdb to build bacteria database.",required=True)
    ap.add_argument("-pDB",'--phageDB',help="Phage database for CRISPRinterLOCAL, you can use PHAinterdb to build phage database. If not provide, use online blastn nr[phage] database")
    ap.add_argument("--pcovs",help="Minimum percent coverage of spacer sequence, default is [0.9].",type=float,default=0.9)
    ap.add_argument("--pident",help="Minimum percent identity of spacer sequence, default is [0.9].",type=float,default=0.9)
    ap.add_argument("--rident",help="Minimum percent identity of repeat sequence, default is [0.9].",type=float,default=0.8)
    ap.add_argument("--MaxProteinNum",help="Maximum number of protein homologs, default is [20]. If the size is too large, NCBI will forbidden your request.",type=int,default=20)
    ap.add_argument("--cctyperDB",help="Input your CRISPRCasTyper DB location")

    return ap

def check_environment():
    ### 检查 minCED 加入 PATH ###
    if os.popen('which minced').read().strip() == '':
        print('minCED not found. Please add minced to PATH.')
        sys.exit(1)
    
    ### 检查 seqkit 加入 PATH ###
    if os.popen('which seqkit').read().strip() == '':
        print('seqkit not found. Please add seqkit to PATH.')
        sys.exit(1)
    
    if os.popen('which cctyper').read().strip() == "" and os.popen('which CRISPRCasTyper').read().strip() == "":
        print("CRISPRCasTyper not found. Please add CRISPRCasTyper to PATH.")
        sys.exit(1)

    if os.popen('which prodigal').read().strip() == '':
        print("Prodigal not found. Please add weblogo to PATH.")
        sys.exit(1)

def check_outdir(outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        print('Output directory %s already exists.' % outdir)
        sys.exit(1)

### 检查互联网连接是否正常
def check_internet_connection():
    ### check internet connection
    ret = os.system('ping baidu.com -n 1')
    if ret == 0:
        return True
    else:
        print('No internet connection. Please check your network.')
        sys.exit(1)

def check_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def check_null_files(infile):
    if os.path.getsize(infile) == 0:
        return True
    else:
        return False

### 检查给定目录下以 .minced 结尾的文件大小是否为0
def check_minced_files(indir):
    acceptable = 0
    for minced_result in os.listdir(indir):
        if minced_result.endswith('.minced'):
            filesize = os.path.getsize(os.path.join(indir,minced_result))
            if filesize == 0:
                os.system("rm "+os.path.join(indir,minced_result)+" "+\
                          os.path.join(indir,minced_result.replace('.minced','.minced.gff'))+" "+\
                          os.path.join(indir,minced_result.replace('.minced','.minced.fa')))
            else:
                acceptable += 1
    if acceptable == 0:
        return False
    else:
        return True

def main():

    check_environment()
    check_internet_connection()


if __name__ == '__main__':
    main()
