#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
CRISPRinterLOCAL - A CRISPR homologous enhancement-based host/phage interaction prediction tool, LOCAL
Author: Author: xintong HUANG, yishu MA, sijia MA
Beijing Normal University - Hong Kong Baptist University United International College, Zhuhai, China
Email: 2030026064@mail.uic.edu.hk; 2030026105@mail.uic.edu.hk; 2030026107@mail.uic.edu.hk
'''
import re

from pyfaidx import Fasta
import math
import sys
import pickle
import os
from Bio import pairwise2
from Bio import SeqIO
import subprocess

def reverse_complete(seq):
    trantab = str.maketrans('ATCGatcg','TAGCtagc')
    return seq.translate(trantab)[::-1]

def run_CRISPRCasTyper(ingen,indb,outdir):
     ccCommand = "cctyper --keep_tmp --db "+indb+" "+ingen+" "+outdir
     os.system(ccCommand) 
     return True

def merge_final_spacer_hit(rawdir,outfile):
    fa = open(outfile,'w')
    for dirs in os.listdir(rawdir):
        if os.path.exists(os.path.join(rawdir,dirs,"SpacerBlast/MergedSpacerBlastOutput.txt")):
            with open(os.path.join(rawdir,dirs,"SpacerBlast/MergedSpacerBlastOutput.txt"),'r') as fb:
                for fbline in fb.readlines():
                    fa.write(fbline)
    fa.close()
    return True

def modify_spacer_file(infile):
    lenidDict = {}
    tempfile = infile+".tmp"
    os.system("cp "+infile+" "+tempfile)
    with open(infile,'w') as fa:
        for records in SeqIO.parse(tempfile,'fasta'):
            #print(records)

            if "#LEVEL-" in str(records.id):
                ### This is a revised spacer, rename the seqid
                crisprIndex = str(records.id).split("_")[-2]
                spacerIndex = str(records.id).split("_")[-1]
                levelRank = "LEVEL-"+str(records.id).split("#LEVEL-")[1].split("_")[0]
                rawContig = str(records.id).split("#LEVEL-")[0]
                spacerID = rawContig+"_"+crisprIndex+"_"+spacerIndex+levelRank
                spacerLen = len(str(records.seq))
                if not spacerID in lenidDict.keys():
                    lenidDict[spacerID] = spacerLen
                fa.write(">"+spacerID+"\n"+str(records.seq)+"\n")
            else:
                spacerID = str(records.id)+"_LEVEL-1"
                spacerLen = len(str(records.seq))
                if not spacerID in lenidDict.keys():
                    lenidDict[spacerID] = spacerLen
                pattern = re.compile(r'^[ATCGatcg]+$')
                # print(type(str(records.seq)), records.seq)
                if bool(pattern.match(str(records.seq))):
                     fa.write(">"+spacerID+"\n"+str(records.seq)+"\n")
    os.system("rm "+tempfile)
    return lenidDict

def get_effector_protein_index(crisprtype,effectordict,genelist,geneindex,contigID):
    searchList = effectordict[crisprtype]
    simplified_gene_list = []
    allCasIDList = []
    for units in genelist:
        simplified_gene_list.append(units.split("_")[0])
    for geneid in searchList:
        if geneid in simplified_gene_list:
            rawindex = simplified_gene_list.index(geneid)
            casindex = geneindex[rawindex]
            casproteinID = contigID+"_"+str(casindex)
            allCasIDList.append([casproteinID,geneid])
    return allCasIDList

def build_temp_protein_fasta(inseqfile,outfile):
    with open(outfile,'w') as fa:
        fa.write(">TempProteinSeq\n"+inseqfile+"\n")
    return True

def run_local_blastp(inseqfile,outfile,indatabase,inflag):
    ### proteinSeq,ProteinBlastResult,blastdb,effectorflag
    targetDB = os.path.join(indatabase,"CasSequenceAccept/",inflag)
    cmd = "blastp -query "+inseqfile+" -db "+targetDB+" -out "+outfile+" -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'"
    os.system(cmd)
    return True

def get_max_length_effector(inlist,refindex):
    maxlen = 0
    maxprot = ""
    maxflag = ""
    for proteininfo in inlist:
        proteinID = proteininfo[0]
        proteinFlag = proteininfo[1]
        proteinseqs = str(refindex[proteinID][::])
        if len(proteinseqs) > maxlen:
            maxprot = proteinseqs
            maxflag = proteinFlag
    return maxprot, maxflag

def make_protein_file(outfile,crisprid,inseq):
    with open(outfile,'w') as ff:
        ff.write(">"+crisprid+"_effectorProt\n"+inseq+"\n")
    return True

def get_CRISPRCas_info(CRISPRCasTyperDir):
    ### This step, get each CRISPR-Cas system info and store in a dictionary
    ### SpacerFile,RepeatSeq,ProteinFile,OutDir,UniqueMode,ProteinBlastResult,ProteinLength,RefMode,ProteinCOV,ProteinIDENT,RepeatIDENT,MaxProteinNum
    ### Dict format: [crisprid]:[ContigID,SpacerFile,RepeatSeq,ProteinSeq,Label]
    ### This step, build by reading the crisprs_all.tab file
    TypeEffectorDict = {'I':['Cas7','Csc2','Cas6','Cas5','Csc1','Cas8','Cse1','Cas11','Cse2'],
                        'III':['Cas7','Csm3','Cas6','Cas5','Csm5'],
                        'IV':['Cas8','Csf1','Cas7','Csf2','Cas6','Cas5','Csf3','Cas10'],
                        'II':['Cas9'],
                        'V':['Cas12'],
                        'VI':['Cas13']}
    CRISPRCasDict = {}
    CasOperonDict = {}
    CRISPROperonDict = {}
    ### First, build cas operon dict, crisprID - cas
    if os.path.exists(os.path.join(CRISPRCasTyperDir,"cas_operons.tab")):
        ### Build protein faa index
        if not os.path.exists(os.path.join(CRISPRCasTyperDir,"proteins.faa.fai")):
            refseqs = Fasta(os.path.join(CRISPRCasTyperDir,"proteins.faa"))
        with open(os.path.join(CRISPRCasTyperDir,"cas_operons.tab"),'r') as fo:
            for foline in fo.readlines():
                if not foline.startswith("Contig\tOperon\tStart"):
                    contigID, operonID = foline.split("\t")[:2]
                    PredictionType = foline.split("\t")[4].split("-")[0]
                    CasGeneList = foline.split("\t")[9].replace("[","").replace("]","").replace("\'","").replace(" ","")
                    if "," in CasGeneList:
                        CasGeneList = CasGeneList.split(",")
                    else:
                        CasGeneList = [CasGeneList]
                    CasGeneRank = foline.split("\t")[10].replace("[","").replace("]","").replace("\'","").replace(" ","")
                    if "," in CasGeneRank:
                        CasGeneRank = CasGeneRank.split(",")
                    else:
                        CasGeneRank = [CasGeneRank]
                    AllEffectorHit = get_effector_protein_index(PredictionType,TypeEffectorDict,CasGeneList,CasGeneRank,contigID)
                    MaxLengthEffector, MaxLengthEffectorFlag = get_max_length_effector(AllEffectorHit,refseqs)
                    if not operonID in CasOperonDict:
                        CasOperonDict[operonID] = [MaxLengthEffector, MaxLengthEffectorFlag]
            fo.close()

    ### This step, merge cas_operons and CRISPR array, based on crisprs_all.tab
    if os.path.exists(os.path.join(CRISPRCasTyperDir,"CRISPR_Cas.tab")):
        with open(os.path.join(CRISPRCasTyperDir,"CRISPR_Cas.tab"),'r') as fc:
            for fcline in fc.readlines():
                if not fcline.startswith("Contig\tOperon\tOperon_Pos"):
                    operonID = fcline.split("\t")[1]
                    CRISPRsList = fcline.split("\t")[4].replace("[","").replace("]","").replace("\'","").replace(" ","")
                    if "," in CRISPRsList:
                        CRISPRsList = CRISPRsList.split(",")
                    else:
                        CRISPRsList = [CRISPRsList]
                    for CRISPRunit in CRISPRsList:
                        CRISPROperonDict[CRISPRunit] = CasOperonDict[operonID]
            fc.close()

    ### CasOperonDict: key: operonID, value: [MaxLengthEffector, MaxLengthEffectorFlag]
    ### CRISPROperonDict: key: crisprID, value: [MaxLengthEffector, MaxLengthEffectorFlag]

    if os.path.exists(os.path.join(CRISPRCasTyperDir, "crisprs_all.tab")):
        CRISPRallFile = os.path.join(CRISPRCasTyperDir, "crisprs_all.tab")
        with open(CRISPRallFile,'r') as fa:
            for faline in fa.readlines():
                if not faline.startswith("Contig\tCRISPR\tStart"):
                    ### Main info, CRISPRCasDict: [crisprid]:[ContigID,SpacerFile,RepeatSeq,ProteinSeq,Label]
                    contigID, crisprID, startPos, EndPos, ConsensusRepeat = faline.split("\t")[:5]
                    spacerFile = os.path.join(CRISPRCasTyperDir,"spacers/"+crisprID+".fa")
                    if crisprID in CRISPROperonDict:
                        ProteinSeq = CRISPROperonDict[crisprID][0]
                        ProteinFlag = CRISPROperonDict[crisprID][1]
                        proteinOutfile = os.path.join(CRISPRCasTyperDir,crisprID+"_effector.faa")
                        make_protein_file(proteinOutfile,crisprID,ProteinSeq)
                        CRISPRCasDict[crisprID] = [contigID, spacerFile, ConsensusRepeat, ProteinSeq, "WithProt", ProteinFlag]
                    else:
                        CRISPRCasDict[crisprID] = [contigID, spacerFile, ConsensusRepeat, None, "NoProt", None]
            fa.close()



    return CRISPRCasDict

def check_CRISPRCas_status(outdir):
    if not os.path.exists(os.path.join(outdir,"spacers/")):
        ### No spacers found
        print("No CRISPR-Cas records on the input genomes, exit.")
        sys.exit(1)
    return True

def near_crispr_array(indict,inid,instart,inend):
    ### indict format:
    ### {[contigid]:[[CRISPR_1_start,CRISPR_1_end],[CRISPR_2_start,CRISPR_2_end]...],...}
    if inid in indict:
        for units in indict[inid]:
            startpos = units[0]
            endpos = units[1]
            ### Inside CRISPR array
            if instart <= endpos and instart >= startpos:
                return True
            if inend <= endpos and instart >= endpos:
                return True
            if abs(inend-endpos) <= 50 or abs(inend-startpos) <= 50:
                return True
            if abs(instart-endpos) <= 50 or abs(instart-startpos) <= 50:
                return True
    return False

def get_local_significant_blastp_hits(blastres,significantres,incov,inident):
    ### 获得significant hits
    ### Significant hits 被定义为evalue <= 1e-5 & qcov >= 0.9 & pident >= 90.0的hits
    ### 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
    ### 只检查前5行最热门的候选
    LevelBox = {}

    with open(blastres,'r') as fa, open(significantres,'w') as fb:
        for records in fa.readlines():
            if records.strip() != "" and not records.startswith("#"):
                Evalue = float(records.split("\t")[10])
                QueryCoverage = float(records.split("\t")[-1])
                QueryIdentity = float(records.split("\t")[2])
                SubjectID = records.split("\t")[1]
                if Evalue <= 1e-30 and QueryCoverage >= 0.98 and QueryIdentity >= 98.0: # if True:
                    LevelBox[SubjectID] = "LEVEL-2"
                if Evalue <= 1e-5 and QueryCoverage >= incov and QueryIdentity >= inident*100:   # if True:
                    fb.write(records)
                    if not SubjectID in LevelBox:
                        LevelBox[SubjectID] = "LEVEL-3"

    fa.close()
    fb.close()
    return LevelBox

def get_significant_blastp_hits(blastres,significantres,proteinlen,incov,inident):
    ### 获得significant hits
    ### Significant hits 被定义为evalue <= 1e-5 & qcov >= 0.9 & pident >= 90.0的hits

    linecount = 0
    ### 只检查前5行最热门的候选
    TopHit = []

    with open(blastres,'r') as fa, open(significantres,'w') as fb:
        for records in fa.readlines():
            if records.strip() != "" and not records.startswith("#"):
                linecount += 1
                Evalue = float(records.split("\t")[10])
                HitLength = int(records.split("\t")[3])
                QueryCoverage = HitLength/proteinlen
                QueryIdentity = float(records.split("\t")[2])
                SubjectID = records.split("\t")[1]

                if Evalue <= 1e-5 and QueryCoverage >= incov and QueryIdentity >= inident*100:   # if True:
                    fb.write(records)
                if Evalue <= 1e-30 and QueryCoverage >= 0.98 and QueryIdentity >= 98.0 and linecount <= 5: # if True:
                    TopHit.append(SubjectID)

    fa.close()
    fb.close()
    return TopHit

def check_tax_ref_level(intaxname,intopbox):
    if intopbox == []:
        return "LEVEL-3"
    else:
        flag = "LEVEL-3"
        for units in intopbox:
            if units == intaxname:
                flag = "LEVEL-1"
            else:
                if units.split()[0] == intaxname.split()[0]:
                    if flag == "LEVEL-3":
                        flag = "LEVEL-2"
                    else:
                        continue
                else:
                    continue
        return flag



def get_seqid2seqtaxid(infodict,tophitbox):
    ### 获得 seqid2seqtaxid
    ### 如果已有tophitID, 则默认选择为no-reference mode; 否则，the tempdict select 将在引用模式下，引用为tophit
    ### If have 100%-100% hit, the top hit is the reference, the phage-protospacer-rank is: species(top)-genus(2nd)-homologous(3rd) --> REFERENCE MODE
    ### If havenot 100%-100% hit, the top hit is non-reference, the phage-protospacer-rank is: original_spacer(top)-homologous(2nd) --> NON-REFERENCE MODE
    
    tempdict = {}
    finaldict = {}

    ### infodict[sbjctID] = [refid,genometaxid,int(orfstart),int(orfend)]   <--  infodict format
    ### tophitbox format: [sbjctID1, sbjctID2, ... , sbjctIDn]
    ### 建 tophitbox_tax
    tophitbox_tax = []
    for wps in tophitbox:
        if wps in infodict:
            tophitbox_tax.append(infodict[wps][1])

    if len(tophitbox_tax) != 0:
        ### Now, use the reference mode, select the tophit
        for wpid, wpinfo in infodict.items():
            genomeid = wpinfo[0]
            taxname = wpinfo[1]
            taxRefLevel = check_tax_ref_level(taxname,tophitbox_tax)
            wpinfo.append(taxRefLevel)
            finaldict[wpid] = wpinfo
    else:
        ### Now, use all mode, select all
        for wpid, wpinfo in infodict.items():
            genomeid = wpinfo[0]
            taxname = wpinfo[1]
            wpinfo.append("LEVEL-3")
            finaldict[wpid] = wpinfo
    return finaldict

def extract_related_seqdump(infile,outdir,refdict):
    refseqs = Fasta(infile)

    ### 基于inhit构建命中位置字典
    duplicatebox = []
    for wpid, wpinfo in refdict.items():
        genomeid = wpinfo[0]
        start = int(wpinfo[2])
        end = int(wpinfo[3])
        levels = wpinfo[4]
        if start >= end:
            start, end = end, start
        if not genomeid in duplicatebox and genomeid in refseqs.keys():
            with open(os.path.join(outdir,genomeid + ".fa"),'a') as fa:
                contiglen = len(str(refseqs[genomeid][::]))
                if start <= 20000:
                    start = 1
                else:
                    start = start - 20000
                if end + 20000 >= contiglen:
                    end = contiglen
                else:
                    end = end + 20000
                fa.write(">"+genomeid+"\n"+str(refseqs[genomeid][start:end])+"\n")
            fa.close()

    return True

def run_minced(indir,outdir):
    for fastafile in os.listdir(indir):
        if fastafile.endswith(".fa"):
            CMD = "minced -spacers "+os.path.join(indir,fastafile)+" "+\
                os.path.join(outdir,fastafile.replace(".fa",".minced"))+" "+\
                os.path.join(outdir,fastafile.replace(".fa",".gff"))
            print("Now running minced for "+fastafile+", command: "+CMD)
            os.system(CMD)
    return True


def select_correct_dr(spacerpkl,tempseq,inrident,casinfopkl,levelbox):
    ### spacerPickle,repeatseq,inrident,CasFullInfoPickle,LevelBox
    ### levelbox: keys: protein_sbjctID/ values: levels
    ### casinfopkl: keys: protein_sbjctID/ values: [flag, assembly, contig, sbjctid, operonid, list(dr), list(crisprid)]
    ### spacerpkl: keys: spacerID/ values: list(spacerSeq)

    ### Read pickle file
    finfo = open(casinfopkl,'rb')
    fspacer = open(spacerpkl,'rb')
    CasInfoDict = pickle.load(finfo)
    SpacerCollection = pickle.load(fspacer)

    finaldict = {}
    for keys, values in levelbox.items():
        if keys in CasInfoDict:
            drseqs = CasInfoDict[keys][0][5]
            crisprID = CasInfoDict[keys][0][6]
            levels = values
            for drindex in range(0,len(drseqs)):
                rptunitseq = drseqs[drindex].upper().replace("U","T")
                topalign = pairwise2.align.globalms(rptunitseq,tempseq.replace("U","T"),10,0,-10,-0.5,one_alignment_only=True)
                formatalign = pairwise2.format_alignment(*topalign[0])
                similarity = formatalign.split("\n")[1].count("|")/len(topalign[0][0].lstrip("-").rstrip("-"))
                if similarity >= inrident:
                    crisprunitID = crisprID[drindex]+".fa"
                    if crisprunitID in SpacerCollection.keys():
                        finaldict[crisprunitID+"#"+levels] = SpacerCollection[crisprunitID]
                    else:
                        print("Error: "+crisprunitID+" not found in spacer collection.")
                else:
                    rcrptunitseq = reverse_complete(rptunitseq)
                    topalign = pairwise2.align.globalms(rcrptunitseq,tempseq.replace("U","T"),10,0,-10,-0.5,one_alignment_only=True)
                    formatalign = pairwise2.format_alignment(*topalign[0])
                    similarity = formatalign.split("\n")[1].count("|")/len(topalign[0][0].lstrip("-").rstrip("-"))
                    if similarity >= inrident:
                        crisprunitID = crisprID[drindex]+".fa"
                        if crisprunitID in SpacerCollection.keys():
                            ### 对于每个间隔符，存储反向补码序列
                            templist = []
                            for units in SpacerCollection[crisprunitID]:
                                templist.append(reverse_complete(units))
                            finaldict[crisprunitID+"#"+levels] = templist
                        else:
                            print("Error: "+crisprunitID+" not found in spacer collection.")
    finfo.close()
    fspacer.close()
    return finaldict

def revise_spacer_file(orispa,indict,inunique,outfile):
    ### CRISPRinterResources.revise_spacer_file(SelectedSpacerFilesDict,inunique,RevisedSpacerFile)
    ### indict: keys: crisprID, values: list(spacerseqs)
    fo = open(outfile,'w')
    uniquebox = []
    if inunique:
        for records in SeqIO.parse(orispa,'fasta'):
            if not str(records.seq) in uniquebox:
                uniquebox.append(str(records.seq))
                fo.write(">"+str(records.id)+"\n"+str(records.seq)+"\n")
        for crisprid, spacerseqlist in indict.items():
            counter = 1
            crisprid = "#".join(crisprid.split("#")[:-1])
            for unitSpacers in spacerseqlist:
                if not unitSpacers in uniquebox:
                    uniquebox.append(unitSpacers)
                    fo.write(">"+crisprid+"_"+str(counter)+"\n"+unitSpacers+"\n")
                    counter += 1
    else:
        for records in SeqIO.parse(orispa,'fasta'):
            uniquebox.append(str(records.seq))
            fo.write(">"+str(records.id)+"\n"+str(records.seq)+"\n")
        for crisprid, spacerseqlist in indict.items():
            counter = 1
            crisprid = "#".join(crisprid.split("#")[:-1])
            for unitSpacers in spacerseqlist:
                if not unitSpacers in uniquebox:
                    uniquebox.append(unitSpacers)
                    fo.write(">"+crisprid+"_"+str(counter)+"\n"+unitSpacers+"\n")
                    counter += 1
    fo.close()
    return True


def split_spacer_file(infile,outdir,seqnum):
    splitpartnum = math.ceil(float(seqnum/30))
    ### Use seqkit to split the file
    CMD = "seqkit split -p "+str(splitpartnum)+" "+infile+" --quiet"
    os.system(CMD)
    splitteddir = os.path.join(infile+".split")
    os.system("mv "+os.path.join(splitteddir,"*")+" "+outdir)
    os.system("rm -r "+splitteddir)
    return True

def rename_spacer_id(infile,inlabel,inorientation):
    finalSpacerBox = {}
    spacercounter = 1
    with open(infile+".tmp",'w') as fa, open(inlabel,'w') as fb:
        for spacers in SeqIO.parse(infile,'fasta'):
            if inorientation == "positive":
                fa.write(">SpacerID"+str(spacercounter)+"\n"+str(spacers.seq)+"\n")
                fb.write(">SpacerID"+str(spacercounter)+"\t"+str(spacers.id)+"\n")
            else:
                fa.write(">SpacerID"+str(spacercounter)+"\n"+reverse_complete(str(spacers.seq))+"\n")
                fb.write(">SpacerID"+str(spacercounter)+"\t"+str(spacers.id)+"\n")
            finalSpacerBox["SpacerID"+str(spacercounter)] = len(str(spacers.seq))
            spacercounter += 1
    fa.close()
    fb.close()
    return finalSpacerBox

def merge_blast_output(indir,outfile):
    with open(outfile,'w') as fa:
        for files in os.listdir(indir):
            with open(os.path.join(indir,files),'r') as fb:
                for fblines in fb.readlines():
                    ### Check blank line
                    if not fblines.strip() == "":
                        fa.write(fblines)
            fb.close()
    fa.close()
    return True

def get_significant_spacer_blast_output(inraw,insig,indict):
    significant_box = {}
    with open(inraw,'r') as fa, open(insig,'w') as fb:
        for faline in fa.readlines():
            if not faline.strip() == "":
                queryid = faline.strip().split("\t")[0]
                mismatches = int(faline.strip().split("\t")[4])
                gapnumbers = int(faline.strip().split("\t")[5])
                evalue = float(faline.strip().split("\t")[-2])
                querystartpos = int(faline.strip().split("\t")[6])
                queryendpos = int(faline.strip().split("\t")[7])
                qlen = indict[queryid]
                qcovs = (queryendpos-querystartpos+1)/qlen
                if evalue <= 1 and qcovs == 1 and mismatches <= 3 and gapnumbers == 0:
                    if not queryid in significant_box:
                        significant_box[queryid] = 0
                    significant_box[queryid] += 1
                    fb.write(faline)
                else:
                    ### This step, check the spacer alignment region, if only start/end 1nt not aligned, then it is still significant. But the mismatch numbers should be less than 2.
                    if querystartpos == 2 or queryendpos == qlen-1:
                        if mismatches <= 2 and gapnumbers == 0:
                            if not queryid in significant_box:
                                significant_box[queryid] = 0
                            significant_box[queryid] += 1
                            fb.write(faline)
    return True

def revise_spacer_length(infile,outfile,inid):
    refseqs = Fasta(infile)
    with open(outfile,'w') as fa:
        for seqid in inid:
            fa.write(">"+seqid+"\n"+str(refseqs[seqid][1:-1])+"\n")
    fa.close()
    return True

def merge_spacer_blast_output(fileA,fileB,outfile):
    ### fileA = significantHit, fileB = revisedSignificantHit
    duplicateBox = []
    with open(fileA,'r') as fa, open(fileB,'r') as fb, open(outfile,'w') as fc:
        for faline in fa.readlines():
            if not faline.strip() == "":
                if int(faline.strip().split("\t")[8]) > int(faline.strip().split("\t")[9]):
                    complexesC = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[9])+1) + "-" + str(int(faline.strip().split("\t")[8])-1)
                    complexesL = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[9])+1) + "-" + str(int(faline.strip().split("\t")[8]))
                    complexesR = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[9])) + "-" + str(int(faline.strip().split("\t")[8])-1)
                else:
                    complexesC = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[8])+1) + "-" + str(int(faline.strip().split("\t")[9])-1)
                    complexesL = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[8])+1) + "-" + str(int(faline.strip().split("\t")[9]))
                    complexesR = faline.strip().split("\t")[0] + "-" + faline.strip().split("\t")[1] + "-" +str(int(faline.strip().split("\t")[8])) + "-" + str(int(faline.strip().split("\t")[9])-1)
                duplicateBox.append(complexesC)
                duplicateBox.append(complexesL)
                duplicateBox.append(complexesR)
                fc.write(faline)
        for fbline in fb.readlines():
            if not fbline.strip() == "":
                if int(fbline.strip().split("\t")[8]) > int(fbline.strip().split("\t")[9]):
                    revisedcomplexes = fbline.strip().split("\t")[0] + "-" + fbline.strip().split("\t")[1] + "-" +str(fbline.strip().split("\t")[8])
                else:
                    revisedcomplexes = fbline.strip().split("\t")[0] + "-" + fbline.strip().split("\t")[1] + "-" +str(fbline.strip().split("\t")[9])
                revisedcomplexes = fbline.strip().split("\t")[0] + "-" + fbline.strip().split("\t")[1] + "-" +str(fbline.strip().split("\t")[8]) + "-" + str(fbline.strip().split("\t")[9])
                if not revisedcomplexes in duplicateBox:
                    fc.write(fbline)
    fa.close()
    fb.close()
    fc.close()
    return True

import subprocess

def local_spacer_blast(infile, outfile, blastdb):
# 构造blastn命令
    blast_cmd = [
    'blastn',
    '-query', infile,
    '-out', outfile + '.temp',
    '-outfmt', '6', # Tabular format
    '-db', blastdb
    ]

    # print(blast_cmd)
# 执行blastn命令
    subprocess.run(blast_cmd)
    return True




def main():
    print("Resources test.")

if __name__ == '__main__':
    main()
