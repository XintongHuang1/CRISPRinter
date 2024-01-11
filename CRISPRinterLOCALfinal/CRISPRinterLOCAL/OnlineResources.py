#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
CRISPRinterLOCAL - A CRISPR homologous enhancement-based host/phage interaction prediction tool, LOCAL
Author: Author: xintong HUANG, yishu MA, sijia MA
Beijing Normal University - Hong Kong Baptist University United International College, Zhuhai, China
Email: 2030026064@mail.uic.edu.hk; 2030026105@mail.uic.edu.hk; 2030026107@mail.uic.edu.hk
'''

import sys
import os
import requests
import time
import pandas as pd
import func_timeout
from func_timeout import func_set_timeout
from urllib.parse import quote
from pyfaidx import Fasta


def run_online_blastp(infile, outfile):
    ### 该在线blastp功能基于NCBI BLAST，用于预测同源蛋白相关的基因组序列。
    ### 搜索细菌数据库, taxid = 2

    ErrorMessage = "CPU usage limit was exceeded"

    UpLoadQuery = quote(infile)
    arguments = "CMD=Put&PROGRAM=blastp&DATABASE=nr&ENTREZ_QUERY=txid2[ORGN]&EXPECT=0.00001&QUERY=" + UpLoadQuery

    r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?" + arguments)
    RID = r.text.split("RID = ")[1].split()[0]
    print("PROTEIN blastp PROGRAM RID: " + RID)

    rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + RID)
    Status = rStatus.text.split("Status=")[1].split()[0]

    while Status != "READY":
        time.sleep(30)
        rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + RID)
        Status = rStatus.text.split("Status=")[1].split()[0]
        print("Refresh in 30s, status: " + Status)

    rResult = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Tabular&RID=" + RID)

    ### 检查rResult错误信息，如果错误信息为“CPU usage limit was exceeded”，则重新运行blastp程序。
    if ErrorMessage in rResult.text:
        print("NCBI blastp CPU usage limit was exceeded, re-run blastp program.")
        run_online_blastp(infile, outfile)

    else:
        with open(outfile + ".temp", 'w') as fot:
            fot.write(rResult.text)
        fot.close()
        flag = 0

        with open(outfile + ".temp", 'r') as fotr, open(outfile, 'w') as foc:
            for fotrecord in fotr.readlines():
                if flag == 1:
                    if fotrecord.strip() != "</PRE>":
                        foc.write(fotrecord)
                if "hits found" in fotrecord:
                    flag = 1
                if fotrecord.strip() == "</PRE>":
                    flag = 0

        os.system("rm " + outfile + ".temp")
        return True

def sbjctID_count(insig):
    df = pd.DataFrame(columns=['sbjctID'])
    with open(insig, 'r') as fa:
        for faline in fa.readlines():
            sbjctID = faline.split("\t")[1]
            new_row = {'subjectid': sbjctID}
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

    value_counts = df["subjectid"].value_counts()  # groupby and count
    df = value_counts.to_frame(name="times")
    df = df[~df.index.str.contains('subjectid')]
    df.reset_index(inplace=True)
    df.rename(columns={df.columns[0]: 'sbjctID'}, inplace=True)
    return df


def sbjctID_2_taxid_host(insig):
    df = pd.DataFrame(columns=['sbjctID', 'taxid', 'host'])
    genomepage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    with open(insig, 'r') as fa:
        for faline in fa.readlines():

            sbjctID = faline.split("\t")[1]
            protein2fasta = genomepage + sbjctID
            # print(protein2fasta)
            retry_counter = 0
            while retry_counter < 3:
                try:
                    proteinfasta = get_protein_sequences(protein2fasta)
                    break
                except:
                    retry_counter += 1
                    print("Retry in: " + sbjctID + ", retry time: " + str(retry_counter))
            taxid = ""
            try:
                start_index = proteinfasta.find('taxname "') + len('taxname "')
                end_index = proteinfasta.find('"', start_index)
                if start_index != -1 and end_index != -1:
                    taxid = proteinfasta[start_index:end_index]
            except:
                ### Into the next loop
                print("No taxname in: " + sbjctID)
                continue

            host = ""
            try:
                start_pattern = 'host,\n              subname "'
                end_pattern = '"'

                start_index = proteinfasta.find(start_pattern)
                end_index = proteinfasta.find(end_pattern, start_index + len('host,\n              subname "'))

                if start_index != -1 and end_index != -1:
                    host = proteinfasta[start_index + len(start_pattern):end_index].strip()
            except:
                continue
            try:
                start_pattern = 'host,\n          name "'
                end_pattern = '"'

                start_index = proteinfasta.find(start_pattern)
                end_index = proteinfasta.find(end_pattern, start_index + len('host,\n          name "'))

                if start_index != -1 and end_index != -1:
                    host = proteinfasta[start_index + len(start_pattern):end_index].strip()
            except:
                continue
            try:
                start_pattern = 'host,\n              name "'
                end_pattern = '"'

                start_index = proteinfasta.find(start_pattern)
                end_index = proteinfasta.find(end_pattern, start_index + len('host,\n              name "'))

                if start_index != -1 and end_index != -1:
                    host = proteinfasta[start_index + len(start_pattern):end_index].strip()
            except:
                ### Into the next loop
                print("No subname in: " + sbjctID)
                continue


            data = [sbjctID, taxid, host]
            #合并表
            new_row = {'sbjctID': data[0], 'taxid': data[1], 'host': data[2]}
            df = pd.concat([df, pd.DataFrame(new_row, index=[0])], ignore_index=True)
            df['host'] = df['host'].str.replace('"', '')
            df['host'] = df['host'].str.replace('"', '')

            to_be_deleted = df.loc[df["host"].str.startswith("y ::= set {")]
            to_be_deleted_taxids = to_be_deleted['taxid'].tolist()
            #print("Final Spacer to be Deleted with Taxid:")
            for taxid in to_be_deleted_taxids:
                print(taxid)
            to_be_deleted = df[~df['taxid'].str.contains(' phage')]
            to_be_deleted_taxids = to_be_deleted['taxid'].tolist()
            for taxid in to_be_deleted_taxids:
                print(taxid)
            #delete
            df = df.loc[~df["host"].str.startswith("y ::= set {")]
            df = df[df['taxid'].str.contains(' phage')]
    return df


def final_prediction(insig, MergeFile, HostCountFile, Host2CountFile):
    # 空FinalSpacerHitFile的情况
    with open(insig, 'r') as file:
        content = file.read()
    if content.strip() == '':
        print("Empty FinalSpacerHitFile. Coudln't make final prediction")
        return

    a = sbjctID_count(insig)

    with open(insig, "r") as file:
        lines = file.readlines()
    # 去重
    lines_unique = []
    seen = set()
    for line in lines:
        columns = line.split()
        if columns[1] not in seen:
            seen.add(columns[1])
            lines_unique.append(line)
    # 重写insig
    with open(insig, "w") as file:
        file.writelines(lines_unique)

    b = sbjctID_2_taxid_host(insig)
    df = pd.merge(a, b)
    if len(df) == 0:
        print("No Reliable Phage Result")
        return
    merge_max_host = df.iloc[0].loc['host'].split(' ')[:2]

    df['times_weight'] = df['times']
    if len(df) >= 2:
        if df.loc[0, 'times'] >= 5 * df.loc[1, 'times']:
            df.loc[0, "times_weight"] *= 20
        elif df.loc[0, 'times'] >= 1.5 * df.loc[1, 'times']:
            df.loc[0, "times_weight"] *= 10
    df.to_csv(MergeFile, sep='\t', index=True)  # Final Prediction Phage

    grouped = df.groupby('host')
    sum_times = grouped['times'].sum()
    sum_times_weight = grouped['times_weight'].sum()
    host_count = pd.DataFrame({'times': sum_times, 'times_weight': sum_times_weight}).reset_index()
    host_count.to_csv(HostCountFile, sep='\t', index=True)

    host_count.loc[:, 'host'] = host_count['host'].str.split().str[:2].apply(' '.join)
    host2_count = host_count.groupby('host').sum().sort_values(by="times", ascending=False).reset_index()
    host2_count.to_csv(Host2CountFile, sep='\t', index=True)

    # find "times" column maxium
    max_row_index = host2_count['times_weight'].idxmax()
    max_row_host = host2_count.loc[max_row_index, 'host']
    max_row_times = host2_count.loc[max_row_index, 'times']
    # sum of "times"
    total_sum = host2_count['times'].sum()
    # Final Prediction and True Positive
    ratio = max_row_times / (total_sum)
    print("Final Prediction Host:", max_row_host)
    print("TP:", ratio)

    # True Positive after delete empty host
    if ratio != 1.0:
        empty_host = host2_count.loc[host2_count['host'] == '']
        if not empty_host.empty:
            empty_times = int(empty_host['times'].iloc[0])
            ratio_deleteEmpty = max_row_times / (total_sum - empty_times)
            print("TP_deleteEmpty:", ratio_deleteEmpty)
        else:
            print("No empty host found.TP_deleteEmpty:", ratio)
    return

def run_spacer_blast(infile, outfile, blastmode):
    ### 这个基于NCBI爆炸的网络囊胚函数被用来预测假定的原太空源。
    ### 搜索病毒数据库
    print(infile)
    UpLoadQuery = quote(open(infile, 'r').read())
    if blastmode == 'relax':
        arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&NUCL_REWARD=1&WORD_SIZE=16&GAPCOSTS=5 2&NUCL_PENALTY=-1&DATABASE=nt&ENTREZ_QUERY=txid10239[ORGN]&QUERY=" + UpLoadQuery
    else:
        arguments = "CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&ENTREZ_QUERY=txid10239[ORGN]&QUERY=" + UpLoadQuery

    r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?" + arguments)
    try:
        RID = ""
        RID = r.text.split("RID = ")[1].split()[0]
        while RID =="RTOE":
            print("RID extraction failed, re-searching...")
            time.sleep(2)
            r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?" + arguments)
            RID = r.text.split("RID = ")[1].split()[0]
    except:
        while RID == "":
            print("RID extraction failed, re-searching...")
            time.sleep(2)
            r = requests.put("https://blast.ncbi.nlm.nih.gov/Blast.cgi?" + arguments)
            RID = r.text.split("RID = ")[1].split()[0]

    print("SPACER BLAST PROGRAM RID: " + RID)

    rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + RID)
    Status = rStatus.text.split("Status=")[1].split()[0]

    while Status != "READY":
        rStatus = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + RID)
        Status = rStatus.text.split("Status=")[1].split()[0]
        print("Refresh in 30s, status: " + Status)
        time.sleep(30)

    rResult = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Tabular&RID=" + RID)

    with open(outfile + ".temp", 'w') as fot:
        fot.write(rResult.text)
    fot.close()
    flag = 0

    with open(outfile + ".temp", 'r') as fotr, open(outfile, 'w') as foc:
        for fotrecord in fotr.readlines():
            if flag == 1:
                if fotrecord.strip() != "</PRE>" and not fotrecord.startswith("#"):
                    foc.write(fotrecord)
            if "hits found" in fotrecord:
                flag = 1
            if fotrecord.strip() == "</PRE>":
                flag = 0

    os.system("rm " + outfile + ".temp")
    r.close()
    rStatus.close()
    rResult.close()
    return True



### 设置get_session_post_text函数的超时时间，超时时间限制为10秒
@func_set_timeout(10)
def get_session_post_text(inid):
    urls = "https://www.ncbi.nlm.nih.gov/ipg/" + inid
    session = requests.Session()
    r5 = session.get(urls)
    url = "https://www.ncbi.nlm.nih.gov/sviewer/ipg/ipg.cgi?"
    post_data = {
        'query': 'omitHeader%3Dfalse%26wt%3Dxml%26indent%3Dtrue%26rows%3D50%26start%3D0%26q%3D' + inid + '%26sort%3Dpriority_i%2520asc'}
    post_data_1 = {'db': 'ipg', 'solr': '1', 'id': inid}
    session.post(url, data=post_data_1).text
    r6 = session.post(url, data=post_data).text
    session.cookies.clear()
    return r6


@func_set_timeout(20)
def get_protein_sequences(inurl):
    proteintext = requests.get(inurl).text
    return proteintext


@func_set_timeout(20)
def get_genome_fasta(inurl):
    genometext = requests.get(inurl).text
    return genometext


def get_seq_dump_file(insig, outdump, inmaxsize):
    infodict = {}
    tempdict = {}
    duplicateBox = []
    genomedupbox = []

    genomepage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    proteinpage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
    with open(insig, 'r') as fa, open(outdump, 'w') as fb:
        for faline in fa.readlines():
            if len(duplicateBox) == inmaxsize:
                break
            sbjctID = faline.split("\t")[1]
            protein2fasta = proteinpage + sbjctID + "&rettype=fasta"
            retry_counter = 0
            while retry_counter < 3:
                try:
                    proteinfasta = get_protein_sequences(protein2fasta)
                    break
                except:
                    retry_counter += 1
                    print("Retry in: " + sbjctID + ", retry time: " + str(retry_counter))
            try:
                proteinsequence = "".join(proteinfasta.split("\n")[1:])
            except:
                ### Into the next loop
                print("Unknown error in: " + sbjctID)
                continue
            if not proteinsequence in duplicateBox:
                URLCorrectFlag = True
                retry_counter = 0
                ### Set timeout for get_session_post_text function, time out limit is 10s
                while retry_counter < 3:
                    try:
                        r6 = get_session_post_text(sbjctID)
                        break
                    except:
                        retry_counter += 1
                        print("Retry in: " + sbjctID + ", retry time: " + str(retry_counter))
                genomerefheader = "<str name=\"cds_accver_s\">"
                genome_tax_header = "<str name=\"org_s\">"
                taxname_header = "<str name=\"org_s\">"
                orfstartheader = "<long name=\"cds_start_l\">"
                orfendheader = "<long name=\"cds_stop_l\">"
                try:
                    refid = r6.split(genomerefheader)[1].split("<")[0]
                    genometaxid = r6.split(genome_tax_header)[1].split("<")[0].split("|")[1]
                    taxname = r6.split(taxname_header)[1].split("<")[0].split("|")[0]
                    taxname = " ".join(taxname.split()[:2])
                    orfstart = r6.split(orfstartheader)[1].split("<")[0]
                    orfend = r6.split(orfendheader)[1].split("<")[0]
                except:
                    URLCorrectFlag = False
                if URLCorrectFlag:
                    duplicateBox.append(proteinsequence)
                    print(
                        "Putative homolog protein : " + sbjctID + ", source: " + refid + ", taxid: " + genometaxid + ", tax name: " + taxname + ", orf start: " + orfstart + ", orf end: " + orfend)
                    if not sbjctID in infodict.keys() and not refid in genomedupbox:
                        infodict[sbjctID] = [refid, taxname, int(orfstart), int(orfend)]
                        genomedupbox.append(refid)
                        tempdict[refid] = sbjctID
                        gi2fasta = genomepage + refid + "&rettype=fasta"
                        ### retry for 3 times
                        retry_counter = 0
                        while retry_counter < 3:
                            try:
                                genomefasta = get_genome_fasta(gi2fasta)
                                if ">" in genomefasta and not "Error" in genomefasta and not "error" in genomefasta:
                                    fb.write(genomefasta[:-1])
                                    fb.flush()
                                    break
                            except:
                                retry_counter += 1
                                print("Time out in: " + sbjctID + ", retry in " + str(retry_counter) + " time(s)...")
                                continue
                else:
                    print("Failed ID convert for: " + sbjctID)
        fa.close()
        fb.close()

        ### Now, check fasta file
        with open(outdump, 'r') as fa, open(outdump + ".tmp", 'w') as fb:
            for faline in fa.readlines():
                if faline.strip() == "" or "error" in faline or "Error" in faline:
                    continue
                else:
                    fb.write(faline)
        fa.close()
        fb.close()

        os.system("mv " + outdump + ".tmp " + outdump)

        refseqs = Fasta(outdump)
        for refseq in refseqs.keys():
            if refseq.split(".")[0] in tempdict.keys():
                targetsbj = tempdict[refseq.split(".")[0]]
                ### Now, revise the infodict refid
                infodict[targetsbj][0] = refseq

    return infodict


@func_set_timeout(10)
def get_genome_request_text(inurl):
    r = requests.get(inurl)
    return r.text


def get_genome_seqdump_files(inhit, outfile):
    seqdumplimit = 500
    seqdumpcounter = 0
    if seqdumpcounter >= 298:
        print("Seqdump numbers is reaching the limit, if you need a larger number, please modify the code.")
    genomepage = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    duplicatebox = []
    with open(inhit, 'r') as fa, open(outfile, 'w') as fb:
        for faline in fa.readlines():
            sbjctid = faline.split()[1]
            if not sbjctid in duplicatebox:
                if seqdumpcounter >= seqdumplimit:
                    break
                duplicatebox.append(sbjctid)
                gi2fasta = genomepage + sbjctid + "&rettype=fasta"
                retry_counter = 0
                ### if false, retry 3 times
                while retry_counter < 3:
                    try:
                        genomefasta = get_genome_request_text(gi2fasta)
                        if ">" in genomefasta and not "Error" in genomefasta and not "error" in genomefasta:
                            seqdumpcounter += 1
                            fb.write(genomefasta[:-1])
                            break
                        else:
                            continue
                    except:
                        retry_counter += 1
                        print("Time out in: " + sbjctid + ", retry for " + str(retry_counter) + " time(s)...")
                        continue
            else:
                continue
    fa.close()
    fb.close()
    return True


def main(infile, outfile, moder):
    if moder == "blastp":
        run_online_blastp(infile, outfile)
    else:
        run_spacer_blast(infile, outfile)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
