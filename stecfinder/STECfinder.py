#!/usr/bin/env python
# Developers: Michael Payne, Thanh Nguyen
# version 0.5
import argparse
import os
import sys
import subprocess
import json
import re
import uuid
import shutil
from time import sleep as sl


def file_type(f, m):
    if m == 'a' and os.path.splitext(f)[-1] == ".fasta":
        return True
    elif m == 'a' and os.path.splitext(f)[-1] == ".fa":
        return True
    elif m == 'a' and os.path.splitext(f)[-1] == ".fna":
        return True
    elif m == 'r':
        if f.endswith('.fastq.gz') or f.endswith('.fastq'):
            return True
    return False


def get_json_data():
    curr_dir = get_currdir()
    json_file = os.path.join(curr_dir, 'resources/genes.json')
    with open(f"{json_file}") as f:
        data = json.load(f)
    return data

def get_stx_data():
    curr_dir = get_currdir()
    json_file = os.path.join(curr_dir, 'resources/stx-subtype-acc.json')
    with open(json_file) as f:
        data = json.load(f)
    return data

def ipaH_detect(genes):
    if 'ipaH' in genes:
        return True
    return False

def stx_detect(genes):
    stx_conv = get_stx_data()
    stx_out = []
    for i in genes:
        if "stx" in i:
            fixi = i.replace(":","_")
            stx_out.append(stx_conv[fixi])

    stx_out = list(set(stx_out))
    if len(stx_out) > 0:
        return ",".join(stx_out)
    else:
        return False

def delete_genes(remove, genes):
    for g in remove:
        del genes[g]
    return genes

def antigen_search(genes):
    h = {}
    o = {}
    for g in genes:
        g = g.replace(" ","")
        if g.startswith("fliC"):
            antigen = g.split("_")[-1]
            h[antigen] = antigen
        elif g.startswith("wz"):
            antigen = g.split("_")[-1]
            o[antigen] = g.split("_")[0] + "_" + g.split("_")[-1]

    return list(o.values()),list(h.values())

def gene_rename(gene):
    return gene.replace(":","_")

def top_ranked_oantigen(genes_set):
    genetypes = ["wzx","wzy","wzt","wzm"]
    tophits = {x:["",0] for x in genetypes}
    for gene in genes_set:
        if gene.startswith("wz"):
            genetype = gene[:3]
            genescore = genes_set[gene]['score']
            if genescore > tophits[genetype][1]:
                tophits[genetype] = [gene,genescore]
    remove = []
    for gene in genes_set:
        if gene.startswith("wz"):
            genetype = gene[:3]
            if gene != tophits[genetype][0]:
                remove.append(gene)
    genes = delete_genes(remove, genes_set)
    return genes


def blastn_cleanup(blast,args):
    blast.remove('')
    genes_set = {}
    outhits = []
    for line in blast:
        info = line.rstrip().split('\t')
        gene = gene_rename(info[0])
        start = int(info[3])
        end = int(info[4])
        perc_identity = float(info[5])
        score = float(info[5])
        len_coverage = 100 * float(info[2]) / float(info[1])
        if start > end:
            save = start
            start = end
            end = save

        if gene in genes_set.keys():
            genes_set[gene]['positions'].extend(list(range(start, end + 1)))
            genes_set[gene]['match'] = len(set(genes_set[gene]['positions']))
            genes_set[gene]['len_coverage'] = 100 * genes_set[gene]['match'] / genes_set[gene]['slength']
            if perc_identity > genes_set[gene]['pident']:
                genes_set[gene]['pident'] = perc_identity
            continue

        if gene == "ipaH":
            cut = args.ipaH_length
        elif gene.startswith("stx"):
            cut = args.stx_length
        elif gene.startswith("STEC"):
            cut = args.length
        elif gene.startswith("wz"):
            cut = args.o_length
        elif gene.startswith("fl"):
            cut = args.h_length
        else:
            cut = 0


        if float(len_coverage) >= cut:
            genes_set[gene] = {}
            genes_set[gene]['slength'] = float(info[1])
            genes_set[gene]['match'] = float(info[2])
            genes_set[gene]['positions'] = list(range(start, end + 1))
            genes_set[gene]['pident'] = perc_identity
            genes_set[gene]['len_coverage'] = 100 * genes_set[gene]['match'] / genes_set[gene]['slength']
            genes_set[gene]['score'] = score
            outhits.append(line)

    return genes_set, outhits

# def h_top_len(genes):
#     remove = []
#     hant = {}
#     for g, l in genes.items():
#         if "fliC" in g:
#             hant[g] = l
#
#     if len(hant) > 0:
#         max_len = hant[min(hant, key=(lambda k: abs(hant[k]['len_coverage'] - 100)))]
#         for a, l in hant.items():
#
#             if l['len_coverage'] != max_len['len_coverage']:
#                 remove.append(a)
#
#     delete_genes(remove, genes)
#     return genes

def top_ranked_hantigen(genes_set):
    tophit = ["",0]
    for gene in genes_set:
        if gene.startswith("fliC"):
            genescore = genes_set[gene]['score']
            if genescore > tophit[1]:
                tophit = [gene,genescore]
    remove = []
    for gene in genes_set:
        if gene.startswith("fliC"):
            if gene != tophit[0]:
                remove.append(gene)
    genes = delete_genes(remove, genes_set)
    return genes

# def h_top_percid(genes):
#     remove = []
#     hant = {}
#     for g, l in genes.items():
#         if "fliC" in g:
#             hant[g] = l
#
#     if len(hant) > 0:
#         max_pident = hant[min(hant, key=(lambda k: abs(hant[k]['pident'] - 100)))]
#         for a, l in hant.items():
#             if l['pident'] != max_pident['pident']:
#                 remove.append(a)
#
#     delete_genes(remove, genes)
#     return genes

def top_ranked_stx(genes_set):

    genetypes = ["stx1", "stx2"]
    tophits = {x: ["", 0] for x in genetypes}
    for gene in genes_set:
        if gene.startswith("stx"):
            genetype = gene[:4]
            genescore = genes_set[gene]['score']
            if genescore > tophits[genetype][1]:
                tophits[genetype] = [gene, genescore]
    remove = []
    for gene in genes_set:
        if gene.startswith("stx"):
            genetype = gene[:4]
            if gene != tophits[genetype][0]:
                remove.append(gene)
    genes = delete_genes(remove, genes_set)

    # remove = []
    # hant = {}
    # for g, l in genes.items():
    #     if "stx1" in g:
    #         hant[g] = l
    #
    # if len(hant) > 0:
    #     max_pident = hant[min(hant, key=(lambda k: abs(hant[k]['pident'] - 100)))]
    #     for a, l in hant.items():
    #
    #         if l['pident'] != max_pident['pident']:
    #             remove.append(a)
    #
    # hant = {}
    # for g, l in genes.items():
    #     if "stx2" in g:
    #         hant[g] = l
    #
    # if len(hant) > 0:
    #     max_pident = hant[min(hant, key=(lambda k: abs(hant[k]['pident'] - 100)))]
    #     for a, l in hant.items():
    #
    #         if l['pident'] != max_pident['pident']:
    #             remove.append(a)
    #
    # delete_genes(remove, genes)
    return genes

def h_duplicate_remove(genes):
    remove = []
    hant = {}
    for g, l in genes.items():
        if "fliC" in g:
            hant[g] = l
            if g in genes.keys():
                continue
            remove.append(g)

    genes.update(hant)
    delete_genes(remove, genes)
    return genes


def map_depth_ratios(use_kma,bam):
    genes_set = []
    depth_cut = mapping_depth_cutoff(use_kma,bam)
    genes_set.append('average 7 HS genes: ' + str(depth_cut))

    if use_kma:
        depthcol = 8
    else:
        depthcol = 6

    for line in bam:
        info = line.split('\t')
        if len(info) > 1 and '#rname' not in line and '#Template' not in line:
            gene = gene_rename(info[0])
            meandepth = float(info[depthcol])
            ratio = 100 * meandepth / depth_cut
            genes_set.append(gene + '\t' + str(ratio))
    return genes_set


def mapping_depth_cutoff(use_kma,bam):
    mlst = ["NC_000913.3:recA", "NC_000913.3:purA", "NC_000913.3:mdh", "NC_000913.3:icd", "NC_000913.3:gyrB",
            "NC_000913.3:fumC", "NC_000913.3:adk"]
    depth = 0

    if use_kma:
        depthcol = 8
    else:
        depthcol = 6

    for line in bam:
        info = line.split('\t')
        if len(info) > 1 and '#rname' not in line and '#Template' not in line:
            gene = info[0]
            meandepth = float(info[depthcol])
            if gene in mlst:
                depth += meandepth

    return depth / 7


def determine_cluster(genes,data):

    cluster_list = {}
    # Check for the all cluster-specific genes found in genes
    for s in data:
        if s != "sporadic":
            if all(item in genes.keys() for item in data[s]["cluster-genes"]):
                cluster_list[s] = data[s]["cluster-genes"]
    big10 = ''
    if len(cluster_list) == 1:
        if list(cluster_list.keys())[0] == "C5":
            c = "C5"
            big10 = "O121:H19"
        elif list(cluster_list.keys())[0] == "C6":
            c = "C6"
            big10 = "O145:H28"
        elif list(cluster_list.keys())[0] == "C7":
            c = "C7"
            big10 = "O91:H14"
        elif list(cluster_list.keys())[0] == "C8":
            c = "C8"
            big10 = "O146:H21"
        elif list(cluster_list.keys())[0] == "C9":
            c = "C9"
            big10 = "O146:H21"
        elif list(cluster_list.keys())[0] == "O103H2":
            c = "Unknown Cluster"
            big10 = "O103:H2"
        elif list(cluster_list.keys())[0] == "O123H2":
            c = "Unknown Cluster"
            big10 = "O123:H2"
        elif list(cluster_list.keys())[0] == "O103H11":
            c = "Unknown Cluster"
            big10 = "O103:H11"
        elif list(cluster_list.keys())[0] == "O118H16":
            c = "Unknown Cluster"
            big10 = "O118:H16"
        elif list(cluster_list.keys())[0] == "O111H8":
            c = "Unknown Cluster"
            big10 = "O111:H8"
        elif list(cluster_list.keys())[0] == "O26H11":
            c = "Unknown Cluster"
            big10 = "O26:H11"
        elif list(cluster_list.keys())[0] in ["O45H2","O45H2-C3","O45H2-AM37"]:
            c = "Unknown Cluster"
            big10 = "O45:H2"
        else:
            c =list(cluster_list.keys())[0]
    elif 'AM6' in cluster_list.keys() and 'O45H2-AM37' in cluster_list.keys():
        c ='AM6'
    elif 'B1M126' in cluster_list.keys() and 'B1M2' in cluster_list.keys():
        c ='B1M126'
    elif 'B1M27' in cluster_list.keys() and 'B1M36' in cluster_list.keys():
        c ='B1M27'
    elif 'B1M28' in cluster_list.keys() and 'AM28' in cluster_list.keys():
        c ='B1M28'
    elif 'B1M30' in cluster_list.keys() and 'B1M78' in cluster_list.keys():
        c ='B1M30'
    elif 'B1M30' in cluster_list.keys() and 'C15' in cluster_list.keys():
        c ='B1M30'
    elif 'B1M57' in cluster_list.keys() and 'B1M33' in cluster_list.keys():
        c ='B1M57'
    elif 'B1M35' in cluster_list.keys() and 'B1M114' in cluster_list.keys():
        c ='B1M35'
    elif 'B1M36' in cluster_list.keys() and 'B1M114' in cluster_list.keys():
        c ='B1M36'
    elif 'B1M40' in cluster_list.keys() and 'C8' in cluster_list.keys():
        c ='B1M40'
    elif 'B1M48' in cluster_list.keys() and 'B1M114' in cluster_list.keys():
        c ='B1M48'
    elif 'B1M47' in cluster_list.keys() and 'B1M46' in cluster_list.keys():
        c ='B1M47'
    elif 'GM2' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c ='GM2'
    elif 'GM2' in cluster_list.keys() and 'GM4' in cluster_list.keys():
        c ='GM2'
    elif 'B1M120' in cluster_list.keys() and 'B1M121' in cluster_list.keys():
        c ='B1M120'
    elif 'C12' in cluster_list.keys() and 'CM4' in cluster_list.keys():
        c ='C12'
    elif 'C13' in cluster_list.keys() and 'CM5' in cluster_list.keys():
        c ='C13'
    elif 'C13' in cluster_list.keys() and 'B1M14' in cluster_list.keys():
        c ='C13'
    elif 'C13' in cluster_list.keys() and 'B1M49' in cluster_list.keys():
        c ='C13'
    elif 'C3' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c ='C3'
    elif 'C3' in cluster_list.keys() and 'B1M77' in cluster_list.keys():
        c ='C3'
    elif 'C3' in cluster_list.keys() and 'CM5' in cluster_list.keys():
        c ='C3'
    elif 'C8' in cluster_list.keys() and 'AM12' in cluster_list.keys():
        c ='C8'
    elif 'C4' in cluster_list.keys() and 'B1M78' in cluster_list.keys():
        c ='C4'
    elif 'C4' in cluster_list.keys() and 'B2M5' in cluster_list.keys():
        c ='C4'
    elif 'C5' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c ='C5'
    elif 'C5' in cluster_list.keys() and 'B1M77' in cluster_list.keys():
        c ='C5'
    elif 'O157H7' in cluster_list.keys() and 'C7' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'AM18' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'B1M25' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'B1M31' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'B2M7' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'DM14' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'O103H2' in cluster_list.keys():
        c ='O157H7'
    elif 'O157H7' in cluster_list.keys() and 'DM8' in cluster_list.keys():
        c ='O157H7'
    elif 'C18' in cluster_list.keys() and 'O157H7' in cluster_list.keys():
        c ='C18'
    elif 'C14' in cluster_list.keys() and 'B2M11' in cluster_list.keys():
        c ='C14'
    elif 'C17' in cluster_list.keys() and 'B1M14' in cluster_list.keys():
        c ='C17'
    elif 'C17' in cluster_list.keys() and 'B1M33' in cluster_list.keys():
        c ='C17'
    #elif 'AM20' in cluster_list.keys() and 'C14' in cluster_list.keys():
        #c ='AM20'
    elif 'B1M50' in cluster_list.keys() and 'C17' in cluster_list.keys():
        c ='B1M50'
    elif 'B1M74' in cluster_list.keys() and 'DM9' in cluster_list.keys():
        c ='B1M74'
    elif 'C7' in cluster_list.keys() and 'B1M67' in cluster_list.keys():
        c ='C7'
    elif 'C7' in cluster_list.keys() and 'B2M11' in cluster_list.keys():
        c ='C7'
    elif 'C7' in cluster_list.keys() and 'C9' in cluster_list.keys():
        c ='C7'
    elif 'B1M8' in cluster_list.keys() and 'B1M9' in cluster_list.keys():
        c ='B1M8'
    elif 'B1M82' in cluster_list.keys() and 'B1M14' in cluster_list.keys():
        c ='B1M82'
    elif 'B1M92' in cluster_list.keys() and 'B1M93' in cluster_list.keys():
        c ='B1M92'
    elif 'B1M94' in cluster_list.keys() and 'EM6' in cluster_list.keys():
        c ='B1M94'
    elif 'B1M95' in cluster_list.keys() and 'B1M94' in cluster_list.keys():
        c ='B1M95'
    elif 'B1M5' in cluster_list.keys() and 'B1M62' in cluster_list.keys():
        c ='B1M5'
    elif 'B1M5' in cluster_list.keys() and 'DM8' in cluster_list.keys():
        c ='B1M5'
    elif 'B1M50' in cluster_list.keys() and 'B1M121' in cluster_list.keys():
        c ='B1M50'
    elif 'B2M6' in cluster_list.keys() and 'B1M27' in cluster_list.keys():
        c ='B2M6'
    elif 'DM13' in cluster_list.keys() and 'DM14' in cluster_list.keys():
        c ='DM13'
    elif 'CM3' in cluster_list.keys() and 'CM1' in cluster_list.keys():
        c ='CM3'
    elif 'EM5' in cluster_list.keys() and 'EM17' in cluster_list.keys():
        c ='EM5'
    elif 'CM2' in cluster_list.keys() and 'CM1' in cluster_list.keys():
        c ='CM2'
    elif 'B1M80' in cluster_list.keys() and 'B1M43' in cluster_list.keys():
        c ='B1M80'
    elif 'B1M41' in cluster_list.keys() and 'B1M43' in cluster_list.keys():
        c ='B1M41'
    elif 'C8' in cluster_list.keys() and 'C7' in cluster_list.keys():
        c ='C8'
    elif 'AM9' in cluster_list.keys() and 'O157H7' in cluster_list.keys():
        c ='AM9'
    elif 'DM2' in cluster_list.keys() and 'O157H7' in cluster_list.keys():
        c ='DM2'
    elif 'B1M41' in cluster_list.keys() and 'C8' in cluster_list.keys():
        c ='B1M41'
    elif 'GM3' in cluster_list.keys() and 'DM8' in cluster_list.keys():
        c ='GM3'
    elif 'B1M44' in cluster_list.keys() and 'B1M32' in cluster_list.keys():
        c ='B1M44'
    elif 'B1M42' in cluster_list.keys() and 'C8' in cluster_list.keys():
        c ='B1M42'
    elif 'GM1' in cluster_list.keys() and 'AM25' in cluster_list.keys():
        c ='GM1'
    elif 'C1' in cluster_list.keys() and 'O103H11' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c = "C1"
        big10 = 'O103:H11'
    elif 'C1' in cluster_list.keys() and 'O118H16' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c = "C1"
        big10 = 'O118/O151:H16'
    elif 'C1' in cluster_list.keys() and 'O111H8' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c = "C1"
        big10 = 'O111:H8'
    elif 'C1' in cluster_list.keys() and 'O26H11' in cluster_list.keys():
        c = "C1"
        big10 = 'O26:H11'
    elif 'C1' in cluster_list.keys() and 'O111H8' in cluster_list.keys():
        c = "C1"
        big10 = 'O111:H8'
    elif 'C3' in cluster_list.keys() and 'O45H2-C3' in cluster_list.keys():
        c = "C3"
        big10 = 'O45:H2'
    elif 'AM37' in cluster_list.keys() and 'O45H2-AM37' in cluster_list.keys():
        c = "AM37"
        big10 = 'O45:H2'
    elif 'C2' in cluster_list.keys() and 'O45H2' in cluster_list.keys() and 'O103H2' in cluster_list.keys():
        c = "C2"
        big10 = 'O45:H2'
    elif 'C2' in cluster_list.keys() and 'O123H2' in cluster_list.keys() and 'O103H2' in cluster_list.keys():
        c = "C2"
        big10 = 'O123:H2'
    elif 'C2' in cluster_list.keys() and 'O103H2' in cluster_list.keys():
        c = "C2"
        big10 = 'O103:H2'

    elif len(cluster_list) > 1:
        c = ','.join(cluster_list.keys())
    else:
        c = "Unknown Cluster"

    return c, big10

def string_result(res, output):
    result = res['sample'] + "\t" + str(res['stx']) + "\t" + res['cluster'] + "\t" + \
             res['big10'] + "\t" + res['serotype'] + "\t" + res['cluster_serotype'] + "\t" + res['oantigens']\
             + "\t" + res['hantigens'] + "\t" + res['ipaH'] + "\t" + res['notes']
    if output:
        result += "\n"

    return result


def get_gene_type(gene):
    gene_type = ""

    mlst = ["NC_000913.3:recA", "NC_000913.3:purA", "NC_000913.3:mdh", "NC_000913.3:icd", "NC_000913.3:gyrB",
            "NC_000913.3:fumC", "NC_000913.3:adk"]
    plasmid = ["acp", "icsA (virG)", "icsA", "icsB", "ipaA", "ipaB", "ipaC", "ipaD", "ipaJ", "ipgA", "ipgB1", "ipgC",
               "ipgD", "ipgE", "ipgF", "mxiA", "mxiC", "mxiD", "mxiE", "mxiG", "mxiH", "mxiI", "mxiJ", "mxiK", "mxiL",
               "mxiM", "mxiN", "spa13", "spa15", "spa24", "spa29", "spa32", "spa33", "spa40", "spa47", "spa9", "virA",
               "virB", "virF"]
    if gene in mlst:
        gene_type = "House Keeping"
    elif gene.startswith("fliC")or gene.startswith("wz"):
        gene_type = "O/H-antigen gene"
    elif gene.startswith("stx"):
        gene_type = "stx toxin gene"
    elif gene.startswith("STEC"):
        gene_type = "Cluster-Specific"
    elif gene == "ipaH":
        gene_type = "ipaH"

    return gene_type

def run_blast(dir, fileA):
    # Get all genes for ipaH & cluster genes
    qry = f'blastn -db "{dir}/resources/genes.fasta" -outfmt "6 sseqid slen length sstart send pident bitscore" ' \
          f'-perc_identity 80 -query "{fileA}"'
    blast_hits = subprocess.check_output(qry, shell=True, stderr=subprocess.STDOUT)
    blast_hits = blast_hits.decode("ascii").split('\n')
    return blast_hits

def run_kma(dir, r1, r2, threads,tmp,strain_id):

    if not os.path.exists(tmp):
        os.mkdir(tmp)
    kma_db = dir + "/resources/genes.fasta"
    kma_cmd = f'kma -mct 0.001 -ipe "{r1}" "{r2}" -t_db "{kma_db}" -t {threads} -ConClave 2 -mrs 0.001 -mrc 0.001 -ID 1 -o "{tmp}/{strain_id}kmatmp_out"'
    subprocess.run(kma_cmd + ">/dev/null 2>&1", shell=True)

def add_duped_genes(genes_set):
    """
    necessary because KMA removes duplicate sequences in its database and
    several genes are in multiple cluster specific gene sets
    """

    curr_dir = get_currdir()
    json_file = os.path.join(curr_dir, 'resources/duplicate_cluster_genes.json')
    genes_set2 = {}
    f = open(f"{json_file}","r")
    dupe_dict = json.load(f)
    for gene in genes_set:
        indupe = False
        for dupeset in dupe_dict:
            dupels = dupe_dict[dupeset]
            if gene in dupels:
                indupe = True
                for g in dupels:
                    genes_set2[g] = genes_set[gene]
                    # print(g,gene)
        if not indupe:
            genes_set2[gene] = genes_set[gene]
    return genes_set2,dupe_dict




def genes_frm_kma_output(strain_id,args):

    outfile = f"{args.tmpdir}/{strain_id}kmatmp_out.res"
    genes_set = {}
    outhits = []
    if os.path.exists(outfile):
        hit_results = open(outfile,'r').read().splitlines()

        depth_cut = mapping_depth_cutoff(args.use_kma, hit_results)
        for l in hit_results[1:]:
            col = l.split("\t")
            gene = col[0].replace(" ","")

            lengthperc = col[5]
            percid = col[4]
            readdepth = col[8]
            slen = col[3]
            score = col[1]

            if gene == "ipaH":
                cut = args.ipaH_length
                dcut = (float(args.ipaH_depth)/100)*depth_cut
            elif gene.startswith("stx"):
                cut = args.stx_length
                dcut = (float(args.stx_depth)/100)*depth_cut
            elif gene.startswith("fl"):
                cut = args.h_length
                dcut = (float(args.h_depth)/100)*depth_cut
            elif gene.startswith("wz"):
                cut = args.o_length
                dcut = (float(args.o_depth)/100)*depth_cut
            elif gene.startswith("STEC"):
                cut = args.length
                dcut = (float(args.cutoff) / 100) * depth_cut
            else:
                cut = 0
                dcut = (float(args.cutoff)/100)*depth_cut

            if float(lengthperc) >= cut and float(readdepth) >= dcut:
                genes_set[gene] = {}
                genes_set[gene]['slength'] = float(slen)
                genes_set[gene]['pident'] = float(percid)
                genes_set[gene]['depth'] = float(readdepth)
                genes_set[gene]['len_coverage'] = float(lengthperc)
                genes_set[gene]['score'] = float(score)
                outhits.append(l)

    else:
        sys.exit(f'KMA output file missing at: {args.tmpdir}/{strain_id}kmatmp_out.res')
    shutil.rmtree(args.tmpdir) #TODO comment to keep tmp KMA output files
    genes_set,dupedict = add_duped_genes(genes_set)
    return genes_set,outhits

def collapse_multiple_o_antigen(oantigens):
    if oantigens == []:
        return oantigens
    genetypes = ["wzx", "wzy", "wzt", "wzm"]
    possibles = []
    """
    cases
    single O antigen - uncertain from second gene - keep single ( O186 and O123/O186 -> O186)
    uncertain O antigen - single second gene - keep single from second ( O123/O186 and O123 -> O123)
    unrelated O antigen ( O126 and O123/O186 -> O126/O123/O186)
    """
    for wz in genetypes:
        for o in oantigens:
            owz = o.split("_")[0]
            if owz == wz:
                olist = o.split("_")[-1].split("/")
                if possibles == []:
                    possibles = [o]
                else:
                    outposs = []
                    for existing in possibles:
                        ewz = existing.split("_")[0]
                        existolist = existing.split("_")[-1].split("/")
                        overlap = list(set(existolist).intersection(set(olist)))
                        if len(overlap) == 0:
                            outposs.append(existing)
                            outposs.append(o)
                        elif len(overlap) > 0:
                            o_out = [owz + "_" + x for x in olist if x in overlap]
                            eo_out = [ewz + "_" + x for x in olist if x in overlap]
                            outposs += o_out
                    possibles = outposs
    oonly = [x.split("_")[-1] for x in possibles]
    dedup_o = list(set(oonly))
    return dedup_o

def serotype_subset(test,existing):
    if "O" in test and "H" in test:
        if test == existing:
            return True
        else:
            existingO = set(existing.split(":")[0].split("/"))
            existingH = set(existing.split(":")[-1].split("/"))
            testO = set(test.split(":")[0].split("/"))
            testH = set(test.split(":")[-1].split("/"))
            if len(existingO.intersection(testO)) > 0 and len(existingH.intersection(testH)) > 0:
                return True
    elif "H" not in test:
        existingO = set(existing.split(":")[0].split("/"))
        testO = set(test.split(":")[0].split("/"))
        if len(existingO.intersection(testO)) > 0:
            return True
    elif "O" not in test:
        existingH = set(existing.split(":")[-1].split("/"))
        testH = set(test.split(":")[-1].split("/"))
        if len(existingH.intersection(testH)) > 0:
            return True
    return False

def cluster_aware_antigen_search(cluster,genes,clustergenes):


    ## get antigen genes
    oantigens, hantigens = antigen_search(genes)

    ## deduplicate O calls from different genes
    if len(oantigens) > 0:
        # print(oantigens)
        collapsed_o_calls = collapse_multiple_o_antigen(oantigens)
        # print(collapsed_o_calls)
        oa = "/".join(collapsed_o_calls)
    else:
        oa = "-"
    if len(hantigens) > 0:
        ha = list(set([x.split("_")[-1] for x in hantigens]))
        ha = "/".join(ha)
    else:
        ha = "-"
    serotype = oa + ":" + ha

    ## TODO cluster based calling
    """
    1 - if antigen not in cluster list add note
    2 - if partial antigen report possible complete antigens
    3 - if perfect match report with no note
    """
    cluster_serotypes = []
    if cluster not in ["Unknown Cluster", "missing_cluster_genes"]:
        cluster_serotypes = [x for y in cluster.split(",") for x in clustergenes[y] if x != "cluster-genes"]
        cluster_limited_serotypes = []
        for expected_serotypes in cluster_serotypes:
            match = serotype_subset(serotype,expected_serotypes)
            if match:
                cluster_limited_serotypes.append(expected_serotypes)

        if len(cluster_limited_serotypes) > 0:
            cluster_limited_serotypes = ",".join(cluster_limited_serotypes)
        else:
            cluster_limited_serotypes = "-"
    else:
        cluster_limited_serotypes = "N/A"

    return serotype, ",".join(oantigens), ",".join(hantigens), cluster_limited_serotypes, cluster_serotypes



def run_typing(dir, files, mode, args):
    result = {}
    if mode == "r":
        name = os.path.basename(files[1])
        result['sample'] = re.search(r'(.*)\_.*\.fastq\.gz', name).group(1)
        run_kma(dir, files[0], files[1], str(args.t), args.tmpdir,result['sample'])
        genes,hit_results = genes_frm_kma_output(result['sample'],args)
    else:
        blast_results = run_blast(dir, files)
        genes,hit_results = blastn_cleanup(blast_results,args)
        result['sample'] = os.path.basename(files).split('.')[0]


    genes = top_ranked_hantigen(genes)
    genes = top_ranked_stx(genes)
    genes = h_duplicate_remove(genes)
    genes = top_ranked_oantigen(genes)

    clustergenes = get_json_data()
    cluster, big10 = determine_cluster(genes,clustergenes)
    result['cluster'] = cluster
    result['big10'] = big10
    """
    ipah    stx	unclustered	outcome
    +	    +	+	        "ipaH+stx+ = Possible EIEC/Shigella, try out other tool shigeifinder!"
    +	    +	-	        "strain in STEC cluster but has ipaH - looks like EIEC/Shigella but clusters with STEC
    +	    -	+	        ipaH+stx- = Possible EIEC/Shigella, try out other tool shigeifinder!"
    +	    -	-	        "Strain with cluster specific genes but no stx detected and ipaH detected - looks like EIEC/Shigella but clusters with STEC"
    -	    +	+	        "STEC not from any major STEC lineages"
    -	    +	-	        CLUSTERED STX
    -	    -	+	        "ipaH-stx- = Non-STEC E.coli"
    -	    -	-	        "Strain with cluster specific genes but no stx detected"
    """
    ipaH = ipaH_detect(genes)
    stx = stx_detect(genes)
    serotype, oantigens, hantigens, cluster_restricted_serotype, cluster_serotypes = cluster_aware_antigen_search(cluster,genes,clustergenes)


    result['serotype'] = serotype
    result['cluster_serotype'] = cluster_restricted_serotype
    result['oantigens'] = oantigens
    result['hantigens'] = hantigens
    result['notes'] = ""
    if not ipaH and not stx and cluster == "Unknown Cluster":
        result['ipaH'] = '-'
        result['stx'] = '-'
        result['cluster'] = "Other_Ecoli"
        result['notes'] += "ipaH-stx- = Non-STEC E.coli"
    elif ipaH and stx and cluster == "Unknown Cluster":
        result['ipaH'] = '+'
        result['stx'] = stx
        result['cluster'] = "EIEC/Shigella"
        result['notes'] += "ipaH+stx+ = Possible EIEC/Shigella, try out other tool shigeifinder!"
    elif not ipaH and stx and cluster == "Unknown Cluster":
        result['ipaH'] = "-"
        result['stx'] = stx
        result['cluster'] = "Unclustered STEC"
        result['notes'] += "STEC not from any major STEC lineages."
    elif not ipaH and stx and cluster != "Unknown Cluster":
        result['ipaH'] = "-"
        result['stx'] = stx
        result['cluster'] = cluster
    elif not ipaH and not stx and cluster != "Unknown Cluster":
        result['ipaH'] = "-"
        result['stx'] = "-"
        result['cluster'] = cluster
        result['notes'] += "Strain with cluster specific genes but no stx detected."
    elif ipaH and not stx and cluster == "Unknown Cluster":
        result['ipaH'] = "+"
        result['stx'] = "-"
        result['cluster'] = cluster
        result['notes'] += "ipaH+stx- = Possible EIEC/Shigella, try out other tool shigeifinder!"
    elif ipaH and not stx and cluster != "Unknown Cluster":
        result['ipaH'] = "+"
        result['stx'] = "-"
        result['cluster'] = cluster
        result['notes'] += "Strain with cluster specific genes but no stx detected and ipaH detected - looks like EIEC/Shigella but clusters with STEC."
    elif ipaH and stx and cluster != "Unknown Cluster":
        result['ipaH'] = '+'
        result['stx'] = stx
        result['cluster'] = cluster
        result['notes'] += "Strain in STEC cluster but has ipaH and stx - looks like EIEC/Shigella but clusters with STEC."
    else:
        result['ipaH'] = "-"
        result['stx'] = "-"
        result['cluster'] = "-"
        result['notes'] += "Unexpected genetic combination, possible bug, please post issue on shigeifinder github repo"

    if big10 != '' and big10 != serotype:
        if len(oantigens.split(",")) == 0 or len(hantigens.split(",")) == 0:
            result['notes'] += " Possible missing antigen in big10 isolate."
        else:
            result['notes'] += " Missmatch between big10 specific genes and antigen genes, possible big10 specific genes false positive."

    if cluster_restricted_serotype == "-":
        result['notes'] += " Detected O and H antigen genes do not match previously known serotypes for this cluster ({})".format(",".join(cluster_serotypes))
    elif ("H" not in serotype and "H" in cluster_restricted_serotype) or ("O" not in serotype and "O" in cluster_restricted_serotype):
        result['notes'] += " Complete serotype inferred from known cluster serotypes"

    # Print result
    if args.output:
        outp = open(args.output, "a+")
        outp.write(string_result(result, args.output))
        outp.close()
    else:
        print(string_result(result, args.output))

    if args.output:
        outp = open(args.output, "a+")
        if args.hits:
            outp.write("--------------- GENE SET ---------------\n")
            outp.write("#gene\tlength_coverage/depth\tgene_type\n")
            for key, item in genes.items():
                outp.write(key + '\t' + str(item) + '\t' + get_gene_type(key) + "\n")
            if mode == "a":
                outp.write("---------- BLAST GENE HITS ----------\n")
                outp.write("#seqid\tslen\tlength\tsstart\tsend\tpident\n")
            elif args.use_kma:
                outp.write("---------- KMA GENE HITS ----------\n")
            else:
                outp.write("---------- BWA/SAMTOOLS GENE HITS ----------\n")
            outp.write('\n'.join(hit_results) + "\n")
            outp.write("----------------------------------------\n")

        if args.dratio:
            outp.write(
                "--------------- RATIOS OF DEPTH SPECFIC GENES TO AVERAGE DEPTH OF 7 HOUSE KEEPING GENES---------------\n")
            for g in map_depth_ratios(args.use_kma,hit_results):
                outp.write(str(g) + "\n")
            outp.write(
                "--------------------------------------------------------------------------------------------------------------\n")
        outp.close()
    else:
        if args.hits:
            print("--------------- GENE SET ---------------")
            print("#gene\tlength_coverage/depth\tgene_type")
            for key, item in genes.items():
                print(key + '\t' + str(item) + '\t' + get_gene_type(key))
            if mode == "a":
                print("---------- BLAST GENE HITS ----------")
                print("#seqid\tslen\tlength\tsstart\tsend\tpident")
            else:
                print("---------- BWA/SAMTOOLS GENE HITS ----------")
            print('\n'.join(hit_results))
            print("----------------------------------------")

        if args.dratio:
            print(
                "--------------- RATIOS OF DEPTH SPECIFIC GENES TO AVERAGE DEPTH OF 7 HOUSE KEEPING GENES---------------")
            for g in map_depth_ratios(args.use_kma,hit_results):
                print(g)
            print(
                "--------------------------------------------------------------------------------------------------------------")


def check_deps(checkonly, args):
    depslist = ["blastn","kma"]
    f = 0
    for dep in depslist:
        rc = subprocess.call(['which', dep], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if rc == 0:
            sys.stderr.write(f'{dep:<10}:{"installed":<10}\n')
        else:
            sys.stderr.write(f'{dep:<10}:{"missing in path, Please install ":<10}{dep}\n')
            f += 1
    if f > 0:
        sys.stderr.write("One or more dependencies are missing.\n")
        sys.exit(1)
    else:
        sys.stderr.write("All dependencies present.\n")
        if checkonly:
            sys.exit(0)
        else:
            return

def get_currdir():
    dir = os.path.dirname(os.path.realpath(__file__))
    if "OneDrive - UNSW" in dir:
        dir = dir.replace("OneDrive - UNSW", "OneDrive")
    return dir


def main():

    debug = False

    if not debug: ### COMMAND LINE MODE
        parser = argparse.ArgumentParser(
            usage='\nSTECFinder.py -i <input_data1> <input_data2> ... OR\nSTECFinder.py -i <directory/*> OR '
                  '\nSTECFinder.py -i <Read1> <Read2> -r [Raw Reads]\n')
        parser.add_argument("-i", nargs="+", help="<string>: path/to/input_data")
        parser.add_argument("-r", action='store_true', help="Add flag if file is raw reads.")
        parser.add_argument("-t", nargs=1, type=int, default='4', help="number of threads. Default 4.")
        parser.add_argument("--hits", action='store_true', help="To show the blast/alignment hits")
        parser.add_argument("--dratio", action='store_true',
                            help="To show the depth ratios of cluster-specific genes to House Keeping genes")
        parser.add_argument("--cutoff", type=float,
                            help="minimum read coverage for gene to be called",default="10.0")
        parser.add_argument("--length", type=float,
                            help="percentage of gene length needed for positive call",default="90.0")
        parser.add_argument("--ipaH_length", type=float,
                            help="percentage of ipaH gene length needed for positive call",default="10.0")
        parser.add_argument("--ipaH_depth", type=float,
                            help="percentage of ipaH gene length needed for positive call",default="10.0")
        parser.add_argument("--stx_length", type=float,
                            help="percentage of ipaH gene length needed for positive call",default="10.0")
        parser.add_argument("--stx_depth", type=float,
                            help="percentage of ipaH gene length needed for positive call",default="10.0")
        parser.add_argument("--o_length", type=float,
                            help="percentage of ipaH gene length needed for positive call", default="10.0")
        parser.add_argument("--o_depth", type=float,
                            help="percentage of ipaH gene length needed for positive call", default="10.0")
        parser.add_argument("--h_length", type=float,
                            help="percentage of ipaH gene length needed for positive call", default="10.0")
        parser.add_argument("--h_depth", type=float,
                            help="percentage of ipaH gene length needed for positive call", default="10.0")
        parser.add_argument("--update_db", action='store_true',
                            help="Add flag if you added new sequences to genes database.")
        parser.add_argument("--output",
                            help="output file to write to (if not used writes to stdout and tmp folder in current dir)")
        parser.add_argument("--check", action='store_true', help="check dependencies are installed")
        args = parser.parse_args()

        args.use_kma = True





    args.runuuid = str(uuid.uuid1())

    if args.output:
        outdir = os.path.dirname(args.output)
        args.tmpdir = outdir + "/" + args.runuuid
    else:
        args.tmpdir = args.runuuid


    # Directory current script is in
    dir = get_currdir()

    # run and get the intermediate files for blast and bwa (makes updating the genes db easier)
    if args.update_db:
        print(dir + "/resources/genes.fasta")
        subprocess.run(
            'makeblastdb -in "' + dir + '/resources/genes.fasta" -parse_seqids -blastdb_version 5 -title '
                                       '"Shigella/EIEC DB" -dbtype nucl -out "' + dir + '/resources/genes.fasta"',
            shell=True)
        print('makeblastdb -in "' + dir + '/resources/genes.fasta" -out "' + dir + '/resources/genes.fasta" -parse_seqids -blastdb_version 5 -title '
                                       '"Shigella/EIEC DB" -dbtype nucl')
        subprocess.run('kma index -i "' + dir + '/resources/genes.fasta" -o "' + dir + '/resources/genes.fasta"', shell=True)
        sys.exit()

    if args.check:
        check_deps(True, args)

    if args.dratio and not args.r:
        parser.error("-dratio requires -r. Only applies for raw reads.")
    if not args.i:
        parser.error("-i is required")

    if not args.check:
        check_deps(False, args)

    if len(sys.argv) == 0:
        os.system("python3.7 " + dir + "/STECFinder.py -h")
    else:
        outheader = "#SAMPLE\tSTX\tcluster\tbig10_serotype\tserotype\tcluster_serotype\toantigens\thantigens\tIPAH\tNotes"
        mode = 'a'
        if args.r:
            mode = 'r'
        for files in args.i:
            if "*" in args.i[0]:
                dir1 = files.replace("*", "")
                if not os.path.isdir(dir1):
                    sys.exit('Invalid Directory Input! Directory: ' + dir1)
                break
            else:
                if not os.path.isfile(files):
                    sys.exit('Invalid Input File(s)! File Not Found:' + files)
                if not file_type(files, mode):
                    sys.exit('Incorrect File Type! File:' + files)
        if mode == 'r':
            # Run Raw Reads version
            # Check that there is 2 Reads inputed
            if len(args.i) < 2 or len(args.i) % 2 != 0:
                if "*" in args.i[0] and len(args.i) == 1:
                    samples = set()
                    for files in os.listdir(dir1):
                        if files.endswith(".fastq.gz"):
                            path = dir1 + files
                            samples.add(path)
                    reads = sorted(samples)
                    if len(reads) % 2 != 0:
                        sys.exit('Missing Input File(s)!!')
                    i = 0
                    if args.output:
                        outp = open(args.output, "w")
                        outp.write(outheader+"\n")
                        outp.close()
                    else:
                        print(outheader)
                    while i < len(reads):
                        f = [reads[i], reads[i + 1]]
                        run_typing(dir, f, mode, args)
                        i += 2
                    sys.exit()
                else:
                    sys.exit('Missing Input File(s)!!')
            elif len(args.i) > 2:
                files = sorted(args.i)
                i = 0
                if args.output:
                    outp = open(args.output, "w")
                    outp.write(outheader+"\n")
                    outp.close()
                else:
                    print(outheader)
                while i < len(args.i):
                    f = [files[i], files[i + 1]]
                    run_typing(dir, f, mode, args)
                    i += 2
                sys.exit()
            if args.output:
                outp = open(args.output, "w")
                outp.write(outheader+"\n")
                outp.close()
            else:
                print(outheader)
            run_typing(dir, args.i, mode, args)
        else:
            # Run assembled genome version
            if args.output:
                outp = open(args.output, "w")
                outp.write(outheader+"\n")
                outp.close()
            else:
                print(outheader)
            if "*" in args.i[0]:
                list_files = os.listdir(dir1)
                for f in list_files:
                    path = dir1 + "/" + f
                    if file_type(path, mode):
                        run_typing(dir, path, mode,args)
            elif len(args.i) > 1:
                for f in args.i:
                    run_typing(dir, f, mode,args)
            else:
                run_typing(dir, args.i[0], mode, args)


if __name__ == '__main__':
    main()