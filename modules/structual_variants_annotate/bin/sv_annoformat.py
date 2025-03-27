#!/usr/bin/env python3

import pandas as pd
import argparse
import re

# Using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-sv", "--sv_table", type=str)
parser.add_argument("-l", "--wgs_pilot_gene_list", type=str)
parser.add_argument("-v", "--vep_tab", type=str)
parser.add_argument("-s", "--sample", type=str)
parser.add_argument("-o", "--outfile", type=str)
args = parser.parse_args()

# manta data
SV_OUT = args.sv_table
SV_TABLE_OUT_data = pd.read_csv(SV_OUT, sep="\t",low_memory=False)

# gene list data
genlist = args.wgs_pilot_gene_list
#genlist = "Genliste_somatisch.csv"
genlist_data = pd.read_csv(genlist, names=["genes"], header=None)
gene_lst = genlist_data["genes"].to_list()

# Check genelist input for " "
for i in range(len(gene_lst)):
    if ' ' in genlist_data.loc[i, "genes"]:
        raise ValueError('Space in gene name! Please correct!')

# group in TRA, INV, others (DEL/DUPS/...)
TRA_index = []
INV_index = []
others_index = []

for idx, sv_type in enumerate(SV_TABLE_OUT_data["SV_type"]):
    if sv_type == "TRA":
        TRA_index.append(idx)
    elif sv_type == "INV":
        INV_index.append(idx)
    else:
        others_index.append(idx)
          
# TRA data
TRA_data = SV_TABLE_OUT_data.loc[TRA_index, :]

# INV data
INV_data = SV_TABLE_OUT_data.loc[INV_index, :]

# other data (del/dup/ins)
others_data = SV_TABLE_OUT_data.loc[others_index, :]

# prepare TRA data - START ----------------------------------------------------
# Reset indices
TRA_data = TRA_data.reset_index(drop="TRUE")

# remove <TRA> values in ALT column
TRA_data = TRA_data[TRA_data["ALT"] != "<TRA>"]

# get full INV data
TRA_data_1_full = TRA_data[TRA_data["Annotation_mode"] == "full"]

# remove nan
TRA_data_2_full = TRA_data_1_full.dropna(subset=["Gene_name"])

# Reset indices
TRA_data_2_full = TRA_data_2_full.reset_index(drop="TRUE")

# filter gene list
gene_filtered_sv_idx = []
for i, genes in enumerate(TRA_data_2_full["Gene_name"]):
    if len(genes.split(";")) > 1:
        start_gene = genes.split(";")[0]
        end_gene = genes.split(";")[-1]
        if start_gene in gene_lst or end_gene in gene_lst:
            gene_filtered_sv_idx.append(i)
    elif len(genes.split(";")) == 1:
        if genes in gene_lst:
            gene_filtered_sv_idx.append(i)

TRA_data_3 = TRA_data_2_full.loc[gene_filtered_sv_idx,:]

# get ids and mate ids
# Reset indices
TRA_data_3 = TRA_data_3.reset_index(drop="TRUE")
TRA_mate_ids_lst = []
for i, mate_ids in enumerate(TRA_data_3["INFO"]):
    mate_id = "MantaBND:" + "".join((re.findall(r'[\d+\:]', mate_ids.split(";")[1]))[1:])
    TRA_mate_ids_lst.append(mate_id)
    
# make id list together
final_TRA_ids = list(TRA_data_3["ID"]) + TRA_mate_ids_lst

# get data from INV all
final_TRAs = TRA_data[TRA_data["ID"].isin(final_TRA_ids) == True]

# get RefSeq from split
final_TRAs_split = final_TRAs[final_TRAs["Annotation_mode"] == "split"]

# make dataframe !!!check Columns before submittung!!!
df_gene_refseq_TRA = pd.DataFrame()

header_col = ["gene", "Refseq"]

for col in header_col:
    df_gene_refseq_TRA[col] = ""
      
# Reset indices
final_TRAs_split = final_TRAs_split.reset_index(drop="TRUE")
    
for i in range(len(final_TRAs_split["Gene_name"])):
    df_gene_refseq_TRA.loc[i,"gene"] = final_TRAs_split.loc[i, "Gene_name"]
    df_gene_refseq_TRA.loc[i,"Refseq"] = str(final_TRAs_split.loc[i, "Tx"]) + "." +\
                                         str(final_TRAs_split.loc[i, "Tx_version"])[0]

df_gene_refseq_TRA_dict = dict(zip(df_gene_refseq_TRA.gene,df_gene_refseq_TRA.Refseq))                                        

# add unknwon to dict
df_gene_refseq_TRA_dict["unknown"] = "-"
                                       
# gt full info
final_TRAs_full = final_TRAs[final_TRAs["Annotation_mode"] == "full"]


# Reset indices
final_TRAs_full = final_TRAs_full.reset_index(drop="TRUE")

# handling NA in gene Name
final_TRAs_full["Gene_name"] = final_TRAs_full["Gene_name"].fillna("unknown")

# END TRA

# prepare INV data - START ----------------------------------------------------
# Reset indices
INV_data = INV_data.reset_index(drop="TRUE")

# remove <INV> values in ALT column
#INV_data_0 = INV_data[INV_data["ALT"] != "<INV>"]

# get full INV data
INV_data_1_full = INV_data[INV_data["Annotation_mode"] == "full"]

# remove nan
INV_data_2_full = INV_data_1_full.dropna(subset=["Gene_name"])

# Reset indices
INV_data_2_full = INV_data_2_full.reset_index(drop="TRUE")

# filter gene list
gene_filtered_sv_idx = []
for i, genes in enumerate(INV_data_2_full["Gene_name"]):
    if len(genes.split(";")) > 1:
        start_gene = genes.split(";")[0]
        end_gene = genes.split(";")[-1]
        if start_gene in gene_lst or end_gene in gene_lst:
            gene_filtered_sv_idx.append(i)
    elif len(genes.split(";")) == 1:
        if genes in gene_lst:
            gene_filtered_sv_idx.append(i)

INV_data_3 = INV_data_2_full.loc[gene_filtered_sv_idx,:]

# get ids and mate ids
# Reset indices
INV_data_3 = INV_data_3.reset_index(drop="TRUE")
INV_mate_ids_lst = []
for i, mate_ids in enumerate(INV_data_3["INFO"]):
    mate_id = "MantaBND:" + "".join((re.findall(r'[\d+\:]', mate_ids.split(";")[1]))[1:])
    INV_mate_ids_lst.append(mate_id)
    
# make id list 
final_INV_ids = list(INV_data_3["ID"]) + INV_mate_ids_lst

# get data from INV all
final_INVs = INV_data[INV_data["ID"].isin(final_INV_ids) == True]

# get RefSeq from split
final_INVs_split = final_INVs[final_INVs["Annotation_mode"] == "split"]

# make dataframe !!!check Columns before submittung!!!
df_gene_refseq_INV = pd.DataFrame()

header_col = ["gene", "Refseq"]

for col in header_col:
    df_gene_refseq_INV[col] = ""
       
# Reset indices
final_INVs_split = final_INVs_split.reset_index(drop="TRUE")
    
for i in range(len(final_INVs_split["Gene_name"])):
    df_gene_refseq_INV.loc[i,"gene"] = final_INVs_split.loc[i, "Gene_name"]
    df_gene_refseq_INV.loc[i,"Refseq"] = str(final_INVs_split.loc[i, "Tx"]) + "." +\
                                         str(final_INVs_split.loc[i, "Tx_version"])[0]

df_gene_refseq_INV_dict = dict(zip\
                         (df_gene_refseq_INV.gene,df_gene_refseq_INV.Refseq))
                                         
# get full info
final_INVs_full = final_INVs[final_INVs["Annotation_mode"] == "full"]

# remove <TRA> values in ALT column
final_INVs_full = final_INVs_full [final_INVs_full["ALT"] != "<INV>"]

# Reset indices
final_INVs_full = final_INVs_full.reset_index(drop="TRUE")

# END INV

# Start other data (del/dup/ins)
# get full info
others_data_full_0 = others_data[others_data["Annotation_mode"] == "full"]

# remove nan
others_data_full = others_data_full_0.dropna(subset=["Gene_name"])

# Reset indices
others_data_full = others_data_full.reset_index(drop="TRUE")

# filter gene list
gene_filtered_sv_idx = []
for i, genes in enumerate(others_data_full["Gene_name"]):
    if len(genes.split(";")) > 1:
        start_gene = genes.split(";")[0]
        end_gene = genes.split(";")[-1]
        if start_gene in gene_lst or end_gene in gene_lst:
            gene_filtered_sv_idx.append(i)
    elif len(genes.split(";")) == 1:
        if genes in gene_lst:
            gene_filtered_sv_idx.append(i)

filtered_others_data_full = others_data_full.loc[gene_filtered_sv_idx,:]

# Reset indices
filtered_others_data_full = filtered_others_data_full.reset_index(drop="TRUE")

# get RefSeq from split
others_data_split = others_data[others_data["Annotation_mode"] == "split"]

# make dataframe !!!check Columns before submittung!!!
df_gene_refseq_others = pd.DataFrame()

header_col = ["gene", "Refseq"]

for col in header_col:
    df_gene_refseq_others[col] = ""
        
# Reset indices
others_data_split = others_data_split.reset_index(drop="TRUE")
    
for i in range(len(others_data_split["Gene_name"])):
    df_gene_refseq_others.loc[i,"gene"] = others_data_split.loc[i, "Gene_name"]
    df_gene_refseq_others.loc[i,"Refseq"] = str(others_data_split.loc[i, "Tx"]) + "." +\
                                         str(others_data_split.loc[i, "Tx_version"])[0]

df_gene_refseq_others_dict = dict(zip\
                        (df_gene_refseq_others.gene,df_gene_refseq_others.Refseq))

# get sequence ontololgy from vep113
VEP_OUT = args.vep_tab
VEP_TABLE_OUT_data = pd.read_csv(VEP_OUT, usecols=["#Uploaded_variation",
                                                 "Consequence"], 
                                 sep="\t",low_memory=False)

VEP_TABLE_OUT_data = VEP_TABLE_OUT_data.rename\
                    (columns={"#Uploaded_variation": "Manta_ID"})

df_vep_SO_dict = dict(zip\
                        (VEP_TABLE_OUT_data.Manta_ID,VEP_TABLE_OUT_data.Consequence))
       
# collect data
df_svs_TRA = pd.DataFrame()
df_svs_INV = pd.DataFrame()
df_svs_others = pd.DataFrame()

header_col = ["Standort", "Sample", "Start-Gen1", \
              "Start-Chr", "Start-Position", "Start-Transkript1", \
              "End-Gen2", "End-Chr", "End-Position", "End-Transkript2", \
              "Konsequenz", "Klasse", "Länge in Basenpaaren"]

for col in header_col:
    df_svs_TRA[col] = ""
    
for col in header_col:
    df_svs_INV[col] = ""
    
for col in header_col:
    df_svs_others[col] = ""    

# others rename sv_type TRA
for i in range(len(final_TRAs_full["SV_type"])):
    if final_TRAs_full.loc[i,"SV_type"] == "TRA":
        final_TRAs_full.loc[i,"SV_type"] = "Fusion"
        
# merge mates
mate1 = []
mate2 = []
for i in range(len(final_TRAs_full["ID"])):
    if int(final_TRAs_full.loc[i,"ID"].split(":")[-1]) == 0:
        mate1.append(i)
    elif int(final_TRAs_full.loc[i,"ID"].split(":")[-1]) == 1:
        mate2.append(i)
        
df_mate1 = final_TRAs_full.loc[mate1,:]
df_mate2 = final_TRAs_full.loc[mate2,:]

# Reset indices
df_mate1  = df_mate1 .reset_index(drop="TRUE")
df_mate2  = df_mate2 .reset_index(drop="TRUE")  

df_mate1["ID_merge"] = ""
df_mate2["ID_merge"] = ""
for i in range(len(df_mate1["ID"])):
    df_mate1.loc[i,"ID_merge"] = ":".join(df_mate1.loc[i,"ID"].split(":")[:-1])
    
for i in range(len(df_mate2["ID"])):
    df_mate2.loc[i,"ID_merge"] = ":".join(df_mate2.loc[i,"ID"].split(":")[:-1])
    
df_tra_merge = pd.merge(df_mate1,df_mate2,
                        left_on = ["ID_merge"],
                        right_on = ["ID_merge"],
                        how = "left")

# Reset indices
df_tra_merge = df_tra_merge.reset_index(drop="TRUE")

# add translocations
for i in range(len(df_tra_merge["ID_x"])):
    df_svs_TRA.loc[i,"Standort"] = "Bonn_Pathologie"
    df_svs_TRA.loc[i,"Sample"] = args.sample
    df_svs_TRA.loc[i,"Start-Gen1"] = df_tra_merge.loc[i,"Gene_name_x"].split(";")[0]
    df_svs_TRA.loc[i,"Start-Chr"] = "chr" + df_tra_merge.loc[i,"SV_chrom_x"]
    df_svs_TRA.loc[i,"Start-Position"] = df_tra_merge.loc[i,"SV_start_x"]
    df_svs_TRA.loc[i,"Start-Transkript1"] = df_gene_refseq_TRA_dict[df_tra_merge.loc\
                                        [i,"Gene_name_x"].split(";")[0]]
    df_svs_TRA.loc[i,"End-Gen2"] = df_tra_merge.loc[i,"Gene_name_y"].split(";")[-1]
    df_svs_TRA.loc[i,"End-Chr"] = "chr" + df_tra_merge.loc[i,"SV_chrom_y"]
    df_svs_TRA.loc[i,"End-Position"] = df_tra_merge.loc[i,"SV_end_y"]
    df_svs_TRA.loc[i,"End-Transkript2"] = df_gene_refseq_TRA_dict[df_tra_merge.loc\
                                        [i,"Gene_name_y"].split(";")[-1]]
    df_svs_TRA.loc[i,"Konsequenz"] = df_vep_SO_dict\
                                                 [df_tra_merge.loc[i,"ID_x"]] + "," + \
                                     df_vep_SO_dict\
                                                 [df_tra_merge.loc[i,"ID_y"]]
    df_svs_TRA.loc[i,"Klasse"] = df_tra_merge.loc[i,"SV_type_y"]
    df_svs_TRA.loc[i,"Länge in Basenpaaren"] = df_tra_merge.loc[i,"SV_length_x"]
    
# others rename sv_type others (del/dup/ins)
for i in range(len(final_INVs_full["SV_type"])):
    if final_INVs_full.loc[i,"SV_type"] == "INV":
        final_INVs_full.loc[i,"SV_type"] = "Inversion"

# add Inversions
for i in range(len(final_INVs_full["ID"])):
    df_svs_INV.loc[i,"Standort"] = "Bonn_Pathologie"
    df_svs_INV.loc[i,"Sample"] = args.sample
    df_svs_INV.loc[i,"Start-Gen1"] = final_INVs_full.loc[i,"Gene_name"].split(";")[0]
    df_svs_INV.loc[i,"Start-Chr"] = "chr" + final_INVs_full.loc[i,"SV_chrom"]
    df_svs_INV.loc[i,"Start-Position"] = final_INVs_full.loc[i,"SV_start"]
    df_svs_INV.loc[i,"Start-Transkript1"] = df_gene_refseq_INV_dict[final_INVs_full.loc\
                                        [i,"Gene_name"].split(";")[0]]
    df_svs_INV.loc[i,"End-Gen2"] = final_INVs_full.loc[i,"Gene_name"].split(";")[-1]
    df_svs_INV.loc[i,"End-Chr"] = "chr" + final_INVs_full.loc[i,"SV_chrom"]
    df_svs_INV.loc[i,"End-Position"] = final_INVs_full.loc[i,"SV_end"]
    df_svs_INV.loc[i,"End-Transkript2"] = df_gene_refseq_INV_dict[final_INVs_full.loc\
                                        [i,"Gene_name"].split(";")[-1]]
    if final_INVs_full.loc[i,"ID"] in df_vep_SO_dict.keys():
        df_svs_INV.loc[i,"Konsequenz"] = df_vep_SO_dict\
                                                 [final_INVs_full.loc[i,"ID"]]
    else:
        df_svs_INV.loc[i,"Konsequenz"] = "unknown"
    df_svs_INV.loc[i,"Klasse"] = final_INVs_full.loc[i,"SV_type"]
    df_svs_INV.loc[i,"Länge in Basenpaaren"] = final_INVs_full.loc[i,"SV_length"]
    
# others rename sv_type others (del/dup/ins)
for i in range(len(filtered_others_data_full["SV_type"])):
    if filtered_others_data_full.loc[i,"SV_type"] == "DEL":
        filtered_others_data_full.loc[i,"SV_type"] = "Deletion"
    elif filtered_others_data_full.loc[i,"SV_type"] == "DUP":
        filtered_others_data_full.loc[i,"SV_type"] = "Duplikation"
    elif filtered_others_data_full.loc[i,"SV_type"] == "INS":
        filtered_others_data_full.loc[i,"SV_type"] = "Insertion"
        
# add others (del/dup/ins)
for i in range(len(filtered_others_data_full["ID"])):
    df_svs_others.loc[i,"Standort"] = "Bonn_Pathologie"
    df_svs_others.loc[i,"Sample"] = args.sample
    df_svs_others.loc[i,"Start-Gen1"] = filtered_others_data_full.loc[i,"Gene_name"].split(";")[0]
    df_svs_others.loc[i,"Start-Chr"] = "chr" + filtered_others_data_full.loc[i,"SV_chrom"]
    df_svs_others.loc[i,"Start-Position"] = filtered_others_data_full.loc[i,"SV_start"]
    df_svs_others.loc[i,"Start-Transkript1"] = df_gene_refseq_others_dict\
                                                     [filtered_others_data_full.loc\
                                                 [i,"Gene_name"].split(";")[0]]
    df_svs_others.loc[i,"End-Gen2"] = filtered_others_data_full.loc[i,"Gene_name"].split\
                                                                      (";")[-1]
    df_svs_others.loc[i,"End-Chr"] = "chr" + filtered_others_data_full.loc[i,"SV_chrom"]
    df_svs_others.loc[i,"End-Position"] = filtered_others_data_full.loc[i,"SV_end"]
    df_svs_others.loc[i,"End-Transkript2"] = df_gene_refseq_others_dict\
                                                          [filtered_others_data_full.loc\
                                                [i,"Gene_name"].split(";")[-1]]
    df_svs_others.loc[i,"Konsequenz"] = df_vep_SO_dict\
                                                 [filtered_others_data_full.loc[i,"ID"]]
    df_svs_others.loc[i,"Klasse"] = filtered_others_data_full.loc[i,"SV_type"]
    df_svs_others.loc[i,"Länge in Basenpaaren"] = filtered_others_data_full.loc[i,"SV_length"]
    

# format concatenated output
frames = [df_svs_TRA, df_svs_INV, df_svs_others]
sample_output = pd.concat(frames)
sample_output["Länge in Basenpaaren"] = sample_output["Länge in Basenpaaren"].apply\
                                        (abs).astype(int)

sample_output.to_csv(args.outfile, sep=",", index=False)


