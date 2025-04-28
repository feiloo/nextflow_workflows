#!/usr/bin/env python

import pandas as pd
import argparse
import variantlist_utils as vu
from transcript_list_08_08_2024 import transcript_list
from datetime import datetime

# Using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--clc", type=str)
parser.add_argument("-v", "--vep", type=str)
parser.add_argument("-t", "--transcripts", type=str)
parser.add_argument("-e", "--encoding", type=str)
parser.add_argument("-D", "--variant_DBi", type=str)
parser.add_argument("-o", "--outfile", type=str)
parser.add_argument("-rv", "--removed_variants", type=str)
args = parser.parse_args()

# Getting current date and time for log information
date_time_now = datetime.now()

# dd/mm/YY H:M:S
dt_string = date_time_now.strftime("%d/%m/%Y %H:%M:%S")
print("Start:", dt_string)

# Processed file
print("Input_CLC_file:", args.clc)
print("Input_VEP_file:", args.vep)
print("Output_file_1:", args.outfile)
print("Output_file_2:", args.removed_variants)

# Script used
print("Script: organize_variant_table.py")

# Get CLC_PAN_data
clc_PAN_file = args.clc
CLC_variant_track_data_PAN = pd.read_csv(clc_PAN_file, delimiter=",", encoding=args.encoding)

# CLC_PAN_data - adjust region_position
CLC_variant_track_data_PAN = vu.adjust_region_position(CLC_variant_track_data_PAN)

# Filtering of variants not equal to "Reference allele"
# Rename column "Reference allele" to "ReferenceAllele"
CLC_variant_track_data_PAN = CLC_variant_track_data_PAN.rename\
                              (columns={"Reference allele": "ReferenceAllele"})
group_reference_allele = CLC_variant_track_data_PAN.groupby\
                                   (CLC_variant_track_data_PAN.ReferenceAllele)
clc_data_ReferenceAllele_NO =  group_reference_allele.get_group("No")

# Column list to filter
columns_to_filter = ["Frequency", "QUAL", "Forward/reverse balance", \
                     "Average quality",\
                     "Read position test probability", \
                     "Read direction test probability"]

def partition_dataframe(dataframe: pd.DataFrame, condition: pd.Series):
    ''' returns dataframe rows that match and that dont match the condition '''
    return dataframe[condition], dataframe[~condition]

def filter_and_collect(dataframe: pd.DataFrame, collector: list, condition: pd.Series) -> tuple[pd.DataFrame, list[pd.DataFrame]]:
    ''' filters out rows of a dataframe and collects the discarded rows '''
    positives, negatives = partition_dataframe(dataframe, condition)
    collector.append(negatives)
    return positives, collector

# Convert string number in float of columns_to_filter
clc_data_ReferenceAllele_NO = vu.convert_coltype_str_to_float\
    (columns_to_filter, clc_data_ReferenceAllele_NO)


# collect filtered out variant rows 
discarded_clc_data : list[pd.DataFrame] = []

# Filter data with QUAL >= 150
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_ReferenceAllele_NO, discarded_clc_data,
    clc_data_ReferenceAllele_NO["QUAL"] >= 150)


# Filter data with Av quality >= 35
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_filtered, discarded_clc_data,
    clc_data_filtered["Average quality"] >= 35)


# Filter data with count >= 2
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_filtered, discarded_clc_data,
    clc_data_filtered["Count"] >= 2)

# Filter data with Frequency >= 5
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_filtered, discarded_clc_data,
    clc_data_filtered["Frequency"] >= 5)


# Filter data with Forward/reverse balance > 0
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_filtered, discarded_clc_data,
    clc_data_filtered["Forward/reverse balance"] > 0)


# Filter data with Read position test probability >= 0.000001
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_filtered, discarded_clc_data,
    clc_data_filtered["Read position test probability"] >= 0.000001)


# Filter data with Read direction test probability >= 0.000001
clc_data_filtered, discarded_clc_data = filter_and_collect(
    clc_data_filtered, discarded_clc_data,
    clc_data_filtered["Read direction test probability"] >= 0.000001)


# filter non-synonymous???

# Generate xlsx file of removed_variants
# concatenate discarded/removed clc_data variants
removed_variants = pd.concat(discarded_clc_data)

# Sort Frequency of removed variants
process_removed_variants = removed_variants.sort_values(by=["Frequency"], ascending=False)

# Extract Deletions, Inserstion and Replacements with AF > 5 and TERT KRAS Variants
biotypes = ["Insertion", "Deletion", "Replacement"]
get_back_biotypes = process_removed_variants[(process_removed_variants\
                                             ["Type"].isin(biotypes)) & \
                                             (process_removed_variants\
                                             ["Frequency"] >= 5) == True]

discareded_remaining_biotypes = process_removed_variants[(process_removed_variants\
                                             ["Type"].isin(biotypes)) & \
                                             (process_removed_variants\
                                             ["Frequency"] >= 5) == False]

extract_genes = ["TERT", "KRAS"]
get_back_genes = process_removed_variants[process_removed_variants\
["gene (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"].isin(extract_genes)]

final_removed_variants = discareded_remaining_biotypes[~discareded_remaining_biotypes\
["gene (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"].isin(extract_genes)]

# write file
final_removed_variants.to_excel(args.removed_variants, \
                            index = False, \
                             engine= None)

# Adding get back biotypes and genes to data
data_to_analyse = [clc_data_filtered, get_back_biotypes, get_back_genes]

# concatenate clc_data variants
clc_data_filtered = pd.concat(data_to_analyse)

# Reindex clc data for further processing necessary
clc_data_filtered = clc_data_filtered.reset_index()

# "Explode" column "coding region change"
# "Explode" necessary since each row contains several transcript:HGVS values
# Result: presentation of one transcript:HGVS per row; a variant shown in a
# variety of rows derived from the total transcript number
clc_data_filtered = clc_data_filtered.rename(columns=\
           {"Coding region change": "Coding_region_change"})

clc_data_filtered = (
    clc_data_filtered.assign(Coding_region_change=clc_data_filtered\
                                  ["Coding_region_change"].str.split("; "))
        .explode("Coding_region_change")
        .reset_index(drop=True)
)

# Get NM/XM and HGVSc derived from transcript:HGVS in "Coding_region_change"
# NM/XM and HGVSc are saved in new columns ("NM_v" and "HGVSc)
for index, NM_tr in enumerate(clc_data_filtered["Coding_region_change"]):

    if pd.isna(NM_tr):
        clc_data_filtered.loc[index, "NM_v"] = clc_data_filtered.loc\
                                               [index, "NM_v"]
        clc_data_filtered.loc[index, "HGVSc"] = clc_data_filtered.loc\
                                                [index, "HGVSc"]

    else:
        clc_data_filtered.loc[index, "NM_v"] = NM_tr.split(":")[0]
        clc_data_filtered.loc[index, "HGVSc"] = NM_tr.split(":")[1]

# Drop Na values from DataFrame (necessary?)
clc_data_filtered_dropNa = clc_data_filtered[clc_data_filtered\
                          ["NM_v"].str.contains("nan") == False]

# Load PANCANCER RefSeq transcripts to list
RefSeq_NM = transcript_list


# Reset indices; necessary for further processing; since index could be not
# in order due to step "Drop Na values from DataFrame" (see above) (necessary?)
clc_data_filtered_dropNa = clc_data_filtered_dropNa.reset_index()

# convert to str (necessary?)
clc_data_filtered_dropNa["NM_v"] = clc_data_filtered_dropNa["NM_v"].astype(str)

# Get index of rows presented in PANCANCER RefSeq list
NM_idx = []
for i in range(len(clc_data_filtered_dropNa["NM_v"])):
    if clc_data_filtered_dropNa["NM_v"][i].split(".")[0] in RefSeq_NM:
        NM_idx.append(i)

# Filter variants accoriing to NM_idx list
# Store result in new variable
pre_final_data = clc_data_filtered_dropNa.loc[NM_idx, :]

# Convert NM + version in NM only for subsequent merging step
# Store NM only in column NM_merge
index_variants = pre_final_data.index.tolist()
for i in range(len(pre_final_data["NM_v"])):
    tmp_NM_name = pre_final_data["NM_v"][index_variants[i]].split(".")
    pre_final_data.loc[[index_variants[i]], "NM_merge"] =  tmp_NM_name[0]

# Generate clinvar Link and add to data
# Add column "Clinvar_Link" and generate link
pre_final_data["Clinvar_Link"] = ""
rs = pre_final_data["name dbsnp_v151_ensembl_hg38_no_alt_analysis_set"]
rs_idx = rs.index.tolist()
link_lst = []

# Get clinvar https list
pre_final_data["Clinvar_Link"] = vu.clinvar_link_list(link_lst,\
                                 rs_idx, rs)

# Log information
print("Process information:")
print("--> Processing CLC_PanCancer data: successful!")

# Process VEP data (obtained via ensembl-vep tool)
vep_file = args.vep
VEP_data = pd.read_csv(vep_file, delimiter="\t")

# Add column "Chromosome" and "Position"
# Get chromosome and position from column "Location" and NM values
VEP_data = vu.adjust_chr_pos_NM_VEP(VEP_data)

# get cosmic id from existing variation!!!
# Result: Ready VEP data for merging with CLC_PAN_data
VEP_data["cosmic_ID"] = ""
for ID_str in range(len(VEP_data["Existing_variation"])):
    if "COSV" in VEP_data["Existing_variation"].iloc[ID_str]:
        tmp_ID_str = list(filter(lambda x: "COSV" in x, VEP_data\
                         ["Existing_variation"].iloc[ID_str].split(",")))[0]
        VEP_data.loc[ID_str, "cosmic_ID"] = tmp_ID_str
    else:
        VEP_data.loc[ID_str, "cosmic_ID"] = "-"

# Merge CLC_PAN_data and VEP data dfs
# Prepare CLC_PAN_data and VEP data for merging
VEP_data["Chromosome"] = VEP_data["Chromosome"].astype(str)
pre_final_data["Chromosome"] = pre_final_data["Chromosome"].astype(str)

VEP_data["Position"] = VEP_data["Position"].astype(int)
pre_final_data["Position"] = pre_final_data["Position"].astype(int)

VEP_data["End Position"] = VEP_data["End Position"].astype(int)
pre_final_data["End Position"] = pre_final_data["End Position"].astype(int)

VEP_data = VEP_data.rename\
                        (columns={"Feature": "NM_merge"})

VEP_data["NM_merge"] = VEP_data["NM_merge"].astype(str)
pre_final_data["NM_merge"] = pre_final_data["NM_merge"].astype(str)

print("--> Processing VEP_Ensembl data: successful!")

# Merge CLC_PAN_data with VEP via columns "Chromosome", "Position", "End Position", "Allele", "NM_merge"
# and add aditional columns from vep
merged = pd.merge(pre_final_data, VEP_data[["Chromosome", "Position", "End Position", \
                 "Allele", "NM_merge", "Feature_type", "Consequence", "HGVSc", "HGVSp", \
                 "SYMBOL", "BIOTYPE", "EXON", "AF", "MAX_AF", \
                 "gnomADe_AF", "gnomADg_AF", "SIFT", "PolyPhen", "CLIN_SIG", "cosmic_ID"]], \
                 on = ["Chromosome", "Position", "End Position", "Allele", "NM_merge"], \
                 how = "left")

# Get HGVSp nomenclature RefSeq
index_variants_merged = merged.index.tolist()
for i in range(len(merged["HGVSp"])):
    if pd.isna(merged["HGVSp"][index_variants_merged[i]]):
        merged.loc[[index_variants_merged[i]], "HGVS_PROTEIN"] = merged\
            ["HGVSp"][index_variants_merged[i]]
    else:
        tmp_p_name = merged["HGVSp"][index_variants_merged[i]].split(":")
        if len(tmp_p_name) >= 2:
            merged.loc[[index_variants_merged[i]], "HGVS_PROTEIN"] =  tmp_p_name[1]

# Merge/join internal variantDB (variantDBi) PAN_CANCER_DATA
# Change to current "Variantenliste" if needed
variantDBi = pd.read_excel(args.variant_DBi)

merged = pd.merge(merged,\
                  variantDBi,\
                  left_on = ["name dbsnp_v151_ensembl_hg38_no_alt_analysis_set"],\
                  right_on = ["name dbsnp_v151_ensembl_hg38_no_alt_analysis_set"],\
                  how = "left")

# getting fields {coding,inRegulatoryElements,
# notInCodingAndNotInRegulatoryElements} in field localization
# required for the MV Oncology report
merged["localization"] = "-"

# field: coding (bool: biotype = protein_coding and HGVSp = xxx.p.xxx / exists)
# Reset indices
merged = merged.reset_index(drop="TRUE")
for i in range(len(merged["BIOTYPE"])):
    if merged.loc[i,"BIOTYPE"] == "protein_coding" and pd.isna\
      (merged.loc[i,"HGVS_PROTEIN"]) == False and merged.loc\
                 [i,"HGVS_PROTEIN"].startswith("p."):
                     merged.loc[i,"localization"] = "coding"

# fields {inRegulatoryElements,notInCodingAndNotInRegulatoryElements} are
# preliminary not catched here in PanCancer, but must be added in WES/WGS

# get field HGNC from clc field
merged = merged.reset_index(drop="TRUE")
merged["HGNC_MV"] = "-"
for hgnc_id in range(len(merged["HGNC (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"])):
    if len(merged.loc[hgnc_id, "HGNC (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"].split(":")) == 3:
        tmp_hgnc_0 = merged.loc[hgnc_id, "HGNC (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"]
        tmp_hgnc_1 =  tmp_hgnc_0.split(":")[1] + ":" +  tmp_hgnc_0.split(":")[2]
        merged.loc[hgnc_id, "HGNC_MV"] = tmp_hgnc_1
    elif len(merged.loc[hgnc_id, "HGNC (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"].split(",")) == 2:
        tmp_hgnc_2 = merged.loc[hgnc_id, "HGNC (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"].split(",")[0]
        tmp_hgnc_2_1 = tmp_hgnc_2.split(":")[1] + ":" + tmp_hgnc_2.split(":")[2]
        tmp_hgnc_3 = merged.loc[hgnc_id, "HGNC (Homo_sapiens_refseq_GRCh38.p14_no_alt_analysis_set_Genes)"].split(",")[1]
        tmp_hgnc_3_1 = tmp_hgnc_3.split(":")[1] + ":" + tmp_hgnc_3.split(":")[2]
        merged.loc[hgnc_id, "HGNC_MV"] = tmp_hgnc_2_1 + ", " + tmp_hgnc_3_1

# Get comprehensive output format
processed_data_final = merged[["Chromosome", "Position", "End Position", \
                       "Reference", "Allele", "Count", "Coverage", \
                        "Frequency", "QUAL", "Forward/reverse balance", \
                        "Average quality","Read position test probability", \
                        "Read direction test probability", "BaseQRankSum",
                        "Homopolymer", "Homopolymer length", \
                        "Homopolymer region", "Repeat region", \
                        "Count (singleton UMI)", "Count (big UMI)", \
                        "Proportion (singleton UMIs)", \
                        "localization", "HGNC_MV", "Feature_type", \
                        "Consequence", "SYMBOL", \
                        "NM_v", "HGVSc_x", "HGVS_PROTEIN", "Non-synonymous", \
                        "EXON", \
                        "func dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "CLNSIG clinvar_20220730_hg38_no_alt_analysis_set", \
                        "CLIN_SIG", \
                        "CLNREVSTAT clinvar_20220730_hg38_no_alt_analysis_set", \
                        "alleles dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "alleleFreqs dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "AF_EXAC clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF_TGP clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF", "MAX_AF", "gnomADe_AF", "gnomADg_AF", "SIFT", \
                        "PolyPhen", "cosmic_ID", "Wertung", "Clinvar_Link"]]

# Round AF to max 2 decimals
processed_data_final.loc[:,"Frequency"]  = processed_data_final\
                                   ["Frequency"].apply(lambda x: round(x,2))

# Sort Frequency
final_variants = processed_data_final.sort_values(by=["Frequency"], ascending=False)

# Rename columns (final names not yet final known from DNPM)
final_variants = final_variants.rename(columns={"Chromosome": "Chromosom",
                                                "SYMBOL": "Gen",
                                                "NM_v": "Transcript_ID",
                                                "HGVSc_x": "cDNA_Change",
                                                "HGVS_PROTEIN": "Amino_Acid_Change",
                                                "EXON": "Exon"})

# Add "submit" and "Interpretation" column
final_variants.insert(loc=43, column="submit", value="-")
final_variants.insert(loc=44, column="Interpretation", value="-")

# Log information
print("--> Combining and formatting CLC and VEP data: successful!")

# Save file
final_variants.to_excel(args.outfile, \
                            index = False, \
                             engine= None)
# Log information
print("--> Writing data to xlsx output files: successful!")

# Getting current date and time for log information
date_time_now = datetime.now()
# dd/mm/YY H:M:S
dt_string = date_time_now.strftime("%d/%m/%Y %H:%M:%S")
print("End:", dt_string)






































