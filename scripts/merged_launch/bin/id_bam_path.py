#!/usr/bin/env python

import os
import pandas as pd
import argparse
import re
from pathlib import Path
import glob

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_dir", type=str)
parser.add_argument("-b", "--id_bam_path", type=str)
parser.add_argument("-p", "--id_purple_path", type=str)
parser.add_argument("-s", "--id_snvs_path", type=str)
args = parser.parse_args()

parent_directory = os.path.abspath(args.target_dir)
target_dirs_wes = [f for f in os.listdir(parent_directory) if re.match(r'wes_', f)]
target_dirs_wgs = [f for f in os.listdir(parent_directory) if re.match(r'wgs_', f)]

if target_dirs_wes == []:
    library_type = "wgs"
    target_dirs = target_dirs_wgs

if target_dirs_wgs == []:
    library_type = "wes"
    target_dirs = target_dirs_wes

# make id_path_df for json generator
size_of_df = len(target_dirs)

id_path_df = pd.DataFrame("", index=list(range(0,size_of_df)),
                          columns = ["patient_id",
                                     "normal_bam",
                                     "normal_bai",
                                     "tumor_bam",
                                     "tumor_bai",
                                     "purple_qc"
                                     "somatic_snvs",
                                     "germline_snvs"])

#print(target_dirs)
for i, target in enumerate(target_dirs):
    print(i, target)

    # get case id and set xlsx_path
    basename = os.path.basename(target)
    patient_id = basename.split(library_type + "_")[-1] # get check !!!
    #print(patient_id)
    id_path_df.loc[i,"patient_id"] = str(patient_id)

    normal_bam = glob.glob(parent_directory + "/" + target + "/alignments/dna/" + "*N*.redux.bam")
    normal_bai = glob.glob(parent_directory + "/" + target + "/alignments/dna/" + "*N*.redux.bam.bai")
    if normal_bam == []:
        normal_bam = glob.glob(parent_directory + "/" + target + "/alignments/dna/" + "*B*.redux.bam")
        normal_bai = glob.glob(parent_directory + "/" + target + "/alignments/dna/" + "*B*.redux.bam.bai")
    tumor_bam = glob.glob(parent_directory + "/" + target + "/alignments/dna/" + "*T*.redux.bam")
    tumor_bai = glob.glob(parent_directory + "/" + target + "/alignments/dna/" + "*T*.redux.bam.bai")
    #print(parent_directory + tumor_bam[0])
    id_path_df.loc[i,"normal_bam"] = Path(normal_bam[0])
    id_path_df.loc[i,"normal_bai"] = Path(normal_bai[0])
    id_path_df.loc[i,"tumor_bam"] = Path(tumor_bam[0])
    id_path_df.loc[i,"tumor_bai"] = Path(tumor_bai[0])

    #somatic_snv = glob.glob(parent_directory + "/" + target + "/purple/" + "*purple.somatic.vcf.gz")
    #somatic_sv = glob.glob(parent_directory + "/" + target + "/purple/" + "*purple.sv.vcf.gz")
    purple_qc = glob.glob(parent_directory + "/" + target + "/purple/" + "*purple.qc")
    #id_path_df.loc[i,"purple_snv"] = Path(somatic_snv[0])
    #id_path_df.loc[i,"purple_sv"] = Path(somatic_sv[0])
    id_path_df.loc[i,"purple_qc"] = Path(purple_qc[0])

    somatic_snvs = glob.glob(parent_directory + "/" + target + "/pave/" + "*pave.somatic.vcf.gz")
    germline_snvs = glob.glob(parent_directory + "/" + target + "/pave/" + "*pave.germline.vcf.gz")
    id_path_df.loc[i,"somatic_snvs"] = Path(somatic_snvs[0])
    id_path_df.loc[i,"germline_snvs"] = Path(germline_snvs[0])

# make samplesheets
bams_df = id_path_df[["patient_id", "normal_bam", "normal_bai", "tumor_bam","tumor_bai"]]
#bams_df = id_path_df[["patient_id","tumor_bam","tumor_bai"]]
bams_df.to_csv(args.id_bam_path, sep=",",index=False)

purple_df = id_path_df[["patient_id","purple_qc"]]
purple_df.to_csv(args.id_purple_path, sep=",",index=False)

cio_vip_somatic_df = id_path_df[["patient_id", "somatic_snvs"]]
cio_vip_somatic_df["patient_id"] = cio_vip_somatic_df["patient_id"].astype(str) + "_Tpavesomatic"
cio_vip_somatic_df = cio_vip_somatic_df.rename(columns={
                                                        "patient_id": "sample",
                                                        "somatic_snvs": "vcf"})
cio_vip_germline_df = id_path_df[["patient_id", "germline_snvs"]]
cio_vip_germline_df["patient_id"] = cio_vip_germline_df["patient_id"].astype(str) + "_Tpavegermline"
cio_vip_germline_df = cio_vip_germline_df.rename(columns={
                                                        "patient_id": "sample",
                                                        "germline_snvs": "vcf"})
cio_vip_frames = [cio_vip_somatic_df, cio_vip_germline_df]
cio_vip_final_df = pd.concat(cio_vip_frames)
cio_vip_final_df  = cio_vip_final_df.reset_index(drop=True)
cio_vip_final_df.to_csv(args.id_snvs_path, sep=",",index=False)
