import argparse
import boto3
import botocore
from pathlib import Path
import tarfile

from botocore import UNSIGNED
from botocore.client import Config


import urllib.request

files = '''
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_alleles_hg38.zip
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_loci_hg38.zip
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/GC_G1000_hg38.zip
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/RT_G1000_hg38.zip
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index/
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/dragmap/
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes/
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/Length/Homo_sapiens_assembly38.len
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz.tbi
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz.tbi
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.known_indels.vcf.gz
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/beta/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/Control-FREEC/out100m2_hg38.gem
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi
'''


sra_testdata = [
	"https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7890824/SRR7890824",
	"https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7890827/SRR7890827"
    ]

vep_cache_version = '113'

#vep_cache = 'https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz'
vep_cache = 'https://ftp.ensembl.org/pub/release-{vep_cache_version}/variation/indexed_vep_cache/homo_sapiens_vep_{vep_cache_version}_GRCh38.tar.gz'

bucket_name = 'ngi-igenomes' # replace with your bucket name
region = 'eu-west-1'
s3 = boto3.client('s3', region_name=region, config=Config(signature_version=UNSIGNED))

def download_file(obj_name, output_dir):
    localf = str(Path(output_dir) / obj_name.lstrip('igenomes/'))
    outdir = Path(localf).absolute().parent
    outdir.mkdir(exist_ok=True,parents=True)
    if Path(localf).exists():
        return
    s3.download_file(bucket_name, obj_name, localf)


def download_recursive(obj_name, output_dir):
    if obj_name.endswith('/'):
        print(f'recursing {obj_name}')
        # note, list_objects only shows the first 1000
        obs = s3.list_objects(Bucket=bucket_name, Prefix=obj_name)
        print(f'got obs: {obs}')
        ks = [x['Key'] for x in obs['Contents']]
        print(f'got ks: {ks}')
        for o in ks:
            download_recursive(o, output_dir)

    else:
        download_file(obj_name, output_dir)



def download_igenomes(output_dir):
    print(f"Beginning downloading igenomes to output directory: {output_dir}")

    filel = files.strip().split('\n')
    print(f"Downloading files: {filel}")

    for f in filel:
        print(f)
        download_recursive(f, output_dir)


def download_testdata_sra(output_dir):
    for srafile in sra_testdata:
        srafile_name = srafile.split('/')[-1]
        urllib.request.urlretrieve(srafile, str(output_dir) + '/testdata/' + srafile_name + '.sra')


def untar_file(tar_file, directory):
    with tarfile.open(tar_file, 'r') as f:
        f.extractall(directory)

def download_vep_cache(output_dir):
    outpath = (Path(output_dir) / 'vep_cache' / vep_cache_version / 'indexed')
    outpath.mkdir(parents=True, exist_ok=True)
    cachename = vep_cache.split('/')[-1]
    outfile = str(outpath / cachename)
    urllib.request.urlretrieve(vep_cache, outfile)


def cli(output_dir, dataset):
    if dataset == 'reference':
        download_igenomes(output_dir)
        download_vep_cache(output_dir)
    elif dataset == 'test':
        download_testdata_sra(output_dir)
    else:
        raise RuntimeError('invalid dataset')

    
def main():
    parser = argparse.ArgumentParser(description="helper tool to download the pipelines reference or test or example")
    parser.add_argument('dataset', choices=['reference', 'test', 'example'], help="which data to download")
    parser.add_argument("output_dir", help="The directory where output will be saved")
    args = parser.parse_args()
    cli(args.output_dir)

if __name__ == "__main__":
    main()
