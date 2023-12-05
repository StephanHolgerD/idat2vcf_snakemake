import pandas as pd
import os
from glob import glob
import yaml


with open('config.yaml') as s:
    d=yaml.safe_load(s)



ref=d['ref']
idatfolder=d['idatfolder']
clusterfile_egt=d['clusterfile_egt']
manifest_bpm=d['manifest_bpm']

copy_rule = []
samplesheet = d.get('samplesheet')
if samplesheet!=None:
    copy_rule = ['../05_vcf/copydone.txt']

print(idatfolder)
run_name = idatfolder.rstrip('/').split('/')[-1]
print(run_name)
print(ref)

bcftools='bin/bcftools/bcftools'
#gtc_out = 

wkd= os.getcwd()
outdir = '/'.join(wkd.split('/')[:-1])

def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

#for a,b in zip(SAMPLES,SAMPLES_PATH):
 #   a=str(a)
 #   os.makedirs("../01_bedgraphs/"+a, exist_ok=True)
 #   check_symlink(b,f'../01_bedgraphs/{a}/{a}.bedgraph')



rule all:
    input:
        '../04_singlesample_vcf/done.txt',
        '../04_singlesample_vcf/done2.txt',
        copy_rule


       # expand('../02_FilteredBedgraphs/{sample}/{sample}_min{filter}.bedgraph',sample=SAMPLES,filter=minAlignments),
       # expand('{outdir}/03_metilene_call_{filter}/dmr_metilene_qval.0.05.bed',filter=minAlignments,outdir=outdir),
       # expand('{outdir}/03_bat_summarize_{filter}/BAT_summarize_summary_case_control.bedgraph.gz',filter=minAlignments,outdir=outdir)
        #f"{outdir}/03_bat/out.txt"
        #'../03_ConcatFilteredBedgraphs/ConcatFilteredBedgraphs.bed',
        #expand('../04_FilterSplitUnionBed/{sample}/{sample}.bed',sample=SAMPLES)


rule CreateGTC:
    input:
        idatfolder = expand({idatfolder},idatfolder=idatfolder) 

    output:
        out = directory(f'../02_gtc_out/{run_name}') 
    threads: 1
    shell:
        "bin/iaap-cli/iaap-cli/iaap-cli gencall {manifest_bpm} {clusterfile_egt} {output.out} --idat-folder {input.idatfolder}  --output-gtc  --gender-estimate-call-rate-threshold -0.1"




rule GTC_2_BCF_VCF:
    input:
        infolder = f'../02_gtc_out/{run_name}'
    
    output:
        multibcf = f'../03_vcf/{run_name}.bcf',
        multivcf = f'../03_vcf/{run_name}.vcf.gz'

    conda:
        'envs/pysam.yaml'
    threads: 12

    shell:
        'export BCFTOOLS_PLUGINS="bin/bcftools/";\
         {bcftools} +gtc2vcf --no-version -Ou --bpm {manifest_bpm} --egt {clusterfile_egt}  --extra test.tsv --gtcs  {input.infolder} --fasta-ref {ref} |  {bcftools} sort -Ou -T ./bcftools. | {bcftools} norm --no-version -Ob -c x -f {ref} > {output.multibcf};\
         {bcftools} view {output.multibcf} | bgzip -c > {output.multivcf}'


rule SplitVFC1:
    input:
        multivcf = f'../03_vcf/{run_name}.vcf.gz'
    threads: 1

    output:
        fin = '../04_singlesample_vcf/done.txt'
    conda:
        'envs/pysam.yaml'
    shell:
        'for sample in $({bcftools} query -l {input.multivcf}); do echo ${{sample}} ; {bcftools} view -c1 -Oz -s ${{sample}} -o ../04_singlesample_vcf/${{sample}}.vcf.gz {input.multivcf} ; done  ; touch {output.fin}'
    

rule SplitVFC2:
    threads: 1

    input:
        multivcf = f'../03_vcf/{run_name}.vcf.gz'
    output:
        fin = '../04_singlesample_vcf/done2.txt'
    conda:
        'envs/pysam.yaml'
    shell:
        'for sample in $({bcftools} query -l {input.multivcf}); do echo ${{sample}} ; {bcftools} view  -Oz -s ${{sample}} -o ../04_singlesample_vcf/${{sample}}.allarraypositions.vcf.gz {input.multivcf} ; done  ; touch {output.fin}'



rule SamplesheetCopy:
    threads: 1

    input:
        samplesheet = {samplesheet},
        splitdone = '../04_singlesample_vcf/done2.txt'
    output:
        fin = '../05_vcf/copydone.txt'
    conda:
        'envs/pysam.yaml'
    shell:
        'python bin/samplesheetcopy/ssheet_copy.py {input.samplesheet} ../04_singlesample_vcf/ ../05_vcf/ && touch {output.fin}'

