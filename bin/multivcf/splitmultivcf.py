import pysam
from tqdm import tqdm
import sys


invcf = sys.argv[1]

outfolder = sys.argv[2]
def get_header(h):
    hh = pysam.VariantHeader()
    hh.add_meta('fileformat', value='VCFv4.3')

    for _ in h.contigs.values():
        hh.add_meta('contig', items=[('ID',_.name)])
    
    for _ in h.formats.values():
        hh.add_meta('FORMAT', items=[('ID',_.name), ('Number',_.number), ('Type',_.type),('Description',_.description)])
    for _ in h.info.values():
        hh.add_meta('INFO', items=[('ID',_.name), ('Number',_.number), ('Type',_.type),('Description',_.description)])
    
    for _ in h.filters.values():
        hh.add_meta('FILTERS', items=[('ID',_.name), ('Number',_.number), ('Type',_.type),('Description',_.description)])
    
    hh.add_line('##vcfsource=conversion from idat to gtc to vcf, vcf split from multisample to single sample ShD')
    
    return hh

with pysam.VariantFile(invcf) as v:
    for r in v:
        s = r.samples.keys()
        h=v.header
        break
        
        
with pysam.VariantFile(invcf) as v:
    for r in v:
        s = r.samples.keys()
        break


        
openf =[]
for samples in s:
    header = get_header(h)
    header.add_sample(samples)
    o = pysam.VariantFile(f'{outfolder}/{samples}.vcf','w',header=header)
    openf.append(o)


with pysam.VariantFile(invcf) as v:
    for r in tqdm(v.fetch()):
        for o,samples in zip(openf,s):
            rr = o.new_record(contig=r.chrom, start=r.start, stop=r.stop,alleles=(r.alleles),info=r.info)
            
            for f in r.format.keys():
                rr.samples[samples][f]=r.samples[samples].get(f)
            
            o.write(rr)
for o in openf:
    o.close()