import os
import shutil
import sys


sample_sheet = sys.argv[1] 
raw_folder = sys.argv[2]
output_folder = sys.argv[3]
suff = ['.allarraypositions.vcf.gz','.vcf.gz']


#raw_folder = '../04_singlesample_vcf/'
#output_folder = '../05_final_vcf/'


with open(sample_sheet)as f:
    cont = f.read()
    x =cont.find('[Data]')
    cont=cont[x:]
    x = cont.find('\n')
    cont=cont[x+1:]
    x = cont.find('\n')
    cont=cont[x+1:].rstrip()
    
    cont=cont.split('\n')
    
    s_name = [f"{x.split(';')[0]}_{x.split(';')[-1]}" for x in cont]
    cell_code_name = [f"{x.split(';')[1]}_{x.split(';')[2]}" for x in cont]

    


for suffix in suff:
    for s,ss in zip(s_name,cell_code_name):
        wanted = f'{raw_folder}/{ss}{suffix}'
        outfile=f'{output_folder}/{s}_{ss}{suffix}'
        if os.path.isfile(wanted):
            shutil.copyfile(wanted,outfile)
    