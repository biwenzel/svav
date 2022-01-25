#strucure variants assembly and validation
import os
import subprocess
import pandas as pd
from clustering_and_creating import write_fasta
from calling_sv import writing_sv
from compare_sv import similar_sv, genotyping
import pysam
from time import time
import argparse

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()

parser.add_argument('sv_call', help='csv file of called structure variants (DEL/INS)')
parser.add_argument('bam', help='bam file corresponding to the variants')
# optional arguments
parser.add_argument('-o', default='svav_run', help='output folder')
# check this again gives out 3 ref right now, as the default is not overwriten
parser.add_argument('-r', action='append', default=['/confidential/tGenVar/ref/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'], help='reference, the first where the sv is called [hg38], the default is on tGenVar')
parser.add_argument('-m', choices=['chr', 'one', 'ref'], default='chr', help='map seperatly to corresponding chromosome[default], map everything together to the reference, map individually to the reference(takes longer)')

args = parser.parse_args()


# create directory and subfolders
path = args.o
directorys = ['','/fasta_files', '/assemble', '/log', '/bam_files', '/result', '/ref','/hap_files','/plot']
for dir in directorys:
    try:
        os.mkdir(path + dir)
    except OSError:
        print ("Creation of the directory %s failed or already exists" % path)


#creating output files with header
h_file = open(path + '/result/primary_call.txt', 'w')
s_file = open(path + '/result/secondary_call.txt', 'w')
header = 'chr,position,type,len,read_name,ref\n'
h_file.write(header)
s_file.write(header)
h_file.close()
s_file.close()

# pre processing sv_call
sv_list = pd.read_csv(args.sv_call)

sv_list['pos1'] = sv_list['pos1'].astype('int')
sv_list['pos2'] = sv_list['pos2'].astype('int')
sv_list['Chromosome'] = sv_list['Chromosome'].astype('str')

idx = True
try:
    sv_list['indicator'] = sv_list['indicator'].astype('str')
except:
    print('no indicator given')
    sv_list['indicator'] = ''


# creating fasta files
# extracting local reads
fasta_time = time()

bam_pb = pysam.AlignmentFile(args.bam, 'rb')
count = 0
for i in range(len(sv_list['pos1'])):
    cont = sv_list['Chromosome'][i]
    start = sv_list['pos1'][i]
    stop = sv_list['pos2'][i]
    idx = sv_list['indicator'][i]

    count += write_fasta(bam_pb, cont, start, stop, path, 100000, idx)


bam_pb.close()
print('completed: writing local fasta files for', count, 'regions from', len(sv_list.pos1))
print(f'it took: {time()-fasta_time:.2f}s')

#assembly
no_asm = 0
asm_time = time()
# second loop for the rest
for i in range(len(sv_list['pos1'])):
    cont = sv_list['Chromosome'][i]
    start = sv_list['pos1'][i]
    stop = sv_list['pos2'][i]
    idx = sv_list['indicator'][i]

    print('_______________________________________________')
    print("looking at region", cont, start, idx)

    ## assembly
    starttime = time()
    # assembly with created fasta files -> check if it is one or two fasta_files
    wtdbg2_command = ['/confidential/FamilyR13/CODE/wtdbg2_2/wtdbg2/wtdbg2.pl', '-t', '16', '-x', 'sq', '-g', str(30000+ stop - start), '--drop-low-cov-edges','-o']  # genome size !not tested!

    if(os.path.isfile(path + '/fasta_files/'+cont+'_'+str(start)+idx+'.fa')):
        ffile = path + '/fasta_files/'+cont+"_"+str(start)+idx+'.fa'
        o_pre = path + '/assemble/asm_'+cont+"_"+str(start)+idx
        log_file = path + '/log/wtdbg2_'+cont+"_"+str(start)+idx+'.log'
        c_file = path + '/assemble/asm_'+cont+"_"+str(start)+idx+'.cns.fa'

        try:
            subprocess.run(wtdbg2_command+[o_pre]+[ffile], check=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
            subprocess.run('sed -i -e "s/ctg/'+cont+'_'+str(start)+idx+'_/g" '+c_file, check=True, shell=True)
        except subprocess.CalledProcessError:
            print('no assembly was possible')
    elif(os.path.isfile(path +'/fasta_files/'+cont+'_'+str(start)+idx+'_h1.fa') and os.path.isfile(path +'/fasta_files/'+cont+'_'+str(start)+idx+'_h2.fa')):
        try:
            for c in ['1','2']:
                ffile = path + '/fasta_files/'+cont+"_"+str(start)+idx+"_h"+c+'.fa'
                o_pre = path + '/assemble/asm_'+cont+"_"+str(start)+idx+"_h"+c
                log_file = path + '/log/wtdbg2_'+cont+"_"+str(start)+idx+"_h"+c+'.log'
                c_file = path + '/assemble/asm_'+cont+"_"+str(start)+idx+'_h'+c+'.cns.fa'

                subprocess.run(wtdbg2_command+[o_pre]+[ffile], check=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
                subprocess.run('sed -i -e "s/ctg/'+cont+'_'+str(start)+idx+'_h'+c+'_/g" '+c_file, check=True, shell=True)
        except subprocess.CalledProcessError:
            print('no haplotype assembly was possible')
            no_asm += 1
    else:
        print(cont, start, 'no fasta files were created, going to analyse next region')
        no_asm += 1
        continue

    print('completed: assembly', f'it took: {time()-starttime:.2f}s')

print('_______________________________________________')
print('completed: total assembly', f'it took: {time()-asm_time:.2f}s')
print('_______________________________________________')


# mapping und sv call
if args.m == 'one':
    ref_count = 1
    for ref in args.r:

        # map consensus seq(s) againts reference (the one from the bam file)
        map_time = time()
        bam_file = path + '/bam_files/ref'+str(ref_count)+'_merged_regions.bam'
        mapping_command = ['minimap2 -t32 -ax map-pb -2 '+ref+' '+ path + '/assemble/asm_*cns.fa | samtools sort -@4 >'+ bam_file]
        log_file = path + '/log/ref'+str(ref_count)+'_minimap2.log'

        subprocess.run(mapping_command, check=True, shell=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
        print('completed: mapping', f'it took: {time() - map_time:.2f}s')

        # create index to bam file
        index_time = time()
        index_command = ['samtools', 'index', '-b', bam_file, bam_file+'.bai']
        subprocess.run(index_command, check=True)
        print('completed: indexing', f'it took: {time() - index_time:.2f}s')

        # calling sv script to get SVs of consensus seq
        call_time = time()
        bam_local = pysam.AlignmentFile(bam_file, 'rb')
        writing_sv(bam_local,path, ref=str(ref_count))
        bam_local.close()
        print('completed: Sv call', f'it took: {time()-call_time:.2f}s')
        print('_______________________________________________')

        ref_count += 1
elif args.m == 'chr':
    for i in range(len(sv_list['pos1'])):
        cont = sv_list['Chromosome'][i]
        start = sv_list['pos1'][i]
        stop = sv_list['pos2'][i]

        #check if reads where found otherwise skip
        if (not os.path.isfile(path + '/fasta_files/'+cont+'_'+str(start)+idx+'.fa')) and (not os.path.isfile(path +'/fasta_files/'+cont+'_'+str(start)+idx+'_h1.fa')):
            print("skip")
            continue
        # create chr references
        ref_count = 1
        for ref in args.r:

            if not (os.path.isfile(path + '/ref/ref'+str(ref_count)+'_'+cont+'.fa')):
                ref_command =['samtools faidx '+ref +' '+ cont + ' > '+ path + '/ref/ref'+str(ref_count)+'_'+cont+'.fa']
                subprocess.run(ref_command, check=True, shell=True)

            # mapping
            map_time = time()
            bam_file = path + '/bam_files/ref'+str(ref_count)+'_'+cont+"_"+str(start)+idx+'.bam'
            mapping_command = ['minimap2 -t32 -ax map-pb -2 '+ path + '/ref/ref'+str(ref_count)+'_'+cont+'.fa '+ path + '/assemble/asm_'+cont+"_"+str(start)+idx+'*cns.fa | samtools sort -@4 >'+ bam_file]
            log_file = path + '/log/minimap2_'+cont+"_"+str(start)+'.log'

            subprocess.run(mapping_command, check=True, shell=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
            print('completed: mapping', f'it took: {time()-map_timetime:.2f}s')

            # create index to bam file
            index_time = time()
            index_command = ['samtools', 'index', '-b', bam_file, bam_file+'.bai']
            subprocess.run(index_command, check=True)
            print('completed: indexing', f'it took: {time()-index_time:.2f}s')

            # calling sv script to get SVs of consensus seq
            call_time = time()
            bam_local = pysam.AlignmentFile(bam_file, 'rb')

            writing_sv(bam_local,path, [cont], str(ref_count))
            bam_local.close()
            print('completed: Sv call', f'it took: {time()-call_time:.2f}s')
            print('_______________________________________________')

            ref_count += 1
elif args.m == 'ref':
    for i in range(len(sv_list['pos1'])):
        cont = sv_list['Chromosome'][i]
        start = sv_list['pos1'][i]
        stop = sv_list['pos2'][i]

        # check if reads where found otherwise skip
        if not os.path.isfile(path + '/fasta_files/'+cont+'_'+str(start)+idx+'.fa') or not os.path.isfile(path +'/fasta_files/'+cont+'_'+str(start)+idx+'_h1.fa'):
            continue
        ref_count = 1

        for ref in args.r:

            # mapping
            map_time = time()

            bam_file = path + '/bam_files/ref'+str(ref_count)+'_'+cont+"_"+str(start)+idx+'.bam'
            mapping_command = ['minimap2 -t32 -ax map-pb -2 '+ref+' '+path + '/assemble/asm_' + cont + "_" + str(start)+idx+'*cns.fa | samtools sort -@4 >'+ bam_file]
            log_file = path + '/log/minimap2_'+cont+"_"+str(start)+idx+'.log'

            subprocess.run(mapping_command, check=True, shell=True, stdout=open(log_file, 'w'), stderr=subprocess.STDOUT)
            print('completed: mapping', f'it took: {time()-map_time:.2f}s')

            # create index to bam file
            index_time = time()
            index_command = ['samtools', 'index', '-b', bam_file, bam_file+'.bai']
            subprocess.run(index_command, check=True)
            print('completed: indexing', f'it took: {time()-index_time:.2f}s')

            # calling sv script to get SVs of consensus seq
            call_time = time()
            bam_local = pysam.AlignmentFile(bam_file, 'rb')
            writing_sv(bam_local, path, [cont], str(ref_count))
            bam_local.close()
            print('completed: Sv call', f'it took: {time()-call_time:.2f}s')
            print('_______________________________________________')

            ref_count += 1
else:
    print('something went wrong with the mapping option, the default was not recognized')

print(no_asm, 'from', len(sv_list.pos1),'regions were not possible to assemble')

# open haplotype and secondary files and add genotypes
p_res = pd.read_csv(path + '/result/primary_call.txt')
s_res = pd.read_csv(path  + '/result/secondary_call.txt')

p_gen = genotyping(p_res)
s_gen = genotyping(s_res)

p_gen.to_csv(path + '/result/primary_call.txt', index=False)
s_gen.to_csv(path  + '/result/secondary_call.txt', index=False)

sv_list['validated'] = 0

for i in range(len(sv_list['pos1'])):
    cont = sv_list['Chromosome'][i]
    start = sv_list['pos1'][i]
    stop = sv_list['pos2'][i]
    idx = sv_list['indicator'][i]
    type = sv_list['type'][i]
    length = sv_list['length'][i]

    val = p_gen[p_gen.ref==1].copy()
    val = val[val.type.str.casefold() == type.casefold()]
    val = val[val.chr.str.casefold() == cont.casefold()]

    for j in val.index:
        len_val = val.loc[j, 'len']
        pos_val = val.loc[j, 'position']

        if similar_sv(length, len_val, start, pos_val, min_pos=1):
            sv_list.loc[i,'validated'] = 1
            continue

not_val = sv_list[sv_list.validated == 0]
did_val = sv_list[sv_list.validated == 1]

not_val.to_csv(path+'/result/unvalidated_svs.csv', index=False)
did_val.to_csv(path+'/result/validated_svs.csv', index=False)

print('It was possible to validate', sum(sv_list.validated), 'from', len(sv_list.validated))
print('that means', sum(sv_list.validated)/len(sv_list.validated)*100, '%')
