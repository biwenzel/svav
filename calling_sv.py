import pysam
import sys
import pandas


#could also make the max size of the variantion dynamic
def calling_sv(cigartouple, start, cont, file, name, ref):
    pos = start
    for cigarCode, cigarLen in cigartouple:
        if cigarCode == 0: # match
            pos += cigarLen
        elif cigarCode == 1:  # insertion
            if cigarLen >= 50:
                file.write(cont+','+str(pos+1)+',Ins,'+str(cigarLen)+','+name+','+ref+'\n')
        elif cigarCode == 2:  # deletion
            if cigarLen >= 50:
                file.write(cont+','+str(pos+1)+',Del,'+str(cigarLen)+','+name+','+ref+'\n')
            pos += cigarLen
        elif cigarCode == 4 or cigarCode == 5:  # clipped
            continue
        else:
            print("unexpected cigarType")

#can I give a differnet message for a homozygous call?
def writing_sv(bam, path, conts = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'], ref='1'):
    # open them before, write header, before mapping loop?
    h_file = open(path + '/result/primary_call.txt', 'a')
    s_file = open(path + '/result/secondary_call.txt', 'a')
    for chr in conts:
        try:
            reads = bam.fetch(chr)
        except ValueError:
            print(chr, "does not exist on reference", ref)
            continue
        # maybe look at everything, but how to get the Chromosome
        for read in reads:
            if read.is_secondary:
                calling_sv(read.cigartuples, read.reference_start, chr, s_file, read.query_name, ref)
            else:
                calling_sv(read.cigartuples, read.reference_start, chr, h_file, read.query_name, ref)
                print(read.query_name, 'mapped to', chr, read.reference_start, 'on reference', ref)
    h_file.close()
    s_file.close()
    return


##how do i keep them in here (use this script seperatly)
#bam_pb = pysam.AlignmentFile(sys.argv[1], 'rb')
#writing_sv(bam_pb, sys.argv[2])
#bam_pb.close()

#need to compare called SVs with given SVs --> are they in the same region?
