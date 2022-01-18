import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
from tools_from_anton import pca_on_read_distance, cluster_reads, plot_with_clustering

def count_back_cigar(cigarCode, cigarLen, pos, goal):
    if cigarCode == 0: # match
        pos += cigarLen
        goal -= cigarLen
    elif cigarCode == 1:  # insertion
        pos += cigarLen
    elif cigarCode == 2:  # deletion
        goal -= cigarLen
    elif cigarCode == 4 or cigarCode == 5:  # clipped
        pos += cigarLen
    else:
        print("unexpected cigar Type \n only expecting: M, I, D, S, H")
    return pos, goal

#buffer size of the SV would be nice
def get_read_SV_coord(read, start, stop, buffer=200):
    #print(read.get_reference_positions(full_length=True))
    r_pos = 0  # for position on read
    sv_pos = start - read.reference_start - buffer  # for reference length
    sv_len = stop - start + buffer #length of SV
    deletion_read_start = None
    deletion_read_stop = None

    for (cigarCode, cigarLen) in read.cigartuples:
        if sv_pos > 0:
            r_pos, sv_pos = count_back_cigar(cigarCode, cigarLen, r_pos, sv_pos)

        elif not deletion_read_start:
            deletion_read_start = r_pos + sv_pos
            r_pos += sv_pos

        elif sv_len > 0  and r_pos < read.query_alignment_end:
            r_pos, sv_len = count_back_cigar(cigarCode, cigarLen, r_pos, sv_len)
        else:
            break
    deletion_read_stop = r_pos #+ sv_len

    if deletion_read_start > deletion_read_stop:
        print('something went wrong here, the intervall ends before it starts, return naive positions instead')
        print(read.query_name, deletion_read_start, deletion_read_stop, "vs",  sv_pos + buffer + read.query_alignment_start, stop - read.reference_start + read.query_alignment_start)
        return sv_pos + buffer + read.query_alignment_start, stop - read.reference_start + read.query_alignment_start
    return deletion_read_start, deletion_read_stop


def write_fasta(bam, cont, start, stop, path, buffer=10000, idx=''):
    seqs_h1 = []
    seqs_h2 = []
    seqs = []
    fasta = []
    no_tag = 0
    sec = 0
    for read in bam.fetch(contig=cont, start=start, stop=stop+1):
        if read.is_secondary:
            sec += 1
            continue

        r_start, r_stop = get_read_SV_coord(read, start, stop, buffer=200)  # get coordinates of SV on read

        f_start = r_start - buffer
        if f_start < read.query_alignment_start:
            f_start = read.query_alignment_start
        f_stop = r_stop + buffer
        if f_stop > read.query_alignment_end:
            f_stop = read.query_alignment_end


        if read.has_tag('HP'):
            if read.get_tag('HP') == 1:
                seqs_h1.append(SeqRecord(Seq(read.query_sequence[f_start:f_stop]), id=read.query_name))
            else:
                seqs_h2.append(SeqRecord(Seq(read.query_sequence[f_start:f_stop]), id=read.query_name))
        else:
            no_tag = no_tag + 1
            #maybe include reads not phased as well for clustering
        try:
            fasta.append(SeqRecord(Seq(read.query_sequence[f_start:f_stop]), id=read.query_name))
            seqs.append(read.query_sequence[r_start:r_stop])
        except TypeError:
            print(cont, start,'for this read, no sequence is found:', read.query_name,)  # this happens a lot, how does it happen? secondary mapping, supplementary? only deletion?
            #print(sys.exc_info())
    if no_tag > 0 or sec > 0:
        print(cont, start, no_tag, 'read(s) did no have the HP tag and', sec, 'where secondary alignments')

    if len(seqs) < 10:
        print(cont, start,'not enough reads for assembly, only:', len(seqs))
        return 0


    if len(seqs_h1) < 10 or len(seqs_h2) < 10 :
#        if(stop-start < 500):
        print(cont, start,'try clustering')
        dots = pca_on_read_distance(seqs) # takes time
        cluster = cluster_reads(dots, eps=0.2, min_samples=10) # changed values from antons default values
        plot_with_clustering(dots, cluster, additional_labels=None)
        #plt.scatter(dots, cluster.labels_)
        plt.savefig(path+'/plot/'+cont+'_'+str(start)+'_read_cluster.png')
        plt.close()

        if sum(np.unique(cluster.labels_)) < 2:
            seqs_cl0 = []
            seqs_cl1 = []
            for i in range(len(cluster.labels_)):
                if cluster.labels_[i] == -1:
                    continue
                if cluster.labels_[i] == 0:
                    seqs_cl0.append(fasta[i])
                else:
                    seqs_cl1.append(fasta[i])
            SeqIO.write(seqs_cl0, path + "/fasta_files/"+cont+"_"+str(start)+idx+'_h1.fa', "fasta")
            SeqIO.write(seqs_cl1, path + "/fasta_files/"+cont+"_"+str(start)+idx+'_h2.fa', "fasta")
            print(cont, start,'two sufficient clusters were found')
            print(cont, start,'HP1', len(seqs_cl0), 'reads and HP2', len(seqs_cl1))
        else:
            SeqIO.write(fasta, path + '/fasta_files/'+cont+'_'+str(start)+idx+'.fa', "fasta")
            print(cont, start, 'not enough information for seperat assembly,', len(fasta), 'reads used')
    else:
        SeqIO.write(seqs_h1, path + '/fasta_files/'+cont+'_'+str(start)+idx+'_h1.fa', "fasta")
        SeqIO.write(seqs_h2, path + '/fasta_files/'+cont+'_'+str(start)+idx+'_h2.fa', "fasta")
        print(cont, start, 'the fasta files were created')
        print(cont, start,'HP1', len(seqs_h1), 'reads and HP2', len(seqs_h2))

    return 1


##how do i keep them in here (use this script seperatly)
#bam_pb = pysam.AlignmentFile(sys.argv[1], 'rb')
#write_fasta(bam_pb, sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
#bam_pb.close()
