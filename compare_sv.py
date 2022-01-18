import pysam
import sys
import pandas as pd
import numpy as np


def similar_sv(len_a, len_b, pos_a, pos_b, min_len=0.1, min_pos=2):
    #diff_len = np.abs(len_a - len_b)
    diff_avg = (len_a + len_b)/2
    ratio_len = np.abs(len_a/len_b - 1)

    if ratio_len <= min_len:  # 10% in size difference is allowed
        diff_pos = np.abs(pos_a - pos_b)
        return diff_pos < diff_avg/min_pos   # must overlap by half
    return False




def genotyping(asm_calls):
    asm_calls['genotype'] = '0/1'

    #calls = asm_calls
    for i in asm_calls.index:
        if i not in asm_calls.index:
            continue
        if 'h1' in asm_calls.loc[i, 'read_name']:
            h2 = asm_calls[asm_calls.loc[i, 'read_name'][:-3]+'2_1' == asm_calls.read_name].copy() #works with the indicator since hi comes after
            h2 = h2[asm_calls.loc[i, 'ref']== h2.ref]
            h2 = h2[asm_calls.loc[i, 'chr'].casefold()== h2.chr.str.casefold()]
            h2 = h2[asm_calls.loc[i, 'type'].casefold()== h2.type.str.casefold()]


            for j in h2.index:
                len_h1 = asm_calls.loc[i, 'len']
                len_h2 = h2.loc[j, 'len']
                pos_h1 = asm_calls.loc[i, 'position']
                pos_h2 = h2.loc[j, 'position']
                if similar_sv(len_h1, len_h2, pos_h1, pos_h2):
                    #change genotype and delete row in asm_calls
                    position_query = '(position != '+str(pos_h2)+') | (read_name != "'+h2.loc[j, 'read_name']+'")'
                    # something does not work here
                    asm_calls.loc[i, 'genotype'] = '1/1'

                    asm_calls = asm_calls.query(position_query)
        elif not 'h2' in asm_calls.loc[i, 'read_name']:
            asm_calls.loc[i, 'genotype'] = '1/1'
    return asm_calls


#def validation()
