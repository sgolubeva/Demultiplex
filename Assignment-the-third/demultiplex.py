#!/usr/bin/env python
"""The purpose of this script is iterate over illumina reads and place 
them into three categories depending on the index reads. One category is 
matched reads which means that index1 and index 2 match each other with index 2 being
a reverse compliment of index 1. Second category is index hopped: index 1 and 2 don't match,
but they are represented in the list of indexes. Third category is unmatched indexes: indexes that 
contain 'N's or indexes not represented in the list of indexes"""
import argparse
import gzip
import itertools
import matplotlib.pyplot as plt

def get_args():
    """sets global variables to the file names for the files to be demultiplexed"""
    parser = argparse.ArgumentParser(description="sets global variables for file names to read1, 2, 3, and 4")
    parser.add_argument("-r1", "--read1", help="holds file name or path to read 1", required=True, type=str)
    parser.add_argument("-r2", "--read2", help="holds file name or path to read 2", required=True, type=str)
    parser.add_argument("-r3", "--read3", help="holds file name or path to read 3", required=True, type=str)
    parser.add_argument("-r4", "--read4", help="holds file name or path to read 4", required=True, type=str)
    return parser.parse_args()


def rev_compliment(index: str) -> str:
    """Takes index as a string and returns reverse compliment of an index"""
    bases = {'A':'T', 'G':'C', 'T': 'A', 'C': 'G', 'N':'N'}
    rev_ind = index[::-1]
    rev_compliment = []
    for base in rev_ind:
        rev_compliment.append(bases[base])
    return ''.join(rev_compliment)
   

def process_record(fh) -> list:
    """Takes file handle object, returns a list containing header, sequence, +,
    quality score line"""
    header = fh.readline().strip()
    if header == '':
        return []
    seq = fh.readline().strip()
    plus = fh.readline().strip()
    q_score = fh.readline().strip()
    return [header, seq, plus, q_score]

def write_file(fh, index1: str, index2: str, r1: list, r4: list):
    """receives a tuple of 2 file handles: one for read 1 and another for read 4,
    index, read 1 record, read 4 record
    writes records into corresponding files"""
    fh_r1 = fh[0]
    fh_r4 = fh[1]
    index_record = ' ' + index1 + '-' + index2 # create index that will be apended to the header in form: 
                                               # header index1-index2
    r1[0] += index_record #add index record to a header of read 1
    r4[0] += index_record #add index record to a header of read 2
    r1 = [f'{x}\n' for x in r1] #add end of line charachters to each string in r1 list
    r4 = [f'{x}\n' for x in r4] #add end of line charachters to each string in r2 list
    r1_record = ''.join(r1) #join contents of the list into a string for read 1 
    r4_record = ''.join(r4) #join contents of the list into a string for read 4
    fh_r1.write(r1_record)
    fh_r4.write(r4_record)



if __name__ == "__main__":
    file_handle_dict = {} #holds file handles for files to write to
    #keys are indexes, values are tuples of file handles for each file
    index_comb_dict =  {} #holds all index combinations
    match_count = {} # hold matched index counts
    index_list = [] #holds a list of indexes
    args = get_args() # get args object
    read1 = args.read1 #read 1 global
    read2 = args.read2 # read 2 global
    read3 = args.read3 # read 3 global
    read4 = args.read4 # read 4 global
    total_rec_count = 0# counts number of total reads
    matched_number = 0# counts number of matched reads
    unkn_count = 0#counts number of unknown reads
    hop_count = 0# counts number of index hopped reads

    # create a dictionary with indexes as keys and file handles as tuples
    with open ('indexes.txt', 'r') as fh: # add indexes as keys 
        for index in fh:
            index = index[:8]
            index_list.append(index)
            match_count[index] = 0
            file_handle_dict[index] = (open(f'R1_{index}_matching.fq', 'w'), open(f'R4_{index}_matching.fq', 'w'))
    file_handle_dict['unknown'] = (open(f'R1_unknown.fq', 'w'), open(f'R4_unknown.fq', 'w'))
    file_handle_dict['index_hopped'] = (open(f'R1_index_hopped.fq', 'w'), open(f'R4_index_hopped.fq', 'w'))

    #populate the dic gandling index hopped counts
    combinations = itertools.permutations(index_list, 2)# get all possible index combinations
    for comb in combinations: # populate dict with all possible combinations
        index_comb_dict[comb] = 0

    #read read1, read2, read3, read4 files
    with gzip.open(read1, 'rt') as fh1, gzip.open(read2, 'rt') as fh2, gzip.open(read3, 'rt') as fh3, gzip.open(read4, 'rt') as fh4:
        while True:
            r1 = process_record(fh1)
            r2 = process_record(fh2)
            r3 = process_record(fh3)
            r4 = process_record(fh4)

            #break the while true loop when process_record returns an empty list
            if r4 == []:
                for index in file_handle_dict:# close file handles for 52 files
                    for fh in file_handle_dict[index]:
                        fh.close()
                break
            rev_R3 = rev_compliment(r3[1])# reverse compliment read 3

            if 'N' in r2[1] or 'N' in rev_R3: #check if either R2 or R3 contain Ns
                write_file(file_handle_dict['unknown'], r2[1], rev_R3, r1, r4)
                unkn_count += 1
                total_rec_count += 1
            elif r2[1] == rev_R3 and r2[1] in file_handle_dict:# check if R2 and R3 are matched
                write_file(file_handle_dict[r2[1]], r2[1], rev_R3, r1, r4)
                total_rec_count += 1
                match_count[r2[1]] += 1
                matched_number += 1
            elif (r2[1] in file_handle_dict) and (rev_R3 in file_handle_dict): # check if index hopped 
                write_file(file_handle_dict['index_hopped'], r2[1], rev_R3, r1, r4)
                hop_count += 1
                total_rec_count += 1
                index_comb_dict[(r2[1], rev_R3)] += 1
            else: #the rest of indexes are unknown
                write_file(file_handle_dict['unknown'], r2[1], rev_R3, r1, r4)
                unkn_count += 1
                total_rec_count += 1

    #write data into tsv files 
    with open('matched_percents.tsv', 'w') as mp, open('general_stats.tsv', 'w') as gs, open('index_hopped_stats.tsv','w') as ih:
        mp.write(f'Index\tCount\tPercent of total\tPercent of matched\n')
        for match in match_count:
            mp.write(f'{match}\t{match_count[match]}\t{(match_count[match]*100)/total_rec_count: .2f}\t{(match_count[match]*100)/matched_number: .2f}\n')

        gs.write(f'\tCount\tPercent\n')
        gs.write(f'Total\t{total_rec_count}\t{(total_rec_count*100)/total_rec_count: .2f}\n')
        gs.write(f'Matched\t{matched_number}\t{(matched_number*100)/total_rec_count: .2f}\n')
        gs.write(f'Index Hopped\t{hop_count}\t{(hop_count*100)/total_rec_count: .2f}\n')
        gs.write(f'Unknown\t{unkn_count}\t{(unkn_count*100)/total_rec_count: .2f}\n')

        ih.write('Combination\tCount\tPercent\n')
        for combination in index_comb_dict:
            ih.write(f'{combination}\t{index_comb_dict[combination]}\t{(index_comb_dict[combination]*100)/hop_count: .2f}\n')
    
