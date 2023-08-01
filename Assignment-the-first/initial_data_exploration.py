"""The purpose of this script is to do initial exploration
on sequencing data it converts quality score characters to phred scores
and generates histograms for average scores per position on a read """

import bioinfo
import matplotlib.pyplot as plt
import gzip
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="sets global variables of file name, read length, histogram name")
    parser.add_argument("-f", "--file", help="set file name for reading", required=True, type=str)
    parser.add_argument("-l", "--length", help="set read length", required=True, type=int)
    parser.add_argument("-his", "--histogram", help="set histogram name", required=True, type=str)
    return parser.parse_args()


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    # YOUR CODE HERE
    lst = [value] * read_length
    return lst



def populate_list(file: str) -> tuple[list, int]:
    """Takes a fastq file name, calculates sum of quality scores at each position, saves sums into a list, returns a list and line count"""
    q_list = init_list([])
   
    with gzip.open(file, 'rt') as fh:
        i=0
        for line in fh:
            
            i+=1
            if i%4 == 0:
                index = 0
                for q_score in line.strip('\n'):
                    if index < len(line):
                        
                        phred_score = bioinfo.convert_phred(q_score)
                        
                        q_list[index]+=phred_score
                        index+=1
                    else:
                        index = 0
    return(q_list, i)

def calc_mean(q_list: list, lc: int):
    """Takes the list of quality scores and the number of lines in the sequence file
    calculates mean quality score at each position"""
    for index, sum_ in enumerate(q_list):
        q_list[index] = sum_ / (lc/ 4)
    return q_list

def plot_q_score_distribution(q_list: list, x_range: int):
    """Takes a list of mean quality scores and an x axis range which is a read lengt
      and plots a histogram of mean q-scores by position """
    plt.bar(range(x_range), q_list)
    plt.title('Quality Score Distribution')
    plt.ylabel('Mean Phred Quality Score')
    plt.xlabel('Position on the Read')
    plt.savefig(hist_name)
    
if __name__ == "__main__":
    args = get_args()
    file: str = args.file#holds file name
    read_length: int = args.length#holds read length
    hist_name: str = args.histogram#holds histogram name

    q_list, lc = populate_list(file)
    calc_mean(q_list, lc)
    plot_q_score_distribution(q_list, len(q_list))