#from bioinfo.py import convert_phred

file = ''#holds path to a file or file name for reading
def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    # YOUR CODE HERE
    lst = [value] * 101
    return lst


def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    phred = (ord(letter) - 33)
    return phred

def populate_list(file: str) -> tuple[list, int]:
    """Takes a fastq file name, calculates sum of quality scores at each position, saves sums into a list, returns a list and line count"""
    q_list = init_list([])
    with open(file) as fh:
        i=0
        for line in fh:
            
            i+=1
            line.strip('\n')
            if i%4 == 0:
                index = 0
                for q_score in line:
                    if index < len(line):
                        phred_score = convert_phred(q_score)
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
    
if __name__ == "__main__":
    q_list, lc = populate_list(file)
    mean = calc_mean(q_list, lc)