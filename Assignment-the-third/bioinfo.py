#!/usr/bin/env python
# Author: Sasha Golubeva sgolubev.uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    phred = (ord(letter) - 33)
    return phred

def qual_score(phred_score: str) -> float:
    '''Takes a line of phred scores and returns an average phred score for this line'''
    phred_score = phred_score.upper()
    phred_sum = 0
    for score in phred_score:
        phred_sum += convert_phred(score)
    return phred_sum/len(phred_score)
              

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    if seq != '':
        return set(seq)<=(RNAbases if RNAflag else DNAbases)
    else:
        return False 

def gc_content(DNA: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    if validate_base_seq(DNA):
        DNA = DNA.upper()         
        Gs = DNA.count("G")       
        Cs = DNA.count("C")       
        return (Gs+Cs)/len(DNA)
    else:
        return 0.0

def oneline_fasta(file_read: str, file_write: str):
    '''takes a path to a multiline fasta file and generates a new one line fasta file by putting
     sequences into the same line '''
    sequence = ''
    header = ''
    with open(file_read, 'r') as fh, open(file_write, 'w') as wh:
        for line in fh:
            line.strip('\n')
            if line.startswith('>'):
                if sequence == '':
                    header = line
                else:
                    wh.write(f'{header}')
                    wh.write(f'{sequence}\n')               
                    sequence = ''
                    header = line
            else:
                sequence += line.strip('\n')
    wh.write(f'{header}')
    wh.write(f'{sequence}\n')

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert validate_base_seq('AATTT') == True
    assert validate_base_seq('aattg') == True
    assert validate_base_seq('ACGU', RNAflag=True) == True
    assert validate_base_seq('') == False
    assert validate_base_seq('furry cat') == False
    print('Validate Nucleic Acid function is working')

    assert gc_content('GCGCGC') == 1
    assert gc_content('cgcgcg') == 1
    assert gc_content('cgat') == 0.5
    assert gc_content('furry cat') == 0.0
    print('GC content is calculating correctly')

    assert qual_score('EEE') == 36
    assert qual_score('eee') == 36
    assert qual_score('ABCD') == 33.5
    print('Average quality scores are calculated correctly')



    