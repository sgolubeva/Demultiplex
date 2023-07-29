#!/usr/bin/env python
# Author: <sgolubev> <optional@email.address>

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
    '''Write your own doc string'''
    pass

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def gc_content():
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    pass

def oneline_fasta():
    '''docstring'''
    pass

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    