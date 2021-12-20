from __future__ import division

import math
import sys
import fractions
import warnings
from collections import Counter

from nltk.translate.bleu_score import corpus_bleu
from nltk.translate.bleu_score import sentence_bleu

try:
    fractions.Fraction(0, 1000, _normalize=False)
    from fractions import Fraction
except TypeError:
    from nltk.compat import Fraction

# def sentence_bleu(references, hypothesis, weights=(0.25, 0.25, 0.25, 0.25),
#                             smoothing_function=None, auto_reweigh=False):
#     return corpus_bleu([references], [hypothesis],
#                        weights, smoothing_function, auto_reweigh)
global local_dir
import os
local_dir = os.path.dirname(os.path.realpath(__file__)).rstrip('/') + '/'

hypothesis_path = sys.argv[1]
references_path = sys.argv[2]

references = []
candidates = []

with open(hypothesis_path, 'r') as f:
    string = ''
    for line in f:
        string += line
    candidates = string.split()
    
print("Print Candidates: \n")
print(candidates)

with open(references_path, 'r') as f:
    string = ''
    for line in f:
        string += line
    references = string.split()

references_wrapper = []
references_wrapper.append(references)

print("\nPrint references: \n")
print(references_wrapper)
score = sentence_bleu(references_wrapper, candidates)

print("The corpus BLEU score is: " + str(score))
