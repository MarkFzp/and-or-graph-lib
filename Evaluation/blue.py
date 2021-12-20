import sys
import os
import subprocess, shlex
from nltk.translate.bleu_score import sentence_bleu

weighted_bleu = True
samples_count = 500

def main():
    if len(sys.argv) != 4:
        print("usage: python blue.py [test_eval] [learned_tree file] [whole data set]")
        exit()

    test_eval = sys.argv[1]
    learned_tree = sys.argv[2]
    dataFile = sys.argv[3]
    name = os.path.splitext(os.path.basename(learned_tree))[0]

    command_sample = './{exe} sample {tree} {count}'.format(exe=test_eval, tree=learned_tree, count=samples_count)
    samples_file = '{}_samples.txt'.format(name)
    log_file = '{}_bleu.txt'.format(name)

    with open(samples_file, "w+") as sample_f:
        p = subprocess.Popen(shlex.split(command_sample), stdout=sample_f)
        p.wait()
    with open(samples_file, 'r') as sample_f, open(dataFile, 'r') as data_f, open(log_file, 'w+') as log_f:
        sum = 0
        total_prob = 0
        references = [line.rstrip().split() for line in data_f]
        for line in sample_f:
            line = line.rstrip().split(': ')
            prob = float(line[0])
            seq = line[1].split()

            len_seq = len(seq)
            bleu = sentence_bleu(references, seq) if len_seq >= 4 else sentence_bleu(references, seq, weights=[1/len_seq for i in range(len_seq)])
            sum += prob * bleu if weighted_bleu else bleu
            total_prob += prob
        
        bleu = sum / total_prob if weighted_bleu else sum / samples_count
        log_f.write('bleu score: {}\n'.format(bleu))
        log_f.write('total likelihood of all sample data and selected parse trees: {}\n'.format(total_prob))

if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Must use Python 3")
    main()
