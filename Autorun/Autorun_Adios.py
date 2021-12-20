"""
autorun script to test recall of benchmark (jAOG in OpenBottle)

The following parameters can be changed:
    train_test_split = 0.7
    iteration = 20
    process_count = 7
"""

import shlex, subprocess, multiprocessing
import random
import sys
import os
import errno
from collections import Counter
from nltk.translate.bleu_score import sentence_bleu
cores_count = multiprocessing.cpu_count()

# parameters can be changed ##########
train_test_split = 0.4357
iteration = 10 # number of train test splits
process_count = 5
weighted_bleu = True
random.seed(a='vcla')

eta = 0.8
alpha = 0.01
window = 5
######################################

# global variables not to be altered
inner_iteration = 3
name = ""
local_dir = ""
adios_dir = ""
dataFile = ""
bleu_flag = False
recall_flag = False
skip_train = False
bleu_scores = None
process_count = process_count if process_count <= (cores_count - 1) else (cores_count - 1)


def readlines_reverse(filename):
    with open(filename) as qfile:
        qfile.seek(0, os.SEEK_END)
        position = qfile.tell()
        line = ''
        while position >= 0:
            qfile.seek(position)
            next_char = qfile.read(1)
            if next_char == "\n":
                yield line[::-1]
                line = ''
            else:
                line += next_char
            position -= 1
        yield line[::-1]


def Generate_Train_Test_Files_And_Convert(dataFile):
    s = set()
    l = []
    l_unique = []
    with open(dataFile, 'r') as f:
        for line in f:
            if line not in s:
                s.add(line)
                l_unique.append(line)
            l.append(line)
    len_unique = len(l_unique)
    print("All data count: {}".format(len(l)))
    print("Unique data: {}".format(len_unique))
    test_set_size = round((1 - train_test_split) * len_unique)
    
    for i in range(1, iteration+1):
        test_set = set(random.sample(l_unique, k=test_set_size))

        try:
            os.makedirs('{}_{}'.format(name, i))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

        f_train = open('{}_{}/train.txt'.format(name, i), 'w+')
        f_test = open('{}_{}/test.txt'.format(name, i), 'w+')

        train_count = 0
        train_tokens = set()
        for data in l:
            if not data in test_set:
                train_count += 1
                f_train.write('* ' + data.rstrip() + ' #\n')
                for token in data.rstrip('\n').split():
                    train_tokens.add(token)
        
        test_count = 0
        for data in l:
            if data in test_set:
                all_tokens_seen = True
                for token in data.rstrip('\n').split():
                    if token not in train_tokens:
                        all_tokens_seen = False
                        break
                if all_tokens_seen:
                    test_count += 1
                    f_test.write('* ' + data.rstrip() + ' #\n')
        
        print("\ttraining size: {}, test size: {}, training set to train-test-combined: {}".format(train_count, test_count, train_count / (train_count + test_count)))
        
        f_train.close()
        f_test.close()


def Run_Command(i):
    if not skip_train:
        working_dir = '{}{}_{}'.format(local_dir, name, i)
        command = 'java -cp {adios_dir} com.ADIOS.Main -a TRAIN -f train -t test -d {local_dir} -E {eta} -S {alpha} -W {window}'.format\
            (adios_dir=adios_dir, local_dir=working_dir, eta=eta, alpha=alpha, window=window)
        for j in range(1, inner_iteration+1):
            logFile = "{}_{}/log_recall_{}".format(name, i, j)
            with open(logFile, "w+") as log_f:
                p = subprocess.Popen(shlex.split(command), stdout=log_f, stderr=log_f)
                p.wait()
            mv_command = 'mv {name}_{i}/train_generate.txt {name}_{i}/train_generate_{j}.txt'.format(name=name, i=i, j=j)
            subprocess.Popen(shlex.split(mv_command)).wait()
            rm_command = 'rm {}_{}/train_formatted.txt'.format(name, i)
            subprocess.Popen(shlex.split(rm_command)).wait()
    
    if bleu_flag:
        sum = 0
        slice = round(1000 / inner_iteration)
        with open(dataFile, 'r') as data_f:
            references = [line.strip().split() for line in data_f]
            samples = Counter()
            for j in range(1, inner_iteration+1):
                samples_file = '{}_{}/train_generate_{}.txt'.format(name, i, j)
                with open(samples_file, 'r') as sample_f:
                    for k, line in enumerate(sample_f, start=0):
                        if k == slice:
                            break
                        seq = line.strip().lstrip('*').rstrip('#').strip().split()
                        samples[tuple(seq)] += 1

            for seq, count in samples.items():
                len_seq = len(seq)
                bleu = sentence_bleu(references, seq) if len_seq >= 4 else sentence_bleu(references, seq, weights=[1/len_seq for i in range(len_seq)])
                sum += bleu * count
        bleu_scores[i-1] = sum / 1000


def main():
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("usage: python Autorun_Benchmark.py [recall/bleu/both] [JAdios directory] [whole data set to be splitted] [(optional) skip_train]")
        exit()

    global local_dir
    local_dir = os.path.dirname(os.path.realpath(__file__)).rstrip() + '/'
    eval_operation = sys.argv[1]
    global adios_dir
    adios_dir = sys.argv[2].rstrip('/') + '/'
    global dataFile
    dataFile = sys.argv[3]
    global name
    name = os.path.splitext(os.path.basename(dataFile))[0]
    try:
        global skip_train
        skip_train = bool(int(sys.argv[4]))
    except IndexError:
        pass

    if eval_operation in ['recall', 'both']:
        global recall_flag
        recall_flag = True
    if eval_operation in ['bleu', 'both']:
        global bleu_flag
        bleu_flag = True

    if not skip_train:
        Generate_Train_Test_Files_And_Convert(dataFile)

    if bleu_flag:
        global bleu_scores
        bleu_scores = multiprocessing.Array('d', [0.0]*iteration)

    # Adios training
    pool = multiprocessing.Pool(process_count)
    for i in range(1, iteration+1):
        pool.apply_async(Run_Command, args=(i, ))
    pool.close()
    pool.join()

    if recall_flag:
        with open('log_recall_summary.txt', 'w+') as log_f:
            sum = 0
            for i in range(1, iteration+1):
                for j in range(1, inner_iteration+1):
                    logFile = '{}_{}/log_recall_{}'.format(name, i, j)
                    for line in readlines_reverse(logFile):
                        if ('Recall = ' in line):
                            break
                    log_f.write(logFile + ' - ')
                    log_f.write(line + '\n')
                    sum += float(line.lstrip('Recall = '))
            average = sum / (iteration * inner_iteration)
            log_f.write('\n')
            log_f.write('sum of all recalls: {}\n'.format(sum))
            log_f.write('average of all recalls: {}\n'.format(average))
    
    if bleu_flag:
        with open('log_bleu_summary.txt', 'w+') as log_f:
            sum = 0
            for i in range(iteration):
                score = bleu_scores[i]
                sum += score
                log_f.write('{}: {}\n'.format(i, score))
            average = sum / iteration
            log_f.write('average: {}\n'.format(average))


if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Must use Python 3")
    main()
