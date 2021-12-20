"""
autorun script to test recall or perplexity or precision

usage: python Autorun_MCTS.py [recall/perplexity/both/bleu/all/precision]
              [test_run] [test_eval] [whole data set to be splitted] 
              [(optional) (0 or 1) do evalution without training]

initial working directionary contains:
    Autorun_MCTS.py
    config.txt
    test_run (exe)
    test_eval (exe)
    [whole data set]
    (optional) learned_tree.txt (groud truth learned tree for precision)

test_run SHOULD accept 2 arguements: 
    directory to train
    name of training file: located in the directory to train

* the script will automatically update the config during offline training for different train test split portion,
* due to possible non uniqueness of data in data set

The following parameters can be changed:
    train_test_split = 0.7 # size of training set
    iteration = 10 # number of train test splits
    process_count = 10
    perpelxity_penalty = '1e-8'
    samples_count = 1000
"""

import shlex, subprocess, multiprocessing
import random
import sys
from shutil import copyfile
import os
import errno
from nltk.translate.bleu_score import sentence_bleu
cores_count = multiprocessing.cpu_count()

# parameters can be changed ##########
train_test_split = 0.2904 # size of training set
iteration = 6 # number of train test splits
process_count = 3
perpelxity_penalty = '1e-8'
samples_count = 10000
# weighted_bleu = True
######################################


# global variables not to be altered
name = ""
test_run = ""
test_eval = ""
dataFile = ""
recall_flag = False
perpelexity_flag = False
precision_flag = False
bleu_flag = False
bleu_scores = None
recall_scores = None
skip_train = False
learned_tree_index = None
tree_count = None
references = None
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


def Generate_Train_Test_Files():
    s = set()
    l = []
    with open(dataFile, 'r') as f:
        for line in f:
            s.add(line)
            l.append(line)
    len_s = len(s)
    print("All data count: {}".format(len(l)))
    print("Unique data: {}".format(len_s))
    
    for i in range(1, iteration+1):
        test_set = set(random.sample(s, k=round((1 - train_test_split) * len_s)))

        try:
            os.makedirs('{}_{}'.format(name, i))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

        f_train = open('{}_{}/train.txt'.format(name, i), 'w+')
        f_test = open('{}_{}/test.txt'.format(name, i), 'w+')
        for data in l:
            if data in test_set:
                f_test.write(data)
            else:
                f_train.write(data)
        f_train.close()
        f_test.close()


def Copy_Config_And_Update():
    config = open('config.txt', 'r')
    lines = config.readlines()
    offline = True
    for line in lines:
        if "OFFLINE=0" in line:
            offline = False
            break
    config.close()
    
    if not offline:
        for i in range(1, iteration+1):
            copyfile('config.txt', '{}_{}/config.txt'.format(name, i))
    else:
        memory_limit_line = -1
        start_expand_line = -1
        for i in range(len(lines)):
            if 'MEMORY_LIMIT' in lines[i]:
                memory_limit_line = i
            elif 'START_EXPAND' in lines[i]:
                start_expand_line = i
            if memory_limit_line != -1 and start_expand_line != -1:
                break
        if (memory_limit_line == -1) or (start_expand_line == -1):
            raise Exception("MEMORY_LIMIT or START_EXPAND is not found in config.txt")
        
        for i in range(1, iteration+1):
            # adapt the MEMORY_LIMIT and START_EXPAND arguments to newly splitted train.txt
            count = 0
            with open('{}_{}/train.txt'.format(name, i), 'r') as f:
                for line in f:
                    count += 1
            
            lines[memory_limit_line] = 'MEMORY_LIMIT={}\n'.format(count)
            lines[start_expand_line] = 'START_EXPAND={}\n'.format(count)

            config = open('{}_{}/config.txt'.format(name, i), 'w+')
            config.write(''.join(lines))
            config.close()


def Find_Tree_Count():
    with open('config.txt', 'r') as config_f:
        for line in config_f:
            if 'SEARCH_RESOURCE' in line:
                global learned_tree_index
                learned_tree_index = int(line.rstrip().lstrip('SEARCH_RESOURCE='))
                global tree_count
                tree_count = learned_tree_index + 1
                break


def bash(command, logFile):
    with open(logFile, 'w+') as log_f:
        p = subprocess.Popen(shlex.split(command), stdout=log_f, stderr=log_f)
        p.wait()


def Run_Command(i):
    if not skip_train:
        logFile_train = '{}_{}/log_train'.format(name, i)
        command_train = './{exe} {name}_{iter}/ train.txt'.format(exe=test_run, name=name, iter=i)
        bash(command_train, logFile_train)

    dir_list = os.listdir('{}_{}'.format(name, i))
    
    if recall_flag:
        logFile_recall = '{}_{}/log_recall'.format(name, i)
        command_recall = './{exe} recall {name}_{iter}/learned_tree.txt {name}_{iter}/test.txt'.format(exe=test_eval, name=name, iter=i)
        bash(command_recall, logFile_recall)
        for j in range(1, tree_count):
            if ('learned_tree_{}.txt'.format(j) in dir_list):
                logFile_recall = '{}_{}/log_recall_{}'.format(name, i, j)
                command_recall = './{exe} recall {name}_{iter}/learned_tree_{j}.txt {name}_{iter}/test.txt'.format(exe=test_eval, name=name, iter=i, j=j)
                bash(command_recall, logFile_recall)

    if perpelexity_flag:
        logFile_perplexity = '{}_{}/log_perplexity'.format(name, i)
        command_perplexity = './{exe} perplexity {name}_{iter}/learned_tree.txt {name}_{iter}/test.txt {penalty}'.format(exe=test_eval, name=name, iter=i, penalty=perpelxity_penalty)
        bash(command_perplexity, logFile_perplexity)
    
    if precision_flag:
        logFile_precision = '{}_{}/log_precision'.format(name, i)
        command_precision = './{exe} precision {name}_{iter}/learned_tree.txt learned_tree.txt'.format(exe=test_eval, name=name, iter=i)
        bash(command_precision, logFile_precision)

    if bleu_flag:
        for j in range(0, tree_count):
            appendix = '' if j == 0 else '_' + str(j)
            if ('learned_tree{}.txt'.format(appendix) in dir_list):
                learned_tree = "{}_{}/learned_tree{}.txt".format(name, i, appendix)
                command_sample = './{exe} sample {tree} {count} 1'.format(exe=test_eval, tree=learned_tree, count=samples_count)
                samples_file = '{}_{}/samples{}.txt'.format(name, i, appendix)
                with open(samples_file, "w+") as sample_f:
                    p = subprocess.Popen(shlex.split(command_sample), stdout=sample_f)
                    p.wait()
                    
                with open(samples_file, 'r') as sample_f:
                    sum = 0
                    total_prob = 0 
                    for line in sample_f:
                        line = line.rstrip().split(': ')
                        prob = float(line[0])
                        seq = line[1].split()

                        len_seq = len(seq)
                        bleu = sentence_bleu(references, seq) if len_seq >= 4 else sentence_bleu(references, seq, weights=[1/len_seq for i in range(len_seq)])
                        sum += prob * bleu
                        total_prob += prob

                    bleu_scores[(i - 1) * (tree_count) + j] = sum / total_prob


def Operation_Summary(operation):
    operation_log = 'log_{}'.format(operation)
    summary_file = '{}_summary.txt'.format(operation_log)
    with open(summary_file, 'w+') as log_f:
        sum = 0
        for i in range(1, iteration+1):
            logFile = '{}_{}/{}'.format(name, i, operation_log)
            for line in readlines_reverse(logFile):
                if (operation in line):
                    break
            log_f.write('{}_{} - '.format(name, i))
            log_f.write(line + '\n')
            sum += float(line.lstrip(operation + ': '))
        average = sum / iteration
        log_f.write('\n')
        log_f.write('average: {}\n'.format(average))


def main():
    if (len(sys.argv) != 5 and len(sys.argv) != 6) or (sys.argv[1] not in ['recall', 'bleu', 'both', 'perplexity', 'all', 'precision']):
        print("usage: python Autorun_MCTS.py [recall/bleu/both/perplexity/all/precision] [test_run] [test_eval] [whole data set to be splitted] [(optional) (0 or 1) do evalution without training]")
        exit()

    eval_operation = sys.argv[1]
    global test_run
    test_run = sys.argv[2]
    global test_eval
    test_eval = sys.argv[3]
    global dataFile
    dataFile = sys.argv[4]
    try:
        global skip_train
        skip_train = bool(int(sys.argv[5]))
    except IndexError:
        pass
    global name
    name = os.path.splitext(os.path.basename(dataFile))[0]

    if eval_operation in ['recall', 'both', 'all']:
        global recall_flag
        recall_flag = True
    if eval_operation in ['perplexity', 'all']:
        global perpelexity_flag
        perpelexity_flag = True
    if eval_operation in ['bleu', 'both', 'all']:
        global bleu_flag
        bleu_flag = True
    if eval_operation == 'precision':
        global precision_flag
        precision_flag = True
        # global skip_train
        skip_train = True
    
    if not skip_train:
        Generate_Train_Test_Files()
        Copy_Config_And_Update()
    
    if recall_flag or bleu_flag:
        Find_Tree_Count()
    
    if bleu_flag:
        global bleu_scores
        bleu_scores = multiprocessing.Array('d', [-1] * (iteration * (tree_count)))
        with open(dataFile, 'r') as data_f:
            global references
            references = [line.rstrip().split() for line in data_f]

    # train, recall, perplexity
    pool = multiprocessing.Pool(process_count)
    for i in range(1, iteration+1):
        pool.apply_async(Run_Command, args=(i,))
    pool.close()
    pool.join()

    # Summary Output
    if recall_flag:
        summary_file = 'log_recall_summary.txt'
        with open(summary_file, 'w+') as log_f:
            global recall_scores
            recall_scores = [-1] * (iteration * tree_count)
            sum = [0.0] * tree_count
            skip_average = [False] * tree_count
            for i in range(1, iteration+1):
                dir_list = os.listdir('{}_{}'.format(name, i))
                for j in range(0, tree_count):
                    appendix = '' if j == 0 else '_' + str(j)
                    if 'log_recall{}'.format(appendix) in dir_list:
                        logFile = '{}_{}/log_recall{}'.format(name, i, appendix)
                        for line in readlines_reverse(logFile):
                            if ('recall' in line):
                                break
                        score = float(line.lstrip('recall: '))
                        log_f.write('shuffle_{}: learned_tree{}: {}\n'.format(i, appendix, score))
                        sum[j] += score
                        recall_scores[(i - 1) * tree_count + j] = score
                    else:
                        skip_average[j] = True
            average = [sum[i] / iteration if not skip_average[i] else -1 for i in range(tree_count)]
            for j in range(tree_count):
                if not skip_average[j]:
                    appendix = ' of CHOSEN learned trees' if j == 0 else '_' + str(j)
                    log_f.write('average{}: {}\n'.format(appendix, average[j]))

    
    if perpelexity_flag:
        summary_file = 'log_perplexity_summary.txt'
        with open(summary_file, 'w+') as log_f:
            sum = 0.0
            for i in range(1, iteration+1):
                logFile = '{}_{}/log_perplexity'.format(name, i)
                for line in readlines_reverse(logFile):
                    if ('perplexity' in line):
                        break
                log_f.write('{}_{} - '.format(name, i))
                log_f.write(line + '\n')
                sum += float(line.lstrip('perplexity: '))
            average = sum / iteration
            log_f.write('\n')
            log_f.write('average: {}\n'.format(average))


    if bleu_flag:
        with open('log_bleu_summary.txt', 'w+') as log_f:
            sum = [0.0 for k in range(tree_count)]
            skip_average = [False] * (tree_count)
            for i in range(iteration * (tree_count)):
                score = bleu_scores[i]
                tree_index = i % (tree_count)
                if score != -1:
                    sum[tree_index] += score
                    folder_index = (i // (tree_count)) + 1
                    log_f.write('shuffle_{}: learned_tree{}: {}\n'.format(folder_index, '' if tree_index == 0 else '_' + str(tree_index), score))
                else:
                    skip_average[tree_index] = True
            average = [sum[i] / iteration if not skip_average[i] else -1 for i in range(tree_count)]
            for j in range(tree_count):
                if not skip_average[j]:
                    appendix = ' of CHOSEN learned trees' if j == 0 else '_' + str(j)
                    log_f.write('average{}: {}\n'.format(appendix, average[j]))
    

    if bleu_flag and recall_flag:
        with open('f_measure_summary.txt', 'w+') as f_measure_f:
            for i in range(iteration * tree_count):
                recall_score = recall_scores[i]
                bleu_score = bleu_scores[i]
                if recall_score != -1 and bleu_score != -1:
                    f_measure_score = 2 * (bleu_score * recall_score / (bleu_score + recall_score))
                    folder_index = (i // (tree_count)) + 1
                    tree_index = i % (tree_count)
                    f_measure_f.write('shuffle_{}: learned_tree{}: {}\n'.format(folder_index, '' if tree_index == 0 else '_' + str(tree_index), f_measure_score))



if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Must use Python 3")
    main()
