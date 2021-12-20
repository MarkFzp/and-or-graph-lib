"""
autorun script to test recall or perplexity or precision

help: python Autorun_MCTS.py --help

test_run SHOULD accept 2 arguements: 
    directory to train
    name of training file: located in the directory to train

* the script will automatically update the config during offline training for different train test split portion,
* due to possible non uniqueness of data in data set

The following parameters can be changed:
    train_test_split = 0.7 # size of training set
    split = 10 # number of train test splits
    perplexity_penalty = '1e-8'
    samples_count = 1000
"""

import argparse
import shlex, subprocess, multiprocessing
import random
import sys
from shutil import copyfile
import os
import errno
# from nltk.translate.bleu_score import sentence_bleu
cores_count = multiprocessing.cpu_count()

# default values can be changed ##########
train_test_split = 0.289351852 # size of training set 1000/3456
split = 6 # number of train test splits
perplexity_penalty = 1e-8
samples_count = 10000
random.seed(a='vcla')
select_tree_count = 10 # incremental
# weighted_bleu = True
######################################


# global variables not to be altered
name = ""
test_run = ""
test_eval = ""
dataFile = ""
working_dir = ""

process_count = None

recall_flag = False
perplexity_flag = False
precision_flag = False
bleu_flag = False

references = None
bleu_scores = None
recall_scores = None

skip_train = False
offline = True

# offline
total_tree = None
top_tree = None
top_tree_with_selected = None
tree_index_map = None

# incremental
folder_tree_count_max = None
folder_tree_counts = None
train_size_max = None
train_sizes = None
select_tree_indicess = None



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


def Check_Offline():
    global offline
    with open('{}config.txt'.format(working_dir), 'r') as config_f:
        for line in config_f:
            if "OFFLINE=0" in line:
                offline = False
                print("offline: {}".format(offline))
                break


def Generate_Train_Test_Files():
    global train_size_max, train_sizes
    if not offline:
        train_sizes = [0] * split
        
    s = set()
    l_unique = []
    l = []
    with open(dataFile, 'r') as f:
        for line in f:
            if line not in s:
                s.add(line)
                l_unique.append(line)
            l.append(line)
    len_unique = len(l_unique)
    print("All data count: {}".format(len(l)))
    print("Unique data: {}\n".format(len_unique))

    test_set_size = round((1 - train_test_split) * len_unique)
    
    for i in range(1, split+1):
        print("split {}:".format(i))
        test_set = set(random.sample(l_unique, k=test_set_size))

        try:
            os.makedirs('{}{}_{}'.format(working_dir, name, i))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise

        f_train = open('{}{}_{}/train.txt'.format(working_dir, name, i), 'w+')
        f_test = open('{}{}_{}/test.txt'.format(working_dir, name, i), 'w+')
        
        train_count = 0
        train_tokens = set()
        for data in l:
            if not data in test_set:
                train_count += 1
                f_train.write(data)
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
                    f_test.write(data)
        
        print("\ttraining size: {}, test size: {}, training set to train-test-combined: {}".format(train_count, test_count, train_count / (train_count + test_count)))
        
        if not offline:
            train_sizes[i-1] = train_count
            if train_size_max is None:
                train_size_max = train_count
            elif train_count > train_size_max:
                train_size_max = train_count

        f_train.close()
        f_test.close()


def Copy_Config_And_Update():
    config = open('{}config.txt'.format(working_dir), 'r')
    lines = config.readlines()
    config.close()
    
    if not offline:
        for i in range(1, split+1):
            copyfile('{}config.txt'.format(working_dir), '{}{}_{}/config.txt'.format(working_dir, name, i))
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
        
        for i in range(1, split+1):
            # adapt the MEMORY_LIMIT and START_EXPAND arguments to newly splitted train.txt
            count = 0
            with open('{}{}_{}/train.txt'.format(working_dir, name, i), 'r') as f:
                for line in f:
                    count += 1
            
            lines[memory_limit_line] = 'MEMORY_LIMIT={}\n'.format(count)
            lines[start_expand_line] = 'START_EXPAND={}\n'.format(count)

            config = open('{}{}_{}/config.txt'.format(working_dir, name, i), 'w+')
            config.write(''.join(lines))
            config.close()


def Find_Tree_Count():
    global folder_tree_counts, folder_tree_count_max
    if offline:
        with open('{}config.txt'.format(working_dir), 'r') as config_f:
            for line in config_f:
                if offline and 'SEARCH_RESOURCE' in line:
                    learned_tree_index = int(line.rstrip().lstrip('SEARCH_RESOURCE='))
                    global total_tree
                    total_tree = learned_tree_index + 1
                    global top_tree
                    if not top_tree or top_tree > learned_tree_index:
                        top_tree = learned_tree_index
                    global top_tree_with_selected
                    top_tree_with_selected = top_tree + 1
                    break
    elif skip_train: # incremental and skip_train
        assert((train_size_max is None) and (train_sizes is None) and (folder_tree_count_max is None) and (folder_tree_counts is None))
        folder_tree_counts = [0] * split
        for i in range(1, split+1):
            tree_index_max = 0
            dir_list = os.listdir('{}{}_{}'.format(working_dir, name, i))
            for item in dir_list:
                if 'learned_tree_' in item:
                    try:
                        tree_index = int(item.lstrip('learned_tree_').rstrip('.txt'))
                    except TypeError:
                        error_msg = "[Error] learned tree files format: 'learned_tree_[0-9]+.txt'"
                        print(error_msg)
                        raise Exception(error_msg)
                    if tree_index > tree_index_max:
                        tree_index_max = tree_index
            folder_tree_counts[i-1] = tree_index_max
            # for j in range(1, folder_tree_count_max+1):
            #     tree_name = 'learned_tree_{}.txt'.format(j)
            #     if tree_name not in dir_list:
            #         error_msg = "[Error] {} is not in folder {}: tree indices need to be consecutive".format(tree_name, str(i))
            #         print(error_msg)
            #         raise Exception(error_msg)
        folder_tree_count_max = max(folder_tree_counts)
    else: # incremental and not skip_train
        with open('{}config.txt'.format(working_dir), 'r') as config_f:
            for line in config_f:
                if 'START_EXPAND' in line:
                    start_expand = int(line.rstrip().lstrip('START_EXPAND='))
                    folder_tree_count_max = train_size_max - start_expand + 1
                    folder_tree_counts = [count - start_expand + 1 for count in train_sizes]
                    break


    if not offline:
        print("folder_tree_count: {}".format(folder_tree_counts))
        global select_tree_indicess
        assert(folder_tree_counts != None)
        for i in range(1, split+1):
            folder_tree_count = folder_tree_counts[i-1]
            if select_tree_count <= folder_tree_count:
                interval = (folder_tree_count - 1) / (select_tree_count - 1)
                select_tree_indices = [round(z * interval) + 1 for z in range(select_tree_count)]
            else:
                select_tree_indices = list(range(1, folder_tree_count+1))
            select_tree_indicess[i-1] = select_tree_indices
    
        print("select_tree_indicess: {}".format(select_tree_indicess))


def bash(command, logFile):
    with open('{}{}'.format(working_dir, logFile), 'w+') as log_f:
        p = subprocess.Popen(shlex.split(command), stdout=log_f, stderr=log_f, cwd=working_dir)
        p.wait()
        log_f.flush()


def Run_Command(i):
    if not skip_train:
        logFile_train = '{}_{}/log_train'.format(name, i)
        command_train = '{exe} {name}_{iter}/ train.txt'.format(exe=test_run, name=name, iter=i)
        bash(command_train, logFile_train)

    dir_list = os.listdir('{}{}_{}'.format(working_dir, name, i)) 

    # preprocessing ##############################################################################################################
    if offline:
        if not perplexity_flag:
            tree_post_m = dict()
            with open('{}{}_{}/all_posteriors.txt'.format(working_dir, name, i), 'r') as log_f:
                # next(log_f)
                for line in log_f:
                    if "Root" in line:
                        continue
                    line = line.split()
                    index = int(line[1])
                    posterior = float(line[3])
                    if index in tree_post_m:
                        raise Exception("same tree index found")
                    tree_post_m[index] = posterior

            # get top_tree number of trees with the largest posteriors in ascending tree index
            tree_post_sorted = [x[0] for x in sorted(tree_post_m.items(), key=lambda x: x[1], reverse=True)[:top_tree]]
            tree_selected_index = tree_post_sorted[0]
            tree_index_sorted = sorted(tree_post_sorted)
            tree_index_map[(i - 1) * top_tree_with_selected] = tree_selected_index
            for j, tree_index in enumerate(tree_index_sorted, start=1):
                tree_index_map[(i - 1) * top_tree_with_selected + j] = tree_index

            print('{}:'.format(i))
            print(tree_post_sorted)
            print(tree_index_sorted)
        
    else:
        # done in Find_Tree_Count()
        global select_tree_indicess
        select_tree_indices = select_tree_indicess[i-1]
        

    ##############################################################################################################################

    if recall_flag:
        if offline:
            for j in range(0, top_tree_with_selected):
                appendix = '' if j == 0 else '_' + str(tree_index_sorted[j - 1])
                if ('learned_tree{}.txt'.format(appendix) in dir_list):
                    logFile_recall = '{}_{}/log_recall{}'.format(name, i, appendix)
                    command_recall = '{exe} recall {name}_{iter}/learned_tree{appendix}.txt {name}_{iter}/test.txt'.format(exe=test_eval, name=name, iter=i, appendix=appendix)
                    bash(command_recall, logFile_recall)
        else:
            for j in select_tree_indices:
                # learned_tree_j should be in the folder, as tested in Find_Tree_Count()
                logFile_recall = '{}_{}/log_recall_{}'.format(name, i, j)
                command_recall = '{exe} recall {name}_{iter}/learned_tree_{j}.txt {name}_{iter}/test.txt'.format(exe=test_eval, name=name, iter=i, j=j)
                bash(command_recall, logFile_recall)

    if perplexity_flag:
        logFile_perplexity = '{}_{}/log_perplexity'.format(name, i)
        command_perplexity = '{exe} perplexity {name}_{iter}/learned_tree.txt {name}_{iter}/test.txt {penalty}'.format(exe=test_eval, name=name, iter=i, penalty=perplexity_penalty)
        bash(command_perplexity, logFile_perplexity)
    
    if precision_flag:
        logFile_precision = '{}_{}/log_precision'.format(name, i)
        command_precision = '{exe} precision {name}_{iter}/learned_tree.txt learned_tree.txt'.format(exe=test_eval, name=name, iter=i)
        bash(command_precision, logFile_precision)

    if bleu_flag:
        if offline:
            # tree_selected_bleu = None
            for j in range(0, top_tree_with_selected):
                appendix = '' if j == 0 else '_' + str(tree_index_sorted[j - 1])
                if ('learned_tree{}.txt'.format(appendix) in dir_list):
                    learned_tree = "{}_{}/learned_tree{}.txt".format(name, i, appendix)
                    command_sample = '{exe} sample {tree} {count} 1'.format(exe=test_eval, tree=learned_tree, count=samples_count)
                    samples_file = '{}{}_{}/samples{}.txt'.format(working_dir, name, i, appendix)
                    with open(samples_file, "w+") as sample_f:
                        p = subprocess.Popen(shlex.split(command_sample), stdout=sample_f, cwd=working_dir)
                        p.wait()
                        
                    command_bleu = '{} bleu {} {}'.format(test_eval, dataFile, samples_file)
                    log_file = '{}_{}/log_bleu{}'.format(name, i, appendix)
                    bash(command_bleu, log_file)
        
        else:
            for j in select_tree_indices:
                learned_tree = "{}_{}/learned_tree_{}.txt".format(name, i, j)
                command_sample = '{exe} sample {tree} {count} 1'.format(exe=test_eval, tree=learned_tree, count=samples_count)
                samples_file = '{}{}_{}/samples_{}.txt'.format(working_dir, name, i, j)
                with open(samples_file, "w+") as sample_f:
                    p = subprocess.Popen(shlex.split(command_sample), stdout=sample_f, cwd=working_dir)
                    p.wait()
                    
                command_bleu = '{} bleu {} {}'.format(test_eval, dataFile, samples_file)
                log_file = '{}_{}/log_bleu_{}'.format(name, i, j)
                bash(command_bleu, log_file)

        # python nltk ####################################################################################
        # with open(samples_file, 'r') as sample_f:
        #     bleu_weighted = None
        #     if j != tree_selected_index:
        #         sum = 0
        #         total_prob = 0 
        #         for line in sample_f:
        #             line = line.rstrip().split(': ')
        #             prob = float(line[0])
        #             seq = line[1].split()
        #             len_seq = len(seq)
        #             bleu = sentence_bleu(references, seq) if len_seq >= 4 else sentence_bleu(references, seq, weights=[1/len_seq for i in range(len_seq)])
        #             sum += prob * bleu
        #             total_prob += prob
        #         bleu_weighted = sum/total_prob
        #     else:
        #         bleu_weighted = tree_selected_bleu
        #     if j == 0:
        #         tree_selected_bleu = bleu_weighted
        #     bleu_scores[(i - 1) * (top_tree_with_selected) + j] = bleu_weighted
        ##################################################################################################


def Operation_Summary(operation):
    operation_log = '{}log_{}'.format(working_dir, operation)
    summary_file = '{}{}_summary.txt'.format(working_dir, operation_log)
    with open(summary_file, 'w+') as log_f:
        sum = 0
        for i in range(1, split+1):
            logFile = '{}{}_{}/{}'.format(working_dir, name, i, operation_log)
            for line in readlines_reverse(logFile):
                if (operation in line):
                    break
            log_f.write('{}_{} - '.format(name, i))
            log_f.write(line + '\n')
            sum += float(line.lstrip(operation + ': '))
        average = sum / split
        log_f.write('\n')
        log_f.write('average: {}\n'.format(average))


def main():

    global test_run, test_eval, dataFile, working_dir, skip_train, name, split, process_count, train_test_split, perplexity_penalty, samples_count, top_tree, select_tree_count
    # for evaluation
    global select_tree_indicess


    parser = argparse.ArgumentParser(description='a python script to run training and evluation for aog_lib')
    parser.add_argument('operation', choices=['recall', 'bleu', 'both', 'perplexity', 'all', 'precision'], \
                        help="evaluations to be performed (both: recall & bleu, all: recall, bleu & perplexity, precision: assume skip training)")
    parser.add_argument('datafile', help='location of the whole data set')
    parser.add_argument('config', help='location of the config file')
    parser.add_argument('--skip_train', action='store_true', default=False, help='a flag to skip training process')
    parser.add_argument('--split', type=int, default=split, help='number of random splits of the data set')
    parser.add_argument('--process_count', type=int, help='max number of subprocesses running in parallel')
    parser.add_argument('--train_test_split', type=float, default=train_test_split, help='ratio of the training set')
    parser.add_argument('--perplexity_penalty', type=float, default=perplexity_penalty)
    parser.add_argument('--samples_count', type=int, default=samples_count, help='number of samples for bleu score')
    parser.add_argument('--select_tree_count', type=int, help='1. [offline] number of trees with top posteriors to be evaluated, 2. [incremental] equally divided all trees indices into continuous parts and select trees with indices at conjuctures')

    args = parser.parse_args()

    eval_operation = args.operation
    aog_lib_dir = os.path.abspath(os.path.abspath(__file__) + '/../..') + '/'
    test_run = '{}test_run'.format(aog_lib_dir)
    test_eval = '{}Evaluation/test_eval'.format(aog_lib_dir)
    dataFile = os.path.abspath(args.datafile)
    print(dataFile)
    name = os.path.splitext(os.path.basename(dataFile))[0]
    working_dir = os.path.dirname(os.path.abspath(args.config)) + '/'
    print('working directory: {}'.format(working_dir))
    skip_train = args.skip_train
    split = args.split
    if not args.process_count:
        process_count = split
    else:
        process_count = args.process_count
    process_count = process_count if process_count <= (cores_count - 1) else (cores_count - 1)
    train_test_split = args.train_test_split
    perplexity_penalty = args.perplexity_penalty
    samples_count = args.samples_count
    if args.select_tree_count:
        top_tree = args.select_tree_count
        select_tree_count = args.select_tree_count
    select_tree_indicess = [list()] * split
    

    if eval_operation in ['recall', 'both', 'all']:
        global recall_flag
        recall_flag = True
    if eval_operation in ['perplexity', 'all']:
        global perplexity_flag
        perplexity_flag = True
    if eval_operation in ['bleu', 'both', 'all']:
        global bleu_flag
        bleu_flag = True
    if eval_operation == 'precision':
        global precision_flag
        precision_flag = True
        # global skip_train
        skip_train = True
    
    Check_Offline()

    if not skip_train:
        Generate_Train_Test_Files()
        Copy_Config_And_Update()
    
    if recall_flag or bleu_flag:
        Find_Tree_Count()
        if offline:
            global tree_index_map
            tree_index_map = multiprocessing.Array('i', [-1] * (split * (top_tree_with_selected)))
        
    # if bleu_flag:
    #     global bleu_scores
    #     bleu_scores = multiprocessing.Array('d', [-1] * (split * top_tree_with_selected))
    #     with open(dataFile, 'r') as data_f:
    #         global references
    #         references = [line.rstrip().split() for line in data_f]

    # train, recall, perplexity
    pool = multiprocessing.Pool(process_count)
    for i in range(1, split+1):
        pool.apply_async(Run_Command, args=(i,))
    pool.close()
    pool.join()

    # Summary Output
    if recall_flag:
        if offline:
            summary_file = '{}log_recall_summary.txt'.format(working_dir)
            with open(summary_file, 'w+') as log_f:
                global recall_scores
                recall_scores = [-1] * (split * top_tree_with_selected)
                sum = [0.0] * top_tree_with_selected
                skip_average = [False] * top_tree_with_selected
                for i in range(1, split+1):
                    dir_list = os.listdir('{}{}_{}'.format(working_dir, name, i))
                    for j in range(0, top_tree_with_selected):
                        appendix = '' if j == 0 else '_' + str(tree_index_map[(i - 1) * top_tree_with_selected + j])
                        if 'log_recall{}'.format(appendix) in dir_list:
                            logFile = '{}{}_{}/log_recall{}'.format(working_dir, name, i, appendix)
                            for line in readlines_reverse(logFile):
                                if ('recall' in line):
                                    break
                            score = float(line.lstrip('recall: '))
                            log_f.write('split_{}: learned_tree{}: {}\n'.format(i, appendix, score))
                            sum[j] += score
                            recall_scores[(i - 1) * top_tree_with_selected + j] = score
                        else:
                            skip_average[j] = True
                    log_f.write('\n')
                    
                average = [sum[i] / split if not skip_average[i] else -1 for i in range(top_tree_with_selected)]
                for j in range(top_tree_with_selected):
                    if not skip_average[j]:
                        appendix = ' of CHOSEN learned trees' if j == 0 else '_' + str(j)
                        log_f.write('average{}: {}\n'.format(appendix, average[j]))
        
        else:
            summary_file = '{}log_recall_summary.txt'.format(working_dir)
            recalls = [([-1] * select_tree_count) for z in range(split)]
            recall_sum = [0] * select_tree_count
            recall_tree_count = [0] * select_tree_count
            average_last_tree_index = 0
            with open(summary_file, 'w+') as log_f:
                for i in range(1, split+1):
                    dir_list = os.listdir('{}{}_{}'.format(working_dir, name, i))
                    select_tree_indices = select_tree_indicess[i-1]
                    assert(len(select_tree_indices) == select_tree_count)
                    average_last_tree_index += select_tree_indices[-1] / split
                    for tree_order, j in enumerate(select_tree_indices):
                        log_recall_file_name = 'log_recall_{}'.format(j)
                        if log_recall_file_name in dir_list:
                            logFile = '{}{}_{}/log_recall_{}'.format(working_dir, name, i, j)
                            for line in readlines_reverse(logFile):
                                if ('recall' in line):
                                    break
                            score = float(line.lstrip('recall: '))
                            log_f.write('split_{}: learned_tree_{}: {}\n'.format(i, j, score))
                            recalls[i-1][tree_order] = score
                            recall_sum[tree_order] += score
                            recall_tree_count[tree_order] += 1
                        else:
                            log_f.write('[Error]: {} is not present\n'.format(log_recall_file_name))
                    log_f.write('\n')
                log_f.write('one epoch: training around {:.2f} data\n'.format((average_last_tree_index - 1) / (select_tree_count - 1)))

                # recall_average used in f_measure_summary
                recall_average = [recall_sum[i] / recall_tree_count[i] for i in range(select_tree_count)]
                for i in range(select_tree_count):
                    log_f.write('epoch {}: {}\n'.format(i, recall_average[i]))

    
    if perplexity_flag:
        summary_file = '{}log_perplexity_summary.txt'.format(working_dir)
        with open(summary_file, 'w+') as log_f:
            sum = 0.0
            for i in range(1, split+1):
                logFile = '{}{}_{}/log_perplexity'.format(working_dir, name, i)
                for line in readlines_reverse(logFile):
                    if ('perplexity' in line):
                        break
                log_f.write('{}_{} - '.format(name, i))
                log_f.write(line + '\n')
                sum += float(line.lstrip('perplexity: '))
            average = sum / split
            log_f.write('\n')
            log_f.write('average: {}\n'.format(average))


    if bleu_flag:
        # python nltk bleu #######################################################################################################
        # summary_file = '{}log_bleu_summary.txt'.format(working_dir)
        # with open(summary_file, 'w+') as log_f:
        #     sum = [0.0] * top_tree_with_selected
        #     skip_average = [False] * top_tree_with_selected
        #     for i in range(split * top_tree_with_selected):
        #         score = bleu_scores[i]
        #         tree_index = 0 if i % top_tree_with_selected == 0 else tree_index_map[i]
        #         if score != -1:
        #             sum[i % top_tree_with_selected] += score
        #             folder_index = (i // top_tree_with_selected) + 1
        #             log_f.write('split_{}: learned_tree{}: {}\n'.format(folder_index, '' if tree_index == 0 else '_' + str(tree_index), score))
        #         else:
        #             skip_average[tree_index] = True
        #         if i % top_tree_with_selected == top_tree_with_selected - 1:
        #             log_f.write('\n')

        #     average = [sum[i] / split if not skip_average[i] else -1 for i in range(top_tree_with_selected)]
        #     for j in range(top_tree_with_selected):
        #         if not skip_average[j]:
        #             appendix = ' of CHOSEN learned trees' if j == 0 else '_' + str(j)
        #             log_f.write('average{}: {}\n'.format(appendix, average[j]))
        ##########################################################################################################################

        summary_file = '{}log_bleu_summary.txt'.format(working_dir)
        with open(summary_file, 'w+') as log_f:
            if offline:
                global bleu_scores
                bleu_scores = [-1] * (split * top_tree_with_selected)
                sum = [0.0] * top_tree_with_selected
                skip_average = [False] * top_tree_with_selected
                for i in range(1, split+1):
                    dir_list = os.listdir('{}{}_{}'.format(working_dir, name, i))
                    for j in range(0, top_tree_with_selected):
                        appendix = '' if j == 0 else '_' + str(tree_index_map[(i - 1) * top_tree_with_selected + j])
                        if 'log_bleu{}'.format(appendix) in dir_list:
                            logFile = '{}{}_{}/log_bleu{}'.format(working_dir, name, i, appendix)
                            for line in readlines_reverse(logFile):
                                if ('bleu' in line):
                                    break
                            score = float(line.lstrip('bleu: '))
                            log_f.write('split_{}: learned_tree{}: {}\n'.format(i, appendix, score))
                            sum[j] += score
                            bleu_scores[(i - 1) * top_tree_with_selected + j] = score
                        else:
                            skip_average[j] = True
                    log_f.write('\n')
                    
                average = [sum[i] / split if not skip_average[i] else -1 for i in range(top_tree_with_selected)]
                for j in range(top_tree_with_selected):
                    if not skip_average[j]:
                        appendix = ' of CHOSEN learned trees' if j == 0 else '_' + str(j)
                        log_f.write('average{}: {}\n'.format(appendix, average[j]))
            
            else:
                summary_file = '{}log_bleu_summary.txt'.format(working_dir)
                bleus = [([-1] * select_tree_count) for z in range(split)]
                bleu_sum = [0] * select_tree_count
                bleu_tree_count = [0] * select_tree_count
                average_last_tree_index = 0
                with open(summary_file, 'w+') as log_f:
                    for i in range(1, split+1):
                        dir_list = os.listdir('{}{}_{}'.format(working_dir, name, i))
                        select_tree_indices = select_tree_indicess[i-1]
                        assert(len(select_tree_indices) == select_tree_count)
                        average_last_tree_index += select_tree_indices[-1] / split
                        for tree_order, j in enumerate(select_tree_indices):
                            log_bleu_file_name = 'log_bleu_{}'.format(j)
                            if log_bleu_file_name in dir_list:
                                logFile = '{}{}_{}/log_bleu_{}'.format(working_dir, name, i, j)
                                for line in readlines_reverse(logFile):
                                    if ('bleu' in line):
                                        break
                                score = float(line.lstrip('bleu: '))
                                log_f.write('split_{}: learned_tree_{}: {}\n'.format(i, j, score))
                                bleus[i-1][tree_order] = score
                                bleu_sum[tree_order] += score
                                bleu_tree_count[tree_order] += 1
                            else:
                                log_f.write('[Error]: {} is not present\n'.format(log_bleu_file_name))
                        log_f.write('\n')
                    log_f.write('one epoch: training around {:.2f} data\n'.format((average_last_tree_index - 1) / (select_tree_count - 1)))

                    # bleu_average used in f_measure_summary
                    bleu_average = [bleu_sum[i] / bleu_tree_count[i] for i in range(select_tree_count)]
                    for i in range(select_tree_count):
                        log_f.write('epoch {}: {}\n'.format(i, bleu_average[i]))
    

    if bleu_flag and recall_flag:
        if offline:
            summary_file = '{}f_measure_summary.txt'.format(working_dir)
            with open(summary_file, 'w+') as f_measure_f:
                for i in range(split * top_tree_with_selected):
                    recall_score = recall_scores[i]
                    bleu_score = bleu_scores[i]
                    if recall_score != -1 and bleu_score != -1:
                        f_measure_score = 2 * (bleu_score * recall_score / (bleu_score + recall_score))
                        folder_index = (i // top_tree_with_selected) + 1
                        f_measure_f.write('split_{}: learned_tree{}: {}\n'.format(folder_index, '' if i % top_tree_with_selected == 0 else '_' + str(tree_index_map[i]), f_measure_score))
                    if i % top_tree_with_selected == top_tree_with_selected - 1:
                        f_measure_f.write('\n')
        else:
            summary_file = '{}f_measure_summary.txt'.format(working_dir)
            with open(summary_file, 'w+') as f_measure_f:
                for i in range(1, split+1):
                    for j in range(select_tree_count):
                        recall_score = recalls[i-1][j]
                        bleu_score = bleus[i-1][j]
                        f_measure = 2 * (recall_score * bleu_score) / (recall_score + bleu_score) if recall_score != -1 or bleu_score != -1 else -1
                        f_measure_f.write('split_{}: learned_tree_{}: {}\n'.format(i, select_tree_indicess[i-1][j], f_measure))
                    f_measure_f.write('\n')



if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Must use Python 3")
    main()
