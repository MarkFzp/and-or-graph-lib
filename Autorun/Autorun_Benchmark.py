"""
autorun script to test recall or perplexity of benchmark (jAOG in OpenBottle)

usage: python Autorun_Benchmark.py [recall/perplexity/both/bleu/all] [OpenBottle directory]
              [test_construct] [test_eval] [whole data set to be splitted] 
              [(optional) (0 or 1) do evalution without training]

The following parameters can be changed:
    train_test_split = 0.7
    iteration = 10 # number of train test splits
    process_count = 10
    perpelxity_penalty = '1e-8'
    samples_count = 1000
    alpha = 0.5
    context_range = 5
"""

import shlex, subprocess, multiprocessing
import random
import sys
import os
import errno
from nltk.translate.bleu_score import sentence_bleu
cores_count = multiprocessing.cpu_count()

# parameters can be changed ##########
train_test_split = 0.7
iteration = 10 # number of train test splits
process_count = 5
perpelxity_penalty = '1e-8'
samples_count = 500
# weighted_bleu = True
random.seed(a='vcla')

alpha = 0.5
context_range = 5
######################################


# global variables not to be altered
name = ""
openbottle_dir = ""
test_construct = ""
test_eval = ""
dataFile = ""
local_dir = ""
recall_flag = False
perpelexity_flag = False
bleu_flag = False
bleu_scores = None
skip_train = False
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


def Generate_Train_Test_Files_And_Mapping():
    s = set()
    l = []
    l_unique = []
    terminal_set = set()
    with open(dataFile, 'r') as f:
        for line in f:
            if line not in s:
                s.add(line)
                l_unique.append(line)
            l.append(line)
            for terminal in line.rstrip().split():
                terminal_set.add(terminal)
    print("All data count: {}".format(len(l)))
    print("Unique data: {}".format(len(l_unique)))

    mapping_contents = ','.join(str(i) for i in range(1, len(terminal_set)+1)) + '\n' + ','.join(sorted(terminal_set)) + '\n'
    test_set_size = round((1 - train_test_split) * len(l_unique))

    for i in range(1, iteration+1):
        print("split {}:".format(i))
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
        
        f_train.close()
        f_test.close()

        with open('{}_{}/mapping_input.txt'.format(name, i), "w+") as f:
            f.write(mapping_contents)


def Convert_Data():
    for i in range (1, iteration+1):
        mapping_file = '{}_{}/mapping_input.txt'.format(name, i)
        train_file = '{}_{}/train.txt'.format(name, i)
        ''' 
        mapping_file:
            Terminal_id1,Terminal_id2,...
            Terminal_Content1,Terminal_Content2,...
        
        train_file:
            data1\n
            data2\n
        '''
        with open(mapping_file, "r+") as mapping_f:
            # print("Reading grammar from {0}".format(sys.argv[2]))
            Terminals = mapping_f.readline().strip('\n').split(',')
            Contents = mapping_f.readline().strip('\n').split(',')
            # print("Terminal #: {0}".format(len(Terminals)))
            # print(Terminals)
            # print("Contents #: {0}".format(len(Contents)))
            # print(Contents)
            
            assert(len(Contents) == len(Terminals))
            
            count = 0
            position_dict = {}
            for id in Terminals:
                position_dict[id] = count
                count += 1

            with open(train_file, "r+") as train_f:
                data = train_f.readlines()
                content_dict = {}
                # terminal_content : terminal_id
                for j in range(len(Contents)):
                    content_dict[Contents[j]] = Terminals[j]

                with open("{}_{}/converted_data.txt".format(name, i), "w+") as output_datafile:
                    output_datafile.write("{0}\n".format(len(Terminals)))
                    for id in Terminals:
                        output_datafile.write("{0}\n".format(id))
                    
                    output_datafile.write("\n")
                    output_datafile.write("{0}\n".format(len(data)))
                    # row: c1a1 c2a2 c2b1 
                    for row in data:
                        # terminals: ["c1a1","c2a2","c2b1"]
                        terminals = row.strip(" \n").split(' ')
                        # print(terminals)
                        count1 = 1.0
                        for terminal in terminals:
                            # print('terminal: {}'.format(terminal))
                            # print(content_dict[terminal])
                            # print(position_dict[content_dict[terminal]])
                            # print()
                            output_datafile.write("{0} [{1} {2} ] ".format(position_dict[content_dict[terminal]], count1, count1))
                            count1 += 1
                        output_datafile.write("\n")

                with open("{}_{}/mapping_output.txt".format(name, i), "w+") as output_mappingfile:
                    for i in range(len(content_dict)):
                        output_mappingfile.write("{0},{1}\n".format(Contents[i], Terminals[i]))


def Run_Command(i):
    converted_data = "{}_{}/converted_data.txt".format(name, i)
    test_file = "{}_{}/test.txt".format(name, i)
    benchmark_grammar = "{}_{}/benchmark_grammar.txt".format(name, i)
    learned_tree = "{}_{}/learned_tree.txt".format(name, i)
    command_benchmark = "java -cp {dir}grammar_induction/jAOG/src:{dir}grammar_induction/jAOG/src/args4j-2.33.jar aog.learn.bc.GrammarLearner \
	-combinerType aog.app.sequence.SequenceCombiner \
	-bcRelationType aog.app.sequence.FollowingRelation \
	-bgmRelationType aog.app.sequence.FollowingRelation \
	-contextType aog.app.sequence.SequenceContext \
    -contextRange {cr} -alpha {alpha} \
    -input {train} -output {output}".format(dir=openbottle_dir, train=local_dir + converted_data, output=local_dir + benchmark_grammar, cr=context_range, alpha=alpha)
    command_construct = "./{exe} {grammar} {mapping} {output}".format(exe=test_construct, grammar=benchmark_grammar, mapping="{}_{}/mapping_output.txt".format(name, i), output=learned_tree)
    
    if not skip_train:
        # benchmark training
        logFile_benchmark = "{}_{}/log_benchmark".format(name, i)
        with open(logFile_benchmark, 'w+') as log_f:
            p = subprocess.Popen(shlex.split(command_benchmark), stdout=log_f, stderr=log_f)
            p.wait()
        
        # construct_tree_from_benchmark
        logFile_construct = "{}_{}/log_construct".format(name, i)
        with open(logFile_construct, "w+") as log_f:
            p = subprocess.Popen(shlex.split(command_construct), stdout=log_f, stderr=log_f)
            p.wait()

    if recall_flag:
        command_recall = './{exe} recall {tree} {test}'.format(exe=test_eval, tree=learned_tree, test=test_file)
        logFile_recall = "{}_{}/log_recall".format(name, i)
        with open(logFile_recall, "w+") as log_f:
            p = subprocess.Popen(shlex.split(command_recall), stdout=log_f, stderr=log_f)
            p.wait()
    
    if perpelexity_flag:
        command_perplexity = './{exe} perplexity {tree} {test} {penalty}'.format(exe=test_eval, tree=learned_tree, test=test_file, penalty=perpelxity_penalty)
        logFile_perplexity = "{}_{}/log_perplexity".format(name, i)
        with open(logFile_perplexity, "w+") as log_f:
            p = subprocess.Popen(shlex.split(command_perplexity), stdout=log_f, stderr=log_f)
            p.wait()
    
    if bleu_flag:
        command_sample = './{exe} sample {tree} {count} 1'.format(exe=test_eval, tree=learned_tree, count=samples_count)
        samples_file = '{}_{}/samples.txt'.format(name, i)
        with open(samples_file, "w+") as sample_f:
            p = subprocess.Popen(shlex.split(command_sample), stdout=sample_f)
            p.wait()
            sample_f.flush()
            
        command_bleu = './{} bleu {} {}'.format(test_eval, dataFile, samples_file)
        logFile_bleu = '{}_{}/log_bleu'.format(name, i)
        with open(logFile_bleu, "w+") as log_f:
            p = subprocess.Popen(shlex.split(command_bleu), stdout=log_f, stderr=log_f)
            p.wait()
            log_f.flush()
        
        for line in readlines_reverse(logFile_bleu):
            if ('bleu' in line):
                break
        score = float(line.lstrip('bleu: '))
        bleu_scores[i-1] = score

        # with open(samples_file, 'r') as sample_f, open(dataFile, 'r') as data_f:
        #     sum = 0
        #     total_prob = 0
        #     references = [line.rstrip().split() for line in data_f]
        #     for line in sample_f:
        #         line = line.rstrip().split(': ')
        #         prob = float(line[0])
        #         seq = line[1].split()

        #         len_seq = len(seq)
        #         bleu = sentence_bleu(references, seq) if len_seq >= 4 else sentence_bleu(references, seq, weights=[1/len_seq for i in range(len_seq)])
        #         sum += prob * bleu
        #         total_prob += prob
        #     bleu_scores[i-1] = sum / total_prob


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
        # log_f.write('sum: {}\n'.format(sum))
        log_f.write('average: {}\n'.format(average))


def main():
    if (len(sys.argv) != 6 and len(sys.argv) != 7) or (sys.argv[1] not in ['recall', 'perplexity', 'both', 'bleu', 'all']):
        print("usage: python Autorun_Benchmark.py [recall/perplexity/both/bleu/all] [OpenBottle directory] [test_construct] [test_eval] [whole data set to be splitted] [(optional) (0 or 1) do evalution without training]")
        exit()

    global local_dir
    local_dir = os.path.dirname(os.path.realpath(__file__)).rstrip('/') + '/'
    eval_operation = sys.argv[1]
    global openbottle_dir
    openbottle_dir = sys.argv[2].rstrip('/') + '/'
    global test_construct
    test_construct = sys.argv[3]
    global test_eval
    test_eval = sys.argv[4]
    global dataFile
    dataFile = sys.argv[5]
    try:
        global skip_train
        skip_train = bool(int(sys.argv[6]))
    except IndexError:
        pass
    global name
    name = os.path.splitext(os.path.basename(dataFile))[0]
    
    if eval_operation in ['recall', 'both', 'all']:
        global recall_flag
        recall_flag = True
    if eval_operation in ['perplexity', 'both', 'all']:
        global perpelexity_flag
        perpelexity_flag = True
    if eval_operation in ['bleu', 'all']:
        global bleu_flag
        bleu_flag = True
   
    if not skip_train:
        Generate_Train_Test_Files_And_Mapping()
        Convert_Data()

    if bleu_flag:
        global bleu_scores
        bleu_scores = multiprocessing.Array('d', [0.0]*iteration)

    # Benchmark training
    pool = multiprocessing.Pool(process_count)
    for i in range(1, iteration+1):
        pool.apply_async(Run_Command, args=(i, ))
    pool.close()
    pool.join()

    # Summary Output
    if recall_flag:
        Operation_Summary('recall')
    
    if perpelexity_flag:
        Operation_Summary('perplexity')
    
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
