import subprocess
import shlex
import itertools as it

split_ratio = [0.2904]
split = [1000]
threshold = [-1.9, -1.1]
alpha = [0.4, 0.6]
res = [5]
stop = [5, 7, 10]
lines = None
with open('config.txt', 'r') as f:
    lines = ''.join(f.readlines())

p = subprocess.Popen(shlex.split('mkdir offline_test'))
p.wait()

for i in it.product(res, split, alpha, threshold,stop):
    r, s, a, thr, ear = i
    dir = 'offline_test/res{}_{}_alpha{}_thres{}_stop{}/'.format(i[0], i[1], i[2], i[3], i[4])
    command = 'mkdir {}'.format(dir)
    p = subprocess.Popen(shlex.split(command))
    p.wait()
    command = 'cp ../aog_lib/test_run ../aog_lib/Evaluation/test_eval {}_mcts.py {}'.format(s, dir)
    p = subprocess.Popen(shlex.split(command))
    p.wait()
    command = 'mv {}/{}_mcts.py {}/Autorun_MCTS.py'.format(dir, s, dir)
    p = subprocess.Popen(shlex.split(command))
    p.wait()
    config = '{}/config.txt'.format(dir)
    with open(config, 'w+') as config_f:
        config_alpha = 'ALPHA_LIKELIHOOD={}\n'.format(a)
        config_thres = 'THRESHOLD={}\n'.format(thr)
        config_res = 'SEARCH_RESOURCE={}\n'.format(r)
        config_stop = 'NEG_GAIN_COUNT={}\n'.format(ear)
        config_f.write(lines + config_alpha + config_thres + config_res + config_stop)
