import random
import sys


def main():
    if len(sys.argv) != 3:
        print("usage: python Split.py [file_to_be_splitted] [test_set_size_ratio]")
        exit()
    fileName = sys.argv[1]
    test_ratio = float(sys.argv[2])
    if test_ratio <= 0 or test_ratio >= 1:
        print("test_set_size_ratio should be in the range (0,1)")
        exit()

    s = set()
    l = []
    with open(fileName, 'r') as f:
        for line in f:
            s.add(line)
            l.append(line)
    len_s = len(s)
    print("All data count: {}".format(len(l)))
    print("Unique data: {}".format(len_s))
    
    for i in range(20):
        test_set = set(random.sample(s, k=round(test_ratio*len_s)))
        f_train = open('{}_train_{}.txt'.format(fileName, i+1), 'w+')
        f_test = open('{}_test_{}.txt'.format(fileName, i+1), 'w+')
        for data in l:
            if data in test_set:
                f_test.write(data)
            else:
                f_train.write(data)
        f_train.close()
        f_test.close()
    

if __name__ == "__main__":
    if sys.version_info[0] < 3:
        print("Must use Python 3")
        exit()
    main()