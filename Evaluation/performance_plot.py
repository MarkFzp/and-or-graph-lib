import os
import sys
import re
import matplotlib.pyplot as plt
import itertools
import argparse

setting_path = ""
num_of_split = 0
res = 0
num_of_plot = 1
sort = False
plot_vertical=False
def plot():
    f_splits = [[] for i in range(num_of_split)] 
    bleu_splits = [[] for i in range(num_of_split)]
    recall_splits = [[] for i in range(num_of_split)]
    posterior_splits = [[] for i in range(num_of_split)]
    prior_splits = [[] for i in range(num_of_split)]

    splt_num = 1
    splt_str = "split_"
    robotic_str = "robotic_commands_"
    #read in f_score
    with open (setting_path + "f_measure_summary.txt") as fscore:
        lines = fscore.readlines()
        if splt_str not in lines[0]:
            splt_str = "shuffle_"
        for line in lines:
            #put the f_score into corresponding split
            if splt_str + str(splt_num) in line:
                f_splits[splt_num - 1].append(float(line.split(':')[-1]))
            elif splt_str + str(splt_num + 1) in line:
                splt_num += 1
                f_splits[splt_num - 1].append(float(line.split(':')[-1])) 

    

    #read in recall score
    splt_num = 1
    with open (setting_path + "log_recall_summary.txt") as recall:
        lines = recall.readlines()
        for line in lines:
            #put the f_score into corresponding split
            if splt_str + str(splt_num) in line:
                recall_splits[splt_num - 1].append(float(line.split(':')[-1]))
            elif splt_str + str(splt_num + 1) in line:
                splt_num += 1
                recall_splits[splt_num - 1].append(float(line.split(':')[-1])) 

           
                

    
    #read in bleu score
    splt_num = 1
    with open (setting_path + "log_bleu_summary.txt") as bleu:
        lines = bleu.readlines()
        for line in lines:
            #put the f_score into corresponding split
            if splt_str + str(splt_num) in line:
                bleu_splits[splt_num - 1].append(float(line.split(':')[-1]))
            elif splt_str + str(splt_num + 1) in line:
                splt_num += 1
                bleu_splits[splt_num - 1].append(float(line.split(':')[-1])) 

           


    #remove the selected learned tree score
    for f,b,r in itertools.zip_longest(f_splits, bleu_splits, recall_splits):
        f.pop(0)
        b.pop(0)
        r.pop(0)

    #reverse f_score, recall and bleu to the learning order
    for i in range(len(f_splits)):
        f_splits[i].reverse()
        bleu_splits[i].reverse()
        recall_splits[i].reverse()
        
    #read in posterior and prior      
    for i in range(num_of_split):
        with open (setting_path + robotic_str + str(i+1) + "/all_posteriors.txt") as learn_file:
            lines = learn_file.readlines()
            lines.pop(0)
            for line in lines:
                temp = line.split(":")
                posterior_splits[i].append( float(temp[2].split(" ")[1]))
                prior_splits[i].append(float(temp[3].split(" ")[1].rstrip("total")))  

    
    # normalize posterior and prior
    if num_of_plot == 1:
        
        for post_split, prior_split in itertools.zip_longest(posterior_splits, prior_splits):
            post_max = max(post_split)
            post_min = min(post_split)
            prior_max = max(prior_split)
            prior_min = min(prior_split)
            for i in range(len(post_split)):
                post_split[i] = (post_split[i] - post_min ) / (post_max-post_min)
            for i in range(len(prior_split)):
                prior_split[i] = (prior_split[i] - prior_min ) / (prior_max-prior_min)

    #sort according to posterior
    if sort:
        for i in range(num_of_split):
            zipped = list(zip(posterior_splits[i],prior_splits[i],f_splits[i],bleu_splits[i],recall_splits[i]))
            zipped.sort()
            a,b,c,d,e = zip(*zipped)
            posterior_splits[i] = a
            prior_splits[i] = b
            f_splits[i] = c
            bleu_splits[i] = d
            recall_splits[i] = e
        

    
    #plot
    x = range(1,res+1)
    for i in range(num_of_split):
        fig = plt.figure(figsize=(20, 10))
        if num_of_plot == 1:
            axs0 = fig.add_subplot(111)

        elif num_of_plot == 2:
            if plot_vertical:
                axs0 = fig.add_subplot(121)
                axs1 = fig.add_subplot(122,sharex=axs0)
            else:
                axs0 = fig.add_subplot(211)
                axs1 = fig.add_subplot(212,sharex=axs0)
            fig.subplots_adjust(hspace=0)
        else:
            if plot_vertical:                
                axs0 = fig.add_subplot(131)
                axs1 = fig.add_subplot(132,sharex=axs0)
                axs2 = fig.add_subplot(133,sharex=axs0)
            else:
                axs0 = fig.add_subplot(311)
                axs1 = fig.add_subplot(322,sharex=axs0)
                axs2 = fig.add_subplot(333,sharex=axs0)
            fig.subplots_adjust(hspace=0)
        x_coord = [i for i in range(res+1)]
        for xc in x_coord:
            axs0.axvline(x = xc,dashes = (1,4))
            if num_of_plot == 2:
                axs1.axvline(x = xc,dashes=(1,4))
            if num_of_plot == 3:
                axs1.axvline(x = xc,dashes=(1,4))
                axs2.axvline(x = xc,dashes=(1,4))
                
        #plot f-score, recall, bleu for each split
        axs0.plot(x,f_splits[i], 'r.-',label='f-score')
        axs0.plot(x,recall_splits[i],'b+--',label='recall')
        axs0.plot(x,bleu_splits[i],'y+--',label='bleu')
        axs0.set_ylim(0,1)
        axs0.legend(loc="lower right")

        if num_of_plot == 1:
            #plot posterior for each split
            axs0.plot(x,posterior_splits[i],'g.-',label = 'posterior')
            axs0.legend(loc="lower right")
            #plot prior for each split
            axs0.plot(x,prior_splits[i],'.-',color="orange",label = 'prior')
            axs0.legend(loc="lower right")
        elif num_of_plot == 2:               
            #plot posterior for each split
            axs1.plot(x,posterior_splits[i],'g.-',label = 'posterior')
            axs1.legend(loc="lower right")
            #plot prior for each split
            axs1.plot(x,prior_splits[i],'.-',color="orange",label = 'prior')
            axs1.legend(loc="lower right")

        else:
            #plot posterior for each split
            axs1.plot(x,posterior_splits[i],'g.-',label = 'posterior')
            axs1.legend(loc="lower right")
            #plot prior for each split
            axs2.plot(x,prior_splits[i],'.-',color="orange",label = 'prior')
            axs2.legend(loc="lower right")
        # fig = plt.figure(figsize=(20, 10))
        # new_plot = fig.add_subplot(111)
        # new_plot.plot(x,f_splits[i], 'ro-',label='f-score')
        # new_plot.plot(x,recall_splits[i],'b+--',label='recall')
        # new_plot.plot(x,bleu_splits[i],'y+--',label='bleu')
        # new_plot.plot(x,posterior_splits[i],'o-',label = 'posterior')
        # new_plot.plot(x,prior_splits[i],'o-',label = 'prior')
        # plt.ylim(0,1)
        # plt.legend(loc="lower right")
        
        
        plt.savefig(setting_path + "split_"+str(i+1)+"_analysis_sort-{}_{}_verti-{}.png".format(str(sort),str(num_of_plot),str(plot_vertical)))



def main():
    global setting_path, num_of_split, res, num_of_plot, sort,plot_vertical

    #argument parsing
    parser = argparse.ArgumentParser(description='Visualization for analysing results')
    parser.add_argument('setting_path', help = 'the relative path to the experiment setting folder')
    parser.add_argument('--sort', action = "store_true", help = 'whether to sort the posterior of each resource within a split before visualization' )
    parser.add_argument('--num_of_plots',type = int, default = 1,choices=[1,2,3], help= "the number of plots want to visualize(1 or 2)")
    parser.add_argument('--vertical',action="store_true", help= "plot horizontally or vertically")
    args = parser.parse_args()
    setting_path = str(os.path.dirname(os.path.abspath(__file__))) + "/"+args.setting_path
    sort = args.sort
    num_of_plot = args.num_of_plots
    plot_vertical = args.vertical
    dirlist = os.listdir(setting_path)

    #calculate the number of splits
    for f in dirlist:
        if re.match("robotic_commands_[1-9]",f) is not None:
            num_of_split += 1

    #find the resources used
    print(setting_path)
    setting_folder = re.search("res[0-9][0-9]*",setting_path)

    if setting_folder is None:
        print("The setting folder does not specify the number of resources!\n")
        exit(1)
    print(setting_folder.group(0))
    res = int(setting_folder.group(0).replace("res",''))
    print("the number of splits is: {}".format(num_of_split))
    print("the number of res is: {}".format(res))

    #plot 
    plot()

if __name__ == "__main__":
    main()
