import matplotlib.pyplot as plt
import sys
import subprocess

def plot_stats(y,ylabel,filename,text=None):
    x = [i for i in range(1,len(y)+1)]
    plt.grid(True)
    plt.plot(x,y)
    plt.ylabel(ylabel)
    plt.xlabel("No. AOF Generated")
    if text is not None:
        plt.title(text)
    plt.savefig(filename)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise ValueError("Arguments needed!\n")
    
    #clear all files in AOG_Sample_Results
    subprocess.run(["rm ../AOG_Sample_Results/*"],shell=True)
    
    #compile test_run
    try:
        subprocess.run(["g++",  "-std=c++11", "-lboost_system", "-lboost_filesystem", "../Tests/test_run.cpp", "-o", "t"],check=True)
    except subprocess.CalledProcessError:
        exit()

    #run test_run to get all stats files
    subprocess.run(["./t" ,sys.argv[1] ,sys.argv[2] ])

    with open("./stats.txt") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            lines[i] = lines[i].rstrip(" \n").split(" ")
    
    
    #parse the filename 
    temp = sys.argv[2].split(".")[0]
    # print(temp)
    temp = temp.split("_")[:-1]
    filename = temp[0] + "_"+temp[1]

    #parse the resource argument
    with open(sys.argv[1]+"config.txt") as f1:
        configs = f1.readlines()
        for config in configs:
            if "SEARCH_RESOURCE" in config:
                resource = config.rstrip("\n").rstrip(" ").split("=")[1]

    store_path = "./{}/{}".format(filename,resource)
    subprocess.run(["mkdir",filename])
    subprocess.run(["mkdir", store_path])

    for i in range(len(lines)):
        lines[i] = list(map(float,lines[i]))
    #plot number of bigrams in each step

    avg_bigram = sum(lines[0]) / len(lines[0])
    bigram_text = "Average number of bigrams: {}".format(avg_bigram)
    plot_stats(lines[0],"Number of Bigrams","{}/{}_bigram_rsc_{}".format(store_path,filename,resource),bigram_text)

    #plot the number of variations
    total_variation = lines[4][0]    
    avg_variation = round(sum(lines[1]) / len(lines[1]),2)
    variation_text = "Total number of variations examined: {}\nAverage number of variations: {}".format(total_variation,avg_variation)
    plot_stats(lines[1],"Number of Variations","{}/{}_variation_rsc_{}".format(store_path,filename,resource),variation_text)

    #plot the best posteriors
    plot_stats(lines[2],"Log Best Posteriors","{}/{}_posterior_rsc_{}".format(store_path,filename,resource))
    
    #calculate KL_divergence
    try:
        subprocess.run(["g++",  "-std=c++11", "../Utils/CalculateKL.cpp", "-o","calKL"],check=True)
    except subprocess.CalledProcessError:
        exit()
    
    
    kl_divs = []
    # truth_path = sys.argv[1]+"../10000_grammar_examples/Output/"+sys.argv[2]    
    truth_path = sys.argv[1]+"../Grammar_Example/Output/"+sys.argv[2]    
    
    for i in range(len(lines[3])):
        name = "sample_" + str(i)
        print(name)
        output = subprocess.run(["./calKL","../AOG_Sample_Results/{}".format(name),truth_path],stdout = subprocess.PIPE)
        kl_divs.append(output.stdout.decode("utf-8").rstrip('\n'))
   

    print(kl_divs)

    for i in range(len(kl_divs)):
        if lines[3][i] == 0:
            kl_divs[i] = kl_divs[i-1]

    kl_divs = list(map(float,kl_divs))
    #plot kl divergence of each round
    print(kl_divs)
    plot_stats(kl_divs,"KL Divergence","{}/{}_kl_divergence_rsc_{}".format(store_path,filename,resource))


    
