######## This script takes 1) the absolute path to phylum superfolder and 2) the absolute path to
######## the parallel-fastq-dump bin folder adn 3) n of cores for parallel-fqd.
######## Conda environment must be activated beforehand.

import sys
import os
import time 

current_time = time.ctime(time.time())
print('Run Started: ' + current_time)

#absolute path to superfolder
superfolder = str(sys.argv[1])
if superfolder[-1] != '/':
    superfolder += '/'

#absolute path to sra toolkit bin directory
pfqd_path = str(sys.argv[2])
if pfqd_path[-1] != '/':
    pfqd_path += '/'

threads = str(sys.argv[3])

# get all phylum folders
phyl_folders = []
for x in range(len(sorted(os.listdir(superfolder)))):
    if os.path.isdir(superfolder + os.listdir(superfolder)[x]):
        phyl_folders.append(os.listdir(superfolder)[x])
phyl_folders = sorted(phyl_folders)

for phyldir in phyl_folders:
    # build path to subfolder
    phylpath = os.path.join(superfolder, phyldir)
    for spdir in sorted(os.listdir(phylpath)):
        if os.path.isdir(os.path.join(phylpath, spdir)):
            os.chdir(os.path.join(phylpath, spdir))
            for accfile in os.listdir(os.path.join(phylpath, spdir)):
                if '.acc' in accfile:
                    for acc in open(os.path.join(phylpath, spdir, accfile)):
                        start = time.time()
                        acc = acc.strip('\n')
                        if str(acc + '_1.fastq.gz') in os.listdir(os.path.join(phylpath, spdir)):
                            print(acc + ' was already downloaded.')
                            continue
                        #download previously partially not downloaded fastq.gz files
                        print(spdir)
                        os.system(pfqd_path + 'prefetch ' + acc + ' --max-size 30G ;' + pfqd_path + 'fasterq-dump ' + acc + " --split-files --threads " + threads + ' ;'
                        + 'rm -r ' + acc + ';' + 'pigz ' + acc + '* ' + '-p' + threads)
                        end = time.time()
                        elapsed = end - start
                        print('Time elapsed with ' + acc + ': ' + str(elapsed/60) + ' minutes.')


current_time = time.ctime(time.time())
print('Run Ended: ' + current_time)
