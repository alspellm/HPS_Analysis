#/bin/env python
#from optparse import OptionParser

import os
import argparse
# Parse the command line options.
#parser = OptionParser()
parser = argparse.ArgumentParser(description="The baseConfig options for hpstr. ")
#input directory (functions being used)
parser.add_argument("-F", "--in_dir", type=str, dest="FuncDir",
            help="Sets the script input directory.", metavar="FuncDir", default="/sdf/group/hps/users/epeets/run/resonance_fitting/functions/")

parser.add_argument("-d", "--dir", type=str, dest="outdir",
    help="Sets the script output directory.", metavar="outdir", default="/sdf/group/hps/users/epeets/run/resonance_fitting/sh/")

#setting min and max window ranges for each script
parser.add_argument("-m", "--win_min_range", type=int, dest="win_min_r",
    help="minimum window value range to iterate on", nargs = "+", metavar="win_min_r",default=(30,105,5))
parser.add_argument("-x", "--win_max_range", type=int, dest="win_max_r",
    help="maximum window value range to iterate on", nargs = "+", metavar="win_max_r",default=(170,215,5))


parser.add_argument("-B", "--number_batches", type=int, dest="batch_n",
    help="sets number of batches to make for each window range", metavar="batch_n",default=1)




options = parser.parse_args()


#using a job_id to make each job somewhat unique
job_id = 1



#looping throught input folder and p 
subJobFile = open("%s/subJob_%i%i%i_%i%i%i_%ib.sh"%(options.outdir, options.win_min_r[0], options.win_min_r[1], options.win_min_r[2], options.win_max_r[0], options.win_max_r[1], options.win_max_r[2], options.batch_n), 'w+')


#looping over number of batches selected for a range 
for K in range(options.batch_n):
    #setting outdirectory for each window range being considered 
    outDir = "%s/sh_%i%i%i_%i%i%i_batch_%i"%(options.outdir, options.win_min_r[0],options.win_min_r[1], options.win_min_r[2],options.win_max_r[0],options.win_max_r[1], options.win_max_r[2],K+1)

    if not os.path.exists(outDir):
        os.mkdir(outDir)

    for fun in os.listdir(options.FuncDir):
        func_name = os.path.splitext(os.path.basename(fun))[0]
        if "_out" in func_name:
            continue
        func_file_name = "%s/%s.sh"%(outDir,func_name)


        tempFile = open(func_file_name, "w+")
        tempFile.write("")
        tempFile.close()

        outputFile = open(func_file_name, "a+")
        outputFile.write("#!/usr/bin/scl enable devtoolset-8 -- /bin/bash\n")
        outputFile.write("#SBATCH --ntasks=1\n")
        outputFile.write("#SBATCH --time=72:00:00\n")
        outputFile.write("#SBATCH --mem=2000M\n")
        outputFile.write("#SBATCH --partition=hps\n")
        outputFile.write("#SBATCH --job-name=%s_fitB_job_%i\n"%(func_name, job_id))
        outputFile.write("#SBATCH --output=/scratch/epeets/log/%s_%i_%i_%i.txt\n"%(func_name,options.win_min_r[0],options.win_min_r[1],job_id))
        outputFile.write('python3 /sdf/group/hps/users/epeets/run/resonance_fitting/global_fit_5.py -i /sdf/group/hps/users/epeets/run/resonance_fitting/functions/%s.txt -P /sdf/group/hps/users/epeets/run/resonance_fitting/parameters/%s.txt -B %i -m %i %i %i -x %i %i %i -R %i -Q 60 -d /sdf/group/hps/users/epeets/run/resonance_fitting/functions/%s_out/%i%i%i_%i%i%i_batch_%i/ -o %s.root' % (func_name, func_name, K, options.win_min_r[0], options.win_min_r[1], options.win_min_r[2], options.win_max_r[0], options.win_max_r[1], options.win_max_r[2], job_id, func_name,options.win_min_r[0], options.win_min_r[1], options.win_min_r[2], options.win_max_r[0], options.win_max_r[1], options.win_max_r[2],K+1, func_name))
                            
        outputFile.close()
        job_id += 1
        subJobFile.write("sbatch %s\n"%(func_file_name))
        print('    Wrote: %s' % (func_file_name))
        pass
subJobFile.close()

