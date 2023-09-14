#!/bin/sh
hps-mc-batch slurm -q hps -W 1 -m 2000 -E /sdf/home/a/alspellm/.bashrc -o -r 1:1000 -d /scratch/alspellm/tritrig_beam/ana -l /scratch/alspellm/tritrig_beam/ana/logs ana /sdf/group/hps/users/alspellm/run/mc/tritrig_beam/kf_ana_05012023/jobs.json -c /sdf/group/hps/users/alspellm/run/mc/tritrig_beam/kf_ana_05012023/.hpsmc 
