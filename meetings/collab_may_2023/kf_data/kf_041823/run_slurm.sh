#!/bin/sh
hps-mc-batch slurm -q hps -W 1 -m 4000 -E /sdf/home/a/alspellm/.bashrc -o -r 1:1958 -d /scratch/alspellm/ana/kf/ -l /scratch/alspellm/ana/kf/logs ana /sdf/group/hps/users/alspellm/run/data/ana/kf_041823/jobs.json -c /sdf/group/hps/users/alspellm/run/data/ana/kf_041823/.hpsmc 
