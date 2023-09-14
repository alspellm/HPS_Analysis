#!/bin/sh
hps-mc-batch slurm -q hps -W 1 -m 2000 -E /sdf/home/a/alspellm/.bashrc -o -r 1:1000 -d /scratch/alspellm/wab_beam/ana -l /scratch/alspellm/wab_beam/ana/logs ana /sdf/group/hps/users/alspellm/run/mc/wab_beam/wab_kf_ana_05022023/jobs.json -c /sdf/group/hps/users/alspellm/run/mc/wab_beam/wab_kf_ana_05022023/.hpsmc 
