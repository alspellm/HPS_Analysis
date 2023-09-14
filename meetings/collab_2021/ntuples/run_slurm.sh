#!/bin/sh
hps-mc-batch slurm -o -q shared -r 1:46 -m 6000 -E /sdf/group/hps/users/alspellm/src/hps-mc/install/bin/hps-mc-env.sh -W 5 -d /scratch/alspellm/test -c /sdf/group/hps/users/alspellm/projects/collab_2021/ntuples/.hpsmc -l /scratch/alspellm/test/logs hpstr /sdf/group/hps/users/alspellm/projects/collab_2021/ntuples/jobs.json
