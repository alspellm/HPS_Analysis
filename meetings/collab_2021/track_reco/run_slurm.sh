#!/bin/sh
hps-mc-batch slurm -o -q shared -r 1:46 -m 5000 -E /sdf/group/hps/users/alspellm/src/hps-mc/install/bin/hps-mc-env.sh -W 5 -d /scratch/alspellm/hps_14552 -c /sdf/group/hps/users/alspellm/projects/collab_2021/track_reco/.hpsmc -l /scratch/alspellm/hps_14552/logs hpstr /sdf/group/hps/users/alspellm/projects/collab_2021/track_reco/jobs.json
