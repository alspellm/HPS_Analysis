#!/bin/bash

#for dir in /u/re/alspellm/work/ana/qual/mc/ap/link/*
for dir in /nfs/slac/g/hps_data2/users/bravo/mc/det16/run16/ap/dis/*
do
    outd=$(basename $dir)
    echo $outd
    rmdnr=$(expr $outd % 10)
    echo $rmndr
    if [ "$rmdnr" -eq 0 ]; then
        python scripts/submit_jobs.py --outdir ~/work/ana/qual/mc/ap/$outd --indir $dir --step=hipster --fileExt root --nevents -1 --queue short --isData 0 --year=2016 --verbose --hpstrFolder ~/work/src/hpstrMaster/ --hpstrCfg ~/work/src/hpstrMaster/processors/config/anaVtxTuple_cfg.py --submit
    fi
done 
