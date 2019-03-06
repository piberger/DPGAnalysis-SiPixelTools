#!/bin/bash
run=`(echo -n $1 |cut -d/ -f 9 && echo -n $1 | cut -d/ -f 10) | tr -d '\n'`
echo ${run}
echo "python p.py /eos/cms/$1 ..."
lumisections=`python p.py /eos/cms/$1 | awk '{print $3}' | tr '\n' ' '`
for i in $lumisections; do
echo "process $1 LS ${i}";
cmsRun runRawToDigi_old_cfg.py $1 $i # > /afs/cern.ch/work/p/piberger/digis_${run}_${i}.log;
done
echo 'done.'

