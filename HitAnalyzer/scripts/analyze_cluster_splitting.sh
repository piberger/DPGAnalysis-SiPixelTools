#!/bin/bash
echo "python p.py /eos/cms/$1 ..."
lumisections=`python p.py /eos/cms/$1 | awk '{print $3}' | tr '\n' ' '`
for i in $lumisections; do
echo "process $1 LS ${i}";
cmsRun runRawToDigi_old_cfg.py $1 $i;
done
echo 'done.'

