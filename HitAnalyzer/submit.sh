#!/bin/bash
RUN=$1
DATASET=`dasgoclient -query="dataset dataset=/ZeroBias/*/RAW run=${RUN}"`
MEMUSAGE="2G"
QUEUE="all.q"
LOGDIR="/t3home/berger_p2/pixel/logs/"
REDIRECTOR="root://cmsxrootd.fnal.gov/"
j=0
for i in `dasgoclient -query="file dataset=$DATASET run=$RUN"`; do j=$((j+1)); qsub -V -cwd -q ${QUEUE} -l h_vmem=${MEMUSAGE} -N RUN$RUN-job$j -j y -o ${LOGDIR}/${RUN}-${j}.txt -pe smp 1 run.sh ${REDIRECTOR}/${i}; done

