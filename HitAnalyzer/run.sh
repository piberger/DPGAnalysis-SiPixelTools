#!/bin/bash
export XRD_NETWORKSTACK=IPv4
VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh 
echo ${CMSSW_BASE}
mkdir /scratch/$USER/clusterNtupleProducer/
hash=`head /dev/urandom | tr -dc A-Za-z0-9 | head -c10`
mkdir /scratch/$USER/clusterNtupleProducer/$hash
cmsenv
cmsRun scripts/runRawToDigi_t3_cfg.py $1 /scratch/$USER/clusterNtupleProducer/$hash
echo "copy: /scratch/$USER/clusterNtupleProducer/$hash/*.root to root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/berger_p2/pixel/clusters/"
xrdcp -f /scratch/$USER/clusterNtupleProducer/$hash/*.root  root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/berger_p2/pixel/clusters/
echo "done."

