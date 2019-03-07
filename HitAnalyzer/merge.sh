#!/bin/bash
RUN=$1
PATHIN=/pnfs/psi.ch/cms/trivcat/store/user/berger_p2/pixel/clusters/
PATHOUT=/pnfs/psi.ch/cms/trivcat/store/user/$USER/pixel/clusterNtuples/
TEMPFILE=/scratch/$USER/Run${RUN}.root
REDIRECTOR=root://t3dcachedb.psi.ch:1094/

for i in `ls -1 ${PATHIN}/*${RUN}*`; do echo "${REDIRECTOR}/${i}"; done | tr '\n' ' ' | xargs hadd -f ${TEMPFILE}; xrdcp ${TEMPFILE} ${REDIRECTOR}/${PATHOUT}

