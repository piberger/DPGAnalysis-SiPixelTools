#!/bin/bash
eval "dasgoclient -query=\"file dataset=/ExpressPhysics/Commissioning2018-Express-v1/FEVT run="$1"\" > files"$1
eval "split -l5 files"$1" files"$1"-"
eval "for f in `ls -1 files"$1"-* | tr '\n' ' '`; do cat batch.template | sed \"s/{file}/\${f}/g\" > batch\${f}.sh; bsub -R \"pool>3000\" -q 8nh -J \${f} < batch\${f}.sh; done"
