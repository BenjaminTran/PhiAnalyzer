#!/bin/bash

rsync -tuPr\
    --exclude=.*\
   /Users/btran/research/PhiAnalyzer/PhiAnalyzer/ btran@lxplus.cern.ch:/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/PhiAnalyzer/PhiAnalyzer/
