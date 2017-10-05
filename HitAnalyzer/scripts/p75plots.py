import ROOT
from plot import *
import sys

runlist = [int(x) for x in sys.argv[-1].split(',')] if len(sys.argv) > 1 else [302328]
print runlist
for run in runlist:
    graphOptions = [
        {'run': run, 'split': 2, 'color': ROOT.kBlue+2, 'size': 1.6, 'style': 33, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]', 'limitpoints': 5000},
    ]
    make_tgraph(graphOptions, "p75_run_%d.png"%run, xrange=[0.0,1.8])
