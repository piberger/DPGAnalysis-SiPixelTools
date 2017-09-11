import ROOT
from plot import *
import sys

runlist = [int(x) for x in sys.argv[-1].split(',')] if len(sys.argv) > 1 else [302328]
print runlist

#graphOptions = [
#{'run': 299593, 'split': 1, 'color': ROOT.kMagenta+2, 'style': 21, 'yrange': [0, 10], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]'},
#{'run': 299478, 'split': 1, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.4},
#{'run': 299479, 'split': 1, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.4},
#{'run': 299380, 'split': 1, 'color': ROOT.kGreen-1, 'style': 23, 'size': 1.4},
#{'run': 300806, 'split': 1, 'color': ROOT.kRed-3, 'style': 33, 'size': 1.4},
#{'run': 300811, 'split': 1, 'color': ROOT.kViolet-3, 'style': 20, 'size': 1.4},
#{'run': 300817, 'split': 1, 'color': ROOT.kTeal-3, 'style': 29, 'size': 1.4},
#{'run': 301298, 'split': 1, 'color': ROOT.kOrange-3, 'style': 21, 'size': 1.4},
#{'run': 301417, 'split': 1, 'color': ROOT.kViolet+3, 'style': 22, 'size': 1.4},
#{'run': 301627, 'split': 1, 'color': ROOT.kBlue+3, 'style': 24, 'size': 1.4},
#{'run': 301664, 'split': 1, 'color': ROOT.kGreen-3, 'style': 20, 'size': 1.4},
#{'run': 301694, 'split': 1, 'color': ROOT.kRed+3, 'style': 23, 'size': 1.4},
#{'run': 302019, 'split': 1, 'color': ROOT.kGreen+2, 'style': 29, 'size': 1.4},
#{'run': 302026, 'split': 1, 'color': ROOT.kMagenta+2, 'style': 33, 'size': 1.4},
#{'run': 302228, 'split': 1, 'color': ROOT.kRed+2, 'style': 21, 'size': 1.4},
#{'run': 302262, 'split': 1, 'color': ROOT.kTeal+2, 'style': 20, 'size': 1.4},
#{'run': 302277, 'split': 1, 'color': ROOT.kCyan+2, 'style': 22, 'size': 1.4},
#]
#runlist = [302328]
#make_tgraph(graphOptions, "p76_all_runs.png", xrange=[0.5,1.7])

for run in runlist: 
    graphOptions = [
        {'run': run, 'split': 1, 'color': ROOT.kBlue+2, 'size': 1.6, 'style': 33, 'yrange': [0, 10], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]', 'limitpoints': 100},
    ]
    make_tgraph(graphOptions, "p76_run_%d.png"%run, xrange=[0.0,1.8])

