from FWCore.PythonUtilities.LumiList   import LumiList
import sys


files = sys.argv[1:]
# put this here after parsing the arguments since ROOT likes to
# grab command line arguments even when it shouldn't.
from DataFormats.FWLite import Lumis, Handle

for f in files:
    #ef = f.replace("/eos/cms","root://eoscms.cern.ch")
    #ef = f.replace("/eos/cms","root://cms-xrd-global.cern.ch")
    ef = f.replace("/eos/cms","root://xrootd-cms.infn.it")
    lumis = Lumis (ef)
    for lum in lumis:
        print f,lum.aux().run(),lum.aux().id().luminosityBlock()
