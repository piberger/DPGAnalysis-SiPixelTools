import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("MyRawToDigi")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# to use no All 
# 2015
#process.GlobalTag.globaltag = 'GR_P_V56' # works for 2469763
#process.GlobalTag.globaltag = 'GR_P_V56' # for 247607
# 2017
process.GlobalTag.globaltag = '92X_dataRun2_Express_v2' # 
# 2016
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v3' # for 266277
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v8' # for 272
# AUTO conditions 
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_design', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2017', '')

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
# accept if 'path_1' succeeds
process.hltfilter = hlt.hltHighLevel.clone(
# Min-Bias	
#    HLTPaths = ['HLT_Physics*'],
#    HLTPaths = ['HLT_Random*'],
    HLTPaths = ['HLT_ZeroBias*'],
#    HLTPaths = ['HLT_L1SingleMuOpen_v*'],
#    HLTPaths = ['HLT_PAZeroBias*'],
#    HLTPaths = ['HLT_PARandom*'],
#    HLTPaths = ['HLT_PAMinBias*'],
# Commissioning:
#    HLTPaths = ['HLT_L1Tech5_BPTX_PlusOnly_v*'],
#    HLTPaths = ['HLT_L1Tech6_BPTX_MinusOnly_v*'],
#    HLTPaths = ['HLT_L1Tech7_NoBPTX_v*'],
#
#    HLTPaths = ['p*'],
#    HLTPaths = ['path_?'],
    andOr = True,  # False = and, True=or
    throw = False
    )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

print "ARGS:",sys.argv
restrictToLumisection = -1
if len(sys.argv) > 3:
    restrictToLumisection = int(sys.argv[3])
    print "-> restrict to lumisection ", restrictToLumisection

process.source = cms.Source("PoolSource",
 # fileNames =  cms.untracked.vstring('file:rawdata.root')
 fileNames =  cms.untracked.vstring(

sys.argv[2],
#"/store/express/Run2017C/ExpressPhysics/FEVT/Express-v1/000/299/480/00000/0047AB9D-6E6D-E711-9A32-02163E0119EE.root",
#"/store/express/Run2017C/ExpressPhysics/FEVT/Express-v1/000/299/479/00000/000AE42E-696D-E711-AC50-02163E0142C8.root",
#"/store/express/Run2017C/ExpressPhysics/FEVT/Express-v1/000/299/478/00000/0455436B-656D-E711-BE85-02163E011F68.root",
#"/store/express/Run2017C/ExpressPhysics/FEVT/Express-v1/000/299/481/00000/000B7C22-9B6D-E711-9223-02163E019DB2.root",
#"/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/219/00000/004B8939-EF55-E711-AF94-02163E01461F.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v2/000/296/643/00000/08070285-5E4F-E711-BB56-02163E01476C.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v2/000/296/664/00000/00F38C8D-A54F-E711-A312-02163E01472F.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v3/000/297/016/00000/0255F076-B552-E711-9E49-02163E011FA9.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v3/000/297/015/00000/0C482900-BA52-E711-A271-02163E014666.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v3/000/297/011/00000/0C3E8DF9-B652-E711-9480-02163E011EB4.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v3/000/297/012/00000/02E7996A-B852-E711-B96F-02163E0143F0.root",
#"/store/express/Run2017A/ExpressPhysics/FEVT/Express-v3/000/297/003/00000/029D5DDE-A252-E711-9166-02163E0146E9.root",
# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/295/439/00000/426B4782-FC43-E711-9328-02163E011A76.root",
# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/295/439/00000/F853E939-FF43-E711-A59D-02163E013479.root",

# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/295/381/00000/00034176-5043-E711-AA6F-02163E019C98.root",

# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/295/318/00000/06085100-B142-E711-A1D0-02163E01A1FA.root",

# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/295/209/00000/005E61C9-D341-E711-BEAE-02163E019C9F.root",

# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/295/128/00000/2E861F6A-1F41-E711-BEEC-02163E01A32B.root",

# "/store/express/Run2017A/ExpressPhysics/FEVT/Express-v1/000/294/929/00000/000E8FD5-1D40-E711-80C9-02163E01467E.root",


# 273725
#"/store/express/Run2016B/ExpressPhysics/FEVT/Express-v2/000/273/725/00000/00030CAB-6C1E-E611-90F9-02163E0137A8.root",
# 273730
#"/store/express/Run2016B/ExpressPhysics/FEVT/Express-v2/000/273/730/00000/00821BA0-381F-E611-BBC5-02163E0142C0.root",


 ),
#    "file:/afs/cern.ch/work/d/dkotlins/public/MC/mu/pt100_71_pre7/raw/raw2.root"
  inputCommands=cms.untracked.vstring(
	  'keep *',
	  'drop CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters__RECO'
  )

)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('273725:83-273725:9999')

# Cabling
#  include "CalibTracker/Configuration/data/Tracker_FakeConditions.cff"
#process.load("CalibTracker.Configuration.SiPixel_FakeConditions_cff")
#process.load("CalibTracker.Configuration.SiPixelCabling.SiPixelCabling_SQLite_cff")
#process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_21X_GLOBALTAG"
#process.GlobalTag.globaltag = "CRAFT_V3P::All"
#process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')
#process.siPixelCabling.connect = 'sqlite_file:cabling.db'
#process.siPixelCabling.toGet = cms.VPSet(cms.PSet(
#    record = cms.string('SiPixelFedCablingMapRcd'),
#    tag = cms.string('SiPixelFedCablingMap_v14')
#))


process.load("EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi")
# for simultaions 
#process.siPixelDigis.InputLabel = 'siPixelRawData'
# for data
#process.siPixelDigis.InputLabel = 'source'
process.siPixelDigis.InputLabel = 'rawDataCollector'
process.siPixelDigis.IncludeErrors = True
process.siPixelDigis.Timing = False 
#process.siPixelDigis.UsePilotBlade = True 
process.siPixelDigis.UsePhase1 = True 

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siPixelDigis'),
    destinations = cms.untracked.vstring('r2d'),
#    r2d = cms.untracked.PSet( threshold = cms.untracked.string('DEBUG'))
    r2d = cms.untracked.PSet( threshold = cms.untracked.string('INFO'))
#    r2d = cms.untracked.PSet( threshold = cms.untracked.string('WARNING'))
)

try:
    runNumber = sys.argv[2].split('/')[-4] + sys.argv[2].split('/')[-3]
except:
    runNumber = None

outputFilename = 'digis_'

if runNumber:
    outputFilename = outputFilename + runNumber + '_'
if (restrictToLumisection > 0):
    outputFilename = outputFilename + '%d_'%restrictToLumisection

outputFilename = outputFilename + sys.argv[2].split('/')[-1].strip()

if outputFilename.endswith('.root'):
    outputFilename = outputFilename[0:-5]
outputFullName = 'file:/afs/cern.ch/work/p/piberger/' + outputFilename + '.root'
print "OUTPUT:", outputFullName

process.out = cms.OutputModule("PoolOutputModule",
    fileName =  cms.untracked.string(outputFullName),
    outputCommands = cms.untracked.vstring("drop *","keep *_siPixelDigis_*_*")
)


process.a = cms.EDAnalyzer("PixDigisTest",
    Verbosity = cms.untracked.bool(False),
    phase1 = cms.untracked.bool(True),
    includezeroadc = cms.untracked.bool(False),
    lumisection = cms.untracked.int32(restrictToLumisection),
# sim in V7
#    src = cms.InputTag("mix"),
# old default
    src = cms.InputTag("siPixelDigis"),
    Select1 = cms.untracked.int32(0),  # select the cut type, o no cut
    Select2 = cms.untracked.int32(0),  # select the cut value   
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFullName)
)

#process.p = cms.Path(process.siPixelDigis)
#process.p = cms.Path(process.siPixelDigis*process.a)
process.p = cms.Path(process.hltfilter*process.siPixelDigis*process.a)
# disable data write to disk 
#process.ep = cms.EndPath(process.out)
