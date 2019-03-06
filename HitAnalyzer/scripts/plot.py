import ROOT
import glob
import array
import math

# all: number of straight long clustersalong rows  of size min to max, complete + broken1 + broken2
# broken1: number of such broken clusters with 1 pixel missing
# ratio1: broken1/all
p76data = {'all': {}, 'broken1': {}, 'broken2': {}, 'ratio1': {}, 'ratio2': {}}

def set_tdr_style():
        tdrStyle = ROOT.TStyle("tdrStyle", "Style for P-TDR")

        tdrStyle.SetCanvasBorderMode(0)
        tdrStyle.SetCanvasColor(ROOT.kWhite)
        tdrStyle.SetCanvasDefH(800)  # Height of canvas
        tdrStyle.SetCanvasDefW(800)  # Width of canvas
        tdrStyle.SetCanvasDefX(0)  # POsition on screen
        tdrStyle.SetCanvasDefY(0)

        tdrStyle.SetPadBorderMode(0)
        tdrStyle.SetPadColor(ROOT.kWhite)
        tdrStyle.SetPadGridX(True)
        tdrStyle.SetPadGridY(True)
        tdrStyle.SetGridColor(0)
        tdrStyle.SetGridStyle(3)
        tdrStyle.SetGridWidth(1)

        tdrStyle.SetFrameBorderMode(0)
        tdrStyle.SetFrameBorderSize(1)
        tdrStyle.SetFrameFillColor(0)
        tdrStyle.SetFrameFillStyle(0)
        tdrStyle.SetFrameLineColor(1)
        tdrStyle.SetFrameLineStyle(1)
        tdrStyle.SetFrameLineWidth(1)

        # tdrStyle.SetHistFillColor(1)
        # tdrStyle.SetHistFillStyle(0)
        tdrStyle.SetHistLineColor(1)
        tdrStyle.SetHistLineStyle(0)
        tdrStyle.SetHistLineWidth(1)
        # tdrStyle.SetLegoInnerR(Float_t rad = 0.5)
        # tdrStyle.SetNumberContours(Int_t number = 20)

        tdrStyle.SetEndErrorSize(2)
        # tdrStyle.SetErrorMarker(20)
        # tdrStyle.SetErrorX(0.)

        tdrStyle.SetMarkerStyle(20)

        # For the fit/function:
        tdrStyle.SetOptFit(1)
        tdrStyle.SetFitFormat("5.4g")
        tdrStyle.SetFuncColor(2)
        tdrStyle.SetFuncStyle(1)
        tdrStyle.SetFuncWidth(1)

        # For the date:
        tdrStyle.SetOptDate(0)
        # tdrStyle.SetDateX(Float_t x = 0.01)
        # tdrStyle.SetDateY(Float_t y = 0.01)

        # For the statistics box:
        tdrStyle.SetOptFile(0)
        tdrStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat("mr")
        tdrStyle.SetStatColor(ROOT.kWhite)
        tdrStyle.SetStatFont(42)
        tdrStyle.SetStatFontSize(0.025)
        tdrStyle.SetStatTextColor(1)
        tdrStyle.SetStatFormat("6.4g")
        tdrStyle.SetStatBorderSize(1)
        tdrStyle.SetStatH(0.1)
        tdrStyle.SetStatW(0.15)
        # tdrStyle.SetStatStyle(Style_t style = 1001)
        # tdrStyle.SetStatX(Float_t x = 0)
        # tdrStyle.SetStatY(Float_t y = 0)

        # Margins:
        tdrStyle.SetPadTopMargin(0.05)
        tdrStyle.SetPadBottomMargin(0.13)
        tdrStyle.SetPadLeftMargin(0.16)
        tdrStyle.SetPadRightMargin(0.06)

        # For the Global title:

        tdrStyle.SetOptTitle(0)
        tdrStyle.SetTitleFont(42)
        tdrStyle.SetTitleColor(1)
        tdrStyle.SetTitleTextColor(1)
        tdrStyle.SetTitleFillColor(10)
        tdrStyle.SetTitleFontSize(0.05)
        # tdrStyle.SetTitleH(0) # Set the height of the title box
        # tdrStyle.SetTitleW(0) # Set the width of the title box
        # tdrStyle.SetTitleX(0) # Set the position of the title box
        # tdrStyle.SetTitleY(0.985) # Set the position of the title box
        # tdrStyle.SetTitleStyle(Style_t style = 1001)
        # tdrStyle.SetTitleBorderSize(2)

        # For the axis titles:

        tdrStyle.SetTitleColor(1, "XYZ")
        tdrStyle.SetTitleFont(42, "XYZ")
        tdrStyle.SetTitleSize(0.06, "XYZ")
        # tdrStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
        # tdrStyle.SetTitleYSize(Float_t size = 0.02)
        tdrStyle.SetTitleXOffset(0.9)
        tdrStyle.SetTitleYOffset(1.85)
        # tdrStyle.SetTitleOffset(1.1, "Y") # Another way to set the Offset

        # For the axis labels:

        tdrStyle.SetLabelColor(1, "XYZ")
        tdrStyle.SetLabelFont(42, "XYZ")
        tdrStyle.SetLabelOffset(0.007, "XYZ")
        tdrStyle.SetLabelSize(0.05, "XYZ")

        # For the axis:

        tdrStyle.SetAxisColor(1, "XYZ")
        tdrStyle.SetStripDecimals(ROOT.kTRUE)
        tdrStyle.SetTickLength(0.03, "XYZ")
        tdrStyle.SetNdivisions(510, "XYZ")
        tdrStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
        tdrStyle.SetPadTickY(1)

        # Change for log plots:
        tdrStyle.SetOptLogx(0)
        tdrStyle.SetOptLogy(0)
        tdrStyle.SetOptLogz(0)

        # Postscript options:
        tdrStyle.SetPaperSize(20., 20.)
        # tdrStyle.SetLineScalePS(Float_t scale = 3)
        # tdrStyle.SetLineStyleString(Int_t i, const char* text)
        # tdrStyle.SetHeaderPS(const char* header)
        # tdrStyle.SetTitlePS(const char* pstitle)

        # tdrStyle.SetBarOffset(Float_t baroff = 0.5)
        # tdrStyle.SetBarWidth(Float_t barwidth = 0.5)
        # tdrStyle.SetPaintTextFormat(const char* format = "g")
        # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
        # tdrStyle.SetTimeOffset(Double_t toffset)
        # tdrStyle.SetHistMinimumZero(ROOT.kTrue)

        tdrStyle.SetHatchesLineWidth(5)
        tdrStyle.SetHatchesSpacing(0.05)

        tdrStyle.cd()

def get_luminosity(runNumber, lumisectionNumber):
    lumiCSVfileName = 'lumi%d.csv'%runNumber
    instaLumi = -1
    with open(lumiCSVfileName, 'r') as lumiCSVfile:
        for line in lumiCSVfile:
            values = [x.strip() for x in line.split(',')]
            try:
                if int(values[0]) == lumisectionNumber:
                    instaLumi = float(values[3]) * 1.0e-4
                    break
            except:
                pass
    return instaLumi


def get_lumisection_data(fileName, options):
    #print "FILE:", fileName
    path = '/'.join(fileName.split('/')[0:-1]).strip()
    runNumber = int(fileName.split('/')[-1].split('_')[1])
    lumisectionNumber = int(fileName.split('/')[-1].split('_')[2])

    # get instantaneous lumi 
    if 'lsplot' in options and options['lsplot']:
        instaLumi = lumisectionNumber
    else:
        instaLumi = get_luminosity(runNumber, lumisectionNumber)
 
    # get cluster splitting probability
    rootFile = ROOT.TFile.Open(fileName)
    if (not rootFile or  rootFile.IsZombie()):
        print "FAILED open file."
        return None
    tree = rootFile.Get('a/tree')
    if not tree:
        print "FAILED getting tree."
        return None

    nClustersComplete = 0
    nClustersBroken = 0

    clusterSplitType = options['split'] if 'split' in options else 2
    clusterSelectionLambda = options['selection'] if 'selection' in options else lambda x: x.size == 7 and '_LYR1_' in x.module 

    global p76data
    
    numEntries = tree.GetEntries()
    for i in range(numEntries):
        tree.GetEntry(i)
        #if tree.size == 7 and '_LYR1_' in tree.module and '_LDR6' in tree.module and ('BpI' in tree.module or 'BmI' in tree.module):
        if clusterSelectionLambda(tree):
            nClustersComplete += tree.complete
            if clusterSplitType == 2:
                nClustersBroken += tree.broken2
            elif clusterSplitType == 1:
                nClustersBroken += tree.broken1
            else:
                print "bad value for split type:", clusterSplitType
                exit()
        #p76 test
        moduleName = tree.module.split('\0')[0]
        if (tree.size == 7):
            if (moduleName not in p76data['all']):
                p76data['all'][moduleName] = []
                p76data['broken1'][moduleName] = []
                p76data['broken2'][moduleName] = []
                p76data['ratio1'][moduleName] = []
                p76data['ratio2'][moduleName] = []
            p76data['ratio1'][moduleName].append(tree.broken1*1.0/(tree.complete+tree.broken1+tree.broken2))
            p76data['ratio2'][moduleName].append(tree.broken2*1.0/(tree.complete+tree.broken1+tree.broken2))
            p76data['broken1'][moduleName].append(tree.broken1)
            p76data['broken2'][moduleName].append(tree.broken2)
            p76data['all'][moduleName].append(tree.broken1+tree.complete+tree.broken2)
    #print "TEST p76:"
    #for module,p76dataPoint in p76data.items():
    #    print module,' ',sum(p76dataPoint)/(1.0*len(p76dataPoint)) if len(p76dataPoint) > 0 else 0
    #raw_input()
    
    nClustersTotal = nClustersComplete + nClustersBroken
    fraction = 1.0 * nClustersBroken / nClustersTotal if nClustersTotal > 0 else -1 
    error = math.sqrt(fraction*(1.0-fraction)/(1.0*nClustersTotal)) if nClustersTotal > 0 else 0
    if 'LSoffset' in options:
        # plot vs LS number (with offset, set offset to 0 for plot vs. LS)
        result = [lumisectionNumber + options['LSoffset'],  fraction*100.0, 0, error*100.0]
    else:
        # plot vs. instantaneous lumi
        result = [instaLumi, fraction*100.0, 0, error*100.0]
    return result
   
def get_bx_data(fileNames, options):
    nBx = 3564 
    bxClustersComplete = [0]*nBx
    bxClustersSplit = [0]*nBx

    for fileName in fileNames:
        fileLs = int(fileName.split('/')[-1].split('_')[2])
        if 'ls' in options and fileLs < options['ls'][0] or fileLs > options['ls'][1]:
            print "SKIP:", fileName
            continue
        rootFile = ROOT.TFile.Open(fileName)
        if (not rootFile or  rootFile.IsZombie()):
            print "FAILED open file.", fileName
            continue
        tree = rootFile.Get('a/bxtree')
        if not tree:
            print "FAILED getting tree.", fileName
            continue

        clusterSplitType = options['split'] if 'split' in options else 2
        clusterSelectionLambda = options['selection'] if 'selection' in options else lambda x: x.size == 7

        numEntries = tree.GetEntries()
        for i in range(numEntries):
            tree.GetEntry(i)
            bx = tree.bx
            if clusterSelectionLambda(tree):
                bxClustersComplete[bx] += tree.complete
                if clusterSplitType == 2:
                    bxClustersSplit[bx] += tree.broken2
                elif clusterSplitType == 1:
                    bxClustersSplit[bx] += tree.broken1
                else:
                    print "bad value for split type:", clusterSplitType
                    exit()
    
    bxClustersRatio = [0]*nBx
    bxClustersRatioError = [0]*nBx
    for bx in range(nBx):
        fraction = 1.0 * bxClustersSplit[bx] / (bxClustersSplit[bx] + bxClustersComplete[bx]) if (bxClustersSplit[bx] + bxClustersComplete[bx]) > 0 else -1 
        error = math.sqrt(fraction*(1.0-fraction)/(1.0*(bxClustersSplit[bx] + bxClustersComplete[bx]))) if (bxClustersSplit[bx] + bxClustersComplete[bx]) > 0 else 0
        bxClustersRatio[bx] = fraction*100.0
        bxClustersRatioError[bx] = error*100.0
    return [bxClustersRatio, bxClustersRatioError]

def get_all_points(options):
    runNumber = options['run']
    result = []
    globMask = ((options['subdir']+'/') if 'subdir' in options else '') + "digis_%d_*_*.root"%runNumber
    clusterTreeFileNames = glob.glob(globMask)
    for clusterTreeFileName in clusterTreeFileNames:
        if 'ls' in options:
            lsRange = options['ls']
            ls = int(clusterTreeFileName.split('/')[-1].split('_')[2])
            if ls < lsRange[0] or ls > lsRange[1]:
                continue
            else:
                print "USE LS:", ls
        result.append(get_lumisection_data(clusterTreeFileName, options))

    '''
    print "TEST p76:"
    
    with open('p76_permodule_%d_all.txt'%runNumber, 'w') as outputFile:
        for module,p76dataPoint in p76data['all'].items():
            outputFile.write('%s %d\n'%(module,sum(p76dataPoint)))
    with open('p76_permodule_%d_broken1.txt'%runNumber, 'w') as outputFile:
        for module,p76dataPoint in p76data['broken1'].items():
            outputFile.write('%s %d\n'%(module,sum(p76dataPoint)))
    with open('p76_permodule_%d_broken2.txt'%runNumber, 'w') as outputFile:
        for module,p76dataPoint in p76data['broken2'].items():
            outputFile.write('%s %d\n'%(module,sum(p76dataPoint)))
    with open('p76_permodule_%d_ratio1.txt'%runNumber, 'w') as outputFile:
        for module,p76dataPoint in p76data['ratio1'].items():
            outputFile.write('%s %f\n'%(module,sum(p76dataPoint)/(1.0*len(p76dataPoint)) if len(p76dataPoint) > 0 else 0))
    with open('p75_permodule_%d_ratio2.txt'%runNumber, 'w') as outputFile:
        for module,p76dataPoint in p76data['ratio2'].items():
            outputFile.write('%s %f\n'%(module,sum(p76dataPoint)/(1.0*len(p76dataPoint)) if len(p76dataPoint) > 0 else 0))
    '''

    return result

def make_graph(options):
    graphPointsListAll = get_all_points(options)
    graphPointsList = [x for x in graphPointsListAll if x and x[0] > 0.001 and x[1] >= 0.1]

    if 'limitpoints' in options and len(graphPointsList) > options['limitpoints']:
        graphPointsList = graphPointsList[0:options['limitpoints']]
        print 'number of points limited to ', options['limitpoints']
    print graphPointsList 
    try:
        selectedPoints = [x for x in graphPointsList if x[0] > 0.001 and x[3]<0.5]
        mean = sum([x[1] for x in selectedPoints]) / len(selectedPoints) if len(selectedPoints) > 0 else 0
    except:
        mean = -1 
    print 'selected:', selectedPoints
    print '-'*80,'\n','RUN:',options['run'],':',mean,'\n','-'*80
    #for p in graphPointsList:
    #    print ('%1.3f'%p[0]).ljust(7),('%1.3f'%p[1]).ljust(7),('%1.3f'%p[3]).ljust(10) 
    #print graphPointsList

    rescaleX = options['scalex'] if 'scalex' in options else 1.0

    x_points = array.array('d', [x[0]*rescaleX for x in graphPointsList])
    y_points = array.array('d', [x[1] for x in graphPointsList])
    x_errors = array.array('d', [x[2]*rescaleX for x in graphPointsList])
    y_errors = array.array('d', [x[3] for x in graphPointsList])
    
    try: 
        tg = ROOT.TGraphErrors(len(graphPointsList), x_points, y_points, x_errors, y_errors)
    except:
        tg = None

    if tg:
        tg.SetMarkerStyle(options['style'] if 'style' in options else 2)
        tg.SetMarkerColor(options['color'] if 'color' in options else ROOT.kBlue+2)
        tg.SetMarkerSize(options['size'] if 'size' in options else 1.5)
        
        if 'xrange' in options:
            tg.GetXaxis().SetRangeUser(options['xrange'][0], options['xrange'][1])
            tg.GetXaxis().SetLimits(options['xrange'][0], options['xrange'][1])
        if 'yrange' in options:
            tg.GetYaxis().SetRangeUser(options['yrange'][0], options['yrange'][1])
        if 'xtitle' in options:
            tg.GetXaxis().SetTitle(options['xtitle']) 
        if 'ytitle' in options:
            tg.GetYaxis().SetTitle(options['ytitle']) 
            tg.GetYaxis().SetTitleOffset(1.05) 
        tg.SetFillStyle(0) 
    #tg.GetXaxis().SetLabelSize(0.03)
    return tg 

def make_tgraph(graphOptions, fileName, legendPos=None, xrange=None):
    first = True
    if legendPos:
        leg = ROOT.TLegend(*legendPos)
    else:
        leg = ROOT.TLegend(0.19,0.65,0.36,0.85)
    for graphOption in graphOptions:
        try:
            graphOption['tgraph'] = make_graph(graphOption)
            if xrange:
                graphOption['tgraph'].GetXaxis().SetLimits(xrange[0],xrange[1])
                graphOption['tgraph'].GetXaxis().SetRangeUser(xrange[0],xrange[1])
            if first:
                graphOption['tgraph'].Draw('AP')
                first = False
            else:
                graphOption['tgraph'].Draw('P;same')
            leg.AddEntry(graphOption['tgraph'], graphOption['name'] if 'name' in graphOption else '%d'%graphOption['run'])
        except Exception as e:
           print 'EXCEPTION:', e
    leg.SetFillColor(ROOT.kWhite)
    leg.Draw()
    c1.SaveAs(fileName)
    fileNameRoot ='.'.join(fileName.split('.')[:-1]) + '.root'
    c1.SaveAs(fileNameRoot)


def make_tgraph_bx(options, legendPos=None, xrange=None):
    runNumber = options['run']
    fileNames = glob.glob("digis_%d_*_*.root"%runNumber)
    bxData = get_bx_data(fileNames, options)
    
    x_points = array.array('d', [i for i, x in enumerate(bxData[0]) if x > 0])
    y_points = array.array('d', [x for i, x in enumerate(bxData[0]) if x > 0])
    x_errors = array.array('d', [0 for i, x in enumerate(bxData[0]) if x > 0])
    y_errors = array.array('d', [bxData[1][i]  for i, x in enumerate(bxData[0]) if x > 0])
    
    tg = ROOT.TGraphErrors(len(x_points), x_points, y_points, x_errors, y_errors)
    tg.SetMarkerStyle(options['style'] if 'style' in options else 2)
    tg.SetMarkerColor(options['color'] if 'color' in options else ROOT.kBlue+2)
    tg.SetMarkerSize(options['size'] if 'size' in options else 1.5)
    
    if 'xrange' in options:
        tg.GetXaxis().SetRangeUser(options['xrange'][0], options['xrange'][1])
        tg.GetXaxis().SetLimits(options['xrange'][0], options['xrange'][1])
    if 'yrange' in options:
        tg.GetYaxis().SetRangeUser(options['yrange'][0], options['yrange'][1])
    if 'xtitle' in options:
        tg.GetXaxis().SetTitle(options['xtitle']) 
    if 'ytitle' in options:
        tg.GetYaxis().SetTitle(options['ytitle']) 
        tg.GetYaxis().SetTitleOffset(1.05) 
    tg.SetFillStyle(0) 
    tg.Draw('AP')	
    tg.GetXaxis().SetTitle('bx')
    tg.GetYaxis().SetTitle('p(7->5) [%]')
    tg.GetYaxis().SetTitleOffset(1.1)
    tg.SetTitle('Run %d'%runNumber)
    fileName = 'cluster_splitting_vs_bx_%d.png'%runNumber
    c1.SaveAs(fileName)
    fileNameRoot ='.'.join(fileName.split('.')[:-1]) + '.root'
    c1.SaveAs(fileNameRoot)

def get_bigcluster_data(fileNames, options):
    bData = []
    for fileName in fileNames:
        tFile = ROOT.TFile.Open(fileName, 'read')
        if tFile and not tFile.IsZombie():
            histo = tFile.Get('a/hNdbCellsLay1')
            print file
            if histo:
                binContents = {x: histo.GetBinContent(histo.GetXaxis().FindBin(x)) for x in [36,37,38,39,40]}
                print binContents
                fileParts = fileName.split('_')
                countBig = sum([v for k,v in binContents.iteritems()])
                fraction = 1.0 * countBig / histo.GetEntries()

                fraction_error = math.sqrt(fraction * (1.0-fraction)/(1.0*histo.GetEntries())) if histo.GetEntries() > 0 else 0
                print "RUN",fileParts[1]," ",fileParts[2]," ",countBig," ",histo.GetEntries()," ",fraction
                bData.append([
                                int(fileParts[2]) if 'LS' in options else get_luminosity(int(fileParts[1]), int(fileParts[2])),
                                fraction,
                                0,
                                fraction_error
                            ])        
            else:
                print "not found!"
    return bData

def get_tgraph_bigclusters(options):
    runNumber = options['run']
    fileNames = glob.glob("digis_%d_*_*.root"%runNumber)
    bData = get_bigcluster_data(fileNames, options)

    bDataClean = [x for x in bData if x[0] > 0]
   
    x_points = array.array('d', [x[0] for x in bDataClean])
    y_points = array.array('d', [x[1] for x in bDataClean])
    x_errors = array.array('d', [x[2] for x in bDataClean]) 
    y_errors = array.array('d', [x[3] for x in bDataClean])
    
    tg = ROOT.TGraphErrors(len(x_points), x_points, y_points, x_errors, y_errors)
    tg.SetMarkerStyle(options['style'] if 'style' in options else 2)
    tg.SetMarkerColor(options['color'] if 'color' in options else ROOT.kBlue+2)
    tg.SetMarkerSize(options['size'] if 'size' in options else 1.5)
    
    if 'xrange' in options:
        tg.GetXaxis().SetRangeUser(options['xrange'][0], options['xrange'][1])
        tg.GetXaxis().SetLimits(options['xrange'][0], options['xrange'][1])
    if 'yrange' in options:
        tg.GetYaxis().SetRangeUser(options['yrange'][0], options['yrange'][1])
    if 'xtitle' in options:
        tg.GetXaxis().SetTitle(options['xtitle']) 
    if 'ytitle' in options:
        tg.GetYaxis().SetTitle(options['ytitle']) 
        tg.GetYaxis().SetTitleOffset(1.05) 
    tg.SetFillStyle(0) 
    tg.Draw('AP')
    tg.GetXaxis().SetTitle('LS' if 'LS' in options else 'luminosity')
    tg.GetYaxis().SetTitle('fraction')
    tg.GetYaxis().SetTitleOffset(1.1)
    tg.SetTitle('Run %d'%runNumber)
    return tg

def make_tgraph_bigclusters(optionsList, fileName, legendPos=None):
    graphs = [ [options, get_tgraph_bigclusters(options)] for options in optionsList]
    first = True
    if legendPos:
        leg = ROOT.TLegend(*legendPos)
    else:
        leg = ROOT.TLegend(0.19,0.65,0.36,0.85)
    for options, graph in graphs:
         if first:
            graph.Draw('AP')
            first = False
         else:
            graph.Draw('P,SAME')
         leg.AddEntry(graph, options['name'] if 'name' in options else '%d'%options['run'])       
    leg.Draw()
    c1.SaveAs(fileName)
    fileNameRoot ='.'.join(fileName.split('.')[:-1]) + '.root'
    c1.SaveAs(fileNameRoot)


set_tdr_style()
c1 = ROOT.TCanvas("c1","c1",800,800)
set_tdr_style()

if __name__ == "__main__":
    graphOptions = [
    {'run': 299593, 'xrange': [0.9,1.8], 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300389, 'color': ROOT.kGreen+2, 'style': 22, 'size': 1.8},
    {'run': 300400, 'color': ROOT.kGreen+4, 'style': 23, 'size': 1.8},
    {'run': 300742, 'color': ROOT.kMagenta-3, 'style': 47, 'size': 1.8},
    {'run': 300806, 'color': ROOT.kRed+2, 'style': 34},
    {'run': 300811, 'color': ROOT.kPink+4, 'style': 29, 'size': 2.0},
    {'run': 300812, 'color': ROOT.kOrange+2, 'style': 33},
    ]

    #make_tgraph(graphOptions, "output.png")

    graphOptions = [
    {'run': 299593, 'color': ROOT.kMagenta+2, 'style': 21, 'xrange': [0.8,1.52], 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 299478, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.8},
    {'run': 299479, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.8},
    {'run': 299380, 'color': ROOT.kOrange+3, 'style': 23, 'size': 1.8},
    {'run': 300811, 'color': ROOT.kGreen+2, 'style': 29, 'size': 2.0},
    ]
    #make_tgraph(graphOptions, "new_timing_compare_resets_p75.png")

    graphOptions = [
    {'run': 299593, 'color': ROOT.kMagenta+2, 'style': 21, 'xrange': [0.8,1.52], 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 299478, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.8},
    {'run': 299479, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.8},
    {'run': 299380, 'color': ROOT.kOrange+3, 'style': 23, 'size': 1.8},
    {'run': 300811, 'color': ROOT.kGreen+2, 'style': 29, 'size': 2.0},
    {'run': 300811, 'name': '300811 (incl. ADC 0)', 'subdir': 'test/adc0included/',  'color': ROOT.kRed, 'style': 29, 'size': 2.0},
    ]
    #make_tgraph(graphOptions, "new_timing_compare_resets_p75_with_adc0included.png")

    graphOptions = [
    {'run': 299593, 'split': 1, 'color': ROOT.kMagenta+2, 'style': 21, 'xrange': [0.8,1.52], 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]'},
    {'run': 299478, 'split': 1, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.8},
    {'run': 299479, 'split': 1, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.8},
    {'run': 299380, 'split': 1, 'color': ROOT.kOrange+3, 'style': 23, 'size': 1.8},
    {'run': 300811, 'split': 1, 'color': ROOT.kGreen+2, 'style': 29, 'size': 2.0},
    {'run': 300811, 'split': 1, 'name': '300811 (incl. ADC 0)', 'subdir': 'test/adc0included/',  'color': ROOT.kRed, 'style': 29, 'size': 2.0},
    ]
    #make_tgraph(graphOptions, "new_timing_compare_resets_p76_with_adc0included.png")


    graphOptions = [
    {'run': 299593, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'color': ROOT.kMagenta+2, 'style': 21, 'xrange': [0.8,1.52], 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'LYR2 cluster splitting p(7->5) [%]'},
    {'run': 299478, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.8},
    {'run': 299479, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.8},
    {'run': 299380, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'color': ROOT.kOrange+3, 'style': 23, 'size': 1.8},
    {'run': 300811, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'color': ROOT.kGreen+2, 'style': 29, 'size': 2.0},
    ]

    #make_tgraph(graphOptions, "new_timing_compare_resets_p75_layer2.png")


    graphOptions = [
    {'run': 299593, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'split': 1, 'color': ROOT.kMagenta+2, 'style': 21, 'xrange': [0.8,1.52], 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'LYR2 cluster splitting p(7->6) [%]'},
    {'run': 299478, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'split': 1, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.8},
    {'run': 299479, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'split': 1, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.8},
    {'run': 299380, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'split': 1, 'color': ROOT.kOrange+3, 'style': 23, 'size': 1.8},
    {'run': 300811, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'split': 1, 'color': ROOT.kGreen+2, 'style': 29, 'size': 2.0},
    ]

    #make_tgraph(graphOptions, "new_timing_compare_resets_p76_layer2.png")

    graphOptions = [
    {'run': 300806, 'color': ROOT.kMagenta+2, 'style': 21, 'xrange': [0.0,1.8], 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300811, 'color': ROOT.kOrange+1, 'style': 22, 'size': 1.8},
    {'run': 300812, 'color': ROOT.kBlue+1, 'style': 34, 'size': 1.8},
    {'run': 301046, 'color': ROOT.kOrange+3, 'style': 23, 'size': 1.8},
    {'run': 301141, 'color': ROOT.kGreen+2, 'style': 29, 'size': 2.0},
    {'run': 301283, 'color': ROOT.kRed, 'style': 3, 'size': 1.3},
    {'run': 301330, 'color': ROOT.kBlack, 'style': 33, 'size': 1.8},
    ]

    #make_tgraph(graphOptions, "cluster_splitting_new_runs.png", legendPos=[0.34,0.65,0.51,0.85])


    graphOptions = [
    {'run': 301281, 'color': ROOT.kBlack, 'style': 21, 'xrange': [0.6,1.0], 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 301283, 'color': ROOT.kGreen+2, 'style': 34, 'size': 1.8},
    {'run': 301179, 'color': ROOT.kBlue+2, 'style': 22, 'size': 1.8},
    {'run': 301180, 'color': ROOT.kOrange+2, 'style': 23, 'size': 1.8},
    ]

    #make_tgraph(graphOptions, "reset_50_70_hz.png", xrange=[0.6,0.9])

    graphOptions = [
    {'run': 301417, 'name': '301417 (1932b)', 'color': ROOT.kBlack, 'style': 21, 'limitpoints': 50, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 299593, 'name': '299593 (2556b)', 'color': ROOT.kGreen+2, 'style': 22, 'size': 1.8, 'limitpoints': 50}, 
    {'run': 301391, 'name': '301391 (1740b)', 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 50} ,
    {'run': 300817, 'name': '300817 (2556b)', 'color': ROOT.kOrange+2, 'style': 34, 'size': 1.8, 'limitpoints': 50}, 
    {'run': 301283, 'name': '301283 (1356b)', 'color': ROOT.kMagenta+2, 'style': 29, 'size': 1.8, 'limitpoints': 50}, 
    {'run': 301298, 'name': '301298 (1548b)', 'color': ROOT.kRed+2, 'style': 33, 'size': 1.8, 'limitpoints': 50}, 
    {'run': 300806, 'name': '300806 (2556b)', 'color': ROOT.kCyan+1, 'style': 21, 'size': 1.5, 'limitpoints': 50},
    {'run': 300811, 'name': '300811 (2556b)', 'color': ROOT.kSpring+5, 'style': 22, 'size': 1.5, 'limitpoints': 50},
    ]

    #make_tgraph(graphOptions, "cluster_splitting_p75_run301417_v2.png", xrange=[0.5,1.75])

    graphOptions = [
    {'run': 301417, 'name': '301417 (1932b)', 'scalex': 1000.0/1919, 'color': ROOT.kBlack, 'style': 21, 'size': 1.5, 'limitpoints': 50, 'yrange': [0, 10], 'xtitle': 'lumi/#colliding_bunches [10^{31} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 299593, 'name': '299593 (2556b)', 'scalex': 1000.0/2544, 'color': ROOT.kGreen+2, 'style': 22, 'size': 1.5, 'limitpoints': 50}, 
    {'run': 301391, 'name': '301391 (1740b)', 'scalex': 1000.0/1727, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.5, 'limitpoints': 50} ,
    {'run': 300817, 'name': '300817 (2556b)', 'scalex': 1000.0/2542, 'color': ROOT.kOrange+2, 'style': 34, 'size': 1.5, 'limitpoints': 50}, 
    {'run': 301283, 'name': '301283 (1356b)', 'scalex': 1000.0/1343, 'color': ROOT.kMagenta+2, 'style': 29, 'size': 1.75, 'limitpoints': 50}, 
    {'run': 301298, 'name': '301298 (1548b)', 'scalex': 1000.0/1535, 'color': ROOT.kRed+2, 'style': 33, 'size': 1.5, 'limitpoints': 50}, 
    {'run': 300806, 'name': '300806 (2556b)', 'scalex': 1000.0/2542, 'color': ROOT.kCyan+1, 'style': 22, 'size': 1.5, 'limitpoints': 50},
    {'run': 300811, 'name': '300811 (2556b)', 'scalex': 1000.0/2542, 'color': ROOT.kSpring+5, 'style': 21, 'size': 1.5, 'limitpoints': 50},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_run301417_nbunchesrescaled_v1.png", xrange=[0.35,0.65])


    graphOptions = [
    {'run': 301417, 'name': '301417 (1932b)', 'split': 1, 'color': ROOT.kBlack, 'style': 21, 'yrange': [0, 15], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]'},
    {'run': 299593, 'name': '299593 (2556b)', 'split': 1, 'color': ROOT.kGreen+2, 'style': 22, 'size': 1.8}, 
    {'run': 301391, 'name': '301391 (1740b)', 'split': 1, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8} ,
    {'run': 300817, 'name': '300817 (2556b)', 'split': 1, 'color': ROOT.kOrange+2, 'style': 34, 'size': 1.8}, 
    {'run': 301283, 'name': '301283 (1356b)', 'split': 1, 'color': ROOT.kMagenta+2, 'style': 29, 'size': 1.8}, 
    {'run': 301298, 'name': '301298 (1548b)', 'split': 1, 'color': ROOT.kRed+2, 'style': 33, 'size': 1.8}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p76_run301417.png", xrange=[0.5,1.2])


    graphOptions = [
    {'run': 299593, 'name': '299593 (2556b 25ns)', 'color': ROOT.kBlack, 'style': 21, 'limitpoints': 500, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 299614, 'name': '299614 (1248b 50ns)', 'color': ROOT.kGreen+2, 'style': 22, 'size': 1.8, 'limitpoints': 50}, 
    {'run': 299616, 'name': '299616 (1248b 50ns)', 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 50} ,
    {'run': 299617, 'name': '299617 (1248b 50ns)', 'color': ROOT.kOrange+2, 'style': 34, 'size': 1.8, 'limitpoints': 50}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_25ns_50ns.png", xrange=[0.0,1.5])


    graphOptions = [
    {'run': 301627, 'color': ROOT.kBlack, 'style': 21, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 301664, 'color': ROOT.kBlue+2, 'style': 22, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 301694, 'color': ROOT.kOrange+2, 'style': 34, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75301627__301664_301694.png", xrange=[0.0,1.5])


    graphOptions = [
    {'run': 300233, 'name': '300233 Vdd=10', 'color': ROOT.kRed+2, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300234, 'name': '300234 Vdd=6', 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 300806, 'name': '300806 Vdd=6', 'color': ROOT.kGreen+2, 'style': 34, 'limitpoints': 50},
    {'run': 300811, 'name': '300811 Vdd=6', 'color': ROOT.kOrange+2, 'style': 29, 'size': 2.0, 'limitpoints': 50},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_300233vdig10_300234vdig6_300806vdig6_300811vdig6.png", xrange=[1.3,1.7])

    graphOptions = [
    {'run': 300233, 'split': 1, 'name': '300233 Vdd=10', 'color': ROOT.kRed+2, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]'},
    {'run': 300234, 'split': 1, 'name': '300234 Vdd=6', 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 300806, 'split': 1, 'name': '300806 Vdd=6', 'color': ROOT.kGreen+2, 'style': 34, 'limitpoints': 50},
    {'run': 300811, 'split': 1, 'name': '300811 Vdd=6', 'color': ROOT.kOrange+2, 'style': 29, 'size': 2.0, 'limitpoints': 50},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p76_300233vdig10_300234vdig6_300806vdig6_300811vdig6.png", xrange=[1.3,1.7])



    graphOptions = [
    {'run': 302019, 'color': ROOT.kRed+2, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 302026, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 301298, 'color': ROOT.kGreen+2, 'style': 34, 'limitpoints': 50},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_302019_302026_301298.png", xrange=[0.85,1.0])


    graphOptions = [
    {'run': 302228, 'name': '302228 (1356b)', 'color': ROOT.kRed+2, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 302262, 'name': '302262 (1550b)', 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 302277, 'name': '302277 (1164b)', 'color': ROOT.kGreen+2, 'style': 34, 'limitpoints': 200},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_302228_302262_302277.png", xrange=[0.5,0.95])

    graphOptions = [
    {'run': 302228, 'name': '302228 (1356b)', 'scalex': 1000.0/1344, 'color': ROOT.kRed+2, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 10], 'xtitle': 'lumi/#colliding_bunches [10^{31} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 302262, 'name': '302262 (1550b)', 'scalex': 1000.0/1538, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 302277, 'name': '302277 (1164b)', 'scalex': 1000.0/1152, 'color': ROOT.kGreen+2, 'style': 34, 'limitpoints': 200},
    ]
    # 'scalex': 1000.0/1919,
    #make_tgraph(graphOptions, "cluster_splitting_p75_302228_302262_302277_rescaled.png", xrange=[0.4,0.7])

    graphOptions = [
    {'run': 300806, 'name': '300806 ', 'color': ROOT.kBlue, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300811, 'name': '300811 ', 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 300817, 'name': '300817 ', 'color': ROOT.kAzure+7, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 303948, 'name': '303948 (8b4e)', 'color': ROOT.kRed+2, 'style': 34, 'limitpoints': 200},
    {'run': 303989, 'name': '303989 (8b4e)', 'color': ROOT.kRed, 'style': 33, 'limitpoints': 200},
    {'run': 304209, 'name': '304209 (8b4e)', 'color': ROOT.kGreen+2, 'style': 21, 'size': 1.1, 'limitpoints': 200},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_300806_300811_300817_303948_303989_304209.png", xrange=[0.0,1.8])

    graphOptions = [
    {'run': 303794, 'split': 1,  'name': '303794', 'color': ROOT.kGray, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [5, 10], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]'},
    {'run': 303794, 'split': 1,  'name': '303794 LS4-16', 'ls':[4,16], 'color': ROOT.kBlue, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 303794, 'split': 1,  'name': '303794 LS17-22', 'ls':[17,22], 'color': ROOT.kRed, 'style': 22, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p76_303794.png", xrange=[0.0298,0.0305])

    graphOptions = [
    {'run': 300155, 'name': '300155', 'color': ROOT.kBlue, 'style': 22, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300233, 'name': '300233', 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 300234, 'name': '300234', 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 300235, 'name': '300235', 'color': ROOT.kOrange+8, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_300155_300233_300234_300235.png", xrange=[0.0,1.8])

    graphOptions = [
    {'run': 304366, 'name': '304366', 'color': ROOT.kGray, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'yrange': [0, 35], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 304366, 'name': '304366 LS1380-1387 ref', 'ls':[1380,1387], 'color': ROOT.kBlack, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 304366, 'name': '304366 LS1392-1396 400V', 'ls':[1392,1396], 'color': ROOT.kBlue, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 304366, 'name': '304366 LS1398-1402 300V', 'ls':[1398,1402], 'color': ROOT.kGreen+2, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 304366, 'name': '304366 LS1403-1407 250V', 'ls':[1403,1407], 'color': ROOT.kRed+2, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 304366, 'name': '304366 LS1408-1412 200V', 'ls':[1408,1412], 'color': ROOT.kMagenta, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 304366, 'name': '304366 LS1413-1417 150V', 'ls':[1413,1417], 'color': ROOT.kOrange+2, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    {'run': 304366, 'name': '304366 LS1418-1422 100V', 'ls':[1418,1422], 'color': ROOT.kCyan+1, 'style': 22, 'size': 1.8, 'limitpoints': 200},
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_304366_hv.png", xrange=[0.775,0.825])
    
    graphOptions = [
    {'run': 304446, 'color': ROOT.kBlue, 'style': 22, 'size': 1.6, 'limitpoints': 200, 'yrange': [5,10], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 304447, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304505, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304506, 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304507, 'color': ROOT.kOrange+8, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304508, 'color': ROOT.kOrange+2, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6272and6275_zoom.png", xrange=[1.2,1.6])
    
    graphOptions = [
    {'run': 304446, 'color': ROOT.kBlue, 'style': 22, 'size': 1.6, 'limitpoints': 200, 'yrange': [0,11], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 304447, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6272.png", xrange=[0,1.8])
    
    graphOptions = [
    {'run': 304505, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'yrange': [0,11], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    {'run': 304506, 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304507, 'color': ROOT.kOrange+8, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304508, 'color': ROOT.kOrange+2, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6275.png", xrange=[0,1.8])
    
    graphOptions = [
    {'run': 304446, 'split':1, 'color': ROOT.kBlue, 'style': 22, 'size': 1.6, 'limitpoints': 200, 'yrange': [4, 12], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->6) [%]'},
    {'run': 304447, 'split':1, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'split':1, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'split':1, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'split':1, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304505, 'split':1, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304506, 'split':1, 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304507, 'split':1, 'color': ROOT.kOrange+8, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304508, 'split':1, 'color': ROOT.kOrange+2, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p76_fill6272and6275.png", xrange=[0.55,1.8])

    graphOptions = [
    {'run': 300806, 'color': ROOT.kGray, 'style': 22, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300811, 'color': ROOT.kGray, 'style': 23, 'size': 1.6, 'limitpoints': 200},
    {'run': 300817, 'color': ROOT.kGray, 'style': 33, 'size': 1.6, 'limitpoints': 200},
    {'run': 304446, 'color': ROOT.kBlue, 'style': 29, 'size': 1.6, 'limitpoints': 200},
    {'run': 304447, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304505, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304506, 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304507, 'color': ROOT.kOrange+8, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304508, 'color': ROOT.kOrange+2, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6061-6272-6275.png", xrange=[0.6,1.8])
    
    graphOptions = [
    {'run': 304446, 'color': ROOT.kBlue, 'style': 29, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 11], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 304447, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304562, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6272-6276.png", xrange=[0.6,1.8])
   
    
    graphOptions = [
    {'run': 304505, 'color': ROOT.kBlue, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 11], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    {'run': 304506, 'color': ROOT.kBlue+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304507, 'color': ROOT.kAzure+8, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304508, 'color': ROOT.kAzure+3, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304562, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6275-6276.png", xrange=[0.6,1.8])
   
    graphOptions = [
    {'run': 300806, 'color': ROOT.kGray, 'style': 22, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 300811, 'color': ROOT.kGray, 'style': 23, 'size': 1.6, 'limitpoints': 200},
    {'run': 300817, 'color': ROOT.kGray, 'style': 33, 'size': 1.6, 'limitpoints': 200},
    {'run': 304446, 'color': ROOT.kBlue, 'style': 29, 'size': 1.6, 'limitpoints': 200},
    {'run': 304447, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304562, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6061-6272-6276.png", xrange=[0.6,1.8])
   
    graphOptions = [
    {'run': 301417, 'name': '301417 (1932b)', 'scalex': 1000.0/1919, 'color': ROOT.kBlack, 'style': 21, 'size': 1.5, 'limitpoints': 50, 'yrange': [0, 30], 'xtitle': 'lumi/#colliding_bunches [10^{31} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    {'run': 299593, 'name': '299593 (2556b)', 'scalex': 1000.0/2544, 'color': ROOT.kCyan+3, 'style': 22, 'size': 1.5, 'limitpoints': 50}, 
    {'run': 301391, 'name': '301391 (1740b)', 'scalex': 1000.0/1727, 'color': ROOT.kGreen+3, 'style': 23, 'size': 1.5, 'limitpoints': 50} ,
    {'run': 301283, 'name': '301283 (1356b)', 'scalex': 1000.0/1343, 'color': ROOT.kTeal+3, 'style': 29, 'size': 1.75, 'limitpoints': 50}, 
    {'run': 300806, 'name': '300806 (2556b)', 'scalex': 1000.0/2542, 'color': ROOT.kGray, 'style': 22, 'size': 1.6, 'limitpoints': 50},
    {'run': 300811, 'name': '300811 (2556b)', 'scalex': 1000.0/2542, 'color': ROOT.kGray, 'style': 23, 'size': 1.6, 'limitpoints': 50},
    {'run': 300817, 'name': '300817 (2556b)', 'scalex': 1000.0/2542, 'color': ROOT.kGray, 'style': 33, 'size': 1.6, 'limitpoints': 50},
    {'run': 304446, 'name': '304446 (1868b)', 'scalex': 1000.0/1866, 'color': ROOT.kBlue, 'style': 29, 'size': 1.6, 'limitpoints': 200},
    {'run': 304447, 'name': '304446 (1868b)', 'scalex': 1000.0/1866, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200}, 
    {'run': 304448, 'name': '304448 (1868b)', 'scalex': 1000.0/1866, 'color': ROOT.kAzure+8, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 304451, 'name': '304451 (1868b)', 'scalex': 1000.0/1866, 'color': ROOT.kAzure+3, 'style': 33, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304452, 'name': '304452 (1868b)', 'scalex': 1000.0/1866, 'color': ROOT.kAzure-7, 'style': 34, 'size': 2.0, 'limitpoints': 200}, 
    {'run': 304562, 'name': '304562 (1868b)', 'scalex': 1000.0/1866, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 50}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_lumi_per_bunch_august_october.png", xrange=[0.0,1.0])

    #make_tgraph_bx({'run': 305063, 'color': ROOT.kRed, 'style': 23, 'size': 1.6, 'limitpoints': 50})
    #make_tgraph_bx({'run': 305059, 'ls': [63, 572], 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6})
    #make_tgraph_bx({'run': 305208, 'ls': [1, 373], 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6})
    #make_tgraph_bx({'run': 304200, 'ls': [1, 322], 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6})
    #make_tgraph_bx({'run': 300806, 'ls': [1,9999], 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6})
    #make_tgraph_bx({'run': 300817, 'ls': [1,9999], 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6})

    graphOptions = [
    {'run': 302479, 'color': ROOT.kBlue, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 11], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    {'run': 303832, 'color': ROOT.kRed, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_302479_303832.png", xrange=[0.6,1.4])

    graphOptions = [
    {'run': 305588, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 11], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    {'run': 305589, 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 305590, 'color': ROOT.kOrange+2, 'style': 33, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6324.png", xrange=[0.4,1.8])

    graphOptions = [
    {'run': 305747, 'color': ROOT.kRed+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 20], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    {'run': 305756, 'color': ROOT.kGreen+2, 'style': 33, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6336.png", xrange=[0.08,0.12])
    
    graphOptions = [
    {'run': 305747, 'color': ROOT.kRed+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'xrange':[0.1, 0.11], 'yrange': [0, 5.0e-6], 'xtitle': 'LS', 'ytitle': 'fraction of >= 36 frames'}, 
    {'run': 305756, 'color': ROOT.kGreen+2, 'style': 33, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph_bigclusters(graphOptions, "fraction_big_clusters_305747_305756.png")

    graphOptions = [
    {'run': 305590, 'color': ROOT.kBlue+2, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'xrange':[0,1.8], 'yrange': [0, 3.0e-5], 'xtitle': 'LS', 'ytitle': 'fraction of >= 36 frames'}, 
    {'run': 305589, 'color': ROOT.kBlue, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 305588, 'color': ROOT.kAzure+8, 'style': 22, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 305747, 'color': ROOT.kRed+2, 'style': 34, 'size': 1.6, 'limitpoints': 200},
    {'run': 305756, 'color': ROOT.kGreen+2, 'style': 33, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph_bigclusters(graphOptions, "fraction_big_clusters_305588_89_90.png")
    
    graphOptions = [
    {'run': 305747,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '1.5 ns', 'color': ROOT.kBlack, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'LSoffset': 40, 'yrange': [0, 6], 'xtitle': 'LS', 'ytitle': 'LYR2 cluster splitting p(5->3) [%]'}, 
    {'run': 305748,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '0 ns', 'color': ROOT.kRed+2, 'style': 29, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 30}, 
    {'run': 305749,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '-2 ns', 'color': ROOT.kOrange+2, 'style': 33, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 20}, 
    {'run': 305750,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '-4 ns', 'color': ROOT.kBlue+2, 'style': 32, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 10}, 
    {'run': 305751,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '-6 ns', 'color': ROOT.kGreen+2, 'style': 34, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 0}, 
    {'run': 305752,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '3 ns', 'color': ROOT.kAzure+8, 'style': 22, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 50}, 
    {'run': 305753,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '5 ns', 'color': ROOT.kMagenta+2, 'style': 23, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 60}, 
    {'run': 305754,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '7 ns', 'color': ROOT.kGray, 'style': 34, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 70}, 
    {'run': 305755,'selection': lambda x: x.size == 5 and '_LYR2_' in x.module, 'name': '9 ns', 'color': ROOT.kTeal+3, 'style': 29, 'size': 1.8, 'limitpoints': 200, 'LSoffset': 80}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p53_timingscan_LYR2.png", xrange=[0,90])
    
    graphOptions = [
    {'run': 305757, 'name': 'HV scan', 'color': ROOT.kBlack, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'LSoffset': 0, 'yrange': [0, 30], 'xtitle': 'LS', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_hvscan.png", xrange=[0,60])

    graphOptions = [
    {'run': 306091, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module, 'color': ROOT.kBlue, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'LYR2 cluster splitting p(7->5) [%]'}, 
    {'run': 306092, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module,  'color': ROOT.kBlue+2, 'style': 29, 'size': 1.8, 'limitpoints': 200}, 
    {'run': 306095, 'selection': lambda x: x.size == 7 and '_LYR2_' in x.module,  'color': ROOT.kAzure+8, 'style': 22, 'size': 1.8, 'limitpoints': 200}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6358_LYR2.png", xrange=[0.4,2.2])

    graphOptions = [
    {'run': 306126, 'color': ROOT.kBlue, 'style': 23, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6360.png", xrange=[0.0,0.5])
    
    graphOptions = [
    {'run': 306432, 'color': ROOT.kBlue+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6370.png", xrange=[0,0.2])

    graphOptions = [
    {'run': 306572, 'color': ROOT.kBlue+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6384.png", xrange=[0,0.16])
    
    graphOptions = [
    {'run': 306896, 'color': ROOT.kBlue+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 3], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6402.png", xrange=[0.001,0.0015])

    #graphOptions = [
    #{'run': 307042, 'color': ROOT.kBlue+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    #]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6413.png", xrange=[0,1.8])

    #graphOptions = [
    #{'run': 307042, 'color': ROOT.kBlue+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'lsplot': True,  'xtitle': 'LS', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    #]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6413_ls.png")

    graphOptions = [
    {'run': 307076, 'color': ROOT.kBlue+2, 'style': 34, 'size': 1.6, 'limitpoints': 200, 'yrange': [0, 30], 'lsplot': True,  'xtitle': 'LS', 'ytitle': 'cluster splitting p(7->5) [%]'}, 
    ]
    #make_tgraph(graphOptions, "cluster_splitting_p75_fill6415_ls.png")

    graphOptions = [
      {'run': 314890, 'color': ROOT.kBlue+2, 'style': 22, 'size': 1.8, 'limitpoints': 200,   'yrange': [0, 10], 'xtitle': 'lumi [10^{34} cm^{-2}s^{-1}]', 'ytitle': 'cluster splitting p(7->5) [%]'},
    ]
    make_tgraph(graphOptions, "cluster_splitting_p75_314890.png", xrange=[0.0,0.3])

