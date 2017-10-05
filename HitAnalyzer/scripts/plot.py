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


def get_lumisection_data(fileName, options):
    #print "FILE:", fileName
    path = '/'.join(fileName.split('/')[0:-1]).strip()
    runNumber = int(fileName.split('/')[-1].split('_')[1])
    lumisectionNumber = int(fileName.split('/')[-1].split('_')[2])

    # get instantaneous lumi 
    lumiCSVfileName = 'lumi%d.csv'%runNumber
    if len(path) > 0:
        lumiCSVfileName = path + '/' + lumiCSVfileName
    instaLumi = 0
    with open(lumiCSVfileName, 'r') as lumiCSVfile:
        for line in lumiCSVfile:
            values = [x.strip() for x in line.split(',')]
            try:
                if int(values[0]) == lumisectionNumber:
                    instaLumi = float(values[3]) * 1.0e-4
                    break
            except:
                pass
    #print "Instantaneous luminosity: %f"%instaLumi
   
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
    return [instaLumi, fraction*100.0, 0, error*100.0]
   
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
   
    tg = ROOT.TGraphErrors(len(graphPointsList), x_points, y_points, x_errors, y_errors)
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
    make_tgraph(graphOptions, "cluster_splitting_p75_300155_300233_300234_300235.png", xrange=[0.0,1.8])

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
