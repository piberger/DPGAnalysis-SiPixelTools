// -*- C++ -*-
//
// Package:    PixDigisTest
// Class:      PixDigisTest
// 
/**\class PixDigisTest PixDigisTest.cc 

 Description: Test pixel digis. 
 Barrel & Forward digis. Uses root histos.
 Adopted for the new simLinks. 
 Added the online detector index. d.k. 11/09
 Works with CMSSW_7
 New detector ID.
 Modified to use "byToken"
 Modify for 4 layers.
*/
//
// Original Author:  d.k.
//         Created:  Jan CET 2006
//
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

// my includes
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelNameUpgrade.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapNameUpgrade.h"

// data formats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
//#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

// For the big pixel recongnition
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

// for simulated Tracker hits
//#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

// For L1
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"

// For HLT
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"


// To use root histos
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


// For ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>

#define HISTOS
#define L1
#define HLT

#include "DataFormats/TCDS/interface/TCDSRecord.h"

using namespace std;

// Enable this to look at simlinks (link simhit->digis)
// Can be used only with simulated data.
//#define USE_SIM_LINKS


//
// class decleration
//

class PixDigisTest : public edm::EDAnalyzer {

public:

  explicit PixDigisTest(const edm::ParameterSet&);
  ~PixDigisTest();
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob(); 

private:
  // ----------member data ---------------------------
  bool PRINT;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> tPixelDigi;
#ifdef USE_SIM_LINKS
  edm::EDGetTokenT< edm::DetSetVector<PixelDigiSimLink> > tPixelDigiSimLink;
#endif
  edm::EDGetTokenT<TCDSRecord> tcdsrecord_; 
  float count0, count1, count2, count3;

#ifdef HISTOS

  TFile* hFile;
  TH1F *hdetunit;
  TH1F *heloss1,*heloss2, *heloss3, *heloss4;
  TH1F *hneloss1,*hneloss2, *hneloss3, *hneloss4;
  TH1F *helossF1,*helossF2,*helossF3;
  TH1F *hpixid,*hpixsubid,*hlayerid,*hshellid,*hsectorid,
    *hladder1id,*hladder2id,*hladder3id,*hladder4id,
    *hz1id,*hz2id,*hz3id,*hz4id;
  TH1F *hcols1,*hcols2,*hcols3,*hrows1,*hrows2,*hrows3,*hcols4,*hrows4;
  TH1F *hcolsF1,*hcolsF2,*hcolsF3,*hrowsF1,*hrowsF2,*hrowsF3;
  TH1F *hdigisPerDet1,*hdigisPerDet2,*hdigisPerDet3,*hdigisPerDet4;
  TH1F *hdigisPerLay1,*hdigisPerLay2,*hdigisPerLay3,*hdigisPerLay4;
  TH1F *hdetsPerLay1,*hdetsPerLay2,*hdetsPerLay3,*hdetsPerLay4;
  TH1F *hdigisPerDetF1,*hdigisPerDetF2,*hdigisPerDetF3;
  TH1F *hdigisPerLayF1,*hdigisPerLayF2,*hdigisPerLayF3;
  TH1F *hdetsPerLayF1,*hdetsPerLayF2,*hdetsPerLayF3;
  TH1F *hdetr, *hdetz, *hdetrF, *hdetzF;
  TH1F *hcolsB,  *hrowsB,  *hcolsF,  *hrowsF;
  TH1F *hcols1big, *hrows1big, *heloss1bigx, *heloss1bigy;
  TH1F *hsimlinks, *hfract;
  TH1F *hblade1, *hblade2, *hblade3;
  TH1F *hClusterSplit1Row1, *hClusterSplit2Row1, *hClusterSplit1Col1, *hClusterSplit2Col1;
  //TH2F *htest, *htest2;
  TH2F *hpdetMap4,*hpdetMap3,*hpdetMap2,*hpdetMap1; 
  TH2F *hdetMap4,*hdetMap3,*hdetMap2,*hdetMap1, 
    *hpixMap1, *hpixMap2, *hpixMap3,*hpixMap4;
  TH2F *hClusterSplit1Pos1;
  TH2F *hClusterSplit2Pos1;

  TH2F *hpdetMaps4,*hpdetMaps3,*hpdetMaps2,*hpdetMaps1; 

  //* hpixMapNoise; 
  TH2F *hxy, *hphiz1, *hphiz2, *hphiz3, *hphiz4; // bpix 
  TH2F *hzr, *hxy11, *hxy12, *hxy21, *hxy22, *hxy31, *hxy32;  // fpix 

  TH1F *hevent, *hlumi, *horbit, *hbx0, *hlumi0, *hlumi1,*hbx1,*hbx2,*hbx3,*hbx4,*hbx5,*hbx6,*hbx;
  TH1F *hdets, *hdigis, *hdigis0, *hdigis1, *hdigis2,*hdigis3,*hdigis4,*hdigis5; 

  TH1F *hNdbCellsLay1;
  TH1F *adcLay1Big, *adcLay1Small;

  TProfile *hadc1ls,*hadc2ls,*hadc3ls,*hadc4ls,*hadc0ls; 
  TProfile *hadc1bx,*hadc2bx,*hadc3bx,*hadc4bx,*hadc0bx; 

#endif

  // test!
  char digiTestModuleName[37];
  int digiTestEvent;
  int digiTestNdigis;
  int digiTestDigiRow[66000];
  int digiTestLS;
  int digiTestDigiCol[66000];
  int digiTestDigiAdc[66000];
  TTree* digiTestTree;
  edm::InputTag src_;  
  bool phase1_;

  // custom clustering
  char clustersModuleName[37];
  int clustersSizeBranch;
  int clustersBxBranch;
  int clustersNcomplete;
  int clustersNbroken1;
  int clustersNbroken2;
  
  // cluster tree -> event info
  unsigned int clusterTree_bx;
  unsigned int clusterTree_event;
  unsigned int clusterTree_lumiblock;
  unsigned int clusterTree_orbit;
  unsigned int clusterTree_orbitOfLastBgo14;
  int clusterTree_tres;


  // cluster tree -> module position
  int clusterTree_shell;
  int clusterTree_sector;
  int clusterTree_ladder;
  int clusterTree_layer;
  int clusterTree_module;
  unsigned int clusterTree_layerC;
  unsigned int clusterTree_ladderC;
  unsigned int clusterTree_zindex;
  unsigned int clusterTree_pos_x;
  unsigned int clusterTree_pos_y;

  // cluster tree -> clusters
  int clusterTree_size;
  int clusterTree_x[480];
  int clusterTree_y[480];
  int clusterTree_adc[480];

  // count clusters for each module < <module, size>, count>
  std::map< std::pair<std::string, int>, int> clustersComplete;
  std::map< std::pair<std::string, int>, int> clustersSplitByOne;
  std::map< std::pair<std::string, int>, int> clustersSplitByTwo;

  // count clusters per BX: < <bx, size>, count>
  std::map< std::pair<int, int>, int> clustersBXComplete;
  std::map< std::pair<int, int>, int> clustersBXSplitByOne;
  std::map< std::pair<int, int>, int> clustersBXSplitByTwo;
  
std::map< std::string, int> moduleNhits;
  std::map< std::string, int> moduleNclusters;
  TTree* clusterTree;
  int restrictToLumisection;
  bool includeZeroAdc;


//  std::vector< std::vector<int> > doubleColumnHitsL1;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PixDigisTest::PixDigisTest(const edm::ParameterSet& iConfig) {

  PRINT = iConfig.getUntrackedParameter<bool>("Verbosity",false);
  src_ =  iConfig.getParameter<edm::InputTag>( "src" );
  tPixelDigi = consumes <edm::DetSetVector<PixelDigi>> (src_);
#ifdef USE_SIM_LINKS
  tPixelDigiSimLink = consumes < edm::DetSetVector<PixelDigiSimLink> > ( src_);
#endif 
  phase1_ =  iConfig.getUntrackedParameter<bool>( "phase1",false);
  restrictToLumisection = iConfig.getUntrackedParameter<int>("lumisection", -1);
  includeZeroAdc = iConfig.getUntrackedParameter<bool>("includezeroadc", true);

  tcdsrecord_ = consumes<TCDSRecord>(edm::InputTag("tcdsRawToDigiProducer","tcdsRecord"));
  cout<<" Construct PixDigisTest "<<endl;
}


PixDigisTest::~PixDigisTest() {
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  cout<<" Destroy PixDigisTest "<<endl;

}

//
// member functions
//
// ------------ method called at the begining   ------------
void PixDigisTest::beginJob() {

   using namespace edm;
   cout << "Initialize PixDigisTest " <<endl;
   count0=count1=count2=count3=0;

#ifdef HISTOS
  edm::Service<TFileService> fs;

  // Histos go to a subdirectory "PixRecHits")
  //TFileDirectory subDir = fs->mkdir( "mySubDirectory" );
  //TFileDirectory subSubDir = subDir.mkdir( "mySubSubDirectory" );

   // put here whatever you want to do at the beginning of the job
    //hFile = new TFile ( "digis_tree.root", "RECREATE" );
   
    clusterTree = fs->make<TTree>("tree","clusters");
    //clusterTree->Branch("module",(void*)clustersModuleName,"string/C",36);
    //clusterTree->Branch("size",&clustersSizeBranch,"size/I");
    //clusterTree->Branch("complete",&clustersNcomplete,"complete/I");
    //clusterTree->Branch("broken1",&clustersNbroken1,"broken1/I");
    //clusterTree->Branch("broken2",&clustersNbroken2,"broken2/I");
    
    // event information
    clusterTree->Branch("bx", &clusterTree_bx, "bx/i");
    clusterTree->Branch("event", &clusterTree_event, "bx/i");
    clusterTree->Branch("lumiblock", &clusterTree_lumiblock, "lumiblock/i");
    clusterTree->Branch("orbit", &clusterTree_orbit, "orbit/i");
    clusterTree->Branch("orbitOfLastBgo14", &clusterTree_orbitOfLastBgo14, "orbitOfLastBgo14/i");
    clusterTree->Branch("tres", &clusterTree_tres, "tres/I");

    // module information
    clusterTree->Branch("shell", &clusterTree_shell, "shell/I");
    clusterTree->Branch("sector", &clusterTree_sector, "sector/I");
    clusterTree->Branch("ladder", &clusterTree_ladder, "ladder/I");
    clusterTree->Branch("layer", &clusterTree_layer, "layer/I");
    clusterTree->Branch("module", &clusterTree_module, "module/I");
    clusterTree->Branch("layerC", &clusterTree_layerC, "layerC/i");
    clusterTree->Branch("ladderC", &clusterTree_ladderC, "ladderC/i");
    clusterTree->Branch("zindex", &clusterTree_zindex, "zindex/i");
    
    // clusters
    clusterTree->Branch("pos_x", &clusterTree_pos_x, "pos_x/i");
    clusterTree->Branch("pos_y", &clusterTree_pos_y, "pos_y/i");
    clusterTree->Branch("size", &clusterTree_size,"size/I");
    
    // pixels
    clusterTree->Branch("x", clusterTree_x, "x[size]/i");
    clusterTree->Branch("y", clusterTree_y, "y[size]/i");
    clusterTree->Branch("adc", clusterTree_adc, "adc[size]/i");

    // < <module_name, cluster_length>, count >

    hdetunit = fs->make<TH1F>( "hdetunit", "Det unit", 1000,
                              302000000.,302300000.);
    hpixid = fs->make<TH1F>( "hpixid", "Pix det id", 10, 0., 10.);
    hpixsubid = fs->make<TH1F>( "hpixsubid", "Pix Barrel id", 10, 0., 10.);
    hlayerid = fs->make<TH1F>( "hlayerid", "Pix layer id", 10, 0., 10.);
    hsectorid = fs->make<TH1F>( "hsectorid", "Pix sector ", 10, 0., 10.);
    hshellid = fs->make<TH1F>( "hshellid", "Shell", 5, 0., 5.);
    hladder1id = fs->make<TH1F>( "hladder1id", "Ladder L1 id", 23, -11.5, 11.5);
    hladder2id = fs->make<TH1F>( "hladder2id", "Ladder L2 id", 35, -17.5, 17.5);
    hladder3id = fs->make<TH1F>( "hladder3id", "Ladder L3 id", 47, -23.5, 23.5);
    hladder4id = fs->make<TH1F>( "hladder4id", "Ladder L4 id", 67, -33.5, 33.5);
    hz1id = fs->make<TH1F>( "hz1id", "Z-index id L1", 11, -5.5, 5.5);
    hz2id = fs->make<TH1F>( "hz2id", "Z-index id L2", 11, -5.5, 5.5);
    hz3id = fs->make<TH1F>( "hz3id", "Z-index id L3", 11, -5.5, 5.5);
    hz4id = fs->make<TH1F>( "hz4id", "Z-index id L4", 11, -5.5, 5.5);
 
    hdigisPerDet1 = fs->make<TH1F>( "hdigisPerDet1", "Digis per det l1", 
			      400, -0.5, 399.5);
    hdigisPerDet2 = fs->make<TH1F>( "hdigisPerDet2", "Digis per det l2", 
			      200, -0.5, 199.5);
    hdigisPerDet3 = fs->make<TH1F>( "hdigisPerDet3", "Digis per det l3", 
			      200, -0.5, 199.5);
    hdigisPerDet4 = fs->make<TH1F>( "hdigisPerDet4", "Digis per det l4", 
			      200, -0.5, 199.5);
    //
    //const float maxSize=199.5;
    const float maxSize=9999.5;
    hdigisPerLay1 = fs->make<TH1F>( "hdigisPerLay1", "Digis per layer l1", 
			      100, -0.5, maxSize);
    hdigisPerLay2 = fs->make<TH1F>( "hdigisPerLay2", "Digis per layer l2", 
			      100, -0.5, maxSize);
    hdigisPerLay3 = fs->make<TH1F>( "hdigisPerLay3", "Digis per layer l3", 
			      100, -0.5, maxSize);
    hdigisPerLay4 = fs->make<TH1F>( "hdigisPerLay4", "Digis per layer l4", 
			      100, -0.5, maxSize);
    hdetsPerLay1 = fs->make<TH1F>( "hdetsPerLay1", "Full dets per layer l1", 
			      161, -0.5, 160.5);
    hdetsPerLay2 = fs->make<TH1F>( "hdetsPerLay2", "Full dets per layer l2", 
			      257, -0.5, 256.5);
    hdetsPerLay3 = fs->make<TH1F>( "hdetsPerLay3", "Full dets per layer l3", 
			      353, -0.5, 352.5);
    hdetsPerLay4 = fs->make<TH1F>( "hdetsPerLay4", "Full dets per layer l4", 
			      513, -0.5, 512.5);

    hdigisPerDetF1 = fs->make<TH1F>( "hdigisPerDetF1", "Digis per det d1", 
			      200, -0.5, 199.5);
    hdigisPerDetF2 = fs->make<TH1F>( "hdigisPerDetF2", "Digis per det d2", 
			      200, -0.5, 199.5);
    hdigisPerDetF3 = fs->make<TH1F>( "hdigisPerDetF3", "Digis per det d3", 
			      200, -0.5, 199.5);
    hdigisPerLayF1 = fs->make<TH1F>( "hdigisPerLayF1", "Digis per layer d1", 
			      2000, -0.5, 1999.5);
    hdigisPerLayF2 = fs->make<TH1F>( "hdigisPerLayF2", "Digis per layer d2", 
			      2000, -0.5, 1999.5);
    hdigisPerLayF3 = fs->make<TH1F>( "hdigisPerLayF3", "Digis per layer d3", 
			      2000, -0.5, 1999.5);
    hdetsPerLayF1 = fs->make<TH1F>( "hdetsPerLayF1", "Full dets per layer d1", 
			      161, -0.5, 160.5);
    hdetsPerLayF2 = fs->make<TH1F>( "hdetsPerLayF2", "Full dets per layer d2", 
			      257, -0.5, 256.5);
    hdetsPerLayF3 = fs->make<TH1F>( "hdetsPerLayF3", "Full dets per layer d3", 
			      257, -0.5, 256.5);

    heloss1 = fs->make<TH1F>( "heloss1", "Pix charge l1", 256, 0., 256.);
    heloss2 = fs->make<TH1F>( "heloss2", "Pix charge l2", 256, 0., 256.);
    heloss3 = fs->make<TH1F>( "heloss3", "Pix charge l3", 256, 0., 256.);
    heloss4 = fs->make<TH1F>( "heloss4", "Pix charge l4", 256, 0., 256.);
    hneloss1 = fs->make<TH1F>( "hneloss1", "Pix noise charge l1", 256, 0., 256.);
    hneloss2 = fs->make<TH1F>( "hneloss2", "Pix noise charge l2", 256, 0., 256.);
    hneloss3 = fs->make<TH1F>( "hneloss3", "Pix noise charge l3", 256, 0., 256.);
    hneloss4 = fs->make<TH1F>( "hneloss4", "Pix noise charge l4", 256, 0., 256.);
    heloss1bigx = fs->make<TH1F>( "heloss1bigx", "L1 big-x pix", 256, 0., 256.);
    heloss1bigy = fs->make<TH1F>( "heloss1bigy", "L1 big-y pix", 256, 0., 256.);

    hcols1 = fs->make<TH1F>( "hcols1", "Layer 1 cols", 500,-1.5,498.5);
    hcols2 = fs->make<TH1F>( "hcols2", "Layer 2 cols", 500,-1.5,498.5);
    hcols3 = fs->make<TH1F>( "hcols3", "Layer 3 cols", 500,-1.5,498.5);
    hcols4 = fs->make<TH1F>( "hcols4", "Layer 4 cols", 500,-1.5,498.5);
    hcols1big = fs->make<TH1F>( "hcols1big", "Layer 1 big cols", 500,-1.5,498.5);
 
    hrows1 = fs->make<TH1F>( "hrows1", "Layer 1 rows", 200,-1.5,198.5);
    hrows2 = fs->make<TH1F>( "hrows2", "Layer 2 rows", 200,-1.5,198.5);
    hrows3 = fs->make<TH1F>( "hrows3", "layer 3 rows", 200,-1.5,198.5);
    hrows4 = fs->make<TH1F>( "hrows4", "layer 4 rows", 200,-1.5,198.5);
    hrows1big = fs->make<TH1F>( "hrows1big", "Layer 1 big rows", 200,-1.5,198.5);
 
    hblade1 = fs->make<TH1F>( "hblade1", "blade num, disk1", 60, 0., 60.);
    hblade2 = fs->make<TH1F>( "hblade2", "blade num, disk2", 60, 0., 60.);
    hblade3 = fs->make<TH1F>( "hblade3", "blade num, disk3", 60, 0., 60.);

    helossF1 = fs->make<TH1F>( "helossF1", "Pix charge d1", 100, 0., 300.);
    helossF2 = fs->make<TH1F>( "helossF2", "Pix charge d2", 100, 0., 300.);
    helossF3 = fs->make<TH1F>( "helossF3", "Pix charge d3", 100, 0., 300.);
    hcolsF1 = fs->make<TH1F>( "hcolsF1", "Disk 1 cols", 500,-1.5,498.5);
    hcolsF2 = fs->make<TH1F>( "hcolsF2", "Disk 2 cols", 500,-1.5,498.5);
    hcolsF3 = fs->make<TH1F>( "hcolsF3", "Disk 3 cols", 500,-1.5,498.5);
    hrowsF1 = fs->make<TH1F>( "hrowsF1", "Disk 1 rows", 200,-1.5,198.5);
    hrowsF2 = fs->make<TH1F>( "hrowsF2", "Disk 2 rows", 200,-1.5,198.5);
    hrowsF3 = fs->make<TH1F>( "hrowsF3", "Disk 3 rows", 200,-1.5,198.5);

    hdetr = fs->make<TH1F>("hdetr","det r",180,0.,18.);
    hdetz = fs->make<TH1F>("hdetz","det z",520,-26.,26.);
    hdetrF = fs->make<TH1F>("hdetrF","det r",150,0.,15.);
    hdetzF = fs->make<TH1F>("hdetzF","det z",700,-70.,70.);

    hcolsB = fs->make<TH1F>("hcolsB","cols per bar det",450,0.,450.);
    hrowsB = fs->make<TH1F>("hrowsB","rows per bar det",200,0.,200.);
    hcolsF = fs->make<TH1F>("hcolsF","cols per for det",450,0.,450.);
    hrowsF = fs->make<TH1F>("hrowsF","rows per for det",200,0.,200.);

    hsimlinks = fs->make<TH1F>("hsimlinks"," track ids",200,0.,200.);
    hfract = fs->make<TH1F>("hfract"," track rractions",100,0.,1.);

    hClusterSplit1Row1 = fs->make<TH1F>("hClusterSplit1Row1","rows per bar det",200,0.,200.);
    hClusterSplit2Row1 = fs->make<TH1F>("hClusterSplit2Row1","rows per bar det",200,0.,200.);
    hClusterSplit1Col1 = fs->make<TH1F>("hClusterSplit1Col1","cols per bar det",450,0.,450.);
    hClusterSplit2Col1 = fs->make<TH1F>("hClusterSplit2Col1","cols per bar det",450,0.,450.);
    //                                             mod      ladder
    hdetMap1 = fs->make<TH2F>("hdetMap1"," ",9,-4.5,4.5,13,-6.5,6.5);
    hdetMap1->SetOption("colz");
    hdetMap2 = fs->make<TH2F>("hdetMap2"," ",9,-4.5,4.5,29,-14.5,14.5);
    hdetMap2->SetOption("colz");
    hdetMap3 = fs->make<TH2F>("hdetMap3"," ",9,-4.5,4.5,45,-22.5,22.5);
    hdetMap3->SetOption("colz");
    hdetMap4 = fs->make<TH2F>("hdetMap4"," ",9,-4.5,4.5,65,-32.5,32.5);
    hdetMap4->SetOption("colz");

    hpdetMap1 = fs->make<TH2F>("hpdetMap1"," ",9,-4.5,4.5,13,-6.5,6.5);
    hpdetMap1->SetOption("colz");
    hpdetMap2 = fs->make<TH2F>("hpdetMap2"," ",9,-4.5,4.5,29,-14.5,14.5);
    hpdetMap2->SetOption("colz");
    hpdetMap3 = fs->make<TH2F>("hpdetMap3"," ",9,-4.5,4.5,45,-22.5,22.5);
    hpdetMap3->SetOption("colz");
    hpdetMap4 = fs->make<TH2F>("hpdetMap4"," ",9,-4.5,4.5,65,-32.5,32.5);
    hpdetMap4->SetOption("colz");

    hpixMap1 = fs->make<TH2F>("hpixMap1"," ",416,0.,416.,160,0.,160.);
    hpixMap1->SetOption("colz");
    hpixMap2 = fs->make<TH2F>("hpixMap2"," ",416,0.,416.,160,0.,160.);
    hpixMap2->SetOption("colz");
    hpixMap3 = fs->make<TH2F>("hpixMap3"," ",416,0.,416.,160,0.,160.);
    hpixMap3->SetOption("colz");
    hpixMap4 = fs->make<TH2F>("hpixMap4"," ",416,0.,416.,160,0.,160.);
    hpixMap4->SetOption("colz");

    hpdetMaps1 = fs->make<TH2F>("hpdetMaps1"," slected hits ",9,-4.5,4.5,13,-6.5,6.5);
    hpdetMaps1->SetOption("colz");
    hpdetMaps2 = fs->make<TH2F>("hpdetMaps2","  slected hits ",9,-4.5,4.5,29,-14.5,14.5);
    hpdetMaps2->SetOption("colz");
    hpdetMaps3 = fs->make<TH2F>("hpdetMaps3","  slected hits ",9,-4.5,4.5,45,-22.5,22.5);
    hpdetMaps3->SetOption("colz");
    hpdetMaps4 = fs->make<TH2F>("hpdetMaps4","  slected hits ",9,-4.5,4.5,65,-32.5,32.5);
    hpdetMaps4->SetOption("colz");

    hClusterSplit1Pos1 = fs->make<TH2F>("hClusterSplit1Pos1"," ",416,0.,416.,160,0.,160.); 
    hClusterSplit1Pos1->SetOption("colz");
    hClusterSplit2Pos1 = fs->make<TH2F>("hClusterSplit2Pos1"," ",416,0.,416.,160,0.,160.); 
    hClusterSplit2Pos1->SetOption("colz");

    //hpixMapNoise = fs->make<TH2F>("hpixMapNoise"," ",416,0.,416.,160,0.,160.);
    //hpixMapNoise->SetOption("colz");

    //htest = fs->make<TH2F>("htest"," ",10,0.,10.,20,0.,20.);
    //htest2 = fs->make<TH2F>("htest2"," ",10,0.,10.,300,0.,300.);
    //htest->SetOption("colz");
    //htest2->SetOption("colz");

  hevent = fs->make<TH1F>("hevent","event",1000,0,10000000.);
  horbit = fs->make<TH1F>("horbit","orbit",100, 0,100000000.);

  hlumi1  = fs->make<TH1F>("hlumi1", "lumi", 2000,0,2000.);
  hlumi0  = fs->make<TH1F>("hlumi0", "lumi", 2000,0,2000.);

  hbx6    = fs->make<TH1F>("hbx6",   "bx",   4000,0,4000.);  
  hbx5    = fs->make<TH1F>("hbx5",   "bx",   4000,0,4000.);  
  hbx4    = fs->make<TH1F>("hbx4",   "bx",   4000,0,4000.);  
  hbx3    = fs->make<TH1F>("hbx3",   "bx",   4000,0,4000.);  
  hbx2    = fs->make<TH1F>("hbx2",   "bx",   4000,0,4000.);  
  hbx1    = fs->make<TH1F>("hbx1",   "bx",   4000,0,4000.);  
  hbx0    = fs->make<TH1F>("hbx0",   "bx",   4000,0,4000.);  
  hbx    = fs->make<TH1F>("hbx",   "bx",   4000,0,4000.);  
  
  hNdbCellsLay1  = fs->make<TH1F>("hNdbCellsLay1", "n DB cells per DC per BX Layer 1", 50, 0, 50.); 
  adcLay1Big  = fs->make<TH1F>("adcLay1Big", "ADC > 36 frames per DC, Layer 1", 256, 0, 256.0); 
  adcLay1Small  = fs->make<TH1F>("adcLay1Small", "ADC <= 36 frames per DC, Layer 1", 256, 0, 256.0); 

  hdets  = fs->make<TH1F>( "hdets",  "Dets with hits", 2000, -0.5, 1999.5);
  const int sizeH=20000;
  const float lowH = -0.5;
  const float highH = 99999.5;
  hdigis  = fs->make<TH1F>( "hdigis", "All Digis", sizeH, lowH, highH);
  hdigis0 = fs->make<TH1F>( "hdigis0", "All Digis zoomed", 2000, lowH, 1999.5);
  hdigis1 = fs->make<TH1F>( "hdigis1", "All Digis for full events", sizeH, lowH, highH);
  hdigis2 = fs->make<TH1F>( "hdigis2", "BPix Digis", sizeH, lowH, highH);
  hdigis3 = fs->make<TH1F>( "hdigis3", "Fpix Digis", sizeH, lowH, highH);
  hdigis4 = fs->make<TH1F>( "hdigis4", "All Digis - on bunch", sizeH, lowH, highH);
  hdigis5 = fs->make<TH1F>( "hdigis5", "All Digis - off bunch ", sizeH, lowH, highH);

  // 2D
  hzr = fs->make<TH2F>("hzr"," ",240,-60.,60.,68,0.,17.);  // x-y plane
  hxy11 = fs->make<TH2F>("hxy11"," ",320,-16.,16.,320,-16.,16.); // x-y pla 
  hxy12 = fs->make<TH2F>("hxy12"," ",320,-16.,16.,320,-16.,16.); // x-y pla
  hxy21 = fs->make<TH2F>("hxy21"," ",320,-16.,16.,320,-16.,16.); // x-y pl
  hxy22 = fs->make<TH2F>("hxy22"," ",320,-16.,16.,320,-16.,16.); // x-y pla
  hxy31 = fs->make<TH2F>("hxy31"," ",320,-16.,16.,320,-16.,16.); // x-y pla
  hxy32 = fs->make<TH2F>("hxy32"," ",320,-16.,16.,320,-16.,16.); // x-y plae

  hxy = fs->make<TH2F>("hxy"," ",340,-17.,17.,340,-17.,17.);  // x-y plane
  hphiz1 = fs->make<TH2F>("hphiz1"," ",108,-27.,27.,140,-3.5,3.5);
  hphiz2 = fs->make<TH2F>("hphiz2"," ",108,-27.,27.,140,-3.5,3.5);
  hphiz3 = fs->make<TH2F>("hphiz3"," ",108,-27.,27.,140,-3.5,3.5);
  hphiz4 = fs->make<TH2F>("hphiz4"," ",108,-27.,27.,140,-3.5,3.5);

  hadc1ls = fs->make<TProfile>("hadc1ls","adc1 vs ls",1000,0,1000,     0.,255.);
  hadc1bx = fs->make<TProfile>("hadc1bx","adc1 vs bx",4000,-0.5,3999.5,0.,255.);
  hadc2ls = fs->make<TProfile>("hadc2ls","adc2 vs ls",1000,0,1000,     0.,255.);
  hadc2bx = fs->make<TProfile>("hadc2bx","adc2 vs bx",4000,-0.5,3999.5,0.,255.);
  hadc3ls = fs->make<TProfile>("hadc3ls","adc3 vs ls",1000,0,1000,     0.,255.);
  hadc3bx = fs->make<TProfile>("hadc3bx","adc3 vs bx",4000,-0.5,3999.5,0.,255.);
  hadc4ls = fs->make<TProfile>("hadc4ls","adc4 vs ls",1000,0,1000,     0.,255.);
  hadc4bx = fs->make<TProfile>("hadc4bx","adc4 vs bx",4000,-0.5,3999.5,0.,255.);
  hadc0ls = fs->make<TProfile>("hadc0ls","adc0 vs ls",1000,0,1000,     0.,255.);
  hadc0bx = fs->make<TProfile>("hadc0bx","adc0 vs bx",4000,-0.5,3999.5,0.,255.);

#endif

    // count number of used buffer cells per BX per DC in layer 1 
    // initialize
    //for (int iDC=0;iDC<16*26;iDC++) {
    //    std::vector<int> emptyDC;
    //    doubleColumnHitsL1.push_back(emptyDC);
    //}
    //cout << "n DC:" << doubleColumnHitsL1.size() << std::endl;
}

// ------------ method called to produce the data  ------------
void PixDigisTest::analyze(const edm::Event& iEvent, 
			   const edm::EventSetup& iSetup) {

  const bool MY_DEBUG = false;
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopo;
  iSetup.get<TrackerTopologyRcd>().get(tTopo);
  const TrackerTopology* tt = tTopo.product();

  using namespace edm;
  if(PRINT) cout<<" Analyze PixDigisTest for phase "<<phase1_<<endl;

  //  int run       = iEvent.id().run();
  int event     = iEvent.id().event();
  int lumiBlock = iEvent.luminosityBlock();
  int bx        = iEvent.bunchCrossing();
  int orbit     = iEvent.orbitNumber();
  edm::Handle<TCDSRecord> tcdsData; 
  iEvent.getByToken( tcdsrecord_, tcdsData ); 
  if( tcdsData.isValid() ) {
      clusterTree_orbitOfLastBgo14 = tcdsData->getOrbitOfLastBgo(14);
      clusterTree_tres = (orbit - clusterTree_orbitOfLastBgo14)*3564 + (bx - 3000);// always sent in BX 3000 
  } else {
      clusterTree_orbitOfLastBgo14 = 0;
      clusterTree_tres = -1;
  }

  //cout << "event " << event << std::endl; 
  digiTestLS = lumiBlock;

  // fill tree
  clusterTree_lumiblock = lumiBlock;
  clusterTree_event     = event;
  clusterTree_bx        = bx;
  clusterTree_orbit     = orbit;

  hbx0->Fill(float(bx));
  hlumi0->Fill(float(lumiBlock));

  // // eliminate bunches with beam
  bool bunch = false;
  // if( bx==410 || bx==460 || bx==510) bunch = true;
  // else if( bx==560 || bx==610 || bx==660 ) bunch = true;
  // else if( bx==1292 || bx==1342 || bx==1392 ) bunch = true;
  // else if( bx==1454 || bx==1504 || bx==1554 ) bunch = true;
  // else if( bx==2501 || bx==2601 ) bunch = true;
  // else if( bx==3080 || bx==3030 || bx == 3180 ) bunch = true;

  // if(bx>=1 && bx<=351) { if( (bx%50) == 1 ) bunch = true; }
  // else if(bx>=892 && bx<=1245) {
  //   if( ((bx-892)%50) == 0 ) bunch = true;
  //   else if( ((bx-895)%50) == 0 ) bunch = true;
  // } else if(bx>=1786 && bx<=2286) { if( ((bx-1786)%50) == 0 ) bunch = true; }
  // else if(bx>=2671 && bx<=3021) { if( ((bx-2671)%50) == 0 ) bunch = true; }

  // if(bunch) {
  //   //cout<<" reject "<<bx<<endl;
  //   hbx2->Fill(float(bx));
  // } else {
  //   if(bx==892) cout<<" something wrong"<<endl;
  //   if(bx==1245) cout<<" something wrong"<<endl;
  //   if(bx==3021) cout<<" something wrong"<<endl;
  //   if(bx==2286) cout<<" something wrong"<<endl;
  // } 

    // Get digis
  edm::Handle< edm::DetSetVector<PixelDigi> > pixelDigis;
  iEvent.getByToken( tPixelDigi , pixelDigis);

#ifdef USE_SIM_LINKS
  // Get simlink data
  edm::Handle< edm::DetSetVector<PixelDigiSimLink> > pixelSimLinks;
  iEvent.getByToken( tPixelDigiSimLink,   pixelSimLinks);
#endif

  // Get event setup (to get global transformation)
  edm::ESHandle<TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);

  int numberOfDetUnits = 0;
  int totalNumOfDigis = 0;

  int numberOfDetUnits1 = 0;
  int totalNumOfDigis1 = 0;
  int numberOfDetUnits2 = 0;
  int totalNumOfDigis2 = 0;
  int numberOfDetUnits3 = 0;
  int totalNumOfDigis3 = 0;
  int numberOfDetUnits4 = 0;
  int totalNumOfDigis4 = 0;
  int numOfDigisPerDet1 = 0;
  int numOfDigisPerDet2 = 0;
  int numOfDigisPerDet3 = 0;
  int numOfDigisPerDet4 = 0;

  int numberOfDetUnitsF1 = 0;
  int totalNumOfDigisF1 = 0;
  int numberOfDetUnitsF2 = 0;
  int totalNumOfDigisF2 = 0;
  int numberOfDetUnitsF3 = 0;
  int totalNumOfDigisF3 = 0;
  int numOfDigisPerDetF1 = 0;
  int numOfDigisPerDetF2 = 0;
  int numOfDigisPerDetF3 = 0;

  count0++; // count events 

  // Iterate on detector units
  edm::DetSetVector<PixelDigi>::const_iterator DSViter;
  for(DSViter = pixelDigis->begin(); DSViter != pixelDigis->end(); DSViter++) {
    bool valid = true;
    unsigned int detid = DSViter->id; // = rawid
    DetId detId(detid);
    //const GeomDetUnit      * geoUnit = geom->idToDetUnit( detId );
    //const PixelGeomDetUnit * pixDet  = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
    unsigned int detType=detId.det(); // det type, tracker=1
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1
    
    if(MY_DEBUG) 
      cout<<"Det: "<<detId.rawId()<<" "<<detId.null()<<" "<<detType<<" "<<subid<<endl;
    
#ifdef HISTOS
    hdetunit->Fill(float(detid));
    hpixid->Fill(float(detType));
    hpixsubid->Fill(float(subid));
#endif
    
    if(detType!=1) continue; // look only at tracker
    ++numberOfDetUnits;
    count1++; // count det untis

    //for (int iDC=0;iDC<16*26;iDC++) {
    //    doubleColumnHitsL1.clear();
    //}
    
    // count number of used buffer cells per BX per DC in layer 1 
    // initialize
    std::vector< std::vector< std::pair< std::pair<int,int>, int> > > doubleColumnHitsL1;
    doubleColumnHitsL1.resize(16*26);
    for(std::vector < std::vector < std::pair< std::pair<int,int>, int>  > >::iterator it
                    =doubleColumnHitsL1.begin();it!=doubleColumnHitsL1.end();++it) {
       it->reserve(100);
    }

    // Get the geom-detector 
    const PixelGeomDetUnit * theGeomDet = 
      dynamic_cast<const PixelGeomDetUnit*> ( theTracker.idToDet(detId) );
    double detX = theGeomDet->surface().position().x();
    double detY = theGeomDet->surface().position().y();
    double detZ = theGeomDet->surface().position().z();
    double detR = theGeomDet->surface().position().perp();
    double detPhi = theGeomDet->surface().position().phi();
    //const BoundPlane plane = theGeomDet->surface(); // does not work
    
//     int cols = theGeomDet->specificTopology().ncolumns();
//     int rows = theGeomDet->specificTopology().nrows();
//     float pitchX = theGeomDet->specificTopology().pitch().first;
//     float pitchY = theGeomDet->specificTopology().pitch().second;
    
    const PixelTopology &topology = theGeomDet->specificTopology(); 
    int cols = topology.ncolumns();
    int rows = topology.nrows();
    float pitchX = topology.pitch().first;
    float pitchY = topology.pitch().second;
     
    unsigned int layerC=0;
    unsigned int ladderC=0;
    unsigned int zindex=0;

    // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
    int shell  = 0; // shell id 
    int sector = 0; // 1-8
    int ladder = 0; // 1-22
    int layer  = 0; // 1-3
    int module = 0; // 1-4
    bool half  = false; // 

    unsigned int disk=0; //1,2,3
    unsigned int blade=0; //1-24
    unsigned int side=0; //size=1 for -z, 2 for +z
    unsigned int panel=0; //panel=1
    
    std::string fullModuleName = "unititialized"; 
   // Subdet it, pix barrel=1, forward=2
    if(subid==2) {   // forward

#ifdef HISTOS
      hdetrF->Fill(detR);
      hdetzF->Fill(detZ);
      hcolsF->Fill(float(cols));
      hrowsF->Fill(float(rows));
#endif

      
      disk=tTopo->pxfDisk(detid); //1,2,3
      blade=tTopo->pxfBlade(detid); //1-24
      zindex=tTopo->pxfModule(detid); //
      side=tTopo->pxfSide(detid); //size=1 for -z, 2 for +z
      panel=tTopo->pxfPanel(detid); //panel=1
      
      if(PRINT) {
	cout<<"FPix (cmssw) disk "<<disk<<", blade "
		    <<blade<<", module "<<zindex<<", side "<<side<<", panel "
		    <<panel<<endl;
	if(MY_DEBUG) cout<<" col/row, pitch "<<cols<<" "<<rows<<" "
			 <<pitchX<<" "<<pitchY<<endl;
      }

      // correct way 
      PixelEndcapName pfn(detid,tt,phase1_);
      if(PRINT) 
	cout<<pfn.diskName()<<" "<<pfn.bladeName()<<" "<<pfn.pannelName()<<" "
	    <<pfn.ringName()<<" "<<pfn.plaquetteName()<<" "<<pfn.name()<<endl;

      // All these produce wrong names for phase1
      // PixelEndcapName pfn1(detid,false);
      // cout<<pfn1.diskName()<<" "<<pfn1.bladeName()<<" "<<pfn1.pannelName()<<" "<<pfn1.ringName()<<" "<<pfn1.plaquetteName()
      // 	  <<" "<<pfn1.name()<<endl;
      // PixelEndcapName pfn2(detid,true);
      // cout<<pfn2.diskName()<<" "<<pfn2.bladeName()<<" "<<pfn2.pannelName()<<" "<<pfn2.ringName()<<" "<<pfn2.plaquetteName()
      // 	  <<" "<<pfn2.name()<<endl;
      // PixelEndcapNameUpgrade pfn3(detid);
      // cout<<pfn3.diskName()<<" "<<pfn3.bladeName()<<" "<<pfn3.pannelName()<<" "<<pfn3.plaquetteName()
      // 	  <<" "<<pfn3.name()<<endl;

      fullModuleName = "forward";

    } else if(subid == 1) { // Barrel 
      
      // Barell layer = 1,2,3
      layerC=tTopo->pxbLayer(detid);
      // Barrel ladder id 1-20,32,44.
      ladderC=tTopo->pxbLadder(detid);
      // Barrel Z-index=1,8
      zindex=tTopo->pxbModule(detid);
      // Convert to online 
      PixelBarrelName pbn(detid,tt,phase1_);
      // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
      PixelBarrelName::Shell sh = pbn.shell(); //enum
      sector = pbn.sectorName();
      ladder = pbn.ladderName();
      layer  = pbn.layerName();
      module = pbn.moduleName();
      half  = pbn.isHalfModule();
      shell = int(sh);
      // change the module sign for z<0
      if(shell==1 || shell==2) module = -module;
      // change ladeer sign for Outer )x<0)
      if(shell==1 || shell==3) ladder = -ladder;
      
      fullModuleName = pbn.name();

      if(PRINT) { 
	cout<<" BPix layer/ladder/module (cmssw) "
	    <<layerC<<" "<<ladderC<<" "<<zindex<<" (online) "
	    <<pbn.name()<<" "<<sh<<"("<<shell<<") "<<sector<<" "
	    <<layer<<" "<<ladder<<" "<<module<<" "<<half<< endl;
	//cout<<" Barrel det, thick "<<detThick<<" "
	//  <<" layer, ladder, module "
	//  <<layer<<" "<<ladder<<" "<<zindex<<endl;
	//cout<<" col/row, pitch "<<cols<<" "<<rows<<" "
	//  <<pitchX<<" "<<pitchY<<endl;

	// wrong for phase1
	//PixelBarrelNameUpgrade pbn1(detid);
	//cout<<pbn1.sectorName()<<" "<<pbn1.ladderName()<<" "<<pbn1.layerName()<<" "<<pbn1.moduleName()<<pbn1.name();



      }      

    } // end fb-bar

#ifdef USE_SIM_LINKS
    // Look at simlink information (simulated data only)

    int numberOfSimLinks = 0;
    edm::DetSetVector<PixelDigiSimLink>::const_iterator 
      isearch = pixelSimLinks->find(detid);

    if(isearch != pixelSimLinks->end()) {      //if it is not empty
      edm::DetSet<PixelDigiSimLink> link_detset = (*pixelSimLinks)[detid];
      edm::DetSet<PixelDigiSimLink>::const_iterator it;
      // Loop over DigisSimLink in this det unit
      for(it = link_detset.data.begin(); 
	  it != link_detset.data.end(); it++) {
	
	numberOfSimLinks++;
	// these methods should be declared const, fixed by M.P.
	// wait for next releasse and then uncomment
	unsigned int chan = it->channel();
	unsigned int simTrack = it->SimTrackId();
	float frac = it->fraction();
#ifdef HISTOS
	hsimlinks->Fill(float(simTrack));
	hfract->Fill(frac);
#endif

	if(PRINT) cout<<" Sim link "<<numberOfSimLinks<<" "<<chan<<" "
		      <<simTrack<<" "<<frac<<endl;
      } // end simlink det loop

    } // end simlink if

#endif  // USE_SIM_LINKS

    unsigned int numberOfDigis = 0;
   
    strncpy(digiTestModuleName, fullModuleName.c_str(), 36);
    digiTestEvent = count0;
    digiTestNdigis = DSViter->data.size(); 

    // simple col/row vector for custom clustering
    std::vector < std::pair<std::pair<int, int>, int> > pixels;

    // fill tree variables
    clusterTree_shell    = shell;
    clusterTree_sector   = sector;
    clusterTree_ladder   = ladder;
    clusterTree_layer    = layer;
    clusterTree_module   = module;
    clusterTree_layerC   = layerC;
    clusterTree_ladderC  = ladderC;
    clusterTree_zindex   = zindex;

      // Look at digis now
    edm::DetSet<PixelDigi>::const_iterator  di;
    int digiTestCounter = 0;
   

      for(di = DSViter->data.begin(); di != DSViter->data.end(); di++) {
	//for(di = begin; di != end; di++) {
	
	numberOfDigis++;
	totalNumOfDigis++;
	count2++;
       int adc = di->adc();    // charge, modifued to unsiged short 
       int col = di->column(); // column 
       int row = di->row();    // row
       //int tof = di->time();    // tof always 0, method deleted

       digiTestDigiCol[digiTestCounter] = col;
       digiTestDigiRow[digiTestCounter] = row;
       digiTestDigiAdc[digiTestCounter] = adc;
       digiTestCounter++;
      
       // count ADC=0 pixels or not? 
       //if (adc > 0 || includeZeroAdc) {
           pixels.push_back(std::make_pair(std::make_pair(col, row), adc));
       //}

       // channel index needed to look for the simlink to simtracks
       int channel = PixelChannelIdentifier::pixelToChannel(row,col);
       if(PRINT) cout <<numberOfDigis<< " Col: " << col << " Row: " << row 
		      << " ADC: " << adc <<" channel = "<<channel<<endl;

       if(col>415) cout<<" Error: column index too large "<<col<<" Barrel layer, ladder, module "
		       <<layer<<" "<<ladder<<" "<<zindex<<endl;
       if(row>159) cout<<" Error: row index too large "<<row<<endl;

#ifdef HISTOS
       bool noise = false;
       if(layer==1) {
	 //noise = (ladder==6) || (module==-2) || (col==364) || (row==1);
	 if(!noise) {		     
	   heloss1->Fill(float(adc));
	   hcols1->Fill(float(col));
	   hrows1->Fill(float(row));
	   hpixMap1->Fill(float(col),float(row));
	   hpdetMap1->Fill(float(module),float(ladder));
	   if(adc<1) hpdetMaps1->Fill(float(module),float(ladder));
	   hadc1ls->Fill(float(lumiBlock),float(adc));
	   hadc1bx->Fill(float(bx),float(adc));

	   totalNumOfDigis1++;
	   //htest2->Fill(float(module),float(adc));
	   numOfDigisPerDet1++;
        
       // count number of used buffer cells per BX per DC in layer 1
       if (row < 80) {
         int dcIndex = int(col / 2);
         int rocRow = row;
         int rocCol = col % 52;
         if (dcIndex < 16*26) {
           doubleColumnHitsL1[dcIndex].push_back(std::make_pair(std::make_pair(rocCol, rocRow), adc));
         } else {
           cout << "invalid DC index:" << dcIndex << " col " << col << std::endl;
         }
       } else {
         int dcIndex = 8*26 + int(col / 2);
         int rocRow = 159 - row;
         int rocCol = col % 52;
         if (dcIndex < 16*26) {
           doubleColumnHitsL1[dcIndex].push_back(std::make_pair(std::make_pair(rocCol, rocRow), adc));
         } else {
           cout << "invalid DC index:" << dcIndex << " col " << col << std::endl;
         }
       }

//old 	   if(RectangularPixelTopology::isItBigPixelInX(row)) {
//new	   if(topology.isItBigPixelInX(row)) { 
// 	     //cout<<" big in X "<<row<<endl;
// 	     heloss1bigx->Fill(float(adc));
// 	     hrows1big->Fill(float(row));
// 	   }
//old	   if(RectangularPixelTopology::isItBigPixelInY(col)) {
//new	   if(topology.isItBigPixelInY(col)) {
// 	     //cout<<" big in Y "<<col<<endl;
// 	     heloss1bigy->Fill(float(adc));
// 	     hcols1big->Fill(float(col));
// 	   }

	 } // noise
       } else if(layer==2) {
	 // look for the noisy pixel
	 //noise = false; // (ladder==6) && (module==-2) && (col==364) && (row==1);
	 if(noise) {
	   //cout<<" noise pixel "<<layer<<" "<<sector<<" "<<shell<<endl;
	   //hpixMapNoise->Fill(float(col),float(row));
	   hneloss2->Fill(float(adc));
	 } else {		     
	   heloss2->Fill(float(adc));
	   hcols2->Fill(float(col));
	   hrows2->Fill(float(row));
	   hpixMap2->Fill(float(col),float(row));
	   hpdetMap2->Fill(float(module),float(ladder));
	   if(adc<1) hpdetMaps2->Fill(float(module),float(ladder));
	   hadc2ls->Fill(float(lumiBlock),float(adc));
	   hadc2bx->Fill(float(bx),float(adc));
	   totalNumOfDigis2++;
	   numOfDigisPerDet2++;
	 } // noise 
       } else if(layer==3) {
	 //noise = false; //(ladder==6) || (module==-2) || (col==364) || (row==1);
	 if(!noise) {		     
	   heloss3->Fill(float(adc));
	   hcols3->Fill(float(col));
	   hrows3->Fill(float(row));
	   //if(ladder==-13 && module==-4) 
	   hpdetMap3->Fill(float(module),float(ladder));
	   if(adc<1) hpdetMaps3->Fill(float(module),float(ladder));
	   hpixMap3->Fill(float(col),float(row));
	   hadc3ls->Fill(float(lumiBlock),float(adc));
	   hadc3bx->Fill(float(bx),float(adc));
	   totalNumOfDigis3++;
	   numOfDigisPerDet3++;
	 } // noise
       } else if(layer==4) {
	 //noise = false; //(ladder==6) || (module==-2) || (col==364) || (row==1);
	 if(!noise) {		     
	   heloss4->Fill(float(adc));
	   hcols4->Fill(float(col));
	   hrows4->Fill(float(row));
	   hpdetMap4->Fill(float(module),float(ladder));
	   if(adc<1) hpdetMaps4->Fill(float(module),float(ladder));
	   hpixMap4->Fill(float(col),float(row));
	   hadc4ls->Fill(float(lumiBlock),float(adc));
	   hadc4bx->Fill(float(bx),float(adc));
	   totalNumOfDigis4++;
	   numOfDigisPerDet4++;
	 } // noise
       } else if(disk==1) {
	 //noise = false; //(ladder==6) || (module==-2) || (col==364) || (row==1);
	 if(!noise) {		     
	   helossF1->Fill(float(adc));
	   hcolsF1->Fill(float(col));
	   hrowsF1->Fill(float(row));
	   hadc0ls->Fill(float(lumiBlock),float(adc));
	   hadc0bx->Fill(float(bx),float(adc));
	   totalNumOfDigisF1++;
	   numOfDigisPerDetF1++;
	 } // noise 

       } else if(disk==2) {
	 // noise = (ladder==6) || (module==-2) || (col==364) || (row==1);
	 if(!noise) {		     
	   helossF2->Fill(float(adc));
	   hcolsF2->Fill(float(col));
	   hrowsF2->Fill(float(row));
	   hadc0ls->Fill(float(lumiBlock),float(adc));
	   hadc0bx->Fill(float(bx),float(adc));
	   totalNumOfDigisF2++;
	   numOfDigisPerDetF2++;
	 } // noise 
       } else if(disk==3) {
	 // noise = (ladder==6) || (module==-2) || (col==364) || (row==1);
	 if(!noise) {		     
	   helossF3->Fill(float(adc));
	   hcolsF3->Fill(float(col));
	   hrowsF3->Fill(float(row));
	   hadc0ls->Fill(float(lumiBlock),float(adc));
	   hadc0bx->Fill(float(bx),float(adc));
	   totalNumOfDigisF3++;
	   numOfDigisPerDetF3++;
	 } // noise 
       } // end if layer

       if(noise) valid = false;

#endif
       
      } // end for digis in detunit

      if (restrictToLumisection < 1 || restrictToLumisection == lumiBlock) {
            //------------------------------------------------------------------------------------------------------------------------------
            // clustering
            //------------------------------------------------------------------------------------------------------------------------------
            //std::vector < std::pair<int, int> > pixels;
            std::vector < std::pair< std::pair<int, int>, int> > pixelsNew;
            std::vector < std::pair< std::pair<int, int>, int> > cluster;
            std::vector < std::vector < std::pair< std::pair<int, int>, int > > > clusters;
            // TODO: define type alias for above stuff...

            // start with all hits from the event

            int nClusters = 0;
            int clusterRadius = 4; //default: 2
            int newPixelsFound = 0;
            
            // as long as pixels left
            while (pixels.size() > 0 || newPixelsFound > 0) {
                // take last pixel to initialize the clusters
                if (cluster.size() < 1) {
                    cluster.push_back(pixels[pixels.size()-1]);
                    pixels.pop_back();
                }

                // loop over all pixels and check if one of them is close to already found cluster and merge pixel with this cluster
                newPixelsFound = 0;
                for (int i=pixels.size()-1;i>=0;i--) {
                    for (unsigned int j=0;j<cluster.size();j++) {
                        if (abs(pixels[i].first.first-cluster[j].first.first) < clusterRadius && abs(pixels[i].first.second-cluster[j].first.second) < clusterRadius) {
                            // add to cluster
                            cluster.push_back(pixels[i]);
                            // remove from list
                            pixels.erase(pixels.begin()+i);
                            newPixelsFound++;
                            break;
                        }
                    }
                }

                if (newPixelsFound < 1) {
                    // complete cluster found!
                    nClusters++;
                    sort(cluster.begin(), cluster.end() );
                    clusters.push_back(cluster);
                    cluster.clear();
                } else {
                    // remove duplicates from cluster
                    sort(cluster.begin(), cluster.end() );
                    cluster.erase(unique( cluster.begin(), cluster.end() ), cluster.end() );
                }
            }
            if (cluster.size() > 0) {
                sort(cluster.begin(), cluster.end() );
                clusters.push_back(cluster);
                cluster.clear();
            }

            //------------------------------------------------------------------------------------------------------------------------------
            // classify the found clusters by shape
            //------------------------------------------------------------------------------------------------------------------------------

            // loop over all clusters
            for (std::vector < std::vector <std::pair < std::pair<int, int>, int> > >::iterator itCluster = clusters.begin(); itCluster != clusters.end(); itCluster++) {

                // find "bottom left" corner (=min row, min col)
                int minRow = 600;
                int minCol = 600;
                for (std::vector <std::pair< std::pair<int, int>, int> >::iterator clusterPix=(*itCluster).begin(); clusterPix != (*itCluster).end(); clusterPix++) {
                    if ((*clusterPix).first.first < minCol) minCol = (*clusterPix).first.first;
                    if ((*clusterPix).first.second < minRow) minRow = (*clusterPix).first.second;
                }
                
                // align cluster with double columns
                minCol -= minCol % 2;

                // create relative cluster shape vector
                std::vector < std::pair< std::pair<int, int>, int> > relativeClusterShape;
                for (std::vector < std::pair < std::pair<int, int>, int> >::iterator clusterPix=(*itCluster).begin(); clusterPix != (*itCluster).end(); clusterPix++) {
                    relativeClusterShape.push_back(std::make_pair(std::make_pair((*clusterPix).first.first - minCol, (*clusterPix).first.second - minRow), (*clusterPix).second));
                }

                // sort it
                std::sort(relativeClusterShape.begin(), relativeClusterShape.end());

                // count clusters of this shape
                // disabled: only look at specific shapes, see below
                //if (clusterShapes.find(relativeClusterShape) != clusterShapes.end()) {
                //    clusterShapes[relativeClusterShape]++;
                //} else {
                //    clusterShapes[relativeClusterShape] = 1;
                //}

                // fill tree
                if (relativeClusterShape.size() < 480) {
                    clusterTree_size = relativeClusterShape.size();
                    for (int iRelPix=0;iRelPix<clusterTree_size;iRelPix++) {
                        clusterTree_x[iRelPix] = relativeClusterShape[iRelPix].first.first;
                        clusterTree_y[iRelPix] = relativeClusterShape[iRelPix].first.second;
                        clusterTree_adc[iRelPix] = relativeClusterShape[iRelPix].second;
                        std::cout << ">" << clusterTree_x[iRelPix] << " " << clusterTree_y[iRelPix] << ",";
                    }
                    clusterTree->Fill();
                }


                // check if the cluster is broken
                std::vector< std::pair< int,int >  > connectedPart;

                // check if the cluster is broken
                // clusters of specific size
                for (unsigned int checkClusterSize = 4; checkClusterSize<12; checkClusterSize++) {

                    // count all complete clusters of this size (includes longer clusters with outside pixels missing...)
                    if (relativeClusterShape.size() == checkClusterSize) {

                        bool specificClusterShapeFound = true;
                        // check all pixels
                        for (unsigned int p=0;p<relativeClusterShape.size();p++) {

                            // only use clusters aligned with ROC rows
                            if (relativeClusterShape[p].first.second != 0) {
                                specificClusterShapeFound = false;
                                break;
                            }

                            // all pixels there
                            if (relativeClusterShape[p].first.first != (int)p) {
                                specificClusterShapeFound = false;
                                break;
                            }
                        }

                        if (specificClusterShapeFound) {
                            // count clusters vs. module ID
                            std::pair<std::string, int> moduleClusterIdentifier = std::make_pair(fullModuleName, checkClusterSize);
                            if (clustersComplete.find(moduleClusterIdentifier) == clustersComplete.end()) {
                                clustersComplete[moduleClusterIdentifier] = 1;
                            } else {
                                clustersComplete[moduleClusterIdentifier]++;
                            }
                            // vs. BX only for layer 1!!!!
                            if (subid == 1 && layer == 1) {
                                std::pair<int, int> moduleClusterIdentifierBX = std::make_pair(bx, checkClusterSize);
                                if (clustersBXComplete.find(moduleClusterIdentifierBX) == clustersBXComplete.end()) {
                                    clustersBXComplete[moduleClusterIdentifierBX] = 1;
                                } else {
                                    clustersBXComplete[moduleClusterIdentifierBX]++;
                                }
                            }

                            break;                           
                        }
                    }

                    // count all broken clusters, which are shorter by 2 pixels
                    if (relativeClusterShape.size() == checkClusterSize - 2) {

                        // check all the n-2 different places where a cluster could break
                        for (unsigned int k=0;k<checkClusterSize-3;k++) {

                            // check all pixels
                            bool specificClusterShapeFound = true;
                            for (unsigned int p=0;p<relativeClusterShape.size();p++) {

                                // only use clusters aligned with ROC rows
                                if (relativeClusterShape[p].first.second != 0) {
                                    specificClusterShapeFound = false;
                                    break;
                                }

                                // in left part of cluster
                                if ((p < k + 1) && relativeClusterShape[p].first.first == (int)p) {
                                    // everything fine
                                // in right part of cluster
                                } else if ((p > k) && relativeClusterShape[p].first.first == (int)(p + 2)) {
                                    // everything fine
                                } else {
                                    specificClusterShapeFound = false;
                                    break;
                                }
                            }

                            if (specificClusterShapeFound) {
                                if (subid == 1 && layer == 1) {
                                    hClusterSplit2Pos1->Fill(minCol+k+1, minRow); 
                                    hClusterSplit2Pos1->Fill(minCol+k+2, minRow); 
                                    hClusterSplit2Row1->Fill(minRow); 
                                    hClusterSplit2Col1->Fill(minCol+k+1); 
                                    hClusterSplit2Col1->Fill(minCol+k+2); 
                                }
                                
                                std::pair<std::string, int> moduleClusterIdentifier = std::make_pair(fullModuleName, checkClusterSize);
                                if (clustersSplitByTwo.find(moduleClusterIdentifier) == clustersSplitByTwo.end()) {
                                    clustersSplitByTwo[moduleClusterIdentifier] = 1;
                                } else {
                                    clustersSplitByTwo[moduleClusterIdentifier]++;
                                }
                                
                                // vs. BX only for layer 1!!!!
                                if (layer ==  1 && subid == 1) {
                                    std::pair<int, int> moduleClusterIdentifierBX = std::make_pair(bx, checkClusterSize);
                                    if (clustersBXSplitByTwo.find(moduleClusterIdentifierBX) == clustersBXSplitByTwo.end()) {
                                        clustersBXSplitByTwo[moduleClusterIdentifierBX] = 1;
                                    } else {
                                        clustersBXSplitByTwo[moduleClusterIdentifierBX]++;
                                    }
                                }
                                //std::cout << "size " << checkClusterSize << " ---- found ---- k = " << k << "\n";
                                //for (int j =0;j<relativeClusterShape.size();j++) {
                                //    std::cout << relativeClusterShape[j].first << "," << relativeClusterShape[j].second << " ";
                                //}
                                //std::cout << "\n";
                                break;
                            }

                        }
                    }

                    // count all broken clusters, which are shorter by 1 pixels
                    if (relativeClusterShape.size() == checkClusterSize - 1) {

                        // check all the n-2 different places where a cluster could break
                        for (unsigned int k=0;k<checkClusterSize-2;k++) {

                            // check all pixels
                            bool specificClusterShapeFound = true;
                            for (unsigned int p=0;p<relativeClusterShape.size();p++) {

                                // only use clusters aligned with ROC rows
                                if (relativeClusterShape[p].first.second != 0) {
                                    specificClusterShapeFound = false;
                                    break;
                                }

                                // in left part of cluster
                                if ((p < k + 1) && relativeClusterShape[p].first.first == (int)p) {
                                    // everything fine
                                // in right part of cluster
                                } else if ((p > k) && relativeClusterShape[p].first.first == (int)(p + 1)) {
                                    // everything fine
                                } else {
                                    specificClusterShapeFound = false;
                                    break;
                                }
                            }

                            if (specificClusterShapeFound) {

                                if (subid == 1 && layer == 1) {
                                    hClusterSplit1Pos1->Fill(minCol+k+1, minRow); 
                                    hClusterSplit1Row1->Fill(minRow); 
                                    hClusterSplit1Col1->Fill(minCol+k+1); 
                                }

                                // vs. module ID
                                std::pair<std::string, int> moduleClusterIdentifier = std::make_pair(fullModuleName, checkClusterSize);
                                if (clustersSplitByOne.find(moduleClusterIdentifier) == clustersSplitByOne.end()) {
                                    clustersSplitByOne[moduleClusterIdentifier] = 1;
                                } else {
                                    clustersSplitByOne[moduleClusterIdentifier]++;
                                }
                                
                                // vs. BX only for layer 1!!!!
                                if (layer ==  1 && subid == 1) {
                                    std::pair<int, int> moduleClusterIdentifierBX = std::make_pair(bx, checkClusterSize);
                                    if (clustersBXSplitByOne.find(moduleClusterIdentifierBX) == clustersBXSplitByOne.end()) {
                                        clustersBXSplitByOne[moduleClusterIdentifierBX] = 1;
                                    } else {
                                        clustersBXSplitByOne[moduleClusterIdentifierBX]++;
                                    }
                                }
                                break;
                            }

                        }
                    }

                } // done looping over the specific cluster lengths to check

                // measure pixel rates
                if (moduleNhits.find(fullModuleName) == moduleNhits.end()) {
                    // add number of pixels in this cluster
                    moduleNhits[fullModuleName] = relativeClusterShape.size(); 
                } else {
                    moduleNhits[fullModuleName] += relativeClusterShape.size(); 
                }

                // measure cluster rates
                if (moduleNclusters.find(fullModuleName) == moduleNclusters.end()) {
                    moduleNclusters[fullModuleName] = 1;
                } else {
                    moduleNclusters[fullModuleName]++;
                }
            } // done looping over all clusters
       } 


       //if(PRINT) 
      //cout<<" for det "<<detid<<" digis = "<<numberOfDigis<<endl;

#ifdef HISTOS
      // Some histos
      if(valid) {  // to igore noise pixels
	if(subid==2) {   // forward
	  
	  hdetrF->Fill(detR);
	  hdetzF->Fill(detZ);
	  hcolsF->Fill(float(cols));
	  hrowsF->Fill(float(rows));
	  hzr->Fill(detZ,detR);
	  
	  if(disk==1) {
	    hblade1->Fill(float(blade));
	    ++numberOfDetUnitsF1;
	    hdigisPerDetF1->Fill(float(numOfDigisPerDetF1));
	    numOfDigisPerDetF1=0;	
	    if(detZ<0.) hxy11->Fill(detX,detY);
	    else        hxy12->Fill(detX,detY);
	  } else if(disk==2) {
	    hblade2->Fill(float(blade));
	    ++numberOfDetUnitsF2;
	    hdigisPerDetF2->Fill(float(numOfDigisPerDetF2));
	    numOfDigisPerDetF2=0;
	    if(detZ<0.) hxy21->Fill(detX,detY);
	    else        hxy22->Fill(detX,detY);

	  } else if(disk==3) {
	    hblade3->Fill(float(blade));
	    ++numberOfDetUnitsF3;
	    hdigisPerDetF3->Fill(float(numOfDigisPerDetF3));
	    numOfDigisPerDetF3=0;
	    if(detZ<0.) hxy31->Fill(detX,detY);
	    else        hxy32->Fill(detX,detY);
	  } // if disk

	} else if (subid==1) { // barrel
	  
	  hdetr->Fill(detR);
	  hdetz->Fill(detZ);
	  hcolsB->Fill(float(cols));
	  hrowsB->Fill(float(rows));
	  hxy->Fill(detX,detY);

	  hlayerid->Fill(float(layer));
	  hshellid->Fill(float(shell));
	  hsectorid->Fill(float(sector));

	  if(layer==1) {
	    
	    hladder1id->Fill(float(ladder));
	    hz1id->Fill(float(module));
	    hdetMap1->Fill(float(module),float(ladder));
	    ++numberOfDetUnits1;
	    hdigisPerDet1->Fill(float(numOfDigisPerDet1));
	    numOfDigisPerDet1=0;
	    hphiz1->Fill(detZ,detPhi);
	    
	  } else if(layer==2) {
	    hladder2id->Fill(float(ladder));
	    hz2id->Fill(float(module));
	    hdetMap2->Fill(float(module),float(ladder));
	    ++numberOfDetUnits2;
	    hdigisPerDet2->Fill(float(numOfDigisPerDet2));
	    numOfDigisPerDet2=0;
	    hphiz2->Fill(detZ,detPhi);
	    
	  } else if(layer==3) {
	    hladder3id->Fill(float(ladder));
	    hz3id->Fill(float(module));
	    hdetMap3->Fill(float(module),float(ladder));
	    ++numberOfDetUnits3;
	    hdigisPerDet3->Fill(float(numOfDigisPerDet3));
	    numOfDigisPerDet3=0;
	    hphiz3->Fill(detZ,detPhi);

	  } else if(layer==4) {
	    hladder4id->Fill(float(ladder));
	    hz4id->Fill(float(module));
	    hdetMap4->Fill(float(module),float(ladder));
	    ++numberOfDetUnits4;
	    hdigisPerDet4->Fill(float(numOfDigisPerDet4));
	    numOfDigisPerDet4=0;
	    hphiz4->Fill(detZ,detPhi);
	    
	  } // layer
	} // if bpix	
      } // if valid
#endif

   if (subid==1 && layer==1) {
     // count number of buffer cells used in layer 1 double columns
     for (int iDC=0;iDC<16*26;iDC++) {
       //std::vector< std::vector<int> > doubleColumnHitsL1;
       int nDbCells = 0;
       if (doubleColumnHitsL1[iDC].size() > 0) {
         std::sort(doubleColumnHitsL1[iDC].begin(), doubleColumnHitsL1[iDC].end()); 
         int lastRow = -2;
         for (size_t iHit=0;iHit<doubleColumnHitsL1[iDC].size();iHit++) {
            int hitRow = doubleColumnHitsL1[iDC][iHit].first.second;
            if (hitRow > lastRow + 1) {
              nDbCells++;
              lastRow = hitRow;
            }        
         }
       }
       if (nDbCells > 35) {
           cout << "BIG CLUSTER of size " << nDbCells << " 2x2 frames in BX " << bx << " " << fullModuleName << " DC " << iDC << "\n";

       }
       hNdbCellsLay1->Fill(nDbCells);
       if (nDbCells > 36) {
           for (size_t iHit=0;iHit<doubleColumnHitsL1[iDC].size();iHit++) {
               adcLay1Big->Fill(doubleColumnHitsL1[iDC][iHit].second);
           }
       } else {
           for (size_t iHit=0;iHit<doubleColumnHitsL1[iDC].size();iHit++) {
               adcLay1Small->Fill(doubleColumnHitsL1[iDC][iHit].second);
           }
       }

     }
   }

  } // end for det-units

  if(PRINT) 
    cout << " Number of full det-units = " <<numberOfDetUnits
	 <<" total digis = "<<totalNumOfDigis<<endl;
  hdets->Fill(float(numberOfDetUnits));
  hdigis->Fill(float(totalNumOfDigis));
  hdigis0->Fill(float(totalNumOfDigis));

  if(numberOfDetUnits>0) {
    count3++; // count full events 
    hevent->Fill(float(event));
    hlumi1->Fill(float(lumiBlock));
    hbx1->Fill(float(bx));
    horbit->Fill(float(orbit));

    hdigis1->Fill(float(totalNumOfDigis));
    float tmp = float(totalNumOfDigis1) + float(totalNumOfDigis2) 
      + float(totalNumOfDigis3); 
    hdigis2->Fill(tmp);
    tmp = float(totalNumOfDigisF1) + float(totalNumOfDigisF2);  
    hdigis3->Fill(tmp);

    if(bunch) hdigis4->Fill(float(totalNumOfDigis));
    else hdigis5->Fill(float(totalNumOfDigis));

    if(numberOfDetUnits>20) hbx3->Fill(float(bx));

    if(totalNumOfDigis>100) hbx4->Fill(float(bx));
    else if(totalNumOfDigis>4) hbx5->Fill(float(bx));
    else hbx6->Fill(float(bx));


  }

  
#ifdef HISTOS
  hdigisPerLay1 ->Fill(float(totalNumOfDigis1));
  hdigisPerLay2 ->Fill(float(totalNumOfDigis2));
  hdigisPerLay3 ->Fill(float(totalNumOfDigis3));
  hdigisPerLay4 ->Fill(float(totalNumOfDigis4));
  if(totalNumOfDigis1>0) hdetsPerLay1->Fill(float(numberOfDetUnits1));
  if(totalNumOfDigis2>0) hdetsPerLay2->Fill(float(numberOfDetUnits2));
  if(totalNumOfDigis3>0) hdetsPerLay3->Fill(float(numberOfDetUnits3));
  if(totalNumOfDigis4>0) hdetsPerLay4->Fill(float(numberOfDetUnits4));

  hdigisPerLayF1 ->Fill(float(totalNumOfDigisF1));
  hdigisPerLayF2 ->Fill(float(totalNumOfDigisF2));
  hdigisPerLayF3 ->Fill(float(totalNumOfDigisF3));
  hdetsPerLayF1 ->Fill(float(numberOfDetUnitsF1));
  hdetsPerLayF2 ->Fill(float(numberOfDetUnitsF2));
  hdetsPerLayF3 ->Fill(float(numberOfDetUnitsF3));
#endif

}
// ------------ method called to at the end of the job  ------------
void PixDigisTest::endJob(){
  cout << " End PixDigisTest " << endl;
  if(count0>0) 
    cout<<" events "<<count0<<" full events "<<count3<<" dets per event "<<count1/count0
	<<" digis per event "<<count2/count0<<endl;
  else 
    cout<<count0<<" "<<count1<<" "<<count2<<" "<<count3<<endl;

  float norm = 1.;
  if(count3>0) {
    norm = 1./float(count3);
  }

  // Rescale all 2D plots
  hdetMap1->Scale(norm);
  hdetMap2->Scale(norm);
  hdetMap3->Scale(norm);
  hdetMap4->Scale(norm);
  hpdetMap1->Scale(norm);
  hpdetMap2->Scale(norm);
  hpdetMap3->Scale(norm);
  hpdetMap4->Scale(norm);

  // print cluster shape dictionary and fill tree
  for ( std::map< std::pair<std::string, int>, int>::iterator clusterDict_it = clustersComplete.begin(); clusterDict_it!=clustersComplete.end(); ++clusterDict_it) {
      std::pair<std::string, int> moduleClusterIdentifier = clusterDict_it->first;
      int clustersCount = clusterDict_it->second;
      int clustersCount1 = 0;
      if (clustersSplitByOne.find(moduleClusterIdentifier) != clustersSplitByOne.end()) {
          clustersCount1 = clustersSplitByOne[moduleClusterIdentifier];
      }
      int clustersCount2 = 0;
      if (clustersSplitByTwo.find(moduleClusterIdentifier) != clustersSplitByTwo.end()) {
          clustersCount2 = clustersSplitByTwo[moduleClusterIdentifier];
      }
     // cout << "module " << moduleClusterIdentifier.first << " size " << moduleClusterIdentifier.second << " count " << clustersCount << " " << clustersCount1 << " " << clustersCount2 << "\n";

      // fill tree
      strncpy(clustersModuleName, moduleClusterIdentifier.first.c_str(), 36);
      clustersSizeBranch = moduleClusterIdentifier.second;
      clustersNcomplete = clustersCount;
      clustersNbroken1 = clustersCount1;
      clustersNbroken2 = clustersCount2;
      // OLD TREE: clusterTree->Fill(); 

  } 

  // print cluster shape dictionary vs. BX (for layer 1 only!!)
  for ( std::map< std::pair<int, int>, int>::iterator clusterDict_it = clustersBXComplete.begin(); clusterDict_it!=clustersBXComplete.end(); ++clusterDict_it) {
      std::pair<int, int> moduleClusterIdentifier = clusterDict_it->first;
      int clustersCount = clusterDict_it->second;
      int clustersCount1 = 0;
      if (clustersBXSplitByOne.find(moduleClusterIdentifier) != clustersBXSplitByOne.end()) {
          clustersCount1 = clustersBXSplitByOne[moduleClusterIdentifier];
      }
      int clustersCount2 = 0;
      if (clustersBXSplitByTwo.find(moduleClusterIdentifier) != clustersBXSplitByTwo.end()) {
          clustersCount2 = clustersBXSplitByTwo[moduleClusterIdentifier];
      }
      // cout << "module " << moduleClusterIdentifier.first << " size " << moduleClusterIdentifier.second << " count " << clustersCount << " " << clustersCount1 << " " << clustersCount2 << "\n";

      // fill tree
      clustersBxBranch = moduleClusterIdentifier.first;
      clustersSizeBranch = moduleClusterIdentifier.second;
      clustersNcomplete = clustersCount;
      clustersNbroken1 = clustersCount1;
      clustersNbroken2 = clustersCount2;

  } 
  for ( std::map< std::string, int>::iterator it = moduleNclusters.begin(); it!=moduleNclusters.end(); ++it) {
    double rate = it->second / (count0 * 2.5e-2 * 16 * 0.656);
    cout << "module " << it->first << ": #clusters = " << it->second << " -> " << rate << "\n";
  }


  for ( std::map< std::string, int>::iterator it = moduleNhits.begin(); it!=moduleNhits.end(); ++it) {
    double rate = it->second / (count0 * 2.5e-2 * 16 * 0.656);
    cout << "module " << it->first << ": #digis = " << it->second << " -> " << rate << "\n";
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(PixDigisTest);
