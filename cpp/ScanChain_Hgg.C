#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "TLorentzVector.h"
#include "math.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Base.h"
#include "../NanoCORE/tqdm.h"
#include "../NanoCORE/ZPrimeTools.cc"
#include "../NanoCORE/XYMETCorrection_withUL17andUL18andUL16.h"
#include "../NanoCORE/Tools/goodrun.h"
#include "../NanoCORE/Tools/dorky.h"
#include "../NanoCORE/Tools/puWeight.h"
#include "../NanoCORE/Tools/bTagEff.h"
#include "../NanoCORE/Tools/btagsf/BTagCalibrationStandalone_v2.h"
#include "../NanoCORE/DiPhotonSelections.h"
#include "../NanoCORE/LeptonSelections.h"
#include "../NanoCORE/DiJetSelections.h"
#include "../NanoCORE/GenPart.h"

#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <fstream>
#include <list>
#include <cmath>

#define H1(name,nbins,low,high,xtitle) TH1D *h_##name = new TH1D(#name,"",nbins,low,high); h_##name->GetXaxis()->SetTitle(xtitle); h_##name->GetYaxis()->SetTitle("Events");

#define Zmass 91.1876

// General flags
bool removeSpikes = true;
bool removeDataDuplicates = false;
bool usePuppiMET = true;

const char* outdir = "temp_data";
int mdir = mkdir(outdir,0755);

using namespace std;
using namespace tas;
using namespace duplicate_removal;
int count_test=0;

ofstream txtout("evtnb.txt", ofstream::app);

TCanvas *canvas = new TCanvas("c", "c", 800, 600);

//histograms for Garrett's cutflow 
TH1F *leadPtHist  = new TH1F("leadPtHist", "Leading Lepton pT", 10, 20, 200);
TH1F *leadEtaHist = new TH1F("leadEtaHist", "Leading Lepton eta", 8, -2.4, 2.4);
TH1F *secPtHist   = new TH1F("secPtHist", "Second Leading Lepton pT", 10, 20, 200);
TH1F *secEtaHist  = new TH1F("secEtaHist", "Second Leading Lepton eta", 8, -2.4, 2.4); 
TH1F *dilepMass   = new TH1F("dilepMass", "Dilepton Mass", 15, 50, 150); //GeV resolution, so use 1 GeV bin width
TH1F *dilepEta    = new TH1F("dilepEta", "Dilepton Eta", 8, -2.4, 2.4); //extend x-axis, labels, y-axis scale
TH1F *dilepP      = new TH1F("dilepP", "Dilepton Momentum", 15, 20, 300);
TH1F *jetUnclean  = new TH1F("jetUnclean", "Number of Jets in Event (No Cleaning)", 6, 0, 6);
TH1F *jetClean    = new TH1F("jetClean", "Number of Jets in Event (w/ Cleaning)", 6, 0, 6);
TH1F *jetCleanEE  = new TH1F("jetCleanEE", "Number of Jets in ee Event (w/ Cleaning)", 6, 0, 6);
TH1F *jetCleanMM  = new TH1F("jetCleanMM", "Number of Jets in mm Event (w/ Cleaning)", 6, 0, 6);
TH1F *dilepType     = new TH1F("nLeps", "DiLepType", 6, -0.5, 5.5);
TH1F *nElectrons  = new TH1F("nElectrons", "Number of Electrons", 6, -0.5, 5.5);
TH1F *nMuons      = new TH1F("nMuons", "Number of Muons", 6, -0.5, 5.5);

/*
//TFile to store histograms in
TFile *ZHists     = new TFile("ZHists.root", "RECREATE" );
*/

unsigned int totalElectrons    = 0;
unsigned int totalMuons        = 0;
unsigned int pTetaCutElectrons = 0;
unsigned int pTetaCutMuons     = 0;
unsigned int rightNumElectrons = 0;
unsigned int rightNumMuons     = 0;
unsigned int oppSignElectrons  = 0;
unsigned int oppSignMuons      = 0;

TLorentzVector p4dilepton;

int ScanChain_Hgg(TChain *ch, double genEventSumw, TString year, TString process, int process_id, const char* outdir="temp_data", int PUWeight=1, int bTagSF=1, int JECUnc=0) {
// Event weights / scale factors:
//  0: Do not apply
//  1: Apply central value
// +2: Apply positive variation
// -2: Apply negative variation

  float factor = 1.0;
  float lumi = 1.0;
  float xsec = 1.0;
  bool isMC = true;
  if ( process == "Data" ) {
    isMC = false;
  }
   
  // Processes and cross-sections (in fb):
  // set this in a different file
  else if ( process == "ttbar" )                              { xsec = 87310.0;                         }
  else if ( process == "DY" )                                 { xsec = 5765400.0;                       }
  else if ( process == "WW" )                                 { xsec = 118700.0;                        }
  else if ( process == "WZ" )                                 { xsec = 47130.0;                         }
  else if ( process == "ZZ" )                                 { xsec = 16523.0;                         }
  else if ( process == "tW" )                                 { xsec = 19550;                           }
  else if ( process == "tbarW" )                              { xsec = 19550;                           }
  else if ( process == "tZq" )                                { xsec = 75.8;                            }
  else if ( process == "TTW" )                                { xsec = 204.3;                           }
  else if ( process == "TTZ" )                                { xsec = 252.9;                           }
  else if ( process == "TTHToNonbb" )                         { xsec = 507.5*(1-0.575);                 }
  else if ( process == "TTHTobb" )                            { xsec = 507.5*0.575;                     }
  else if ( process == "TTGG" )                               { xsec = 0.01687 * 1000;                  }
  else if ( process == "TTGJets" )                            { xsec = 4.078 * 1000;                    }
  else if ( process == "TTJets" )                             { xsec = 831.76 * 1000;                   }
  else if ( process == "VBFH_M125" )                          { xsec = 0.00858514 *1000;                }
  else if ( process == "VH_M125" )                            { xsec = 0.00512 *1000;                   }
  else if ( process == "ggHToDiPhoM125" )                     { xsec = 0.1118429*1000 ;                 }
  else if ( process == "ttH_M125" )                           { xsec = 0.5071 * 1000 * 0.00227;         }
  else if ( process == "GJets_HT-40To100" )                   { xsec = 23100*1000;                      }
  else if ( process == "GJets_HT-100To200" )                  { xsec = 8631.0*1000;                     }
  else if ( process == "GJets_HT-200To400" )                  { xsec = 2280.0*1000;                     }
  else if ( process == "GJets_HT-400To600" )                  { xsec = 273*1000;                        }
  else if ( process == "GJets_HT-600ToInf" )                  { xsec = 1*1000;                          }
  else if ( process == "DiPhoton" )                           { xsec = 84.4*1000 ;                      }
  else if ( process == "DiPhotonLow" )                        { xsec = 303.2*1000 ;                     }
  else if ( process == "HHbbgg" )                             { xsec = 0.03105*1000*0.00262230;         }
  else if ( process == "WG" )                                 { xsec = 191.4*1000 ;                     }
  else if ( process == "ZG" )                                 { xsec = 55.6*1000 ;                      }
  else if ( process == "WW" )                                 { xsec = 75.8*1000 ;                      }
  else if ( process == "WZ" )                                 { xsec = 27.6*1000 ;                      }
  else if ( process == "ZZ" )                                 { xsec = 12.14*1000 ;                     }
  else if ( process == "Data" )                               { xsec = 1 ;                              }
  else if ( process.Contains("NMSSM_XToYHTo2G2B") )           { xsec = 1 ; process.ReplaceAll("-","_"); }
  else {
    cout<<"Non-valid process: Exiting!"<<endl;
    return 1;
  }

  // Configuration setup: NanoCORE/Config.{h,cc}
  gconf.nanoAOD_ver = 9;
  gconf.GetConfigs(year.Atoi());
  lumi = gconf.lumi;
  if (year == "2018") lumi = 59.8; //synchronizing with HiggsDNA

  // Golden JSON files
  if ( !isMC ) {
    if ( year == "2016APV" )
      set_goodrun_file_json("../utils/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
    if ( year == "2016nonAPV" )
      set_goodrun_file_json("../utils/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
    if ( year == "2017" )
      set_goodrun_file_json("../utils/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt");
    if ( year == "2018" )
      set_goodrun_file_json("../utils/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt");
  }

  if ( isMC )
    factor = xsec*lumi/genEventSumw;


  //all Modify the name of the output file to include arguments of ScanChain function (i.e. process, year, etc.)
  int mdir = mkdir(outdir,0755);
  TString oDir(outdir);
  TFile* fout = new TFile(oDir+"/output_"+process+"_"+year+".root", "RECREATE");
  TTree* tout = new TTree("tout","Tree with photon variables");


  // define histograms, to be put in a different file TODO
  H1(LeadPhoton_sieie, 20, 0, 0.05, "");
  H1(LeadPhoton_pfPhoIso03, 20, 0, 10, "");
  H1(LeadPhoton_chargedHadronIso, 20, 0, 10, "");
  //H1(LeadPhoton_trkSumPtHollowConeDR03, 20, 0, 10, ""); // To be readded
  H1(SubleadPhoton_sieie, 20, 0, 0.05, "");
  H1(SubleadPhoton_pfPhoIso03, 20, 0, 10, "");
  H1(SubleadPhoton_chargedHadronIso, 20, 0, 10, "");
  //H1(SubleadPhoton_trkSumPtHollowConeDR03, 20, 0, 10, ""); // To be readded


  // Variables for output branches
  //float xcand_pt, xcand_eta, xcand_phi, xcand_mass;

  //float LeadPhoton_pt, LeadPhoton_eta, LeadPhoton_phi, LeadPhoton_mvaID; // LeadPhoton_mass 
  //bool LeadPhoton_pixelSeed, SubleadPhoton_pixelSeed;
  //float SubleadPhoton_pt, SubleadPhoton_eta, SubleadPhoton_phi, SubleadPhoton_mvaID; // SubleadPhoton_mass
  //float Diphoton_pt, Diphoton_eta, Diphoton_phi, Diphoton_mass, Diphoton_pt_mgg, Diphoton_dR;
  //float LeadPhoton_sieie, LeadPhoton_pfPhoIso03, LeadPhoton_chargedHadronIso, LeadPhoton_r9; // LeadPhoton_trkSumPtHollowConeDR03 
  //float SubleadPhoton_sieie, SubleadPhoton_pfPhoIso03, SubleadPhoton_chargedHadronIso, SubleadPhoton_r9; // SubleadPhoton_trkSumPtHollowConeDR03

  int n_electrons, n_muons, electron_1_id, electron_2_id, electron_3_id, muon_1_id, muon_2_id, muon_3_id;

  float electron_1_pt, electron_1_eta, electron_1_phi, electron_1_mass;
  float electron_2_pt, electron_2_eta, electron_2_phi, electron_2_mass;
  float electron_3_pt, electron_3_eta, electron_3_phi, electron_3_mass;

  float muon_1_pt, muon_1_eta, muon_1_phi, muon_1_mass;
  float muon_2_pt, muon_2_eta, muon_2_phi, muon_2_mass;
  float muon_3_pt, muon_3_eta, muon_3_phi, muon_3_mass;

  int n_jets, jet_1_jetId, jet_2_jetId, jet_3_jetId;
  float jet_1_pt, jet_1_eta, jet_1_phi, jet_1_mass, jet_1_btagDeepFlavB;
  float jet_2_pt, jet_2_eta, jet_2_phi, jet_2_mass, jet_2_btagDeepFlavB;
  float jet_3_pt, jet_3_eta, jet_3_phi, jet_3_mass, jet_3_btagDeepFlavB;
  //float dijet_lead_pt, dijet_lead_eta, dijet_lead_phi, dijet_lead_mass, dijet_lead_btagDeepFlavB;
  //float dijet_sublead_pt, dijet_sublead_eta, dijet_sublead_phi, dijet_sublead_mass, dijet_sublead_btagDeepFlavB;
  //float dijet_pt, dijet_eta, dijet_phi, dijet_mass, dijet_dR;
  float pfmet_pt, puppimet_pt;
  int year_out, eventNum;
  float weight_central, weight_central_initial, weight_central_no_lumi;

  //unsigned int LeadPhoton_genPartFlav = 0, SubleadPhoton_genPartFlav = 0;
  //int n_gen_matched_jets = 0, n_gen_matched_in_dijet = 0;
  //bool dijet_lead_gen_match=false, dijet_sublead_gen_match=false;
  //float GenHiggs_pt=-999, GenHiggs_eta=-999, GenHiggs_phi=-999, GenHiggs_mass=-999, GenHiggs_dR=-999;
  //float GenY_pt=-999, GenY_eta=-999, GenY_phi=-999, GenY_mass=-999, GenY_dR=-999;
  //float GenX_pt=-999, GenX_eta=-999, GenX_phi=-999, GenX_mass=-999, GenX_dR=-999;
  //float GenBFromHiggs_1_pt=-999, GenBFromHiggs_1_eta=-999, GenBFromHiggs_1_phi=-999, GenBFromHiggs_1_mass=-999;
  //float GenBFromHiggs_2_pt=-999, GenBFromHiggs_2_eta=-999, GenBFromHiggs_2_phi=-999, GenBFromHiggs_2_mass=-999;


  // Branch booking
  //tout->Branch("xcand_pt", &xcand_pt, "xcand_pt/F");
  //tout->Branch("xcand_eta", &xcand_eta, "xcand_eta/F");
  //tout->Branch("xcand_phi", &xcand_phi, "xcand_phi/F");
  //tout->Branch("xcand_mass", &xcand_mass, "xcand_mass/F");

  //tout->Branch("LeadPhoton_pt",&LeadPhoton_pt,"LeadPhoton_pt/F");
  //tout->Branch("LeadPhoton_eta",&LeadPhoton_eta,"LeadPhoton_eta/F");
  //tout->Branch("LeadPhoton_phi",&LeadPhoton_phi,"LeadPhoton_phi/F");
  ////tout->Branch("LeadPhoton_mass",&LeadPhoton_mass,"LeadPhoton_mass/F");
  //tout->Branch("LeadPhoton_pixelSeed",&LeadPhoton_pixelSeed,"LeadPhoton_pixelSeed/B");
  //tout->Branch("LeadPhoton_r9",&LeadPhoton_r9,"LeadPhoton_r9/F");
  //tout->Branch("LeadPhoton_sieie",&LeadPhoton_sieie,"LeadPhoton_sieie/F");
  //tout->Branch("LeadPhoton_pfPhoIso03",&LeadPhoton_pfPhoIso03,"LeadPhoton_pfPhoIso03/F");
  ////tout->Branch("LeadPhoton_trkSumPtHollowConeDR03",&LeadPhoton_trkSumPtHollowConeDR03,"LeadPhoton_trkSumPtHollowConeDR03/F"); // To be readded
  //tout->Branch("LeadPhoton_chargedHadronIso",&LeadPhoton_chargedHadronIso,"LeadPhoton_chargedHadronIso/F");
  //tout->Branch("LeadPhoton_mvaID",&LeadPhoton_mvaID,"LeadPhoton_mvaID/F");

  //tout->Branch("SubleadPhoton_pt",&SubleadPhoton_pt,"SubleadPhoton_pt/F");
  //tout->Branch("SubleadPhoton_eta",&SubleadPhoton_eta,"SubleadPhoton_eta/F");
  //tout->Branch("SubleadPhoton_phi",&SubleadPhoton_phi,"SubleadPhoton_phi/F");
  ////tout->Branch("SubleadPhoton_mass",&SubleadPhoton_mass,"SubleadPhoton_mass/F");
  //tout->Branch("SubleadPhoton_pixelSeed",&SubleadPhoton_pixelSeed,"SubleadPhoton_pixelSeed/B");  
  //tout->Branch("SubleadPhoton_r9",&SubleadPhoton_r9,"SubleadPhoton_r9/F");
  //tout->Branch("SubleadPhoton_sieie",&SubleadPhoton_sieie,"SubleadPhoton_sieie/F");
  //tout->Branch("SubleadPhoton_pfPhoIso03",&SubleadPhoton_pfPhoIso03,"SubleadPhoton_pfPhoIso03/F");
  ////tout->Branch("SubleadPhoton_trkSumPtHollowConeDR03",&SubleadPhoton_trkSumPtHollowConeDR03,"SubleadPhoton_trkSumPtHollowConeDR03/F"); // To be readded
  //tout->Branch("SubleadPhoton_chargedHadronIso",&SubleadPhoton_chargedHadronIso,"SubleadPhoton_chargedHadronIso/F");
  //tout->Branch("SubleadPhoton_mvaID",&SubleadPhoton_mvaID,"SubleadPhoton_mvaID/F");

  //tout->Branch("Diphoton_pt",&Diphoton_pt,"Diphoton_pt/F");
  //tout->Branch("Diphoton_eta",&Diphoton_eta,"Diphoton_eta/F");
  //tout->Branch("Diphoton_phi",&Diphoton_phi,"Diphoton_phi/F");
  //tout->Branch("Diphoton_mass",&Diphoton_mass,"Diphoton_mass/F");
  //tout->Branch("Diphoton_pt_mgg",&Diphoton_pt_mgg,"Diphoton_pt_mgg/F");
  //tout->Branch("Diphoton_dR",&Diphoton_dR,"Diphoton_dR/F");

  // Electron vars
  tout->Branch("n_electrons",&n_electrons,"n_electrons/I");  

  tout->Branch("electron_1_id",&electron_1_id,"electron_1_id/I");  
  tout->Branch("electron_2_id",&electron_2_id,"electron_2_id/I");  
  tout->Branch("electron_3_id",&electron_3_id,"electron_3_id/I");  

  tout->Branch("electron_1_pt",&electron_1_pt,"electron_1_pt/F");
  tout->Branch("electron_1_eta",&electron_1_eta,"electron_1_eta/F");
  tout->Branch("electron_1_phi",&electron_1_phi,"electron_1_phi/F");
  tout->Branch("electron_1_mass",&electron_1_mass,"electron_1_mass/F");
  
  tout->Branch("electron_2_pt",&electron_2_pt,"electron_2_pt/F");
  tout->Branch("electron_2_eta",&electron_2_eta,"electron_2_eta/F");
  tout->Branch("electron_2_phi",&electron_2_phi,"electron_2_phi/F");
  tout->Branch("electron_2_mass",&electron_2_mass,"electron_2_mass/F");
 
  tout->Branch("electron_3_pt",&electron_3_pt,"electron_3_pt/F");
  tout->Branch("electron_3_eta",&electron_3_eta,"electron_3_eta/F");
  tout->Branch("electron_3_phi",&electron_3_phi,"electron_3_phi/F");
  tout->Branch("electron_3_mass",&electron_3_mass,"electron_3_mass/F");

  // Muon vars
  tout->Branch("n_muons",&n_muons,"n_muons/I");  

  tout->Branch("muon_1_id",&muon_1_id,"muon_1_id/I");  
  tout->Branch("muon_2_id",&muon_2_id,"muon_2_id/I");  
  tout->Branch("muon_3_id",&muon_3_id,"muon_3_id/I");  

  tout->Branch("muon_1_pt",&muon_1_pt,"muon_1_pt/F");
  tout->Branch("muon_1_eta",&muon_1_eta,"muon_1_eta/F");
  tout->Branch("muon_1_phi",&muon_1_phi,"muon_1_phi/F");
  tout->Branch("muon_1_mass",&muon_1_mass,"muon_1_mass/F");
  
  tout->Branch("muon_2_pt",&muon_2_pt,"muon_2_pt/F");
  tout->Branch("muon_2_eta",&muon_2_eta,"muon_2_eta/F");
  tout->Branch("muon_2_phi",&muon_2_phi,"muon_2_phi/F");
  tout->Branch("muon_2_mass",&muon_2_mass,"muon_2_mass/F");
 
  tout->Branch("muon_3_pt",&muon_3_pt,"muon_3_pt/F");
  tout->Branch("muon_3_eta",&muon_3_eta,"muon_3_eta/F");
  tout->Branch("muon_3_phi",&muon_3_phi,"muon_3_phi/F");
  tout->Branch("muon_3_mass",&muon_3_mass,"muon_3_mass/F");

  // Jet vars
  tout->Branch("n_jets",&n_jets,"n_jets/I");  

  tout->Branch("jet_1_jetId",&jet_1_jetId,"jet_1_jetId/I");  
  tout->Branch("jet_2_jetId",&jet_2_jetId,"jet_2_jetId/I");  
  tout->Branch("jet_3_jetId",&jet_3_jetId,"jet_3_jetId/I");  

  tout->Branch("jet_1_pt",&jet_1_pt,"jet_1_pt/F");
  tout->Branch("jet_1_eta",&jet_1_eta,"jet_1_eta/F");
  tout->Branch("jet_1_phi",&jet_1_phi,"jet_1_phi/F");
  tout->Branch("jet_1_mass",&jet_1_mass,"jet_1_mass/F");
  tout->Branch("jet_1_btagDeepFlavB",&jet_1_btagDeepFlavB,"jet_1_btagDeepFlavB/F");
  
  tout->Branch("jet_2_pt",&jet_2_pt,"jet_2_pt/F");
  tout->Branch("jet_2_eta",&jet_2_eta,"jet_2_eta/F");
  tout->Branch("jet_2_phi",&jet_2_phi,"jet_2_phi/F");
  tout->Branch("jet_2_mass",&jet_2_mass,"jet_2_mass/F");
  tout->Branch("jet_2_btagDeepFlavB",&jet_2_btagDeepFlavB,"jet_2_btagDeepFlavB/F");
 
  tout->Branch("jet_3_pt",&jet_3_pt,"jet_3_pt/F");
  tout->Branch("jet_3_eta",&jet_3_eta,"jet_3_eta/F");
  tout->Branch("jet_3_phi",&jet_3_phi,"jet_3_phi/F");
  tout->Branch("jet_3_mass",&jet_3_mass,"jet_3_mass/F");
  tout->Branch("jet_3_btagDeepFlavB",&jet_3_btagDeepFlavB,"jet_3_btagDeepFlavB/F");

  //tout->Branch("dijet_lead_pt",&dijet_lead_pt,"dijet_lead_pt/F");
  //tout->Branch("dijet_lead_eta",&dijet_lead_eta,"dijet_lead_eta/F");
  //tout->Branch("dijet_lead_phi",&dijet_lead_phi,"dijet_lead_phi/F");
  //tout->Branch("dijet_lead_mass",&dijet_lead_mass,"dijet_lead_mass/F");
  //tout->Branch("dijet_lead_btagDeepFlavB",&dijet_lead_btagDeepFlavB,"dijet_lead_btagDeepFlavB/F");

  //tout->Branch("dijet_sublead_pt",&dijet_sublead_pt,"dijet_sublead_pt/F");
  //tout->Branch("dijet_sublead_eta",&dijet_sublead_eta,"dijet_sublead_eta/F");
  //tout->Branch("dijet_sublead_phi",&dijet_sublead_phi,"dijet_sublead_phi/F");
  //tout->Branch("dijet_sublead_mass",&dijet_sublead_mass,"dijet_sublead_mass/F");
  //tout->Branch("dijet_sublead_btagDeepFlavB",&dijet_sublead_btagDeepFlavB,"dijet_sublead_btagDeepFlavB/F");

  //tout->Branch("dijet_pt",&dijet_pt,"dijet_pt/F");
  //tout->Branch("dijet_eta",&dijet_eta,"dijet_eta/F");  
  //tout->Branch("dijet_phi",&dijet_phi,"dijet_phi/F");  
  //tout->Branch("dijet_mass",&dijet_mass,"dijet_mass/F");  
  //tout->Branch("dijet_dR",&dijet_dR,"dijet_dR/F"); 

  tout->Branch("pfmet_pt",&pfmet_pt,"pfmet_pt/F"); 
  tout->Branch("puppimet_pt",&puppimet_pt,"puppimet_pt/F"); 

  tout->Branch("year",&year_out,"year/I");
  tout->Branch("weight_central",&weight_central,"weight_central/F");
  tout->Branch("weight_central_initial",&weight_central_initial,"weight_central_initial/F");
  tout->Branch("weight_central_no_lumi",&weight_central_no_lumi,"weight_central_no_lumi/F");
  tout->Branch("event",&eventNum,"event/I");
  tout->Branch("process_id",&process_id,"process_id/I");

  if (year=="2016nonAPV" || year=="2016APV") year_out = 2016;
  else if (year=="2017") year_out = 2017;
  else if (year=="2018") year_out = 2018;
  else year_out = 0;
  
  //if (isMC) {
  //  tout->Branch("LeadPhoton_genPartFlav",&LeadPhoton_genPartFlav,"LeadPhoton_genPartFlav/C");
  //  tout->Branch("SubleadPhoton_genPartFlav",&SubleadPhoton_genPartFlav,"SubleadPhoton_genPartFlav/C");
  //  tout->Branch("n_gen_matched_jets",&n_gen_matched_jets,"n_gen_matched_jets/I");
  //  tout->Branch("n_gen_matched_in_dijet",&n_gen_matched_in_dijet,"n_gen_matched_in_dijet/I");
  //  tout->Branch("dijet_lead_gen_match",&dijet_lead_gen_match,"dijet_lead_gen_match/B");
  //  tout->Branch("dijet_sublead_gen_match",&dijet_sublead_gen_match,"dijet_sublead_gen_match/B");
  //  tout->Branch("GenHiggs_pt",&GenHiggs_pt,"GenHiggs_pt/F");
  //  tout->Branch("GenHiggs_eta",&GenHiggs_eta,"GenHiggs_eta/F");
  //  tout->Branch("GenHiggs_phi",&GenHiggs_phi,"GenHiggs_phi/F");
  //  tout->Branch("GenHiggs_mass",&GenHiggs_mass,"GenHiggs_mass/F");
  //  tout->Branch("GenHiggs_dR",&GenHiggs_dR,"GenHiggs_dR/F");
  //  tout->Branch("GenY_pt",&GenY_pt,"GenY_pt/F");
  //  tout->Branch("GenY_eta",&GenY_eta,"GenY_eta/F");
  //  tout->Branch("GenY_phi",&GenY_phi,"GenY_phi/F");
  //  tout->Branch("GenY_mass",&GenY_mass,"GenY_mass/F");
  //  tout->Branch("GenY_dR",&GenY_dR,"GenY_dR/F");
  //  tout->Branch("GenX_pt",&GenX_pt,"GenX_pt/F");
  //  tout->Branch("GenX_eta",&GenX_eta,"GenX_eta/F");
  //  tout->Branch("GenX_phi",&GenX_phi,"GenX_phi/F");
  //  tout->Branch("GenX_mass",&GenX_mass,"GenX_mass/F");
  //  tout->Branch("GenX_dR",&GenX_dR,"GenX_dR/F");
  //  tout->Branch("GenBFromHiggs_1_pt",&GenBFromHiggs_1_pt,"GenBFromHiggs_1_pt/F");
  //  tout->Branch("GenBFromHiggs_1_eta",&GenBFromHiggs_1_eta,"GenBFromHiggs_1_eta/F");
  //  tout->Branch("GenBFromHiggs_1_phi",&GenBFromHiggs_1_phi,"GenBFromHiggs_1_phi/F");
  //  tout->Branch("GenBFromHiggs_1_mass",&GenBFromHiggs_1_mass,"GenBFromHiggs_1_mass/F");
  //  tout->Branch("GenBFromHiggs_2_pt",&GenBFromHiggs_2_pt,"GenBFromHiggs_2_pt/F");
  //  tout->Branch("GenBFromHiggs_2_eta",&GenBFromHiggs_2_eta,"GenBFromHiggs_2_eta/F");
  //  tout->Branch("GenBFromHiggs_2_phi",&GenBFromHiggs_2_phi,"GenBFromHiggs_2_phi/F");
  //  tout->Branch("GenBFromHiggs_2_mass",&GenBFromHiggs_2_mass,"GenBFromHiggs_2_mass/F");
  //}


  // Define histo info maps
  map<TString, int> nbins { };
  map<TString, float> low { };
  map<TString, float> high { };
  map<TString, vector<float>> binsx { };
  map<TString, TString> title { };

  // Define histos
  H1(cutflow,20,0,20,"");
  H1(weight,1,0,1,"");
  H1(weight_full,1,0,1,"");


  //if ( PUWeight!=0 ) set_puWeights(); // FIXME to be enabled later?

  int nEventsTotal = 0;
  int nDuplicates = 0;
  int nEventsChain = ch->GetEntries();
  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);
  tqdm bar;

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TFile *file = TFile::Open( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    TString filename(currentFile->GetTitle());

    tree->SetCacheSize(128*1024*1024);
    tree->SetCacheLearnEntries(100);

    nt.Init(tree);
  
    
	for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {
      nt.GetEntry(event);
      tree->LoadTree(event);

      nEventsTotal++;
      bar.progress(nEventsTotal, nEventsChain);

      float weight = 1.0;
      if ( isMC ) {
        weight = nt.genWeight();
        //if(removeSpikes && weight*factor>1e2) continue; //comment out for synchronizing

        // Apply PU reweight // FIXME to be enabled later?
/*        
        if ( PUWeight!=0 ) {
          unsigned int nTrueInt = nt.Pileup_nTrueInt();
          TString whichPUWeight = "central";
          if ( PUWeight==2 ) whichPUWeight = "up";
          else if ( PUWeight==-2 ) whichPUWeight = "down";
          TString puyear = year;
          if ( useOnlyRun2018B ) puyear = "2018B";
          weight *= get_puWeight(nTrueInt, puyear, whichPUWeight);
          if (event <  100) cout << "genWeight after PU reWeight is: " << weight << endl;
        }
*/        
      }

      h_weight_full->Fill(0.5, weight*factor);

      unsigned int runnb = nt.run();
      unsigned int lumiblock = nt.luminosityBlock();
      unsigned long int evtnb = nt.event();
      int npv = nt.PV_npvs();


      // Apply Golden JSON
      if ( !isMC ) {
        if ( !(goodrun(runnb, lumiblock)) )
          continue;
        if ( removeDataDuplicates ) {
          DorkyEventIdentifier id(runnb, evtnb, lumiblock);
          if ( is_duplicate(id) ) {
            ++nDuplicates;
            continue;
          }
        }
      }


      // HLT selection
      // HLT selection not applied in MC for synchronizing (same as HiggsDNA) - FIXME: To be applied later
      if (!isMC){
        if ( (year=="2016nonAPV" || year=="2016APV") &&
            !( (tree->GetBranch("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55") ? nt.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() : 0)
              || (tree->GetBranch("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55") ? nt.HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() : 0) ) ) continue;
        if ( (year=="2017") &&
            !( (tree->GetBranch("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55") ? nt.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55() : 0)  )  ) continue;
        if ( (year=="2018") &&
            !( (tree->GetBranch("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto") ? nt.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto() : 0) ) ) continue;
      }


      // Object selection
      //Photons photons = getPhotons();
      //if (isMC) {
      //  for (auto pho : photons)
      //    pho.setGenPartFlav(pho.idx());
      //}
      //DiPhotons diphotons = DiPhotonPreselection(photons);

      //if (diphotons.size() == 0 ) continue; 

      //DiPhoton selectedDiPhoton = diphotons[0];
      //Photons selectedPhotons={selectedDiPhoton.leadPho, selectedDiPhoton.subleadPho};

      //LeadPhoton_mvaID = selectedDiPhoton.leadPho.mvaID();
      //SubleadPhoton_mvaID = selectedDiPhoton.subleadPho.mvaID();
      //if ( !(SubleadPhoton_mvaID>-0.7 && LeadPhoton_mvaID>-0.7) ) continue;

      Electrons electrons = getElectrons();
      Muons muons = getMuons();
      nElectrons -> Fill(electrons.size());
	  nMuons -> Fill(muons.size());
	  std::cout << "num els = " << electrons.size() << "num mus = " << muons.size() << std::endl;
	  
	  /*
	  
	  BEGIN CUTLFOW:
	  We scan over the leptons in each event. We first throw out leptons whose kinematics do not fall within our desired bounds. Second, we
	  check that the event had exactly two same-flavor, opposite-sign leptons, and throw out events that were not.
	  
	  */
	  
	  Electrons candElectrons;
	  totalElectrons += electrons.size();	//counts number of electrons in all events, no cuts applied yet
	  
	  //apply kinematic cuts to electrons
	  for (unsigned int n = 0; n < electrons.size(); n++){ 
		Electron e = electrons[n];
		if ( !(abs(e.eta()) < 2.4) ) continue;
		if ( (e.pt() < 20)) continue;
		candElectrons.push_back(e); //add electron to candElectrons if it passes all cuts
		pTetaCutElectrons++; //counts the number of electrons that pass kinematic cuts
		}
	  
	  //apply kinematic cuts to muons
	  Muons candMuons;
	  totalMuons += muons.size();
	  for (unsigned int n = 0; n < muons.size(); n++){ 
		Muon mu = muons[n];
		if ( !(abs(mu.eta()) < 2.4) ) continue;
		if ( (mu.pt() < 20)) continue;
		candMuons.push_back(mu);
		pTetaCutMuons++;
	    }
	  
	  // Only allow events with 2 leptons 
      if ((candElectrons.size() + candMuons.size()) != 2) continue; 
	  //count the number of leptons (wrt flavor) in events that passed the previous cuts and now the 2 lepton requirement 
	  rightNumElectrons += candElectrons.size();
	  rightNumMuons += candMuons.size();
	  
	  //dilep type: ee = 0, emu = 1, mumu = 2
	  dilepType -> Fill(candMuons.size());
	  
	  //Only allow opposite signed leptons
	  //First check if the 2 leptons are opposite sign electrons
	  if (candElectrons.size() == 2 && candElectrons[0].id()*candElectrons[1].id() < 0){ 
		  leadPtHist  -> Fill(candElectrons[0].pt());	//fill histograms for leading lepton and second-leading leptons pt and eta 
		  leadEtaHist -> Fill(candElectrons[0].eta());
		  secPtHist   -> Fill(candElectrons[1].pt());
		  secEtaHist  -> Fill(candElectrons[1].eta());
		  p4dilepton = candElectrons[0].p4() + candElectrons[1].p4();	//calculate dilepton 4-vector
		  oppSignElectrons += 2;	//count number of electrons that passed previous cuts and now the opposite sign cut
	  }
	  //Next check if the 2 leptons are opposite sign muons
	  else if (candMuons.size() == 2 && candMuons[0].id()*candMuons[1].id() < 0){
		  leadPtHist  -> Fill(candMuons[0].pt());
		  leadEtaHist -> Fill(candMuons[0].eta());
		  secPtHist   -> Fill(candMuons[1].pt());
		  secEtaHist  -> Fill(candMuons[1].eta());
		  p4dilepton = candMuons[0].p4() + candMuons[1].p4();
		  oppSignMuons += 2;
	  }
	  //if neither of the condions above are met, skip event (the leptons are not opposite-sign, same-flavor)
	  else continue; 
	  
	  
	  //fill dilepton histograms
	  dilepMass -> Fill(p4dilepton.Mag());
	  dilepEta  -> Fill(p4dilepton.Eta());
	  dilepP    -> Fill(p4dilepton.P());
       
	  Jets jets = getJets();
	  
	  Jets nocleanJets;		//jets before cleaning
	  Jets cleanJets;		//jets after cleaning
	  Jets cleanJetsEE;	
	  Jets cleanJetsMM;
	  
	  //loop over all jets in event
	  for (unsigned int n = 0; n < jets.size(); n++){ 
		Jet j = jets[n];
		//kinematic cuts for our jets
		if ( abs(j.eta()) > 2.4 ) continue;
		if ( (j.pt() < 30)) continue;
		nocleanJets.push_back(j); //if jet passes our kinematic cuts, push_back to nocleanjets
		
		//For each jet, loop over the leptons in the events. If the jet is within a delta R (dR) of 0.4 with any lepton, set JetIsLep = True
		bool JetIsLep = false;
		for (unsigned int k = 0; k < candElectrons.size(); k++){
			double dR = sqrt(pow(acos(cos(j.phi() - candElectrons[k].phi())), 2) + pow(j.eta() - candElectrons[k].eta(), 2)); //calculate dR
			if (dR < 0.4){ 
				JetIsLep = true;	//if the "jet" is within 0.4 of a lepton, it is likely that the "jet" is actually that lepton being doubly reconstructed
				goto cnt;			//as a both a lepton and a jet
			}
		}
	    
		for (unsigned int z = 0; z < candMuons.size(); z++){
			double dR = sqrt(pow(acos(cos(j.phi() - candMuons[z].phi())), 2) + pow(j.eta() - candMuons[z].eta(), 2));
			if (dR < 0.4){ 
				JetIsLep = true;
				goto cnt;
		    }
	    }
	    
		//if the jet is outside dR <= 0.4 for every lepton, push back to cleanJets
		if (JetIsLep == false) {
			cleanJets.push_back(j);
		}
		
		cnt:;
	  }
	  //finally, fill nJet histograms 
      jetUnclean -> Fill(nocleanJets.size());
	  jetClean -> Fill(cleanJets.size());
	  if (candElectrons.size() == 2) jetCleanEE -> Fill(cleanJets.size());
	  if (candMuons.size() == 2) jetCleanMM -> Fill(cleanJets.size());
	  //if (jets.size() < 2) continue; 

      //DiJets dijets = DiJetPreselection(jets);
      //DiJet selectedDiJet = dijets[0];

      //if (dijets[0].p4.M()<50) continue;


      // Setting output variables
      //TLorentzVector x_cand = selectedDiPhoton.p4 + selectedDiJet.p4;
      //xcand_pt = x_cand.Pt();
      //xcand_eta = x_cand.Eta();
      //xcand_phi = x_cand.Phi();
      //xcand_mass = x_cand.M();

      //LeadPhoton_pt = selectedDiPhoton.leadPho.pt();
      //LeadPhoton_eta = selectedDiPhoton.leadPho.eta();
      //LeadPhoton_phi = selectedDiPhoton.leadPho.phi();
      ////LeadPhoton_mass = selectedDiPhoton.leadPho.mass();
      //LeadPhoton_r9 = selectedDiPhoton.leadPho.r9();
      //LeadPhoton_pixelSeed = selectedDiPhoton.leadPho.pixelSeed();
      //LeadPhoton_sieie = selectedDiPhoton.leadPho.sieie();
      //LeadPhoton_pfPhoIso03 = selectedDiPhoton.leadPho.phoIso();
      ////LeadPhoton_trkSumPtHollowConeDR03 = selectedDiPhoton.leadPho.trkIso(); // To be readded
      //LeadPhoton_chargedHadronIso = selectedDiPhoton.leadPho.chargedHadIso();

      //SubleadPhoton_pt = selectedDiPhoton.subleadPho.pt();
      //SubleadPhoton_eta = selectedDiPhoton.subleadPho.eta();
      //SubleadPhoton_phi = selectedDiPhoton.subleadPho.phi();
      ////SubleadPhoton_mass = selectedDiPhoton.subleadPho.mass();
      //SubleadPhoton_r9 = selectedDiPhoton.subleadPho.r9();
      //SubleadPhoton_pixelSeed = selectedDiPhoton.subleadPho.pixelSeed();
      //SubleadPhoton_sieie = selectedDiPhoton.subleadPho.sieie();
      //SubleadPhoton_pfPhoIso03 = selectedDiPhoton.subleadPho.phoIso();
      ////SubleadPhoton_trkSumPtHollowConeDR03 = selectedDiPhoton.subleadPho.trkIso(); // To be readded
      //SubleadPhoton_chargedHadronIso = selectedDiPhoton.subleadPho.chargedHadIso();

      //Diphoton_pt = selectedDiPhoton.p4.Pt();
      //Diphoton_eta = selectedDiPhoton.p4.Eta();
      //Diphoton_phi = selectedDiPhoton.p4.Phi();
      //Diphoton_mass = selectedDiPhoton.p4.M();
      //Diphoton_pt_mgg = Diphoton_pt/Diphoton_mass;
      //Diphoton_dR = selectedDiPhoton.dR;
	  
      n_electrons = electrons.size();
      n_muons = muons.size();

      int want_electrons = 3;
      if (n_electrons < want_electrons) { 
        Electron dummy_electron = Electron(999);
        for (int i = 0; i < (want_electrons - n_electrons); i++) {
          electrons.push_back(dummy_electron);
        }
      }
      electron_1_pt = electrons[0].p4().Pt();
      electron_1_eta = electrons[0].p4().Eta();
      electron_1_phi = electrons[0].p4().Phi();
      electron_1_mass = electrons[0].p4().M();
      electron_1_id = electrons[0].id();

      electron_2_pt = electrons[1].p4().Pt();
      electron_2_eta = electrons[1].p4().Eta();
      electron_2_phi = electrons[1].p4().Phi();
      electron_2_mass = electrons[1].p4().M();
      electron_2_id = electrons[1].id();

      electron_3_pt = electrons[2].p4().Pt();
      electron_3_eta = electrons[2].p4().Eta();
      electron_3_phi = electrons[2].p4().Phi();
      electron_3_mass = electrons[2].p4().M();
      electron_3_id = electrons[2].id();

      int want_muons = 3;
      if (n_muons < want_muons) { 
        Muon dummy_muon = Muon(999);
        for (int i = 0; i < (want_muons - n_muons); i++) {
          muons.push_back(dummy_muon);
        }
      }
      muon_1_pt = muons[0].p4().Pt();
      muon_1_eta = muons[0].p4().Eta();
      muon_1_phi = muons[0].p4().Phi();
      muon_1_mass = muons[0].p4().M();
      muon_1_id = muons[0].id();

      muon_2_pt = muons[1].p4().Pt();
      muon_2_eta = muons[1].p4().Eta();
      muon_2_phi = muons[1].p4().Phi();
      muon_2_mass = muons[1].p4().M();
      muon_2_id = muons[1].id();

      muon_3_pt = muons[2].p4().Pt();
      muon_3_eta = muons[2].p4().Eta();
      muon_3_phi = muons[2].p4().Phi();
      muon_3_mass = muons[2].p4().M();
      muon_3_id = muons[2].id();

      n_jets = jets.size();
      int want_jets = 3;
      if (n_jets < want_jets) { 
        Jet dummy_jet = Jet(999);
        for (int i = 0; i < (want_jets - n_jets); i++) {
          jets.push_back(dummy_jet);
        }
      }
      jet_1_pt = jets[0].p4().Pt();
      jet_1_eta = jets[0].p4().Eta();
      jet_1_phi = jets[0].p4().Phi();
      jet_1_mass = jets[0].p4().M();
      jet_1_jetId = jets[0].jetId();
      jet_1_btagDeepFlavB = jets[0].btagDeepFlavB();

      jet_2_pt = jets[1].p4().Pt();
      jet_2_eta = jets[1].p4().Eta();
      jet_2_phi = jets[1].p4().Phi();
      jet_2_mass = jets[1].p4().M();
      jet_2_jetId = jets[1].jetId();
      jet_2_btagDeepFlavB = jets[1].btagDeepFlavB();

      jet_3_pt = jets[2].p4().Pt();
      jet_3_eta = jets[2].p4().Eta();
      jet_3_phi = jets[2].p4().Phi();
      jet_3_mass = jets[2].p4().M();
      jet_3_jetId = jets[2].jetId();
      jet_3_btagDeepFlavB = jets[2].btagDeepFlavB();

      //dijet_lead_pt = selectedDiJet.leadJet.pt();
      //dijet_lead_eta = selectedDiJet.leadJet.eta();
      //dijet_lead_phi = selectedDiJet.leadJet.phi();
      //dijet_lead_mass = selectedDiJet.leadJet.mass();
      //dijet_lead_btagDeepFlavB = selectedDiJet.leadJet.btagDeepFlavB();
      //dijet_sublead_pt = selectedDiJet.subleadJet.pt();
      //dijet_sublead_eta = selectedDiJet.subleadJet.eta();
      //dijet_sublead_phi = selectedDiJet.subleadJet.phi();
      //dijet_sublead_mass = selectedDiJet.subleadJet.mass();
      //dijet_sublead_btagDeepFlavB = selectedDiJet.subleadJet.btagDeepFlavB();
      //dijet_pt = selectedDiJet.p4.Pt();
      //dijet_eta = selectedDiJet.p4.Eta();
      //dijet_phi = selectedDiJet.p4.Phi();
      //dijet_mass = selectedDiJet.p4.M();
      //dijet_dR = selectedDiJet.dR;
      pfmet_pt = isMC ? nt.MET_T1_pt() : nt.MET_pt();
      puppimet_pt = nt.PuppiMET_pt();

      weight_central = weight*factor;
      weight_central_initial = weight;
      weight_central_no_lumi = weight*factor/lumi;
      eventNum = evtnb;

      count_test++;

      // Gen info - FIXME: Doesn't work because all if-statements are trivially false
      //if (isMC) {
      //  LeadPhoton_genPartFlav = selectedDiPhoton.leadPho.genPartFlav();
      //  SubleadPhoton_genPartFlav = selectedDiPhoton.subleadPho.genPartFlav();

      //  TLorentzVector GenHiggs, GenX, GenY;
      //  GenParts genparts = getGenParts();
      //  vector<TLorentzVector> gen_child_xyh, gen_child_ygg, gen_child_hbb;
      //  for (int igenpart=0; igenpart<genparts.size(); igenpart++) {
      //    if (genparts[igenpart].isxyh()) {GenX = genparts[igenpart].mother_p4(); gen_child_xyh.push_back(genparts[igenpart].p4());}
      //    if (genparts[igenpart].isygg()) {GenY = genparts[igenpart].mother_p4(); gen_child_ygg.push_back(genparts[igenpart].p4());}
      //    if (genparts[igenpart].ishbb()) {GenHiggs = genparts[igenpart].mother_p4(); gen_child_hbb.push_back(genparts[igenpart].p4());}
      //  }
      //  TLorentzVector sort_GenPart;
      //  //if (gen_child_xyh[0].Pt()<gen_child_xyh[1].Pt()) {sort_GenPart = gen_child_xyh[0]; gen_child_xyh[0]= gen_child_xyh[1]; gen_child_xyh[1] = sort_GenPart;}
      //  //if (gen_child_ygg[0].Pt()<gen_child_ygg[1].Pt()) {sort_GenPart = gen_child_ygg[0]; gen_child_ygg[0]= gen_child_ygg[1]; gen_child_ygg[1] = sort_GenPart;}
      //  if (!abs(GenX_pt+999)<0.0001) {
      //    GenX_pt = GenX.Pt();
      //    GenX_eta = GenX.Eta();
      //    GenX_phi = GenX.Phi();
      //    GenX_mass = GenX.M();
      //    GenX_dR = gen_child_xyh[0].DeltaR(gen_child_xyh[1]);
      //  }
      //  if (!abs(GenY_pt+999)<0.0001) {
      //    GenY_pt = GenY.Pt();
      //    GenY_eta = GenY.Eta();
      //    GenY_phi = GenY.Phi();
      //    GenY_mass = GenY.M();
      //    GenY_dR = gen_child_ygg[0].DeltaR(gen_child_ygg[1]);
      //  }
      //  if (!abs(GenHiggs_pt+999)<0.0001) {
      //    GenHiggs_pt = GenHiggs.Pt();
      //    GenHiggs_eta = GenHiggs.Eta();
      //    GenHiggs_phi = GenHiggs.Phi();
      //    GenHiggs_mass = GenHiggs.M();
      //    GenHiggs_dR = gen_child_hbb[0].DeltaR(gen_child_hbb[1]);
      //    if (gen_child_hbb[0].Pt()<gen_child_hbb[1].Pt()) {sort_GenPart = gen_child_hbb[0]; gen_child_hbb[0]= gen_child_hbb[1]; gen_child_hbb[1] = sort_GenPart;}
      //    GenBFromHiggs_1_pt = gen_child_hbb[0].Pt();
      //    GenBFromHiggs_1_eta = gen_child_hbb[0].Eta();
      //    GenBFromHiggs_1_phi = gen_child_hbb[0].Phi();
      //    GenBFromHiggs_1_mass = gen_child_hbb[0].M();
      //    GenBFromHiggs_2_pt = gen_child_hbb[1].Pt();
      //    GenBFromHiggs_2_eta = gen_child_hbb[1].Eta();
      //    GenBFromHiggs_2_phi = gen_child_hbb[1].Phi();
      //    GenBFromHiggs_2_mass = gen_child_hbb[1].M();
      //    if (selectedDiJet.leadJet.p4().DeltaR(gen_child_hbb[0])<=0.4 || selectedDiJet.leadJet.p4().DeltaR(gen_child_hbb[1])<=0.4) {dijet_lead_gen_match=true; n_gen_matched_in_dijet++;}
      //    if (selectedDiJet.subleadJet.p4().DeltaR(gen_child_hbb[0])<=0.4 || selectedDiJet.subleadJet.p4().DeltaR(gen_child_hbb[1])<=0.4) {dijet_sublead_gen_match=true; n_gen_matched_in_dijet++;}
      //    for (int ijet=0; ijet<jets.size(); ijet++) {
      //      if (jets[ijet].p4().DeltaR(gen_child_hbb[0])<=0.4 || jets[ijet].p4().DeltaR(gen_child_hbb[1])<=0.4)
      //        n_gen_matched_jets++;
      //    }
      //  }
      //}


      // Histo filling
      //h_LeadPhoton_sieie->Fill(LeadPhoton_sieie);
      //h_LeadPhoton_pfPhoIso03->Fill(LeadPhoton_pfPhoIso03);
      //h_LeadPhoton_chargedHadronIso->Fill(LeadPhoton_chargedHadronIso);
      ////h_LeadPhoton_trkSumPtHollowConeDR03->Fill(LeadPhoton_trkSumPtHollowConeDR03); // To be readded
      //h_SubleadPhoton_sieie->Fill(SubleadPhoton_sieie);
      //h_SubleadPhoton_pfPhoIso03->Fill(SubleadPhoton_pfPhoIso03);
      //h_SubleadPhoton_chargedHadronIso->Fill(SubleadPhoton_chargedHadronIso);
      //h_SubleadPhoton_trkSumPtHollowConeDR03->Fill(SubleadPhoton_trkSumPtHollowConeDR03); // To be readded
      h_weight->Fill(0.5, weight*factor);

      tout->Fill();
    } // Event loop
    
	delete file;
  } // File loop
  
  //Save histograms
  canvas -> Draw();
  leadPtHist -> Draw();
  canvas -> SaveAs("leadPtHist.png");
  canvas -> Clear();
  leadEtaHist -> Draw();
  canvas -> SaveAs("leadEtaHist.png");
  canvas -> Clear();
  secPtHist -> Draw();
  canvas -> SaveAs("secPtHist.png");
  canvas -> Clear();
  secEtaHist -> Draw();
  canvas -> SaveAs("secEtaHist.png");
  canvas -> Clear();
  dilepMass -> Draw();
  canvas -> SaveAs("dilepMass.png");
  canvas -> Clear();
  dilepEta -> Draw();
  canvas -> SaveAs("dilepEta.png");
  canvas -> Clear();
  dilepP -> Draw();
  canvas -> SaveAs("dilepP.png");
  canvas -> Clear();
  jetUnclean -> Draw();
  canvas -> SaveAs("jetUnclean.png");
  canvas -> Clear();
  jetClean -> Draw();
  canvas -> SaveAs("jetClean.png");
  canvas -> Clear();
  jetCleanEE -> Draw();
  canvas -> SaveAs("jetCleanEE.png");
  canvas -> Clear();
  jetCleanMM -> Draw();
  canvas -> SaveAs("jetCleanMM.png");
  
  /*
  ZHists -> WriteObject(&leadPtHist, "leadPtHist");
  ZHists -> WriteObject(&leadEtaHist, "leadEtaHist");
  ZHists -> WriteObject(&secPtHist, "secPtHist");
  ZHists -> WriteObject(&secEtaHist, "secEtaHist");
  ZHists -> WriteObject(&dilepMass, "dilepMass");
  ZHists -> WriteObject(&dilepEta, "dilepEta");
  ZHists -> WriteObject(&dilepP, "dilepP");
  ZHists -> WriteObject(&jetUnclean, "jetUnclean");
  ZHists -> WriteObject(&jetClean, "jetClean");
  ZHists -> WriteObject(&jetCleanEE, "jetCleanEE");
  ZHists -> WriteObject(&jetCleanMM, "jetCleanMM");
  */
  
  
  /*
  // Print number of leptons that pass our selections in total
  //commented out for now because it was causing problems with the progress bar.
  std::cout << "Total electrons: " << ::totalElectrons << " Total Muons: " << ::totalMuons << "\n";
  std::cout << "Kinematic cut \n";
  std::cout << "Electrons: " << pTetaCutElectrons << " Muons: " << pTetaCutMuons << "\n";
  std::cout << "2 leptons cut \n";
  std::cout << "Electrons: " << rightNumElectrons << " Muons: " << rightNumMuons << "\n";
  std::cout << "Opposite sign cut \n";
  std::cout << "Electrons: " << oppSignElectrons << " Muons: " << oppSignMuons << "\n";
  */
  
  bar.finish();
  cout << "nTotal: " << h_weight_full->GetBinContent(1) << ", nPass: " << h_weight->GetBinContent(1) << ", eff: " << h_weight->GetBinContent(1)/h_weight_full->GetBinContent(1) << endl;
  cout << endl;

  if ( removeDataDuplicates )
    cout << "Number of duplicates found: " << nDuplicates << endl;

  txtout.close();
  fout->Write();
  fout->Close();

  return 0;
}
