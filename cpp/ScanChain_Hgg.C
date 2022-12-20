#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
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
#include "../NanoCORE/Tools/muonRecoSF.h"
#include "../NanoCORE/Tools/muonIDSF.h"
#include "../NanoCORE/Tools/muonIsoSF.h"
#include "../NanoCORE/Tools/muonTriggerSF.h"
#include "../NanoCORE/Tools/bTagEff.h"
#include "../NanoCORE/Tools/btagsf/BTagCalibrationStandalone_v2.h"
#include "../NanoCORE/Tools/jetcorr/JetCorrectionUncertainty.h"
#include "../NanoCORE/DiPhotonSelections.h"
#include "../NanoCORE/LeptonSelections.h"
#include "../NanoCORE/DiJetSelections.h"

#include "configuration_Zp.h"

#include <iostream>
#include <iomanip>
#include <sys/stat.h>

#define SUM(vec) std::accumulate((vec).begin(), (vec).end(), 0);
#define SUM_GT(vec,num) std::accumulate((vec).begin(), (vec).end(), 0, [](float x,float y) { return ((y > (num)) ? x+y : x); });
#define COUNT_GT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x > (num); });
#define COUNT_LT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x < (num); });

#define H1(name,nbins,low,high,xtitle) TH1F *h_##name = new TH1F(#name,"",nbins,low,high); h_##name->GetXaxis()->SetTitle(xtitle); h_##name->GetYaxis()->SetTitle("Events");
#define HTemp(name,nbins,low,high,xtitle) TH1F *h_temp = new TH1F(name,"",nbins,low,high); h_temp->GetXaxis()->SetTitle(xtitle); h_temp->GetYaxis()->SetTitle("Events");

// #define DEBUG
#define Zmass 91.1876

// For testing purposes only
bool useOnlyRun2018B = true;

// Looper setup flags
bool muonDebug = false;
bool doMllBins = false;
bool doNbTagBins = true;
bool doTTEnriched = false;
bool doDYEnriched = false;
bool doMuDetRegionBins = false;

// General flags
bool removeSpikes = true;
bool removeDataDuplicates = false;
bool useTuneP = true;
bool usePuppiMET = true;
bool fillRooDataSet = true;
//

const char* outdir = "temp_data";
int mdir = mkdir(outdir,0755);

struct debugger { template<typename T> debugger& operator , (const T& v) { cerr<<v<<" "; return *this; } } dbg;
#ifdef DEBUG
#define debug(args...) do {cerr << #args << ": "; dbg,args; cerr << endl;} while(0)
#else
#define debug(args...)
#endif

using namespace std;
using namespace tas;
using namespace duplicate_removal;
using namespace RooFit;
int count_test=0;

int ScanChain_Hgg(TChain *ch, double genEventSumw, TString year, TString process, int topPtWeight=1, int PUWeight=1, int muonSF=1, int triggerSF=1, int bTagSF=1, int JECUnc=0) {
// Event weights / scale factors:
//  0: Do not apply
//  1: Apply central value
// +2: Apply positive variation
// -2: Apply negative variation

  float factor = 1.0;
  float lumi = 1.0;
  float xsec = 1.0;
  bool isMC = true;

  cout << "Process " << process << endl;

  if ( process.Contains("data") ) {
    isMC = false;
  }
  // SM processes and cross-sections:
  // set this in a different file
  else if ( process == "ttbar" )             xsec = 87310.0; // fb
  else if ( process == "DY" )                xsec = 5765400.0; // fb
  else if ( process == "WW" )                xsec = 118700.0; // fb 
  else if ( process == "WZ" )                xsec = 47130.0; // fb
  else if ( process == "ZZ" )                xsec = 16523.0; // fb
  else if ( process == "tW" )                xsec = 19550; // fb
  else if ( process == "tbarW" )             xsec = 19550; // fb
  else if ( process == "tZq" )               xsec = 75.8; // fb
  else if ( process == "TTW" )               xsec = 204.3; // fb
  else if ( process == "TTZ" )               xsec = 252.9; // fb
  else if ( process == "TTHToNonbb" )        xsec = 507.5*(1-0.575); // fb
  else if ( process == "TTHTobb" )           xsec = 507.5*0.575; // fb
  else if ( process == "ggHToDiPhoM125" )    xsec = 33.93 ; // fb
                                                               // TODO change xs
  else if ( process == "diPhoton" )    xsec = 82.51 ; // 
  else if ( process == "HHggtautau" )    xsec = 0.02669 ; // 
//  111.8429 ; // fb
  // Signal processes and cross-sections:
  else
    {
      cout<<"Non-valid process: Exiting!"<<endl;
      return 1;
    }

  // Configuration setup: NanoCORE/Config.{h,cc}
  gconf.nanoAOD_ver = 9;
  gconf.GetConfigs(year.Atoi());
  lumi = gconf.lumi;

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

  // Modify the name of the output file to include arguments of ScanChain function (i.e. process, year, etc.)
  TFile* fout = new TFile("temp_data/output_"+process+"_"+year+".root", "RECREATE");
  TTree* tout = new TTree("tout","Tree with photon variables");

  // define histograms, to be put in a different file TODO
  H1(leadPho_sieie, 20, 0, 0.05, "");
  H1(leadPho_phoIso, 20, 0, 10, "");
  H1(leadPho_chgIso, 20, 0, 10, "");
  H1(leadPho_trkIso, 20, 0, 10, "");
  H1(subleadPho_sieie, 20, 0, 0.05, "");
  H1(subleadPho_phoIso, 20, 0, 10, "");
  H1(subleadPho_chgIso, 20, 0, 10, "");
  H1(subleadPho_trkIso, 20, 0, 10, "");


  Float_t leadPho_pt, leadPho_eta, leadPho_phi;
  Float_t subleadPho_pt, subleadPho_eta, subleadPho_phi;
  Float_t diPho_pt, diPho_eta, diPho_phi, diPho_mass;
  Float_t leadPho_sieie, leadPho_phoIso, leadPho_trkIso, leadPho_chgIso, leadPho_mvaID;
  Float_t subleadPho_sieie, subleadPho_phoIso, subleadPho_trkIso, subleadPho_chgIso, subleadPho_mvaID;
  Int_t eventNum;

  tout->Branch("leadPho_pt",&leadPho_pt,"leadPho_pt/F");
  tout->Branch("leadPho_eta",&leadPho_eta,"leadPho_eta/F");
  tout->Branch("leadPho_phi",&leadPho_phi,"leadPho_phi/F");
  tout->Branch("leadPho_sieie",&leadPho_sieie,"leadPho_sieie/F");
  tout->Branch("leadPho_phoIso",&leadPho_phoIso,"leadPho_phoIso/F");
  tout->Branch("leadPho_trkIso",&leadPho_trkIso,"leadPho_trkIso/F");
  tout->Branch("leadPho_chgIso",&leadPho_chgIso,"leadPho_chgIso/F");
  tout->Branch("subleadPho_pt",&subleadPho_pt,"subleadPho_pt/F");
  tout->Branch("subleadPho_eta",&subleadPho_eta,"subleadPho_eta/F");
  tout->Branch("subleadPho_phi",&subleadPho_phi,"subleadPho_phi/F");
  tout->Branch("subleadPho_sieie",&subleadPho_sieie,"subleadPho_sieie/F");
  tout->Branch("subleadPho_phoIso",&subleadPho_phoIso,"subleadPho_phoIso/F");
  tout->Branch("subleadPho_trkIso",&subleadPho_trkIso,"subleadPho_trkIso/F");
  tout->Branch("subleadPho_chgIso",&subleadPho_chgIso,"subleadPho_chgIso/F");
  tout->Branch("diPho_pt",&diPho_pt,"diPho_pt/F");
  tout->Branch("diPho_eta",&diPho_eta,"diPho_eta/F");
  tout->Branch("diPho_phi",&diPho_phi,"diPho_phi/F");
  tout->Branch("diPho_mass",&diPho_mass,"diPho_mass/F");
  tout->Branch("eventNum",&eventNum,"eventNum/I");

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
  histoDefinition(nbins, low, high, binsx, title);

  // Define RooDataSet's for fit
  /*
  RooRealVar mfit("mfit", "mfit", 150.0, 3000.0);
  RooRealVar roow("roow", "roow", 1.0);
  map<TString, RooDataSet> roods;
  for ( unsigned int imll=0; imll < mllbin.size(); imll++ ) {
    for ( unsigned int inb=0; inb < nbtag.size(); inb++ ) {
      for (unsigned int iMuDet = 0; iMuDet < MuDetRegion.size(); iMuDet++){
        TString dname = TString("d_") + mllbin[imll] + TString("_") + nbtag[inb] + TString("_") + MuDetRegion[iMuDet];
	TString slice = mllbin[imll] + TString("_") + nbtag[inb] + TString("_") + MuDetRegion[iMuDet];
	if ( fillRooDataSet )
	  roods.insert({slice, RooDataSet(dname,dname,RooArgSet(mfit,roow),WeightVar(roow))});
      }
    }
  }
  */

  // Define selection
  cout << "define selection" << endl;
  vector<TString> selection = { };
  selection.push_back("sel0"); // Skimming + HLT + Good PV
  selection.push_back("sel1"); // 2 high-pT ID muons
  selection.push_back("sel2"); // pT > 53 GeV && |eta| < 2.4 muons

  vector<TString> plot_names = { };
  plot_names.push_back("pfmet_pt");

  map<TString, TH1F*> histos;

  //if ( PUWeight!=0 ) set_puWeights();

  // Setting up JEC uncertainties
  JetCorrectionUncertainty* jec_unc = new JetCorrectionUncertainty(
    "../NanoCORE/Tools/jetcorr/data/"
    + gconf.jecEraMC 
    + "/"
    + gconf.jecEraMC
    + "_Uncertainty_AK4PFchs.txt"
  );

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

    // Before any cuts
    int icutflow = 0;
    TString label = "Total";
    TString slicedlabel = label;

    for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

      nt.GetEntry(event);
      tree->LoadTree(event);

      nEventsTotal++;
      bar.progress(nEventsTotal, nEventsChain);

      float weight = 1.0;
      if ( isMC ) {
	weight = nt.genWeight();
	if(removeSpikes && weight*factor>1e2) continue;

	// Apply PU reweight
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
	  if ( is_duplicate(id) ){
	    ++nDuplicates;
	    continue;
	  }
	}
      }
/*
      // HLT selection
      if ( (year=="2016nonAPV" || year=="2016APV") &&
             !( (tree->GetBranch("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55") ? nt.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() : 0)
               || (tree->GetBranch("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55") ? nt.HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() : 0) ) ) continue;
      if ( (year=="2017") &&
             !( (tree->GetBranch("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55") ? nt.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55() : 0)  )  ) continue;
      if ( (year=="2018") &&
              !( (tree->GetBranch("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto") ? nt.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto() : 0) ) ) continue;
*/
      // do diphoton selection here
      Photons photons = getPhotons(); //sort by pt
      DiPhotons diphotons = DiPhotonPreselection(photons);

      h_weight_full->Fill(weight*factor);

      if (diphotons.size() == 0 ) continue; 

      DiPhoton selectedDiPhoton = diphotons[0];
      Photons vector_photons={};
      vector_photons.push_back(selectedDiPhoton.leadPho);
      vector_photons.push_back(selectedDiPhoton.subleadPho);
      leadPho_mvaID = selectedDiPhoton.leadPho.mvaID();
      subleadPho_mvaID = selectedDiPhoton.subleadPho.mvaID();
//      if (!(subleadPho_mvaID>-0.7 && leadPho_mvaID>-0.7)) continue;
/*
      Electrons electrons = getElectrons(photons);
      Muons muons = getMuons(photons);
      if (electrons.size() != 0 ) continue; 
      if (muons.size() != 0 ) continue; 

      Jets jets = getJets(photons); //sort by b score
      DiJets dijets = DiJetPreselection(jets);
      if (jets.size() < 2) continue; 
      if (dijets.size() == 0 ) continue; 
*/
      leadPho_pt = selectedDiPhoton.leadPho.pt();
      leadPho_eta = selectedDiPhoton.leadPho.eta();
      leadPho_phi = selectedDiPhoton.leadPho.phi();
      leadPho_sieie = selectedDiPhoton.leadPho.sieie();
      leadPho_phoIso = selectedDiPhoton.leadPho.phoIso();
      leadPho_trkIso = selectedDiPhoton.leadPho.trkIso();
      leadPho_chgIso = selectedDiPhoton.leadPho.chargedHadIso();
      subleadPho_pt = selectedDiPhoton.subleadPho.pt();
      subleadPho_eta = selectedDiPhoton.subleadPho.eta();
      subleadPho_phi = selectedDiPhoton.subleadPho.phi();
      subleadPho_sieie = selectedDiPhoton.subleadPho.sieie();
      subleadPho_phoIso = selectedDiPhoton.subleadPho.phoIso();
      subleadPho_trkIso = selectedDiPhoton.subleadPho.trkIso();
      subleadPho_chgIso = selectedDiPhoton.subleadPho.chargedHadIso();
      diPho_pt = selectedDiPhoton.p4.Pt();
      diPho_eta = selectedDiPhoton.p4.Eta();
      diPho_phi = selectedDiPhoton.p4.Phi();
      diPho_mass = selectedDiPhoton.p4.M();
      count_test++;
      if (count_test>65025) cout<<diPho_mass<<" "<<diPho_pt<<" "<<diPho_eta<<" "<<diPho_phi<<" "<<photons.size()<<endl;

      h_leadPho_sieie->Fill(leadPho_sieie);
      h_leadPho_phoIso->Fill(leadPho_phoIso);
      h_leadPho_chgIso->Fill(leadPho_chgIso);
      h_leadPho_trkIso->Fill(leadPho_trkIso);
      h_subleadPho_sieie->Fill(subleadPho_sieie);
      h_subleadPho_phoIso->Fill(subleadPho_phoIso);
      h_subleadPho_chgIso->Fill(subleadPho_chgIso);
      h_subleadPho_trkIso->Fill(subleadPho_trkIso);

      tout->Fill();
      h_weight->Fill(weight*factor);
    } // Event loop

    delete file;

  } // File loop
    
  bar.finish();
  cout << "nTotal: " << h_weight_full->GetBinContent(1) << ", nPass: " << h_weight->GetBinContent(1) << ", eff: " << h_weight->GetBinContent(1)/h_weight_full->GetBinContent(1) << endl;
  cout << endl;

  if ( muonSF!=0 ) {
    reset_muonRecoSF();
    reset_muonIDSF();
    reset_muonIsoSF();
  }
  if ( triggerSF!=0 ) reset_triggerSF();
  if ( bTagSF!=0 ) reset_bTagEff();

  if ( removeDataDuplicates )
    cout << "Number of duplicates found: " << nDuplicates << endl;

  // Avoid histograms with unphysical negative bin content (due to negative GEN weights)
  map<TString, TH1F*>::iterator it;
  for ( it = histos.begin(); it != histos.end(); it++ ) {
    for ( unsigned int b=1; b<(it->second)->GetNbinsX()+1; b++ ) { 
      if ( (it->second)->GetBinContent(b)<0.0) {
	(it->second)->SetBinContent(b,0.0);
	(it->second)->SetBinError(b,0.0);
      }
    }
  }

  //tout->Write();
  fout->Write();
  //if ( fillRooDataSet ) {
  //  fout->cd();
  //  for (const auto& d : roods )
  //    d.second.Write();
  //}
  fout->Close();

  return 0;
}
