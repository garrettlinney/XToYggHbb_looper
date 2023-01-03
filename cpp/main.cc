#include "TROOT.h"
#include "ScanChain_Hgg.C"
#include <fstream>

double getSumOfGenEventSumw(TChain *chaux, bool isMC)
{
  double genEventSumw, sumOfGenEventSumw=0.0;
  if (isMC) {
    chaux->SetBranchAddress("genEventSumw",&genEventSumw);
    for (unsigned int run = 0; run < chaux->GetEntriesFast(); run++)
      {
        chaux->GetEntry(run);
        sumOfGenEventSumw += genEventSumw;
      }
  }
  else{
    sumOfGenEventSumw=1;
  }
  return sumOfGenEventSumw;
}

int main() {
  // Event weights / scale factors:
  //  0: Do not apply
  //  1: Apply central value
  // +2: Apply positive variation
  // -2: Apply negative variation

  vector<TString> samples = { };
  map<TString,TString> sample_names = { };
  map<TString,map<TString,vector<TString>>> sample_prod = { };
  map<TString,map<TString,int>> sample_nfiles = { };

  // TODO is there a way to read formatted config (.json, .yaml. .config) using C++?
/*
  // ggH 
  samples.push_back("ggHToDiPhoM125");
  sample_names.insert({"ggHToDiPhoM125","GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
  sample_nfiles.insert({"ggHToDiPhoM125", {{"2018", 11 } } });
  sample_prod.insert({"ggHToDiPhoM125", { { "2018",       { "GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2_TESTS/220225_213121/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // diPhoton 
  samples.push_back("diPhoton");
  sample_names.insert({"diPhoton","DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa"});
  sample_nfiles.insert({"diPhoton", {{"2018", 39 } } });
  sample_prod.insert({"diPhoton", { { "2018",       { "DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_RunIISummer19UL18MiniAODv2/skimNano-TestUL_DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM_final_TESTS/220225_204254/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
//  "test_sync_samples/"
*/
/*
   // DY
  samples.push_back("DY");
  sample_names.insert({"DY","DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX"});
  sample_nfiles.insert({"DY", {{"2018", 121 } } });
  sample_prod.insert({"DY", { { "2018",       { "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX_2018/skimNano-TestUL_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX_2018_TESTS/220712_183424/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
*/
/*
  // HHggtautau 
  samples.push_back("HHggtautau");
  sample_names.insert({"HHggtautau","GluGluToHHTo2G2Tau_node_cHHH1_TuneCP5_13TeV-powheg-pythia8"});
  sample_nfiles.insert({"HHggtautau", {{"2018", 5 } } });
  sample_prod.insert({"HHggtautau", { { "2018",       { "GluGluToHHTo2G2Tau_node_cHHH1_TuneCP5_13TeV-powheg-pythia8_2018_final/skimNano-TestUL_GluGluToHHTo2G2Tau_node_cHHH1_TuneCP5_13TeV-powheg-pythia8_2018_final_TESTS/220225_214005/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
*/

  // partial 2018D data 
  samples.push_back("EGamma_Run2018D");
  sample_names.insert({"EGamma_Run2018D","EGamma_Run2018D-UL2018"});
  sample_nfiles.insert({"EGamma_Run2018D", {{"2018", 498 } } });
  sample_prod.insert({"EGamma_Run2018D", { { "2018",       { "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220403_112326/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  //DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_RunIISummer19UL18MiniAODv2/skimNano-TestUL_DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM_final_TESTS/220225_204254/0000/
  //GluGluToHHTo2G2Tau_node_cHHH1_TuneCP5_13TeV-powheg-pythia8_2018_final/skimNano-TestUL_GluGluToHHTo2G2Tau_node_cHHH1_TuneCP5_13TeV-powheg-pythia8_2018_final_TESTS/220225_214005/0000/
  //EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220403_112326/0000/
  //EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220403_112326/0000/tree_


  TString year = "2018";
//  TString basedir = "/ceph/cms/store/user/yagu/Run2018/";
  TString basedir = "/ceph/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/";
  
  int topPtWeight=1;
  int PUWeight=1;
  int muonSF=1;
  int triggerSF=1;
  int bTagSF=1;
  int JECUnc=0; // No central value, set to +/-2 to get

  for ( int isample=0; isample<samples.size(); isample++ )
  {

      TString sample = samples[isample];

//      if (sample == "EGamma_Run2018D") continue; //TODO need to check if file exist, if not increase counter by 1
      //if (sample != "HHggtautau") continue;
                                                 //
      int nFiles = sample_nfiles[sample][year];
      TString prod = sample_prod[sample][year][0];

      TChain *ch_temp = new TChain("Events");
      TChain *chaux_temp = new TChain("Runs");

      cout << "nFiles: " << nFiles << endl;
      ifstream i_file;
      int filenameshifter = 0;
      for ( int ifile=1; ifile<=nFiles; ifile++ )
      {

          TString filename = Form(basedir + prod + "tree_%d.root", ifile+filenameshifter);

          ifstream i_file;
          i_file.open(filename);
          while ( !i_file )
          {
              filenameshifter++;
              filename = Form(basedir + prod + "tree_%d.root", ifile+filenameshifter);
              i_file.open(filename);
          }
//          cout<<filename<<endl;
          ch_temp->Add(filename);
          chaux_temp->Add(filename);

      }
      bool isMC=1;
      if (sample == "EGamma_Run2018D") isMC=0;
      cout<<1<<endl;
      ScanChain_Hgg(ch_temp,getSumOfGenEventSumw(chaux_temp, isMC),year,sample,topPtWeight,PUWeight,muonSF,triggerSF,bTagSF,JECUnc);
      cout<<2<<endl;
  }

}
