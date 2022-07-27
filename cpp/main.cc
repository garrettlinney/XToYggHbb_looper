#include "TROOT.h"
#include "ScanChain_Hgg.C"

double getSumOfGenEventSumw(TChain *chaux)
{
  double genEventSumw, sumOfGenEventSumw=0.0;
  chaux->SetBranchAddress("genEventSumw",&genEventSumw);
  for (unsigned int run = 0; run < chaux->GetEntriesFast(); run++)
    {
      chaux->GetEntry(run);
      sumOfGenEventSumw += genEventSumw;
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

  // ggH 
  samples.push_back("ggHToDiPhoM125");
  sample_names.insert({"ggHToDiPhoM125","GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
  sample_prod.insert({"ggHToDiPhoM125", { { "2018",       { "RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                 { "2017",       { "RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1" } },
                                 { "2016APV",    { "RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                 { "2016nonAPV", { "RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1" } } } });

  TString year = "2018";
  
  int topPtWeight=1;
  int PUWeight=1;
  int muonSF=1;
  int triggerSF=1;
  int bTagSF=1;
  int JECUnc=0; // No central value, set to +/-2 to get

  TChain *ch_temp = new TChain("Events");
  TChain *chaux_temp = new TChain("Runs");

  ch_temp->Add("/ceph/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2_TESTS/220225_213121/0000/tree_11.root");
  chaux_temp->Add("/ceph/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2_TESTS/220225_213121/0000/tree_11.root");

  TString sample = samples[0];
  ScanChain(ch_temp,getSumOfGenEventSumw(chaux_temp),year,sample,topPtWeight,PUWeight,muonSF,triggerSF,bTagSF,JECUnc);

}
