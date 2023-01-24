#include "TROOT.h"
#include "ScanChain_Hgg.C"
#include <fstream>

double getSumOfGenEventSumw(TChain *chaux, bool isMC)
{
  double genEventSumw, sumOfGenEventSumw=0.0;
  if (isMC) {
    chaux->SetBranchAddress("genEventSumw",&genEventSumw);
    for (unsigned int run = 0; run < chaux->GetEntriesFast(); run++) {
      chaux->GetEntry(run);
      sumOfGenEventSumw += genEventSumw;
    }
  }
  else {
    sumOfGenEventSumw=1.0;
  }
  return sumOfGenEventSumw;
}

int main() {
  // Map definitions
  vector<TString> samples = { };
  map<TString,TString> sample_names = { };
  map<TString,map<TString,vector<TString>>> sample_prod = { };

  // TODO is there a way to read formatted config (.json, .yaml. .config) using C++?
/*
  // 2018D data
  samples.push_back("EGamma_Run2018D");
  sample_names.insert({"EGamma_Run2018D","EGamma_Run2018D-UL2018"});
  sample_prod.insert({"EGamma_Run2018D", { { "2018",       { "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220308_013616/0000/",
                                                             "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220308_013616/0001/",
                                                             "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220308_013616/0002/"
                                                            } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });


  samples.push_back("EGamma_Run2018D_more");
  sample_names.insert({"EGamma_Run2018D_more","EGamma_Run2018D-UL2018"});
  sample_prod.insert({"EGamma_Run2018D_more", { { "2018",       { "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220330_223617/0000/",
                                                                  "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220330_223617/0001/",
                                                                  "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220330_223617/0002/",
                                                                  "EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final/skimNano-TestUL_EGamma_Run2018D-UL2018_MiniAODv2-v2_MINIAOD_final_TESTS/220403_112326/0000/"
                                                                } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
                       
  samples.push_back("EGamma_Run2018C");
  sample_names.insert({"EGamma_Run2018C","EGamma_Run2018C-UL2018"});
  sample_prod.insert({"EGamma_Run2018C", { { "2018",       { "EGamma_Run2018C-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018C-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220308_013506/0000/",
                                                             "EGamma_Run2018C-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018C-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220330_223512/0000/",
                                                             "EGamma_Run2018C-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018C-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220403_112220/0000/"
                                                            } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  samples.push_back("EGamma_Run2018B");
  sample_names.insert({"EGamma_Run2018B","EGamma_Run2018B-UL2018"});
  sample_prod.insert({"EGamma_Run2018B", { { "2018",       { "EGamma_Run2018B-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018B-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220308_013359/0000/",
                                                             "EGamma_Run2018B-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018B-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220330_223407/0000/",
                                                             "EGamma_Run2018B-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018B-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220403_112116/0000/"
                                                            } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
             
  samples.push_back("EGamma_Run2018A");
  sample_names.insert({"EGamma_Run2018A","EGamma_Run2018A-UL2018"});
  sample_prod.insert({"EGamma_Run2018A", { { "2018",       { "EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220308_013253/0000/",
                                                             "EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220308_013253/0001/",
                                                             "EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220330_223303/0000/",
                                                             "EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220330_223303/0001/",
                                                             "EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_EGamma_Run2018A-UL2018_MiniAODv2-v1_MINIAOD_final_TESTS/220403_112010/0000/"
                                                            } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // diPhoton 
  samples.push_back("diPhoton");
  sample_names.insert({"diPhoton","DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa"});
  sample_prod.insert({"diPhoton", { { "2018",       { "DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_RunIISummer19UL18MiniAODv2/skimNano-TestUL_DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_RunIISummer19UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM_final_TESTS/220225_204254/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // TTGG
  samples.push_back("TTGG");
  sample_names.insert({"TTGG","TTGG_0Jets_TuneCP5_13TeV-amcatnlo-madspin-pythia8"});
  sample_prod.insert({"TTGG", { { "2018",       { "TTGG_0Jets_TuneCP5_13TeV-amcatnlo-madspin-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_TTGG_0Jets_TuneCP5_13TeV-amcatnlo-madspin-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_222646/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // TTGJets
  samples.push_back("TTGJets");
  sample_names.insert({"TTGJets","TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
  sample_prod.insert({"TTGJets", { { "2018",       { "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18MiniAODv2/skimNano-TestUL_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18MiniAODv2_TESTS/220225_223107/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // TTJets                               
  samples.push_back("TTJets");
  sample_names.insert({"TTJets","TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
  sample_prod.insert({"TTJets", { { "2018",       { "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_223214/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // VH_M125
  samples.push_back("VH_M125");
  sample_names.insert({"VH_M125","VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
  sample_prod.insert({"VH_M125", { { "2018",       { "VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18MiniAODv2/skimNano-TestUL_VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18MiniAODv2_TESTS/220225_224148/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });  

  // VBFH_M125
  samples.push_back("VBFH_M125");
  sample_names.insert({"VBFH_M125","VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8"});
  sample_prod.insert({"VBFH_M125", { { "2018",       { "VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8_storeWeights_RunIISummer19UL18MiniAODv2/skimNano-TestUL_VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8_storeWeights_RunIISummer19UL18MiniAODv2_TESTS/220225_223810/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // ttH_M125
  samples.push_back("ttH_M125");
  sample_names.insert({"ttH_M125","ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
  sample_prod.insert({"ttH_M125", { { "2018",       { "ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18MiniAODv2/skimNano-TestUL_ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8_RunIISummer20UL18MiniAODv2_TESTS/220225_231014/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // ggH 
  samples.push_back("ggHToDiPhoM125");
  sample_names.insert({"ggHToDiPhoM125","GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
  sample_prod.insert({"ggHToDiPhoM125", { { "2018",       { "GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights_RunIISummer19UL18MiniAODv2_TESTS/220225_213121/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // GJets_HT-40To100
  samples.push_back("GJets_HT-40To100");
  sample_names.insert({"GJets_HT-40To100","GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8"});
  sample_prod.insert({"GJets_HT-40To100", { { "2018",       { "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_211817/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
              
  // GJets_HT-100To200
  samples.push_back("GJets_HT-100To200");
  sample_names.insert({"GJets_HT-100To200","GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8"});
  sample_prod.insert({"GJets_HT-100To200", { { "2018",       { "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_210621/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // GJets_HT-200To400
  samples.push_back("GJets_HT-200To400");
  sample_names.insert({"GJets_HT-200To400","GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8"});
  sample_prod.insert({"GJets_HT-200To400", { { "2018",       { "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_211040/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // GJets_HT-400To600
  samples.push_back("GJets_HT-400To600");
  sample_names.insert({"GJets_HT-400To600","GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8"});
  sample_prod.insert({"GJets_HT-400To600", { { "2018",       { "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_211501/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // GJets_HT-600ToInf
  samples.push_back("GJets_HT-600ToInf");
  sample_names.insert({"GJets_HT-600ToInf","GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8"});
  sample_prod.insert({"GJets_HT-600ToInf", { { "2018",       { "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2/skimNano-TestUL_GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer19UL18MiniAODv2_TESTS/220225_212236/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

   // DY
  samples.push_back("DY");
  sample_names.insert({"DY","DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX"});
  sample_prod.insert({"DY", { { "2018",       { "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX_2018/skimNano-TestUL_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX_2018_TESTS/220712_183424/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

   // WGamma
  samples.push_back("WGamma");
  sample_names.insert({"WGamma","WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
  sample_prod.insert({"WGamma", { { "2018",       { "WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2/skimNano-TestUL_WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_TESTS/220225_224602/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

   // ZGamma
  samples.push_back("ZGamma");
  sample_names.insert({"ZGamma","ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
  sample_prod.insert({"ZGamma", { { "2018",       { "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2/skimNano-TestUL_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_TESTS/220225_230245/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
*/

  // HHbbgg
  samples.push_back("HHbbgg");
  sample_names.insert({"HHbbgg","GluGluToHHTo2B2G_node_cHHH1_TuneCP5_13TeV-powheg-pythia8"});
  sample_prod.insert({"HHbbgg", { { "2018",       { "GluGluToHHTo2B2G_node_cHHH1_TuneCP5_13TeV-powheg-pythia8_2018_v0/skimNano-TestUL_GluGluToHHTo2B2G_node_cHHH1_TuneCP5_13TeV-powheg-pythia8_2018_v0_TESTS/220729_201156/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
/*
  // NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90
  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90","NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_600_MY_90_2018_v0_070722_TESTS/220711_225445/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95
  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95","NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_600_MY_95_2018_v0_070722_TESTS/220711_225902/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });


  // NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100
  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100","NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_600_MY_100_2018_v0_070722_TESTS/220711_225030/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  // NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100
  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_650_MY_90");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_650_MY_90","NMSSM_XYH_Y_gg_H_bb_MX_650_MY_90"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_650_MY_90", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_650_MY_90_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_650_MY_90_2018_v0_070722_TESTS/220711_230731/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_650_MY_95");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_650_MY_95","NMSSM_XYH_Y_gg_H_bb_MX_650_MY_95"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_650_MY_95", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_650_MY_95_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_650_MY_95_2018_v0_070722_TESTS/220711_231148/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100","NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_650_MY_100_2018_v0_070722_TESTS/220711_230317/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_700_MY_90");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_700_MY_90","NMSSM_XYH_Y_gg_H_bb_MX_700_MY_90"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_700_MY_90", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_700_MY_90_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_700_MY_90_2018_v0_070722_TESTS/220711_232022/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_700_MY_95");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_700_MY_95","NMSSM_XYH_Y_gg_H_bb_MX_700_MY_95"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_700_MY_95", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_700_MY_95_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_700_MY_95_2018_v0_070722_TESTS/220711_232439/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });

  samples.push_back("NMSSM_XYH_Y_gg_H_bb_MX_700_MY_100");
  sample_names.insert({"NMSSM_XYH_Y_gg_H_bb_MX_700_MY_100","NMSSM_XYH_Y_gg_H_bb_MX_700_MY_100"});
  sample_prod.insert({"NMSSM_XYH_Y_gg_H_bb_MX_700_MY_100", { { "2018",       { "NMSSM_XYH_Y_gg_H_bb_MX_700_MY_100_2018_v0_070722/skimNano-TestUL_NMSSM_XYH_Y_gg_H_bb_MX_700_MY_100_2018_v0_070722_TESTS/220711_231606/0000/" } },
                                 { "2017",       { "" } },
                                 { "2016APV",    { "" } },
                                 { "2016nonAPV", { "" } } } });
*/
  TString year = "2018";
  TString basedir = "/ceph/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/";
  
  int PUWeight=1;
  int bTagSF=1;
  int JECUnc=0; // No central value, set to +/-2 to get

  for ( int isample=0; isample<samples.size(); isample++ ) {
    TString sample = samples[isample];

    for (int i_subfile=0; i_subfile<sample_prod[sample][year].size(); i_subfile++){ // FIXME: To be simplified with simpler directory structure
      TString prod = sample_prod[sample][year][i_subfile];

      TChain *ch_temp = new TChain("Events");
      TChain *chaux_temp = new TChain("Runs");

      TString filename = Form(basedir + prod + "tree_*.root");
      ch_temp->Add(filename);
      chaux_temp->Add(filename);

      bool isMC=1;
      if ( sample.Contains("_Run201") ) isMC=0;
      ScanChain_Hgg(ch_temp,getSumOfGenEventSumw(chaux_temp, isMC),year,sample,PUWeight,bTagSF,JECUnc);
    }
  }
}
