import sys,os,copy
import math
import ROOT

def print_header(fout):
    fout.write('#ifndef MUONRECOSF_H\n')  
    fout.write('#define MUONRECOSF_H\n')  
    fout.write('\n')
    fout.write('#include <map>\n')  
    fout.write('#include <vector>\n')  
    fout.write('#include <iostream>\n')  
    fout.write('#include <TString.h>\n')  
    fout.write('\n')
    fout.write('using namespace std;\n')
    fout.write('\n')
    fout.write('map<TString,map<TString,map<TString,float>>> recosfhpt = { };\n')
    fout.write('map<TString,map<TString,map<TString,float>>> recosfunchpt = { };\n')
    fout.write('map<TString,map<TString,map<TString,float>>> recosfmpt = { };\n')
    fout.write('map<TString,map<TString,map<TString,float>>> recosfuncmpt = { };\n')
    fout.write('\n')
    fout.write('vector<float> thresholds_p_hpt;\n')
    fout.write('vector<TString> pbins_hpt;\n')
    fout.write('inline void set_pbins_reco_hpt() {\n')
    tf = ROOT.TFile(indir+"/muon_highpt_recosf_2018.root")
    th   = tf.Get("hHighPtMuonRecoSF").Clone("thist")
    for b in range(1, th.GetNbinsY()+1):
        fout.write('\tthresholds_p_hpt.push_back('+str(th.GetYaxis().GetBinLowEdge(b))+');\n')
        fout.write('\tpbins_hpt.push_back("p'+str(b-1)+'");\n')
    tf.Close()
    fout.write('}\n')
    fout.write('\n')
    fout.write('vector<float> thresholds_eta_hpt;\n')
    fout.write('vector<TString> etabins_hpt;\n')
    fout.write('inline void set_etabins_reco_hpt() {\n')
    tf = ROOT.TFile(indir+"/muon_highpt_recosf_2018.root")
    th   = tf.Get("hHighPtMuonRecoSF").Clone("thist")
    for b in range(1, th.GetNbinsX()+1):
        fout.write('\tthresholds_eta_hpt.push_back('+str(th.GetXaxis().GetBinLowEdge(b))+');\n')
        fout.write('\tetabins_hpt.push_back("eta'+str(b-1)+'");\n')
    tf.Close()
    fout.write('}\n')
    fout.write('\n')
    fout.write('vector<float> thresholds_pt_mpt;\n')
    fout.write('vector<TString> ptbins_mpt;\n')
    fout.write('inline void set_ptbins_reco_mpt() {\n')
    tf = ROOT.TFile(indir+"/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root")
    th   = tf.Get("NUM_TrackerMuons_DEN_genTracks").Clone("thist")
    for b in range(1, th.GetNbinsY()+1):
        fout.write('\tthresholds_pt_mpt.push_back('+str(th.GetYaxis().GetBinLowEdge(b))+');\n')
        fout.write('\tptbins_mpt.push_back("pt'+str(b-1)+'");\n')
    tf.Close()
    fout.write('}\n')
    fout.write('\n')
    fout.write('vector<float> thresholds_eta_mpt;\n')
    fout.write('vector<TString> etabins_mpt;\n')
    fout.write('inline void set_etabins_reco_mpt() {\n')
    tf = ROOT.TFile(indir+"/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root")
    th   = tf.Get("NUM_TrackerMuons_DEN_genTracks").Clone("thist")
    for b in range(1, th.GetNbinsX()+1):
        fout.write('\tthresholds_eta_mpt.push_back('+str(th.GetXaxis().GetBinLowEdge(b))+');\n')
        fout.write('\tetabins_mpt.push_back("eta'+str(b-1)+'");\n')
    tf.Close()
    fout.write('}\n')
    fout.write('\n')

def print_function_h(fout):
    fout.write('TString get_muonRecoSFPBin_hpt(const float p);\n')
    fout.write('\n')
    fout.write('TString get_muonRecoSFEtaBin_hpt(const float abseta);\n')
    fout.write('\n')
    fout.write('TString get_muonRecoSFPtBin_mpt(const float pt);\n')
    fout.write('\n')
    fout.write('TString get_muonRecoSFEtaBin_mpt(const float abseta);\n')
    fout.write('\n')
    fout.write('float get_muonRecoSF( const float p, const float pt, const float eta, const TString year, const TString variation="central" );\n')
    fout.write('\n')
    fout.write('inline void reset_muonRecoSF() {\n')
    fout.write('\n')
    fout.write('\trecosfhpt.clear();\n')
    fout.write('\trecosfunchpt.clear();\n')
    fout.write('\trecosfmpt.clear();\n')
    fout.write('\trecosfuncmpt.clear();\n')
    fout.write('\n')
    fout.write('}\n')
    fout.write('\n')
    fout.write('inline void set_muonRecoSF() {\n')
    fout.write('\tset_pbins_reco_hpt();\n')
    fout.write('\tset_etabins_reco_hpt();\n')
    fout.write('\tset_ptbins_reco_mpt();\n')
    fout.write('\tset_etabins_reco_mpt();\n')
    fout.write('\n')

def print_footer(fout):
    fout.write('}\n')
    fout.write('\n')
    fout.write('#endif\n')
    
def get_muonRecoSF(fout, year = "2018"):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    tf = ROOT.TFile(indir+"/muon_highpt_recosf_"+year+".root")
    th   = tf.Get("hHighPtMuonRecoSF").Clone("thist")
    sfname = "recosfhpt"
    uncname = "recosfunchpt"
    fout.write('\t'+sfname +'.insert({"'+year+'", { }});\n')
    fout.write('\t'+uncname+'.insert({"'+year+'", { }});\n')
    for bx in range(1, th.GetNbinsX()+1):
        fout.write('\t'+sfname +'["'+year+'"].insert({"eta'+str(bx-1)+'", { }});\n')
        fout.write('\t'+uncname+'["'+year+'"].insert({"eta'+str(bx-1)+'", { }});\n')
        for by in range(1, th.GetNbinsY()+1):
            tsf  = th.GetBinContent(bx, by)
            tunc = th.GetBinError  (bx, by)
            fout.write('\t'+sfname +'["'+year+'"]["eta'+str(bx-1)+'"].insert({"p'+str(by-1)+'", '+str(tsf) +'});\n')
            fout.write('\t'+uncname+'["'+year+'"]["eta'+str(bx-1)+'"].insert({"p'+str(by-1)+'", '+str(tunc)+'});\n')
    tf.Close()
    fout.write('\n')

    tf = ROOT.TFile(indir+"/Efficiency_muon_generalTracks_Run"+year+"_UL_trackerMuon.root")
    th   = tf.Get("NUM_TrackerMuons_DEN_genTracks").Clone("thist")
    sfname = "recosfmpt"
    uncname = "recosfuncmpt"
    fout.write('\t'+sfname +'.insert({"'+year+'", { }});\n')
    fout.write('\t'+uncname+'.insert({"'+year+'", { }});\n')
    for bx in range(1, th.GetNbinsX()+1):
        fout.write('\t'+sfname +'["'+year+'"].insert({"eta'+str(bx-1)+'", { }});\n')
        fout.write('\t'+uncname+'["'+year+'"].insert({"eta'+str(bx-1)+'", { }});\n')
        for by in range(1, th.GetNbinsY()+1):
            tsf  = th.GetBinContent(bx, by)
            tunc = th.GetBinError  (bx, by)
            fout.write('\t'+sfname +'["'+year+'"]["eta'+str(bx-1)+'"].insert({"pt'+str(by-1)+'", '+str(tsf) +'});\n')
            fout.write('\t'+uncname+'["'+year+'"]["eta'+str(bx-1)+'"].insert({"pt'+str(by-1)+'", '+str(tunc)+'});\n')
    tf.Close()
    fout.write('\n')

def print_cc(fout,fname):
    fout.write('#include "%s.h"\n'%fname)
    fout.write('\n')
    fout.write('TString get_muonRecoSFPBin_hpt(const float p) {\n')
    fout.write('\tfor ( unsigned int b=thresholds_p_hpt.size()-1; b>=1; b-- ) {\n')
    fout.write('\t\tif ( p > thresholds_p_hpt.at(b) ) {\n')
    fout.write('\t\t\treturn pbins_hpt.at(b);\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('\treturn pbins_hpt.at(0);\n')
    fout.write('}\n')
    fout.write('\n')
    fout.write('TString get_muonRecoSFEtaBin_hpt(const float abseta) {\n')
    fout.write('\tfor ( unsigned int b=thresholds_eta_hpt.size()-1; b>=1; b-- ) {\n')
    fout.write('\t\tif ( abseta > thresholds_eta_hpt.at(b) ) {\n')
    fout.write('\t\t\treturn etabins_hpt.at(b);\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('\treturn etabins_hpt.at(0);\n')
    fout.write('}\n')
    fout.write('\n')
    fout.write('TString get_muonRecoSFPtBin_mpt(const float pt) {\n')
    fout.write('\tfor ( unsigned int b=thresholds_pt_mpt.size()-1; b>=1; b-- ) {\n')
    fout.write('\t\tif ( pt > thresholds_pt_mpt.at(b) ) {\n')
    fout.write('\t\t\treturn ptbins_mpt.at(b);\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('\treturn ptbins_mpt.at(0);\n')
    fout.write('}\n')
    fout.write('\n')
    fout.write('TString get_muonRecoSFEtaBin_mpt(const float abseta) {\n')
    fout.write('\tfor ( unsigned int b=thresholds_eta_mpt.size()-1; b>=1; b-- ) {\n')
    fout.write('\t\tif ( abseta > thresholds_eta_mpt.at(b) ) {\n')
    fout.write('\t\t\treturn etabins_mpt.at(b);\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('\treturn etabins_mpt.at(0);\n')
    fout.write('}\n')
    fout.write('\n')
    fout.write('float get_muonRecoSF(const float p, const float pt, const float eta, const TString year, const TString variation ) {\n')
    fout.write('\tif ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {\n')
    fout.write('\t\tstd::cout << "WARNING: unknown year, returning unity muon RECO SF!" << std::endl;\n')
    fout.write('\t\treturn 1.0;\n')
    fout.write('\t}\n')
    fout.write('\tif ( pt < 10.0 || fabs(eta) > 2.4 ) return 1.0;\n')
    fout.write('\tif ( p > 50.0 && ( fabs(eta) < 1.6 || p > 100.0 ) ) {\n')
    fout.write('\t\tTString etabin = get_muonRecoSFEtaBin_hpt(fabs(eta));\n')
    fout.write('\t\tTString pbin = get_muonRecoSFPBin_hpt(p);\n')
    fout.write('\t\tif ( variation == "central" ) return recosfhpt[year][etabin][pbin];\n')
    fout.write('\t\telse if ( variation == "up" ) return recosfhpt[year][etabin][pbin]+recosfunchpt[year][etabin][pbin];\n')
    fout.write('\t\telse if ( variation == "down" ) return recosfhpt[year][etabin][pbin]-recosfunchpt[year][etabin][pbin];\n')
    fout.write('\t\telse {\n')
    fout.write('\t\t\tstd::cout << "WARNING: unknown variation (central, up, down), returning unity scale factor!" << std::endl;\n')
    fout.write('\t\t\treturn 1.0;\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('\telse {\n')
    fout.write('\t\tTString etabin = get_muonRecoSFEtaBin_mpt(fabs(eta));\n')
    fout.write('\t\tTString ptbin = get_muonRecoSFPtBin_mpt(pt);\n')
    fout.write('\t\tif ( variation == "central" ) return recosfmpt[year][etabin][ptbin];\n')
    fout.write('\t\telse if ( variation == "up" ) return recosfmpt[year][etabin][ptbin]+recosfuncmpt[year][etabin][ptbin];\n')
    fout.write('\t\telse if ( variation == "down" ) return recosfmpt[year][etabin][ptbin]-recosfuncmpt[year][etabin][ptbin];\n')
    fout.write('\t\telse {\n')
    fout.write('\t\t\tstd::cout << "WARNING: unknown variation (central, up, down), returning unity scale factor!" << std::endl;\n')
    fout.write('\t\t\treturn 1.0;\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('}\n')

indir = "./data/"
outdir="./NanoCORE/Tools/"
fname = "muonRecoSF"

fcc = open(outdir+"/"+fname+".cc",'w')
print_cc(fcc, fname)

fheader = open(outdir+"/"+fname+".h",'w')
print_header(fheader)
print_function_h(fheader)

for year in ["2016nonAPV","2016APV","2017","2018"]:
    get_muonRecoSF(fheader, year)

print_footer(fheader)
