import sys,os,copy
import math
import ROOT

def print_header(fout):
    fout.write('#ifndef PHOMVAIDWP90SF_H\n')  
    fout.write('#define PHOMVAIDWP90SF_H\n')  
    fout.write('\n')
    fout.write('#include <map>\n')  
    fout.write('#include <iostream>\n')  
    fout.write('#include <TString.h>\n')  
    fout.write('\n')
    fout.write('using namespace std;\n')
    fout.write('\n')
    fout.write('namespace phoMVAIDWP90SF {\n')
    fout.write('\n')
    fout.write('\tmap<TString,map<TString,map<TString,float>>> idsf = { };\n')
    fout.write('\tmap<TString,map<TString,map<TString,float>>> idsfunc = { };\n')
    fout.write('\n')

def print_function_h(fout):
    fout.write('\tTString get_phoMVAIDWP90SFPtBin(const float pt);\n')
    fout.write('\n')
    fout.write('\tTString get_phoMVAIDWP90SFEtaBin(const float eta);\n')
    fout.write('\n')
    fout.write('\tfloat get_phoMVAIDWP90SF( const float pt, const float eta, const TString year, const TString variation );\n')
    fout.write('\n')
    fout.write('\tinline void reset_phoMVAIDWP90SF() {\n')
    fout.write('\t\tidsf.clear();\n')
    fout.write('\t\tidsfunc.clear();\n')
    fout.write('\t}\n')
    fout.write('\n')

    tf = ROOT.TFile(indir+"EGM2D_PHO_MVA90_UL2018.root")
    th = tf.Get("EGamma_SF2D").Clone("thist")
    fout.write('\tconstexpr int nptbins = '+str(th.GetNbinsY())+';\n')
    fout.write('\tconstexpr int netabins = '+str(th.GetNbinsX())+';\n')
    fout.write('\tfloat *thresholds_pt = new float[nptbins];\n')
    fout.write('\tTString *ptbins = new TString[nptbins];\n')
    fout.write('\tfloat *thresholds_eta = new float[netabins];\n')
    fout.write('\tTString *etabins = new TString[netabins];\n')
    fout.write('\n')
    fout.write('\tinline void set_ranges() {\n')
    for b in range(1, th.GetNbinsY()+1):
        fout.write('\t\tthresholds_pt['+str(b-1)+'] = '+str(th.GetYaxis().GetBinLowEdge(b))+';\n')
        fout.write('\t\tptbins['+str(b-1)+'] = "pt'+str(b-1)+'";\n')
    fout.write('\n')
    for b in range(1, th.GetNbinsX()+1):
        fout.write('\t\tthresholds_eta['+str(b-1)+'] = '+str(th.GetXaxis().GetBinLowEdge(b))+';\n')
        fout.write('\t\tetabins['+str(b-1)+'] = "eta'+str(b-1)+'";\n')
    tf.Close()
    fout.write('\t}\n')
    fout.write('\n')

def print_footer(fout):
    fout.write('}')
    fout.write('\n')
    fout.write('#endif\n')
    
def get_phoMVAIDWP90SF(fout):
    fout.write('\tinline void set_phoMVAIDWP90SF() {\n')
    fout.write('\t\tset_ranges();\n')
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for year in ["2016nonAPV","2016APV","2017","2018"]:
        tf = ROOT.TFile(indir+"EGM2D_PHO_MVA90_UL"+year+".root")
        th   = tf.Get("EGamma_SF2D").Clone("thist")
        sfname = "idsf"
        uncname = "idsfunc"
        fout.write('\t\t'+sfname +'.insert({"'+year+'", { }});\n')
        fout.write('\t\t'+uncname+'.insert({"'+year+'", { }});\n')
        for bx in range(1, th.GetNbinsX()+1):
            fout.write('\t\t'+sfname +'["'+year+'"].insert({"eta'+str(bx-1)+'", { }});\n')
            fout.write('\t\t'+uncname+'["'+year+'"].insert({"eta'+str(bx-1)+'", { }});\n')
            for by in range(1, th.GetNbinsY()+1):
                tsf  = th.GetBinContent(bx, by)
                tunc = th.GetBinError  (bx, by)
                fout.write('\t\t'+sfname +'["'+year+'"]["eta'+str(bx-1)+'"].insert({"pt'+str(by-1)+'", '+str(tsf) +'});\n')
                fout.write('\t\t'+uncname+'["'+year+'"]["eta'+str(bx-1)+'"].insert({"pt'+str(by-1)+'", '+str(tunc)+'});\n')
        tf.Close()
        fout.write('\n')
    fout.write('\t}\n')
    fout.write('\n')

def print_cc(fout,fname):
    fout.write('#include "%s.h"\n'%fname)
    fout.write('\n')
    fout.write('namespace phoMVAIDWP90SF {')
    fout.write('\n')
    fout.write('\tTString get_phoMVAIDWP90SFPtBin(const float pt) {\n')
    fout.write('\t\tfor ( unsigned int b=nptbins-1; b>=1; b-- ) {\n')
    fout.write('\t\t\tif ( pt > thresholds_pt[b] ) {\n')
    fout.write('\t\t\t\treturn ptbins[b];\n')
    fout.write('\t\t\t}\n')
    fout.write('\t\t}\n')
    fout.write('\t\treturn ptbins[0];\n')
    fout.write('\t}\n')
    fout.write('\n')
    fout.write('\tTString get_phoMVAIDWP90SFEtaBin(const float eta) {\n')
    fout.write('\t\tfor ( unsigned int b=netabins-1; b>=1; b-- ) {\n')
    fout.write('\t\t\tif ( eta > thresholds_eta[b] ) {\n')
    fout.write('\t\t\t\treturn etabins[b];\n')
    fout.write('\t\t\t}\n')
    fout.write('\t\t}\n')
    fout.write('\t\treturn etabins[0];\n')
    fout.write('\t}\n')
    fout.write('\n')
    fout.write('\tfloat get_phoMVAIDWP90SF( const float pt, const float eta, const TString year, const TString variation ) {\n')
    fout.write('\t\tif ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {\n')
    fout.write('\t\t\tstd::cout << "WARNING: unknown year, returning unity b-tagging efficiency!" << std::endl;\n')
    fout.write('\t\t\treturn 1.0;\n')
    fout.write('\t\t}\n')
    fout.write('\t\tif ( pt < 20.0 || fabs(eta) > 2.5 ) return 1.0;\n')
    fout.write('\t\tTString etabin = get_phoMVAIDWP90SFEtaBin(eta);\n')
    fout.write('\t\tTString ptbin = get_phoMVAIDWP90SFPtBin(pt);\n')
    fout.write('\t\tif ( variation == "central" ) return idsf[year][etabin][ptbin];\n')
    fout.write('\t\telse if ( variation == "up" ) return idsf[year][etabin][ptbin]+idsfunc[year][etabin][ptbin];\n')
    fout.write('\t\telse if ( variation == "down" ) return idsf[year][etabin][ptbin]-idsfunc[year][etabin][ptbin];\n')
    fout.write('\t\telse {\n')
    fout.write('\t\t\tstd::cout << "WARNING: unknown variation (central, up, down), returning unity efficiency!" << std::endl;\n')
    fout.write('\t\t\treturn 1.0;\n')
    fout.write('\t\t}\n')
    fout.write('\t}\n')
    fout.write('}\n')

indir = "./data/"
outdir="./NanoCORE/Tools/"
fname = "phoMVAIDWP90SF"

fcc = open(outdir+"/"+fname+".cc",'w')
print_cc(fcc, fname)

fheader = open(outdir+"/"+fname+".h",'w')
print_header(fheader)
print_function_h(fheader)

get_phoMVAIDWP90SF(fheader)

print_footer(fheader)

