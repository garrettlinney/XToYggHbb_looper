#ifndef MATERIALUNC_H
#define MATERIALUNC_H

#include <map>

using namespace std;

namespace materialUnc {
  map<TString,map<TString,map<TString,float>>> idsf = { };
  map<TString,map<TString,map<TString,float>>> idsfunc = { };

  constexpr int netabins = 3;

  float *thresholds_eta = new float[netabins];
  TString *etabins = new TString[netabins];
  inline void set_etabins() {
    thresholds_eta[0] = 0.0;
    etabins[0] = "eta0";
    thresholds_eta[1] = 1.0;
    etabins[1] = "eta1";
    thresholds_eta[2] = 1.5;
    etabins[2] = "eta2";
  }

  TString get_etaBin(const float abseta) {
    for ( unsigned int b=netabins-1; b>=1; b-- ) {
      if ( abseta > thresholds_eta[b] ) {
        return etabins[b];
      }
    }
    return etabins[0];
  }

  constexpr int nr9bins = 2;

  float *thresholds_r9 = new float[nr9bins];
  TString *r9bins = new TString[nr9bins];
  inline void set_r9bins() {
    thresholds_r9[0] = 0.0;
    r9bins[0] = "r90";
    thresholds_r9[1] = 0.94;
    r9bins[1] = "r91";
  }

  TString get_r9Bin(const float r9) {
    for ( unsigned int b=nr9bins-1; b>=1; b-- ) {
      if ( r9 > thresholds_r9[b] ) {
        return r9bins[b];
      }
    }
    return r9bins[0];
  }

  inline void reset_materialUnc() {
    idsf.clear();
    idsfunc.clear();
  }

  inline void set_materialUnc() {
    set_etabins();
    set_r9bins();

    // Copied from flashgg:
    // https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2016_cfi.py#L369-L380
    idsf.insert({"2016", { }});
    idsfunc.insert({"2016", { }});
      idsf["2016"].insert({"eta0", { }});
      idsfunc["2016"].insert({"eta0", { }});
        idsf["2016"]["eta0"].insert({"r90", 1.0});
        idsfunc["2016"]["eta0"].insert({"r90", 0.000455});
        idsf["2016"]["eta0"].insert({"r91", 1.0});
        idsfunc["2016"]["eta0"].insert({"r91", 0.000233});
      idsf["2016"].insert({"eta1", { }});
      idsfunc["2016"].insert({"eta1", { }});
        idsf["2016"]["eta1"].insert({"r90", 1.0});
        idsfunc["2016"]["eta1"].insert({"r90", 0.002089});
        idsf["2016"]["eta1"].insert({"r91", 1.0});
        idsfunc["2016"]["eta1"].insert({"r91", 0.002089});
      idsf["2016"].insert({"eta2", { }});
      idsfunc["2016"].insert({"eta2", { }});
        idsf["2016"]["eta2"].insert({"r90", 1.0});
        idsfunc["2016"]["eta2"].insert({"r90", 0.001090});
        idsf["2016"]["eta2"].insert({"r91", 1.0});
        idsfunc["2016"]["eta2"].insert({"r91", 0.002377});

    // Copied from flashgg:
    // https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L377-L399
    idsf.insert({"2017", { }});
    idsfunc.insert({"2017", { }});
      idsf["2017"].insert({"eta0", { }});
      idsfunc["2017"].insert({"eta0", { }});
        idsf["2017"]["eta0"].insert({"r90", 1.0});
        idsfunc["2017"]["eta0"].insert({"r90", 0.000455});
        idsf["2017"]["eta0"].insert({"r91", 1.0});
        idsfunc["2017"]["eta0"].insert({"r91", 0.000233});
      idsf["2017"].insert({"eta1", { }});
      idsfunc["2017"].insert({"eta1", { }});
        idsf["2017"]["eta1"].insert({"r90", 1.0});
        idsfunc["2017"]["eta1"].insert({"r90", 0.002089});
        idsf["2017"]["eta1"].insert({"r91", 1.0});
        idsfunc["2017"]["eta1"].insert({"r91", 0.002089});
      idsf["2017"].insert({"eta2", { }});
      idsfunc["2017"].insert({"eta2", { }});
        idsf["2017"]["eta2"].insert({"r90", 1.0});
        idsfunc["2017"]["eta2"].insert({"r90", 0.001090});
        idsf["2017"]["eta2"].insert({"r91", 1.0});
        idsfunc["2017"]["eta2"].insert({"r91", 0.002377});

    // Copied from flashgg:
    // https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2018_Legacy_cfi.py#L360-L382
    idsf.insert({"2018", { }});
    idsfunc.insert({"2018", { }});
      idsf["2018"].insert({"eta0", { }});
      idsfunc["2018"].insert({"eta0", { }});
        idsf["2018"]["eta0"].insert({"r90", 1.0});
        idsfunc["2018"]["eta0"].insert({"r90", 0.000455});
        idsf["2018"]["eta0"].insert({"r91", 1.0});
        idsfunc["2018"]["eta0"].insert({"r91", 0.000233});
      idsf["2018"].insert({"eta1", { }});
      idsfunc["2018"].insert({"eta1", { }});
        idsf["2018"]["eta1"].insert({"r90", 1.0});
        idsfunc["2018"]["eta1"].insert({"r90", 0.002089});
        idsf["2018"]["eta1"].insert({"r91", 1.0});
        idsfunc["2018"]["eta1"].insert({"r91", 0.002089});
      idsf["2018"].insert({"eta2", { }});
      idsfunc["2018"].insert({"eta2", { }});
        idsf["2018"]["eta2"].insert({"r90", 1.0});
        idsfunc["2018"]["eta2"].insert({"r90", 0.001090});
        idsf["2018"]["eta2"].insert({"r91", 1.0});
        idsfunc["2018"]["eta2"].insert({"r91", 0.002377});
  }

  inline float get_materialUnc(const float eta, const float r9, const TString year, const TString variation ) {
    if ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {
      std::cout << "WARNING: unknown year, returning unity electron veto SF!" << std::endl;
      return 1.0;
    }
    TString actualYear;
    if (year=="2016nonAPV" || year=="2016APV") actualYear = "2016";
    else actualYear = year;

    float abseta = fabs(eta);
    TString etabin = get_etaBin(abseta);
    TString r9bin = get_r9Bin(r9);

    if ( variation == "central" ) return idsf[actualYear][etabin][r9bin];
    else if ( variation == "up" ) return idsf[actualYear][etabin][r9bin]+idsfunc[actualYear][etabin][r9bin];
    else if ( variation == "down" ) return idsf[actualYear][etabin][r9bin]-idsfunc[actualYear][etabin][r9bin];
    else {
      std::cout << "WARNING: unknown variation (central, up, down), returning unity scale factor!" << std::endl;
      return 1.0;
    }
  }
}

#endif
