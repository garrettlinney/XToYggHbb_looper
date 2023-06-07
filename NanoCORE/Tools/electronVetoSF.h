#ifndef ELECTRONVETOSF_H
#define ELECTRONVETOSF_H

#include <map>

using namespace std;

namespace electronVetoSF {
  map<TString,map<TString,map<TString,float>>> idsf = { };
  map<TString,map<TString,map<TString,float>>> idsfunc = { };

  constexpr int netabins = 2;

  float *thresholds_eta = new float[netabins];
  TString *etabins = new TString[netabins];
  inline void set_etabins() {
    thresholds_eta[0] = 0.0;
    etabins[0] = "eta0";
    thresholds_eta[1] = 1.5;
    etabins[1] = "eta1";
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
  inline void set_r9bins(const float abseta) {
    thresholds_r9[0] = 0.0;
    r9bins[0] = "r90";
    thresholds_r9[1] = ( abseta <= 1.5 ? 0.85 : 0.9 );
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

  inline void reset_electronVetoSF() {
    idsf.clear();
    idsfunc.clear();
  }

  inline void set_electronVetoSF() {
    // Copied from flashgg:
    // https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2016_cfi.py#L25-L34 
    idsf.insert({"2016", { }});
    idsfunc.insert({"2016", { }});
      idsf["2016"].insert({"eta0", { }});
      idsfunc["2016"].insert({"eta0", { }});
        idsf["2016"]["eta0"].insert({"r90", 0.9973});
        idsfunc["2016"]["eta0"].insert({"r90", 0.0024});
        idsf["2016"]["eta0"].insert({"r91", 0.9963});
        idsfunc["2016"]["eta0"].insert({"r91", 0.0006});
      idsf["2016"].insert({"eta1", { }});
      idsfunc["2016"].insert({"eta1", { }});
        idsf["2016"]["eta1"].insert({"r90", 0.9700});
        idsfunc["2016"]["eta1"].insert({"r90", 0.0075});
        idsf["2016"]["eta1"].insert({"r91", 0.9909});
        idsfunc["2016"]["eta1"].insert({"r91", 0.0018});

    // Copied from flashgg:
    // https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2017_Legacy_cfi.py#L26-L35
    idsf.insert({"2017", { }});
    idsfunc.insert({"2017", { }});
      idsf["2017"].insert({"eta0", { }});
      idsfunc["2017"].insert({"eta0", { }});
        idsf["2017"]["eta0"].insert({"r90", 0.9838});
        idsfunc["2017"]["eta0"].insert({"r90", 0.0024});
        idsf["2017"]["eta0"].insert({"r91", 0.9913});
        idsfunc["2017"]["eta0"].insert({"r91", 0.0009});
      idsf["2017"].insert({"eta1", { }});
      idsfunc["2017"].insert({"eta1", { }});
        idsf["2017"]["eta1"].insert({"r90", 0.9777});
        idsfunc["2017"]["eta1"].insert({"r90", 0.0180});
        idsf["2017"]["eta1"].insert({"r91", 0.9784});
        idsfunc["2017"]["eta1"].insert({"r91", 0.0026});

    // Copied from flashgg:
    // https://github.com/cms-analysis/flashgg/blob/58c243d9d1f794d7dca8a94fdcd390aed91cb49c/Systematics/python/flashggDiPhotonSystematics2018_Legacy_cfi.py#L25-L34
    idsf.insert({"2018", { }});
    idsfunc.insert({"2018", { }});
      idsf["2018"].insert({"eta0", { }});
      idsfunc["2018"].insert({"eta0", { }});
        idsf["2018"]["eta0"].insert({"r90", 0.9819,});
        idsfunc["2018"]["eta0"].insert({"r90", 0.0020});
        idsf["2018"]["eta0"].insert({"r91", 0.9931});
        idsfunc["2018"]["eta0"].insert({"r91", 0.0005});
      idsf["2018"].insert({"eta1", { }});
      idsfunc["2018"].insert({"eta1", { }});
        idsf["2018"]["eta1"].insert({"r90", 0.9680});
        idsfunc["2018"]["eta1"].insert({"r90", 0.0069});
        idsf["2018"]["eta1"].insert({"r91", 0.9785});
        idsfunc["2018"]["eta1"].insert({"r91", 0.0018});
  }

  inline float get_electronVetoSF(const float eta, const float r9, const TString year, const TString variation ) {
    if ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {
      std::cout << "WARNING: unknown year, returning unity electron veto SF!" << std::endl;
      return 1.0;
    }
    TString actualYear;
    if (year=="2016nonAPV" || year=="2016APV") actualYear = "2016";
    else actualYear = year;

    float abseta = fabs(eta);
    set_etabins();
    set_r9bins(abseta);
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
