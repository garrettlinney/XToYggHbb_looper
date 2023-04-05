#ifndef LOWMASSHGGPRESELSF_H
#define LOWMASSHGGPRESELSF_H

#include <map>

using namespace std;

namespace lowMassHggPreselSF {
  map<TString,map<TString,map<TString,float>>> idsf = { };
  map<TString,map<TString,map<TString,float>>> idsfunc = { };

  constexpr int netabins = 3;

  float *thresholds_eta = new float[netabins];
  TString *etabins = new TString[netabins];
  inline void set_etabins() {
    thresholds_eta[0] = 0.0;
    etabins[0] = "eta0";
    thresholds_eta[1] = 1.5;
    etabins[1] = "eta1";
    thresholds_eta[2] = 6.0;
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

  inline void reset_lowMassHggPreselSF() {
    idsf.clear();
    idsfunc.clear();
  }

  inline void set_lowMassHggPreselSF() {
    idsf.insert({"2016APV", { }});
    idsfunc.insert({"2016APV", { }});
      idsf["2016APV"].insert({"eta0", { }});
      idsfunc["2016APV"].insert({"eta0", { }});
        idsf["2016APV"]["eta0"].insert({"r90", 0.999});
        idsfunc["2016APV"]["eta0"].insert({"r90", 0.020});
        idsf["2016APV"]["eta0"].insert({"r91", 1.007});
        idsfunc["2016APV"]["eta0"].insert({"r91", 0.023});
      idsf["2016APV"].insert({"eta1", { }});
      idsfunc["2016APV"].insert({"eta1", { }});
        idsf["2016APV"]["eta1"].insert({"r90", 1.000});
        idsfunc["2016APV"]["eta1"].insert({"r90", 0.000});
        idsf["2016APV"]["eta1"].insert({"r91", 1.015});
        idsfunc["2016APV"]["eta1"].insert({"r91", 0.010});
      idsf["2016APV"].insert({"eta2", { }});
      idsfunc["2016APV"].insert({"eta2", { }});
        idsf["2016APV"]["eta2"].insert({"r90", 1.000});
        idsfunc["2016APV"]["eta2"].insert({"r90", 0.000});
        idsf["2016APV"]["eta2"].insert({"r91", 1.000});
        idsfunc["2016APV"]["eta2"].insert({"r91", 0.000});

    idsf.insert({"2016nonAPV", { }});
    idsfunc.insert({"2016nonAPV", { }});
      idsf["2016nonAPV"].insert({"eta0", { }});
      idsfunc["2016nonAPV"].insert({"eta0", { }});
        idsf["2016nonAPV"]["eta0"].insert({"r90", 1.001});
        idsfunc["2016nonAPV"]["eta0"].insert({"r90", 0.020});
        idsf["2016nonAPV"]["eta0"].insert({"r91", 1.006});
        idsfunc["2016nonAPV"]["eta0"].insert({"r91", 0.019});
      idsf["2016nonAPV"].insert({"eta1", { }});
      idsfunc["2016nonAPV"].insert({"eta1", { }});
        idsf["2016nonAPV"]["eta1"].insert({"r90", 1.000});
        idsfunc["2016nonAPV"]["eta1"].insert({"r90", 0.0000});
        idsf["2016nonAPV"]["eta1"].insert({"r91", 0.980});
        idsfunc["2016nonAPV"]["eta1"].insert({"r91", 0.015});
      idsf["2016nonAPV"].insert({"eta2", { }});
      idsfunc["2016nonAPV"].insert({"eta2", { }});
        idsf["2016nonAPV"]["eta2"].insert({"r90", 1.000});
        idsfunc["2016nonAPV"]["eta2"].insert({"r90", 0.000});
        idsf["2016nonAPV"]["eta2"].insert({"r91", 1.000});
        idsfunc["2016nonAPV"]["eta2"].insert({"r91", 0.000});

    idsf.insert({"2017", { }});
    idsfunc.insert({"2017", { }});
      idsf["2017"].insert({"eta0", { }});
      idsfunc["2017"].insert({"eta0", { }});
        idsf["2017"]["eta0"].insert({"r90", 1.007});
        idsfunc["2017"]["eta0"].insert({"r90", 0.019});
        idsf["2017"]["eta0"].insert({"r91", 1.017});
        idsfunc["2017"]["eta0"].insert({"r91", 0.016});
      idsf["2017"].insert({"eta1", { }});
      idsfunc["2017"].insert({"eta1", { }});
        idsf["2017"]["eta1"].insert({"r90", 1.000});
        idsfunc["2017"]["eta1"].insert({"r90", 0.000});
        idsf["2017"]["eta1"].insert({"r91", 1.037});
        idsfunc["2017"]["eta1"].insert({"r91", 0.007});
      idsf["2017"].insert({"eta2", { }});
      idsfunc["2017"].insert({"eta2", { }});
        idsf["2017"]["eta2"].insert({"r90", 1.000});
        idsfunc["2017"]["eta2"].insert({"r90", 0.000});
        idsf["2017"]["eta2"].insert({"r91", 1.000});
        idsfunc["2017"]["eta2"].insert({"r91", 0.000});

    idsf.insert({"2018", { }});
    idsfunc.insert({"2018", { }});
      idsf["2018"].insert({"eta0", { }});
      idsfunc["2018"].insert({"eta0", { }});
        idsf["2018"]["eta0"].insert({"r90", 0.983});
        idsfunc["2018"]["eta0"].insert({"r90", 0.025});
        idsf["2018"]["eta0"].insert({"r91", 0.998});
        idsfunc["2018"]["eta0"].insert({"r91", 0.024});
      idsf["2018"].insert({"eta1", { }});
      idsfunc["2018"].insert({"eta1", { }});
        idsf["2018"]["eta1"].insert({"r90", 1.000});
        idsfunc["2018"]["eta1"].insert({"r90", 0.000});
        idsf["2018"]["eta1"].insert({"r91", 0.957});
        idsfunc["2018"]["eta1"].insert({"r91", 0.015});
      idsf["2018"].insert({"eta2", { }});
      idsfunc["2018"].insert({"eta2", { }});
        idsf["2018"]["eta2"].insert({"r90", 1.000});
        idsfunc["2018"]["eta2"].insert({"r90", 0.000});
        idsf["2018"]["eta2"].insert({"r91", 1.000});
        idsfunc["2018"]["eta2"].insert({"r91", 0.000});
  }

  inline float get_lowMassHggPreselSF(const float eta, const float r9, const TString year, const TString variation ) {
    if ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {
      std::cout << "WARNING: unknown year, returning unity electron veto SF!" << std::endl;
      return 1.0;
    }

    float abseta = fabs(eta);
    set_etabins();
    set_r9bins(abseta);
    TString etabin = get_etaBin(abseta);
    TString r9bin = get_r9Bin(r9);

    if ( variation == "central" ) return idsf[year][etabin][r9bin];
    else if ( variation == "up" ) return idsf[year][etabin][r9bin]+idsfunc[year][etabin][r9bin];
    else if ( variation == "down" ) return idsf[year][etabin][r9bin]-idsfunc[year][etabin][r9bin];
    else {
      std::cout << "WARNING: unknown variation (central, up, down), returning unity scale factor!" << std::endl;
      return 1.0;
    }
  }
}

#endif
