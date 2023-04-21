#include "phoMVAIDWP90SF.h"

namespace phoMVAIDWP90SF {
	TString get_phoMVAIDWP90SFPtBin(const float pt) {
		for ( unsigned int b=nptbins-1; b>=1; b-- ) {
			if ( pt > thresholds_pt[b] ) {
				return ptbins[b];
			}
		}
		return ptbins[0];
	}

	TString get_phoMVAIDWP90SFEtaBin(const float eta) {
		for ( unsigned int b=netabins-1; b>=1; b-- ) {
			if ( eta > thresholds_eta[b] ) {
				return etabins[b];
			}
		}
		return etabins[0];
	}

	float get_phoMVAIDWP90SF( const float pt, const float eta, const TString year, const TString variation ) {
		if ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {
			std::cout << "WARNING: unknown year, returning unity b-tagging efficiency!" << std::endl;
			return 1.0;
		}
		if ( pt < 20.0 || fabs(eta) > 2.5 ) return 1.0;
		TString etabin = get_phoMVAIDWP90SFEtaBin(eta);
		TString ptbin = get_phoMVAIDWP90SFPtBin(pt);
		if ( variation == "central" ) return idsf[year][etabin][ptbin];
		else if ( variation == "up" ) return idsf[year][etabin][ptbin]+idsfunc[year][etabin][ptbin];
		else if ( variation == "down" ) return idsf[year][etabin][ptbin]-idsfunc[year][etabin][ptbin];
		else {
			std::cout << "WARNING: unknown variation (central, up, down), returning unity efficiency!" << std::endl;
			return 1.0;
		}
	}
}
