#include "DiJetSelections.h"

using namespace tas;

Jets getJets(Photons photons, const int JESUnc, const int JERUnc) {
    Jets jets;
    for (unsigned int ijet = 0; ijet < nt.nJet(); ijet++) {
        Jet cand_jet = Jet(ijet);

        // Apply JES unc
        if ( JESUnc!=0 ) {
          // JESUnc is only an uncertainty => JESUnc==1 does not change the nominal value
          if ( JESUnc==2  ) cand_jet.setPt( cand_jet.pt_jesTotalUp() );
          if ( JESUnc==-2 ) cand_jet.setPt( cand_jet.pt_jesTotalDown() );
        }
        // Apply JER unc
        if ( JERUnc!=0 ) {
          // JERUnc is only an uncertainty => JERUnc==1 does not change the nominal value
          if ( JERUnc==2  ) cand_jet.setPt( cand_jet.pt_jerUp() );
          if ( JERUnc==-2 ) cand_jet.setPt( cand_jet.pt_jerDown() );
        }

        if ( !(abs(cand_jet.eta()) < 2.4) ) continue;
        if ( !(cand_jet.pt() > 25) ) continue;
        if ( !(cand_jet.jetId()>=1) ) continue;
        bool clean_with_photon = true;
        for (unsigned int iphoton = 0; iphoton < photons.size(); iphoton++){
            if (cand_jet.p4().DeltaR(photons.at(iphoton).p4())<0.4) clean_with_photon = false;
        }
        if (clean_with_photon == false) continue;
        jets.push_back(cand_jet);
    }
    sort(jets.begin(), jets.end(), sortBybscore);

    return jets;
}

DiJets DiJetPreselection(Jets &jets) {
    DiJets dijets;   
    if (jets.size() > 1) {
        DiJet dijet = DiJet(jets[0], jets[1]);
        dijets.push_back(dijet);
    }
    return dijets;    
}
