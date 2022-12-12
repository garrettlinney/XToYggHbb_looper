#include "DiJetSelections.h"
#include "Config.h"
#include "DiPhotonSelections.h"
// too complicated, haven't used that two lepton selections yet
//#include "ElectronSelections.h"
//#include "MuonSelections.h"

using namespace tas;

Jets getJets(Photons photons) {
    Jets jets;
    for (unsigned int ijet = 0; ijet < nt.nJet(); ijet++) {
        Jet cand_jet = Jet(ijet);
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
    if (jets.size()>1){
        DiJet dijet = DiJet(jets[0], jets[1]);
        if (dijet.p4.M() > 50) dijets.push_back(dijet);
    }
    // thre could be either 1 or 0 dijet in the output
    return dijets;    
}
