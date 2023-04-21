#include "LeptonSelections.h"

using namespace tas;

Electrons getElectrons(Photons photons) {
    Electrons electrons;
    for (unsigned int ielec = 0; ielec < nt.nElectron(); ielec++) {
        Electron elec = Electron(ielec);
        if ( !(abs(elec.eta()) < 2.5) ) continue;
        if ( !(abs(elec.eta()) < 1.4442 || abs(elec.eta()) > 1.566) ) continue;
        if ( !(elec.pt() > 10) ) continue;
        if ( !(abs(elec.dxy()) < 0.045) ) continue;
        if ( !(abs(elec.dz()) < 0.2) ) continue;
        if ( !(elec.mvaFall17V2Iso_WP90() == true )) continue;
        bool clean_with_photon = true;
        for (unsigned int iphoton = 0; iphoton < photons.size(); iphoton++){
            if (elec.p4().DeltaR(photons.at(iphoton).p4())<0.2) clean_with_photon = false;
        }
        if (clean_with_photon == false) continue;
        electrons.push_back(elec);
    }

    return electrons;
}

Muons getMuons(Photons photons) {
    Muons muons;
    for (unsigned int imu = 0; imu < nt.nMuon(); imu++) {
        Muon mu = Muon(imu);
        if ( !(abs(mu.eta()) < 2.4) ) continue;
        if ( !(mu.pt() > 15) ) continue;
        if ( !(abs(mu.dxy()) < 0.045) ) continue;
        if ( !(abs(mu.dz()) < 0.2) ) continue;
        if ( !(mu.mediumId()) ) continue;
        if ( !(mu.pfRelIso03_all()<0.3) ) continue;
        if ( !(mu.isGlobal())) continue;
        bool clean_with_photon = true;
        for (unsigned int iphoton = 0; iphoton < photons.size(); iphoton++){
            if (mu.p4().DeltaR(photons.at(iphoton).p4())<0.2) clean_with_photon = false;
        }
        if (clean_with_photon == false) continue;
        muons.push_back(mu);
    }

    return muons;
}
