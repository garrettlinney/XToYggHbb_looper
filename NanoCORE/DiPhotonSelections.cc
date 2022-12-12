#include "DiPhotonSelections.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"

using namespace tas;

bool UseLowR9Photon(Photon pho, bool isEB) {
    bool useThisPhoton = false;
    if (isEB) {
        if ( !(pho.sieie() < 0.015) ) return useThisPhoton;       
        if ( !(pho.trkIso() < 6.0) ) return useThisPhoton;       
        if ( !(pho.phoIso() - 0.16544*pho.perEvtRho() < 4.0) ) return useThisPhoton;       
    } else {
        if ( !(pho.sieie() < 0.035) ) return useThisPhoton;       
        if ( !(pho.trkIso() < 6.0) ) return useThisPhoton;       
        if ( !(pho.phoIso() - 0.13212*pho.perEvtRho() < 4.0) ) return useThisPhoton;       

    }
    // 0.16544 and 0.13212 are copied from flashggPreselectedDiPhotons_cfi.py
    useThisPhoton = true;
    return useThisPhoton;
}

Photons getPhotons() {
    Photons photons;
    for (unsigned int ipho = 0; ipho < nt.nPhoton(); ipho++) {
        Photon pho = Photon(ipho);
        //cout << "eta: " << pho.eta() << ", r9: " << pho.r9() << endl;
        if (pho.pt()<18) continue;
        if ( !(abs(pho.eta()<2.5)) ) continue;
        if ( !(abs(pho.eta())<1.4442 || abs(pho.eta())>1.566) ) continue;
        if ( !(pho.hoe()<0.08) ) continue;
        if (pho.pixelSeed() > 0.5) continue; // this is not standard photon selections, but sam used this to suppress electrons on the DY peak
        if (pho.eveto() < 0.5) continue;
        if ( !(pho.egPhoId() > -0.7) ) continue;
        bool pho_EB_highR9 = abs(pho.eta()) < 1.5 && pho.r9() > 0.85; 
        bool pho_EE_highR9 = abs(pho.eta()) > 1.5 && pho.r9() > 0.90; 
        bool pho_EB_lowR9 = abs(pho.eta()) < 1.5 && pho.r9() < 0.85 && pho.r9() > 0.50 && UseLowR9Photon(pho, true);
        bool pho_EE_lowR9 = abs(pho.eta()) > 1.5 && pho.r9() < 0.90 && pho.r9() > 0.80 && UseLowR9Photon(pho, false);
        if ( !(pho_EB_highR9 || pho_EE_highR9 || pho_EB_lowR9 || pho_EE_lowR9) ) continue;
        photons.push_back(pho);
    }

    sort(photons.begin(), photons.end(), sortByPt);

    return photons;
}

DiPhotons DiPhotonPreselection(Photons &photons/*, bool verbose=false*/) {
    DiPhotons diphotons; 
    float maxSumDiphoPt = 0;    
    for (unsigned int i1 = 0; i1 < photons.size(); i1++) {
        for (unsigned int i2 = i1+1; i2 < photons.size(); i2++) {
            Photon pho1 = photons[i1];
            Photon pho2 = photons[i2];
            DiPhoton dipho = DiPhoton(pho1, pho2);
            if ( !(dipho.leadPho.pt() > 30.0 && dipho.subleadPho.pt() > 18.0) ) continue;
            //if ( !((dipho.leadPho.r9() > 0.8 && dipho.leadPho.chargedHadIso() < 20) || dipho.leadPho.chargedHadIso()/dipho.leadPho.pt() < 0.3) ) continue;
            //if ( !((dipho.subleadPho.r9() > 0.8 && dipho.subleadPho.chargedHadIso() < 20) || dipho.subleadPho.chargedHadIso()/dipho.subleadPho.pt() < 0.3) ) continue;
            if ( !(dipho.leadPho.r9() > 0.8 || dipho.leadPho.chargedHadIso() < 20 || dipho.leadPho.chargedHadIso()/dipho.leadPho.pt() < 0.3) ) continue;
            if ( !(dipho.subleadPho.r9() > 0.8 || dipho.subleadPho.chargedHadIso() < 20 || dipho.subleadPho.chargedHadIso()/dipho.subleadPho.pt() < 0.3) ) continue;

            
            if ( !(dipho.leadPho.pt()/dipho.p4.M() > 0.33 && dipho.subleadPho.pt()/dipho.p4.M() > 0.25) ) continue;
            if (dipho.p4.M() < 55 || dipho.p4.M() > 999999) continue;

            float sumDiPhoPt = dipho.leadPho.pt() + dipho.subleadPho.pt();

            if ( diphotons.size() == 0 )
            {
                diphotons.push_back(dipho);
                maxSumDiphoPt = sumDiPhoPt;
            } 
            else if (sumDiPhoPt > maxSumDiphoPt) 
            {
                diphotons.pop_back();
                diphotons.push_back(dipho);
                maxSumDiphoPt = sumDiPhoPt;
            }

        }
    } 

    // thre could be either 1 or 0 diphoton in the output
    return diphotons;    
}
