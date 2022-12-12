#ifndef LeptonSELECTIONS_H
#define LeptonSELECTIONS_H
#include "Nano.h"
#include "Base.h"
#include "TLorentzVector.h"
#include "DiPhotonSelections.h"

struct Electron {
    Electron(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.Electron_pt()[idx_];
        eta_ = nt.Electron_eta()[idx_];
        dxy_ = nt.Electron_dxy()[idx_];
        dz_ = nt.Electron_dz()[idx_];
        mvaFall17V2Iso_WP90_ = nt.Electron_mvaFall17V2Iso_WP90()[idx_];
        phi_ = nt.Electron_phi()[idx_];
        id_ = nt.Electron_pdgId()[idx_];
        p4_.SetPtEtaPhiM(nt.Electron_pt()[idx_], nt.Electron_eta()[idx_], nt.Electron_phi()[idx_], nt.Electron_mass()[idx_]);
    }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    int id() { return id_; }
    unsigned int idx() { return idx_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float eta() { return eta_; }
    float dxy() { return dxy_; }
    float dz() { return dz_; }
    bool mvaFall17V2Iso_WP90() { return mvaFall17V2Iso_WP90_; }
    float phi() { return phi_; }

  private:
    int id_;
    float pt_ = 0.;
    float eta_ = 0.;
    float dxy_ = 0.;
    float dz_ = 0.;
    TLorentzVector p4_;
    bool mvaFall17V2Iso_WP90_ = false;
    float phi_ = 0.;
    unsigned int idx_;
};

vector<Electron> getElectrons(Photons photons);
typedef std::vector<Electron> Electrons;

struct Muon {
    Muon(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.Muon_pt()[idx_];
        eta_ = nt.Muon_eta()[idx_];
        dxy_ = nt.Muon_dxy()[idx_];
        dz_ = nt.Muon_dz()[idx_];
        phi_ = nt.Muon_phi()[idx_];
        id_ = nt.Muon_pdgId()[idx_];
        mediumId_ = nt.Muon_mediumId()[idx_];
        pfRelIso03_all_ = nt.Muon_pfRelIso03_all()[idx_];
        p4_.SetPtEtaPhiM(nt.Muon_pt()[idx_], nt.Muon_eta()[idx_], nt.Muon_phi()[idx_], nt.Muon_mass()[idx_]);
    }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    int id() { return id_; }
    unsigned int idx() { return idx_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float eta() { return eta_; }
    float dxy() { return dxy_; }
    float dz() { return dz_; }
    float phi() { return phi_; }
    bool mediumId() { return mediumId_; }
    float pfRelIso03_all() {return pfRelIso03_all_; }

  private:
    int id_;
    float pt_ = 0.;
    float eta_ = 0.;
    float dxy_ = 0.;
    float dz_ = 0.;
    TLorentzVector p4_;
    float phi_ = 0.;
    bool mediumId_ = false;
    unsigned int idx_;
    float pfRelIso03_all_ = 0.;
};

vector<Muon> getMuons(Photons photons);
typedef std::vector<Muon> Muons;

#endif
