#ifndef GenPartSELECTIONS_H
#define GenPartSELECTIONS_H
#include "Nano.h"
#include "Base.h"
#include "TLorentzVector.h"

struct GenPart {
    GenPart(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.GenPart_pt()[idx_];
        mass_ = nt.GenPart_mass()[idx_];
        eta_ = nt.GenPart_eta()[idx_];
        phi_ = nt.GenPart_phi()[idx_];
        pdgId_ = nt.GenPart_pdgId()[idx_];
        p4_.SetPtEtaPhiM(pt_, nt.GenPart_eta()[idx_], nt.GenPart_phi()[idx_], nt.GenPart_mass()[idx_]);
        motherIdx_ = nt.GenPart_genPartIdxMother()[idx_];
        isxyh_ = (nt.GenPart_pdgId()[motherIdx_] == 45 && (abs(pdgId_) == 35 || abs(pdgId_) == 25));
        isygg_ = (nt.GenPart_pdgId()[motherIdx_] == 35 && (abs(pdgId_) == 22 || abs(pdgId_) == 22));
        ishbb_ = (nt.GenPart_pdgId()[motherIdx_] == 25 && (abs(pdgId_) == 5 || abs(pdgId_) == 5));
    }

    int id() { return id_; }
    unsigned int idx() { return idx_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float mass() { return mass_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    int pdgId() { return pdgId_; }
    bool isxyh() { return isxyh_; }
    bool isygg() { return isygg_; }
    bool ishbb() { return ishbb_; }
    int motherIdx() { return motherIdx_; }
    struct GenPart* gen_mother;

  private:
    int id_;
    float pt_ = 0.;
    float eta_ = 0.;
    float mass_ = 0;
    TLorentzVector p4_;
    float phi_ = 0.;
    int pdgId_ = 0;
    unsigned int idx_;
    bool isxyh_ = false;
    bool isygg_ = false;
    bool ishbb_ = false;
    int motherIdx_ = 0;
};

vector<GenPart> getGenParts();
typedef std::vector<GenPart> GenParts;

#endif
