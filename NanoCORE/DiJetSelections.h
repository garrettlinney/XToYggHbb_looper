#ifndef DiJetSELECTIONS_H
#define DiJetSELECTIONS_H
#include "Nano.h"
#include "Base.h"
#include "TLorentzVector.h"
#include "DiPhotonSelections.h"

struct Jet {
    Jet(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.Jet_pt()[idx_];
        eta_ = nt.Jet_eta()[idx_];
        jetId_ = nt.Jet_jetId()[idx_];
        p4_.SetPtEtaPhiM(nt.Jet_pt()[idx_], nt.Jet_eta()[idx_], nt.Jet_phi()[idx_], nt.Jet_mass()[idx_]);
        btagDeepFlavB_ = nt.Jet_btagDeepFlavB()[idx_];
    }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    int id() { return id_; }
    unsigned int idx() { return idx_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    int jetId() { return jetId_; }
    float btagDeepFlavB() { return btagDeepFlavB_; }

  private:
    int id_;
    float pt_ = 0.;
    float eta_ = 0.;
    TLorentzVector p4_;
    float phi_ = 0.;
    int jetId_ = 0;
    unsigned int idx_;
    float btagDeepFlavB_ = 0;
};

vector<Jet> getJets(Photons photons);
typedef std::vector<Jet> Jets;

struct DiJet{
    Jet leadJet;
    Jet subleadJet;
    TLorentzVector p4;
    DiJet(Jet p1, Jet p2)
    {
        leadJet = p1;
        subleadJet = p2;
        TLorentzVector leadJetp4 = TLorentzVector();
        TLorentzVector subleadJetp4 = TLorentzVector();
        leadJetp4=p1.p4();
        subleadJetp4=p2.p4();
        p4 = leadJetp4 + subleadJetp4;
    }
};

typedef std::vector<DiJet> DiJets;

inline bool sortBybscore(Jet &p1, Jet &p2)
{
        return p1.btagDeepFlavB() > p2.btagDeepFlavB();    
}

DiJets DiJetPreselection(Jets &jets); 

#endif
