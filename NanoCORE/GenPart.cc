#include "GenPart.h"
#include "Config.h"

using namespace tas;

GenParts getGenParts() {
    GenParts genParts;
    for (unsigned int iGenPart = 0; iGenPart < nt.nGenPart(); iGenPart++) {
        GenPart cand_GenPart = GenPart(iGenPart);
        genParts.push_back(cand_GenPart);
    }
    return genParts;
}

