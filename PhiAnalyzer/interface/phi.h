#ifndef PHIANALYZER__PHIMESON_H
#define PHIANALYZER__PHIMESON_H

#include "PhiAnalyzer/PhiAnalyzer/interface/particle.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"

class PhiMeson : public Particle{
    private:
        std::vector<kaon> KaonDau_;

    public:
        PhiMeson();
        PhiMeson(double mass, TVector3 momentum, TLorentzVector PtEtaPhiE, bool isMatched = false);

        PhiMeson BuildPhi(kaon Pkp, kaon Pkm, bool isMatched = false);
        std::vector<PhiMeson> EventCombinatorialPhi(std::vector<kaon> PKp_, std::vector<kaon> PKm_);
        void addKaonDau(kaon dau);

        //Getter
        kaon getKaonDau(int dauID);


};

#endif
