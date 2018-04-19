#ifndef PHIANALYZER__PHIMESON_H
#define PHIANALYZER__PHIMESON_H
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"

class PhiMeson{
    private:
        double mass_;
        double eta_;
        double rapidity_;
        double phi_;
        double pt_;
        double p_;
        bool isMatched_;
        std::vector<kaon> KaonDau_;

    public:
        PhiMeson();
        PhiMeson(double mass, double pt, double eta, double phi, double rapidity, double p = 0, bool isMatched = false);
        PhiMeson BuildPhi(kaon Pkp, kaon Pkm, bool isMatched = false);
        std::vector<PhiMeson> EventCombinatorialPhi(std::vector<kaon> PKp_, std::vector<kaon> PKm_);
        double getMass();
        double getEta();
        double getRapidity();
        double getPhi();
        double getPt();
        double getP();
        bool getIsGen();
        kaon getKaonDau();

        void setMass(double mass);
        void setEta(double eta);
        void setRapidity(double rapidity);
        void setPhi(double phi);
        void setPt(double pt);
        void setP(double p);
        void setIsGen(bool isMatched);

        void addKaonDau(kaon dau);

};

#endif
