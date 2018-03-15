#ifndef PHIANALYZER__PHIGENMATCH_H
#define PHIANALYZER__PHIGENMATCH_H

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

class PhiGenMatch : public edm::EDAnalyzer {
    public:
        explicit PhiGenMatch(const edm::ParameterSet&);
        ~PhiGenMatch();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::GenParticleCollection> _genCollection;

        TH1D* h_nEvt;
        TH1D* h_phi_yield_rap_1;
        TH1D* h_phi_yield_norap;

        edm::Service<TFileService> fs;


};

#endif
