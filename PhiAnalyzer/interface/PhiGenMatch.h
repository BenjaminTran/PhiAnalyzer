#ifndef PHIANALYZER__PHIGENMATCH_H
#define PHIANALYZER__PHIGENMATCH_H

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "DataFormats/Candidate/interface/Candidate.h"

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
        edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
        edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;

        TH1D* h_nEvt;
        TH1D* h_phi_yield_rap_1;
        TH1D* h_phi_yield_norap;
        TH1D* h_phid1_mass;
        TH1D* h_phid2_mass;

        edm::Service<TFileService> fs;


};

#endif
