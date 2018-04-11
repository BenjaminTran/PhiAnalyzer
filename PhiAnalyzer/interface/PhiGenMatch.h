#ifndef PHIANALYZER__PHIGENMATCH_H
#define PHIANALYZER__PHIGENMATCH_H

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"

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
        TH1D* h_momid;

        edm::Service<TFileService> fs;

        TTree* Signal;
        TTree* Background;

        struct kaon_particle{
            float momentum;
            float pt;
            float ptError;
            float energy;
            float dedx;
            float charge;
            float dz;
            float dzError;
            float dxy;
            float dxyError;
            float eta;
            float rapidity;
            float phi;
            float vx;
            float vy;
            float vz;
            float px;
            float py;
            float pz;
            float vzFlip;
            float chi2;
            float chi2norm;
            float ndof;
            float nhits;
        };

        kaon_particle sigTrack_particle_;
        kaon_particle bckTrack_particle_;

};

#endif
