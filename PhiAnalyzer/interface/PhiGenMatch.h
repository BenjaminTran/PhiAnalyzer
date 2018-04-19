#ifndef PHIANALYZER__PHIGENMATCH_H
#define PHIANALYZER__PHIGENMATCH_H

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/phi.h"

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
        edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > _Dedx_Harmonic2;

        struct tree_particle{
            float mass;
            float momentum_1;
            float pt_1;
            float ptError_1;
            float energy_1;
            float dedx_1;
            float charge_1;
            float dz_1;
            float dzError_1;
            float dxy_1;
            float dxyError_1;
            float eta_1;
            float rapidity_1;
            float phi_1;
            float vx_1;
            float vy_1;
            float vz_1;
            float px_1;
            float py_1;
            float pz_1;
            float vzFlip_1;
            float chi2_1;
            float chi2norm_1;
            float ndof_1;
            float nhits_1;
            float momentum_2;
            float pt_2;
            float ptError_2;
            float energy_2;
            float dedx_2;
            float charge_2;
            float dz_2;
            float dzError_2;
            float dxy_2;
            float dxyError_2;
            float eta_2;
            float rapidity_2;
            float phi_2;
            float vx_2;
            float vy_2;
            float vz_2;
            float px_2;
            float py_2;
            float pz_2;
            float vzFlip_2;
            float chi2_2;
            float chi2norm_2;
            float ndof_2;
            float nhits_2;
        };

        TH1D* h_nEvt;
        TH1D* h_phi_yield_rap_1;
        TH1D* h_phi_yield_norap;
        TH1D* h_phid1_mass;
        TH1D* h_phid2_mass;
        TH1D* h_momid;

        edm::Service<TFileService> fs;

        TTree* Signal;
        TTree* Background;

        tree_particle signalStruct;
        tree_particle backgroundStruct;

        std::string dedxConstraint_;

        void FillTreeStruct(tree_particle& treeStruct, PhiMeson phi);

};

#endif
