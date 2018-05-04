#ifndef PHIANALYZER__PHIKAONTREE_H
#define PHIANALYZER__PHIKAONTREE_H

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

class PhiKaonTree : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit PhiKaonTree(const edm::ParameterSet&);
        ~PhiKaonTree();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
        edm::EDGetTokenT<reco::VertexCollection> _vtxSrc;
        edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > _Dedx_Harmonic2;

        TH1D* h_nEvt;

        TTree* phiKaonTree;
};

#endif
