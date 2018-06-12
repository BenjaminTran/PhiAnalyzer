#ifndef PHIANALYZER__PHIBDTAPP_H
#define PHIANALYZER__PHIBDTAPP_H

#include <iostream>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include "TMVA/Reader.h"

// user include files

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/phi.h"

class BDTApp : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit BDTApp(const edm::ParameterSet&);
        ~BDTApp();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
        edm::EDGetTokenT<reco::VertexCollection> _vertexCollName;
        edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > _Dedx_Harmonic2;

        TH1D* h_nEvt;
        TH1D* h_BDTresponse;
        TH2D* h_masspt;

        double mass;
        double BDTresponse;
        double pt;

        std::string weightFileName;

        TTree* phiKaonTree;
        //TTree* BDTPhiTree;

        utility::tree_particle phiKaonCandidate;
};

#endif
