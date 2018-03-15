#ifndef UTIL_H
#define UTIL_H
// system include files
#include <memory>
#include <TH1D.h>
#include <TH2D.h>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


namespace utility
{
    extern double const kaonMass;

    struct myVertex{
        double bestvx = -999;
        double bestvy = -999;
        double bestvz = -999;
        double bestvxError = -999;
        double bestvyError = -999;
        double bestvzError = -999;

        const reco::Vertex &vtx;
        math::XYZPoint bestvtx;

        myVertex(double vx, double vy, double vz, double vxErr, double vyErr, double vzErr, const reco::Vertex &vtx_, math::XYZPoint bestvtx_) :
            bestvx(vx), bestvy(vy), bestvz(vz), bestvxError(vxErr), bestvyError(vyErr), bestvzError(vzErr), vtx(vtx_), bestvtx(bestvtx_) {}
    };

    struct track_combo{
        reco::TrackCollection::const_iterator track;
        reco::TrackRef track_ref;

        track_combo(reco::TrackCollection::const_iterator &track_, reco::TrackRef track_ref_) :
            track(track_), track_ref(track_ref_) {}
    };

    myVertex MyVertexBuild(edm::Handle<reco::VertexCollection> vertices);

    bool isTrackGood(reco::TrackCollection::const_iterator &track, myVertex myVtx, bool trackPtCut);

    bool SelectionCut(reco::TrackCollection::const_iterator &track, myVertex myVtx, bool trackPtCut, double dzdca, double dxydca, double eta, double ptCutVal, int nhits);

    int Multiplicity(edm::Handle<reco::TrackCollection> tracks,
            edm::Handle<reco::VertexCollection> vertices);

    bool AcceptTrackDeDx(track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, std::string constraint);

    bool GetCollection(edm::EDGetTokenT<T> const& tag, edm::Handle<T>& result);

    double getDeDx(utility::track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack);
}

#endif
