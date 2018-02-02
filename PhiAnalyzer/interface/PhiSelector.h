#ifndef PHIANALYZER__PHISELECTOR_H
#define PHIANALYZER__PHISELECTOR_H
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
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PhiSelector : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PhiSelector(const edm::ParameterSet&);
      ~PhiSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      struct kaon{
          kaon(double p_, double dedx_, double charge_) :
              p(p_), dedx(dedx_), charge(charge_) {}

          double p;
          double dedx;
          double charge;
      };
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void DeDxFiller(reco::TrackCollection::const_iterator &track, edm::ValueMap<reco::DeDxData> DeDxTrack, TH2D* dedx_p);
      void getDeDx(reco::TrackCollection::const_iterator &track, edm::ValueMap<reco::DeDxData> DeDxTrack);
      void FillKaonContainer(reco::TrackCollection::const_iterator &track, edm::ValueMap<reco::DeDxData> DeDxTrack, std::vector<kaon> &pkp, std::vector<kaon> &pkm);

      const double kaonMass = 0.493677;

      // ----------member data ---------------------------
       edm::EDGetTokenT<reco::TrackCollection> _trkSrc;
       edm::EDGetTokenT<reco::VertexCollection> _vtxSrc;
       edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > _Dedx_Harmonic2;
       edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > _Dedx_Trunc40;


       //Vectors to hold kaons to perform combinatorial mass reconstruction. Following PEN naming scheme
       std::vector<PhiSelector::kaon> PKp_Harm; //Positive charged
       std::vector<PhiSelector::kaon> PKm_Harm; //Negative charged
       std::vector<PhiSelector::kaon> PKp_Trun; //Positive charged
       std::vector<PhiSelector::kaon> PKm_Trun; //Negative charged

       int multMin_;
       int multMax_;

       /* If you want to create a TTree
        */
       //float H2dedx;
       //float T2dedx;

       TH1D* h_nEvt;
       TH1D* h_mult;
       TH2D* h_Dedx_p_Harm;
       TH2D* h_Dedx_p_Trun;
};

#endif
