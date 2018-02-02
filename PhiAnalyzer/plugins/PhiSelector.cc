// -*- C++ -*-
//
// Package:    PhiSelector/PhiSelector
// Class:      PhiSelector
// 
/**\class PhiSelector PhiSelector.cc PhiSelector/PhiSelector/plugins/PhiSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Tran
//         Created:  Tue, 30 Jan 2018 19:09:22 GMT
//
//

#include "/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/PhiAnalyzer/PhiAnalyzer/interface/PhiSelector.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"


//
// constructors and destructor
//
PhiSelector::PhiSelector(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

    multMin_ = iConfig.getUntrackedParameter<int>("multMin");
    multMax_ = iConfig.getUntrackedParameter<int>("multMax");
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vtxSrc = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc"));
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("Dedx_Harmonic2"));
    _Dedx_Trunc40 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("Dedx_Trunc40"));

}


PhiSelector::~PhiSelector()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void
PhiSelector::DeDxFiller(reco::TrackCollection::const_iterator &track, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, TH2D* dedx_p)
{
        double dedx     = -999.9;
        double momentum = track->p();
        dedx = getDeDx(track,DeDxTrack);
        dedx_p->Fill(dedx,momentum);
}

double
PhiSelector::getDeDx(reco::TrackCollection::const_iterator &track, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack)
{
    double dedx_ = -999;
    if(DeDxTrack.isValid())
    {
        const edm::ValueMap<reco::DeDxData> dedxTrack = *DeDxTrack.product();
        dedx_ = dedxTrack[track].dEdx();
    }

    return dedx_;
}

void
PhiSelector::FillKaonContainer(reco::TrackCollection::const_iterator &track, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, std::vector<kaon> &pkp, std::vector<kaon> &pkm)
{
        //positive kaons
        if(track->charge() == 1)
        {
            kaon pk(track->p(), getDeDx(track,DeDxTrack), track->charge());
            pkp.push_back(pk);
        }

        //negative kaons
        if(track->charge() == -1)
        {
            kaon pk(track->p(), getDeDx(track,DeDxTrack), track->charge());
            pkm.push_back(pk);
        }
}


// ------------ method called for each event  ------------
void
PhiSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   h_nEvt->Fill(1);

   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByToken(_trkSrc,tracks);

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(_vtxSrc,vertices);

   edm::Handle<edm::ValueMap<reco::DeDxData> > DeDx_Harm;
   iEvent.getByToken(_Dedx_Harmonic2,DeDx_Harm);

   edm::Handle<edm::ValueMap<reco::DeDxData> > DeDx_Trun;
   iEvent.getByToken(_Dedx_Trunc40,DeDx_Trun);

   int mult = utility::trackFilter(tracks,vertices);
   h_mult->Fill(mult);

   for(reco::TrackCollection::const_iterator it = tracks.begin();
           it != tracks.end();
           ++it)
   {
       DeDxFiller(it,DeDx_Harm,h_Dedx_p_Harm);
       DeDxFiller(it,DeDx_Trun,h_Dedx_p_Trun);
       FillKaonContainer(it,DeDx_Harm,PKp_Harm,PKm_Harm);
       FillKaonContainer(it,DeDx_Trun,PKp_Trun,PKm_Trun);
   }
}


// ------------ method called once each job just before starting event loop  ------------
void
PhiSelector::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    h_nEvt = fs->make<TH1D>("nEvt","",10,0,10);
    h_mult = fs->make<TH1D>("mult","",400,0,400);
    h_Dedx_p_Harm = fs->make<TH2D>("Dedx_harm","",200,0,20,1500,0,15);
    h_Dedx_p_Trun = fs->make<TH2D>("Dedx_Trun","",200,0,20,1500,0,15);
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhiSelector::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhiSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
  desc.addUntracked<edm::InputTag>("vtxSrc",edm::InputTag("offlinePrimaryVertices"));
  desc.addUntracked<edm::InputTag>("Dedx_Harmonic2",edm::InputTag("dedxHarmonic2"));
  desc.addUntracked<edm::InputTag>("Dedx_Trunc40",edm::InputTag("dedxTruncated40"));
  desc.addUntracked<int>("multMin",0);
  desc.addUntracked<int>("multMax",999);
  descriptions.add("PhiSelector",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiSelector);
