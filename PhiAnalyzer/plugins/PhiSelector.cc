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

using namespace std;

//
// constructors and destructor
//
PhiSelector::PhiSelector(const edm::ParameterSet& iConfig)
{
    multMin_ = iConfig.getUntrackedParameter<int>("multMin");
    multMax_ = iConfig.getUntrackedParameter<int>("multMax");
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vtxSrc = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc"));
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
}


PhiSelector::~PhiSelector()
{
}


//
// member functions
//

void
PhiSelector::DeDxFiller(track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, TH2D* dedx_p)
{
        double dedx     = -999.9;
        double momentum = track_combo_.track->p();
        dedx = getDeDx(track_combo_,DeDxTrack);
        dedx_p->Fill(momentum,dedx);
}

double
PhiSelector::getDeDx(track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack)
{
    double dedx_ = 1;
    const edm::ValueMap<reco::DeDxData> dedxTrack = *DeDxTrack.product();
    dedx_ = dedxTrack[track_combo_.track_ref].dEdx();

    return dedx_;
}

void
PhiSelector::FillKaonContainer(track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, std::vector<kaon> &pkp, std::vector<kaon> &pkm)
{
    double energy = sqrt(TMath::Power(kaonMass,2) + TMath::Power(track_combo_.track->p(),2));
    kaon pk(track_combo_.track->p(), getDeDx(track_combo_,DeDxTrack), energy, track_combo_.track->charge());

    //positive kaons
    if(track_combo_.track->charge() == 1)
        pkp.push_back(pk);

    //negative kaons
    if(track_combo_.track->charge() == -1)
        pkm.push_back(pk);
}

void
PhiSelector::CombinatorialMass(std::vector<PhiSelector::kaon> PKp, std::vector<PhiSelector::kaon> PKm, TH1D* h_mass_)
{
    for(PhiSelector::kaon Pkp : PKp)
    {
        for(PhiSelector::kaon Pkm : PKm)
        {
            double mass = sqrt(2*TMath::Power(kaonMass,2) + 2*(Pkp.energy*Pkm.energy - Pkp.p*Pkm.p));
            h_mass_->Fill(mass);
        }
    }
}


// ------------ method called for each event  ------------
void
PhiSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Vectors to hold kaons to perform combinatorial mass reconstruction. Following PEN naming scheme
   std::vector<PhiSelector::kaon> PKp_Harm; //Positive charged
   std::vector<PhiSelector::kaon> PKm_Harm; //Negative charged

   h_nEvt->Fill(1);

   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByToken(_trkSrc,tracks);

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(_vtxSrc,vertices);

   edm::Handle<edm::ValueMap<reco::DeDxData> > DeDx_Harm;
   iEvent.getByToken(_Dedx_Harmonic2,DeDx_Harm);
   if(!DeDx_Harm.isValid())
   {
       cout << "Bad DeDx collection" << endl;
       return;
   }

   // Multiplicity selection
   int mult = utility::trackFilter(tracks,vertices);
   h_mult->Fill(mult);

   for(reco::TrackCollection::const_iterator it = tracks->begin();
           it != tracks->end();
           ++it)
   {
       //Use only high purity tracks
       if(!it->quality(reco::TrackBase::highPurity)) continue;

       reco::TrackRef track_ref = reco::TrackRef(tracks,it - tracks->begin());
       track_combo track_bundle(it,track_ref);

       //Fill in the dEdx histograms
       DeDxFiller(track_bundle,DeDx_Harm,h_Dedx_p_Harm);

       // Make the vector of kaons to calculate invariant mass at the end
       FillKaonContainer(track_bundle,DeDx_Harm,PKp_Harm,PKm_Harm);
   }
   CombinatorialMass(PKp_Harm,PKm_Harm,h_mass_Harm);
}


// ------------ method called once each job just before starting event loop  ------------
void
PhiSelector::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    h_nEvt = fs->make<TH1D>("nEvt","",10,0,10);
    h_mult = fs->make<TH1D>("mult","",400,0,400);
    h_mass_Harm = fs->make<TH1D>("mass_harm","",300,0.945,1.095);
    h_Dedx_p_Harm = fs->make<TH2D>("Dedx_harm","",200,0,20,1500,0,15);
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
  desc.addUntracked<int>("multMin",0);
  desc.addUntracked<int>("multMax",999);
  descriptions.add("PhiSelector",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiSelector);
