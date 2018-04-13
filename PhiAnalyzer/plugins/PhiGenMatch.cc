#include "PhiAnalyzer/PhiAnalyzer/interface/PhiGenMatch.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/phi.h"

using namespace std;

//Constructors and destructor
PhiGenMatch::PhiGenMatch(const edm::ParameterSet& iConfig)
{
    _genCollection = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genCollection"));
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexCollName"));
    dedxConstraint_ = iConfig.getUntrackedParameter<string>("dedxConstraint");
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
}

PhiGenMatch::~PhiGenMatch()
{}

// Member functions

// ------------ method called for each event  ------------
void
PhiGenMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    h_nEvt->Fill(1);

    std::vector<kaon> genDauKaons;

    edm::Handle<reco::GenParticleCollection> gencand;
    edm::Handle<reco::TrackCollection> trkSrc;
    edm::Handle<reco::VertexCollection> vertices;
    edm::Handle<edm::ValueMap<reco::DeDxData> > DeDx_Harm;
    try{
        utility::GetCollection<reco::GenParticleCollection>(iEvent,_genCollection,gencand);
        utility::GetCollection<reco::TrackCollection>(iEvent,_trkSrc,trkSrc);
        utility::GetCollection<reco::VertexCollection>(iEvent,_vertexCollName,vertices);
        utility::GetCollection<edm::ValueMap<reco::DeDxData> >(iEvent,_Dedx_Harmonic2,DeDx_Harm);
    }
    catch(const std::invalid_argument& e){
        std::cerr << e.what();
        return;
    }

    utility::myVertex vertex = utility::MyVertexBuild(vertices);


    for(reco::GenParticleCollection::const_iterator gncand = gencand->begin();
            gncand != gencand->end();
            ++gncand)
    {
        int id = gncand->pdgId();
        int mid = 0;
        double rap = gncand->rapidity();
        double eta = gncand->eta();

        if(gncand->numberOfDaughters() > 2) continue;

        if(gncand->numberOfMothers() == 1)
        {
            const reco::Candidate *mom = gncand->mother();
            mid = mom->pdgId();
            //if(mom->numberOfMothers() == 1)
            //{
                //const reco::Candidate * mom1 = mom->mother();
                //mid = mom1->pdgId();
            //}
            h_momid->Fill(mid);
        }

        if(fabs(id) == 333 && fabs(rap) < 1.0)
        {
            h_phi_yield_rap_1->Fill(gncand->mass());
        }

        if(fabs(id) == 333 && fabs(eta) < 2.4)
        {
            h_phi_yield_norap->Fill(gncand->mass());
            const reco::Candidate *d1 = gncand->daughter(0);
            const reco::Candidate *d2 = gncand->daughter(1);

            if(fabs(d1->pdgId()) != 321 && fabs(d2->pdgId()) != 321) continue;

            h_phid1_mass->Fill(d1->mass());
            h_phid2_mass->Fill(d2->mass());

            genDauKaons.push_back(kaon(TVector3(d1->px(), d1->py(), d1->pz()), d1->eta(), d1->phi(), d1->charge(), true));
            genDauKaons.push_back(kaon(TVector3(d2->px(), d2->py(), d2->pz()), d2->eta(), d2->phi(), d2->charge(), true));
        }
    }

    //Perform Matching
    for(reco::TrackCollection::const_iterator trk = trkSrc->begin();
            trk != trkSrc->end();
            ++trk)
    {
        //if(!utility::SelectionCut(trk,vertex,false,1.0,1.0,2.4,0,5))
        reco::TrackRef track_ref = reco::TrackRef(trkSrc,trk - trkSrc->begin());
        utility::track_combo track_bundle(trk, track_ref);
        kaon K(TVector3(trk->px(), trk->py(), trk->pz()), trk->eta(), trk->phi(), trk->charge(), false);

        for(kaon genK : genDauKaons)
        {
            sigTrack_particle_.momentum = trk->p();
            sigTrack_particle_.px       = trk->px();
            sigTrack_particle_.py       = trk->py();
            sigTrack_particle_.pz       = trk->pz();
            sigTrack_particle_.pt       = trk->pt();
            sigTrack_particle_.ptError  = trk->ptError();
            sigTrack_particle_.energy   = K.getEnergy();
            sigTrack_particle_.dedx     = utility::getDeDx(track_bundle,DeDx_Harm);
            sigTrack_particle_.charge   = trk->charge();
            sigTrack_particle_.dz       = trk->dz(vertex.bestvtx);
            sigTrack_particle_.dzError  = sqrt(TMath::Power(trk->dzError(),2) + TMath::Power(vertex.bestvzError,2));
            sigTrack_particle_.dxy      = trk->dxy(vertex.bestvtx);
            sigTrack_particle_.dxyError = sqrt(TMath::Power(trk->d0Error(),2) + vertex.bestvxError*vertex.bestvyError);
            sigTrack_particle_.eta      = trk->eta();
            sigTrack_particle_.rapidity = K.getRapidity();
            sigTrack_particle_.phi      = trk->phi();
            sigTrack_particle_.ndof     = trk->ndof();
            sigTrack_particle_.vx       = trk->vx();
            sigTrack_particle_.vy       = trk->vy();
            sigTrack_particle_.vz       = trk->vz();
            sigTrack_particle_.vzFlip   = -(trk->vz());
            sigTrack_particle_.chi2     = trk->chi2();
            sigTrack_particle_.chi2norm = trk->normalizedChi2();
            sigTrack_particle_.nhits    = trk->numberOfValidHits();
            try
            {
                if(K.matched(genK))
                {
                    //Enter into signal tree
                    //sigTrack_particle_.momentum = trk->p();
                    //sigTrack_particle_.px       = trk->px();
                    //sigTrack_particle_.py       = trk->py();
                    //sigTrack_particle_.pz       = trk->pz();
                    //sigTrack_particle_.pt       = trk->pt();
                    //sigTrack_particle_.ptError  = trk->ptError();
                    //sigTrack_particle_.energy   = K.getEnergy();
                    //sigTrack_particle_.dedx     = utility::getDeDx(track_bundle,DeDx_Harm);
                    //sigTrack_particle_.charge   = trk->charge();
                    //sigTrack_particle_.dz       = trk->dz(vertex.bestvtx);
                    //sigTrack_particle_.dzError  = sqrt(TMath::Power(trk->dzError(),2) + TMath::Power(vertex.bestvzError,2));
                    //sigTrack_particle_.dxy      = trk->dxy(vertex.bestvtx);
                    //sigTrack_particle_.dxyError = sqrt(TMath::Power(trk->d0Error(),2) + vertex.bestvxError*vertex.bestvyError);
                    //sigTrack_particle_.eta      = trk->eta();
                    //sigTrack_particle_.rapidity = K.getRapidity();
                    //sigTrack_particle_.phi      = trk->phi();
                    //sigTrack_particle_.ndof     = trk->ndof();
                    //sigTrack_particle_.vx       = trk->vx();
                    //sigTrack_particle_.vy       = trk->vy();
                    //sigTrack_particle_.vz       = trk->vz();
                    //sigTrack_particle_.vzFlip   = -(trk->vz());
                    //sigTrack_particle_.chi2     = trk->chi2();
                    //sigTrack_particle_.chi2norm = trk->normalizedChi2();
                    //sigTrack_particle_.nhits    = trk->numberOfValidHits();
                    Signal->Fill();
                    break;
                }
                else
                {
                    //Enter into background tree
                    //bckTrack_particle_.momentum = trk->p();
                    //bckTrack_particle_.px       = trk->px();
                    //bckTrack_particle_.py       = trk->py();
                    //bckTrack_particle_.pz       = trk->pz();
                    //bckTrack_particle_.pt       = trk->pt();
                    //bckTrack_particle_.ptError  = trk->ptError();
                    //bckTrack_particle_.energy   = K.getEnergy();
                    //bckTrack_particle_.dedx     = utility::getDeDx(track_bundle,DeDx_Harm);
                    //bckTrack_particle_.charge   = trk->charge();
                    //bckTrack_particle_.dz       = trk->dz(vertex.bestvtx);
                    //bckTrack_particle_.dzError  = sqrt(TMath::Power(trk->dzError(),2) + TMath::Power(vertex.bestvzError,2));
                    //bckTrack_particle_.dxy      = trk->dxy(vertex.bestvtx);
                    //bckTrack_particle_.dxyError = sqrt(TMath::Power(trk->d0Error(),2) + vertex.bestvxError*vertex.bestvyError);
                    //bckTrack_particle_.eta      = trk->eta();
                    //bckTrack_particle_.rapidity = K.getRapidity();
                    //bckTrack_particle_.phi      = trk->phi();
                    //bckTrack_particle_.ndof     = trk->ndof();
                    //bckTrack_particle_.vx       = trk->vx();
                    //bckTrack_particle_.vy       = trk->vy();
                    //bckTrack_particle_.vz       = trk->vz();
                    //bckTrack_particle_.vzFlip   = -(trk->vz());
                    //bckTrack_particle_.chi2     = trk->chi2();
                    //bckTrack_particle_.chi2norm = trk->normalizedChi2();
                    //bckTrack_particle_.nhits    = trk->numberOfValidHits();
                    Background->Fill();
                }
            }
            catch(const std::invalid_argument& e)
            {
                std::cerr << e.what();
            }
        }
    }


}

void
PhiGenMatch::beginJob()
{
    TH1::SetDefaultSumw2();

    h_nEvt = fs->make<TH1D>("h_nEvt","Events",10,0,10);
    h_phi_yield_rap_1 = fs->make<TH1D>("h_phi_yield_rap_1","Gen Phi Yield rap < 1",50,1,1.05);
    h_phi_yield_norap = fs->make<TH1D>("h_phi_yield_norap","Gen Phi Yield no rap",50,1,1.05);
    h_phid1_mass = fs->make<TH1D>("h_phid1_mass","Dau1 mass",100,0.45,0.55);
    h_phid2_mass = fs->make<TH1D>("h_phid2_mass","Dau2 mass",100,0.45,0.55);
    h_momid = fs->make<TH1D>("h_momid","mother ids",1000000,0,1000000);

    Signal = fs->make<TTree>("SignalTree","SignalTree");
    Background = fs->make<TTree>("BackgroundTree","BackgroundTree");

    Signal->Branch("momentum" , &sigTrack_particle_.momentum);
    Signal->Branch("px"       , &sigTrack_particle_.px      );
    Signal->Branch("py"       , &sigTrack_particle_.py      );
    Signal->Branch("pz"       , &sigTrack_particle_.pz      );
    Signal->Branch("pt"       , &sigTrack_particle_.pt      );
    Signal->Branch("ptError"  , &sigTrack_particle_.ptError );
    Signal->Branch("energy"   , &sigTrack_particle_.energy  );
    Signal->Branch("dedx"     , &sigTrack_particle_.dedx    );
    Signal->Branch("dz"       , &sigTrack_particle_.dz      );
    Signal->Branch("dzError"  , &sigTrack_particle_.dzError );
    Signal->Branch("dxy"      , &sigTrack_particle_.dxy     );
    Signal->Branch("dxyError" , &sigTrack_particle_.dxyError);
    Signal->Branch("eta"      , &sigTrack_particle_.eta     );
    Signal->Branch("rapidity" , &sigTrack_particle_.rapidity);
    Signal->Branch("phi"      , &sigTrack_particle_.phi     );
    Signal->Branch("vx"       , &sigTrack_particle_.vx      );
    Signal->Branch("vy"       , &sigTrack_particle_.vy      );
    Signal->Branch("vz"       , &sigTrack_particle_.vz      );
    Signal->Branch("vzFlip"   , &sigTrack_particle_.vzFlip  );
    Signal->Branch("chi2"     , &sigTrack_particle_.chi2    );
    Signal->Branch("chi2norm" , &sigTrack_particle_.chi2norm);
    Signal->Branch("ndof"     , &sigTrack_particle_.ndof    );
    Signal->Branch("nhits"    , &sigTrack_particle_.nhits   );
    Signal->Branch("charge"   , &sigTrack_particle_.charge  );

    Background->Branch("momentum" , &sigTrack_particle_.momentum);
    Background->Branch("px"       , &sigTrack_particle_.px      );
    Background->Branch("py"       , &sigTrack_particle_.py      );
    Background->Branch("pz"       , &sigTrack_particle_.pz      );
    Background->Branch("pt"       , &sigTrack_particle_.pt      );
    Background->Branch("ptError"  , &sigTrack_particle_.ptError );
    Background->Branch("energy"   , &sigTrack_particle_.energy  );
    Background->Branch("dedx"     , &sigTrack_particle_.dedx    );
    Background->Branch("dz"       , &sigTrack_particle_.dz      );
    Background->Branch("dzError"  , &sigTrack_particle_.dzError );
    Background->Branch("dxy"      , &sigTrack_particle_.dxy     );
    Background->Branch("dxyError" , &sigTrack_particle_.dxyError);
    Background->Branch("eta"      , &sigTrack_particle_.eta     );
    Background->Branch("rapidity" , &sigTrack_particle_.rapidity);
    Background->Branch("phi"      , &sigTrack_particle_.phi     );
    Background->Branch("vx"       , &sigTrack_particle_.vx      );
    Background->Branch("vy"       , &sigTrack_particle_.vy      );
    Background->Branch("vz"       , &sigTrack_particle_.vz      );
    Background->Branch("vzFlip"   , &sigTrack_particle_.vzFlip  );
    Background->Branch("chi2"     , &sigTrack_particle_.chi2    );
    Background->Branch("chi2norm" , &sigTrack_particle_.chi2norm);
    Background->Branch("ndof"     , &sigTrack_particle_.ndof    );
    Background->Branch("nhits"    , &sigTrack_particle_.nhits   );
    Background->Branch("charge"   , &sigTrack_particle_.charge  );

}

void
PhiGenMatch::endJob()
{
}

void
PhiGenMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.addUntracked<edm::InputTag>("genCollection",edm::InputTag("genParticles"));
    desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
    desc.addUntracked<edm::InputTag>("vertexCollName",edm::InputTag("offlinePrimaryVertices"));
    desc.addUntracked<string>("dedxConstraint","loose");
    descriptions.add("PhiGenMatch",desc);
}

DEFINE_FWK_MODULE(PhiGenMatch);
