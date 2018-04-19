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

    std::vector<std::vector<kaon> > genDauKaons;
    std::vector<PhiMeson> SignalPhis;
    std::vector<PhiMeson> BackgroundPhis;
    /*
     * For background kaons that did not match. Will be used to construct background phis
     */
    std::vector<kaon> bkgPKp;
    std::vector<kaon> bkgPKm;

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

            std::vector<kaon> tmpKaonDau;
            tmpKaonDau.push_back(kaon(TVector3(d1->px(), d1->py(), d1->pz()), d1->eta(), d1->phi(), d1->charge(), true));
            tmpKaonDau.push_back(kaon(TVector3(d2->px(), d2->py(), d2->pz()), d2->eta(), d2->phi(), d2->charge(), true));

            genDauKaons.push_back(tmpKaonDau);

            //genDauKaons.push_back(kaon(TVector3(d1->px(), d1->py(), d1->pz()), d1->eta(), d1->phi(), d1->charge(), true));
            //genDauKaons.push_back(kaon(TVector3(d2->px(), d2->py(), d2->pz()), d2->eta(), d2->phi(), d2->charge(), true));
        }
    }

    int NumberOfGenPhis = genDauKaons.size();

    std::array<std::vector<kaon>, NumberOfGenPhis> trackKaonPairs;

    //Perform Matching
    for(reco::TrackCollection::const_iterator trk = trkSrc->begin();
            trk != trkSrc->end();
            ++trk)
    {
        //if(!utility::SelectionCut(trk,vertex,false,1.0,1.0,2.4,0,5))
        reco::TrackRef track_ref = reco::TrackRef(trkSrc,trk - trkSrc->begin());
        utility::track_combo track_bundle(trk, track_ref);

        kaon::cutVariables kaonCutVariables;
        kaonCutVariables_.ptError = trk->ptError();
        kaonCutVariables_.dz = trk->dz(vertex.bestvtx);
        kaonCutVariables_.dzError = sqrt(TMath::Power(trk->dzError(),2) + TMath::Power(vertex.bestvzError,2));
        kaonCutVariables_.dxy = trk->dxy(vertex.bestvtx);
        kaonCutVariables_.dxyError = sqrt(TMath::Power(trk->d0Error(),2) + vertex.bestvxError*vertex.bestvyError);
        kaonCutVariables_.nhits = trk->numberOfValidHits();
        kaonCutVariables_.chi2 = trk->chi2();
        kaonCutVariables_.chi2norm = trk->normalizedChi2();
        kaonCutVariables_.vx = trk->vx();
        kaonCutVariables_.vy = trk->vy();
        kaonCutVariables_.vz = trk->vz();
        kaonCutVariables_.ndof = trk->ndof();

        kaon K(TVector3(trk->px(), trk->py(), trk->pz()), trk->eta(), trk->phi(), kaonCutVariables, trk->charge(), utility::getDeDx(track_bundle, DeDx_Harm), false);

        bool kaonMatched = false;

        for(int i=0; i<genDauKaons.size(); i++)
        {
            std::vector<kaon> genKaonPair = genDauKaons[i];
            for(int j=0; j<genKaonPair.size(); j++)
            {
                try
                {
                    if(K.matched(genKaonPair[j]))
                    {
                        trackKaonPairs.at(i).push_back(K);
                        kaonMatched = true;
                        break;
                    }
                }
                catch(const std::invalid_argument& e)
                {
                    std::cerr << e.what();
                }
            }
            if(kaonMatched)
                break;
        }
        if(!kaonMatched)
        {
            if(K.getCharge() == 1)
                bkgPKp.push_back(K);
            else if(K.getCharge() == -1)
                bkgPKm.push_back(K);
        }
    }

    //Build signal Phis
    for(std::vector<kaon> kaonPair : trackKaonPairs)
    {
        //Check if there are two kaons in each container. If not then skip
        if(kaonPair.size() < 2) continue;

        PhiMeson phi = BuildPhi(kaonPair[0],kaonPair[1],true);

        SignalPhis.push_back(phi);
    }

    //Build background Phis
    BackgroundPhis = PhiMeson::EventCombinatorialPhi(bkgPKp, bkgPKm);

    for(PhiMeson phi : SignalPhis)
    {
        kaon dau1 = phi.getKaonDau(0);
        kaon dau2 = phi.getKaonDau(1);
        particle_.mass = phi.getMass();
        particle_.momentum_1 = dau1.getP();
        particle_.pt_1 = dau1.getPt();
        particle_.ptError_1 = dau1.
        particle_.energy_1 = dau1.getEnergy();
        particle_.dedx_1 = dau1.getDedx();
        particle_.charge_1 = dau1.getCharge();
        particle_.dz_1 = dau1.
        particle_.dzError_1 = dau1.
        particle_.dxy_1 = dau1.
        particle_.dxyError_1 = dau1.
        particle_.eta_1 = dau1.getEta();
        particle_.rapidity_1 = dau1.getRapidity();
        particle_.phi_1 = phi.getPhi();
        particle_.vx_1 = dau1.
        particle_.vy_1 = dau1.
        particle_.vz_1 = dau1.
        particle_.px_1 = dau1.getPx();
        particle_.py_1 = dau1.getPy();
        particle_.pz_1 = dau1.getPz();
        particle_.vzFlip_1 = dau1.
        particle_.chi2_1 = dau1.
        particle_.chi2norm_1 = dau1.
        particle_.ndof_1 = dau1.
        particle_.momentum_2 = phi.
        particle_.pt_2 = phi.
        particle_.ptError_2 = phi.
        particle_.energy_2 = phi.
        particle_.dedx_2 = phi.
        particle_.charge_2 = phi.
        particle_.dz_2 = phi.
        particle_.dzError_2 = phi.
        particle_.dxy_2 = phi.
        particle_.dxyError_2 = phi.
        particle_.eta_2 = phi.
        particle_.rapidity_2 = phi.
        particle_.phi_2 = phi.
        particle_.vx_2 = phi.
        particle_.vy_2 = phi.
        particle_.vz_2 = phi.
        particle_.px_2 = phi.
        particle_.py_2 = phi.
        particle_.pz_2 = phi.
        particle_.vzFlip_2 = phi.
        particle_.chi2_2 = phi.
        particle_.chi2norm_2 = phi.
        particle_.ndof_2 = phi.
        particle_.nhits_2 = phi.
        particle_.nhits_2 = phi.
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

    Signal->Branch("momentum" , &Track_particle_.momentum);
    Signal->Branch("px"       , &Track_particle_.px      );
    Signal->Branch("py"       , &Track_particle_.py      );
    Signal->Branch("pz"       , &Track_particle_.pz      );
    Signal->Branch("pt"       , &Track_particle_.pt      );
    Signal->Branch("ptError"  , &Track_particle_.ptError );
    Signal->Branch("energy"   , &Track_particle_.energy  );
    Signal->Branch("dedx"     , &Track_particle_.dedx    );
    Signal->Branch("dz"       , &Track_particle_.dz      );
    Signal->Branch("dzError"  , &Track_particle_.dzError );
    Signal->Branch("dxy"      , &Track_particle_.dxy     );
    Signal->Branch("dxyError" , &Track_particle_.dxyError);
    Signal->Branch("eta"      , &Track_particle_.eta     );
    Signal->Branch("rapidity" , &Track_particle_.rapidity);
    Signal->Branch("phi"      , &Track_particle_.phi     );
    Signal->Branch("vx"       , &Track_particle_.vx      );
    Signal->Branch("vy"       , &Track_particle_.vy      );
    Signal->Branch("vz"       , &Track_particle_.vz      );
    Signal->Branch("vzFlip"   , &Track_particle_.vzFlip  );
    Signal->Branch("chi2"     , &Track_particle_.chi2    );
    Signal->Branch("chi2norm" , &Track_particle_.chi2norm);
    Signal->Branch("ndof"     , &Track_particle_.ndof    );
    Signal->Branch("nhits"    , &Track_particle_.nhits   );
    Signal->Branch("charge"   , &Track_particle_.charge  );

    Background->Branch("momentum" , &Track_particle_.momentum);
    Background->Branch("px"       , &Track_particle_.px      );
    Background->Branch("py"       , &Track_particle_.py      );
    Background->Branch("pz"       , &Track_particle_.pz      );
    Background->Branch("pt"       , &Track_particle_.pt      );
    Background->Branch("ptError"  , &Track_particle_.ptError );
    Background->Branch("energy"   , &Track_particle_.energy  );
    Background->Branch("dedx"     , &Track_particle_.dedx    );
    Background->Branch("dz"       , &Track_particle_.dz      );
    Background->Branch("dzError"  , &Track_particle_.dzError );
    Background->Branch("dxy"      , &Track_particle_.dxy     );
    Background->Branch("dxyError" , &Track_particle_.dxyError);
    Background->Branch("eta"      , &Track_particle_.eta     );
    Background->Branch("rapidity" , &Track_particle_.rapidity);
    Background->Branch("phi"      , &Track_particle_.phi     );
    Background->Branch("vx"       , &Track_particle_.vx      );
    Background->Branch("vy"       , &Track_particle_.vy      );
    Background->Branch("vz"       , &Track_particle_.vz      );
    Background->Branch("vzFlip"   , &Track_particle_.vzFlip  );
    Background->Branch("chi2"     , &Track_particle_.chi2    );
    Background->Branch("chi2norm" , &Track_particle_.chi2norm);
    Background->Branch("ndof"     , &Track_particle_.ndof    );
    Background->Branch("nhits"    , &Track_particle_.nhits   );
    Background->Branch("charge"   , &Track_particle_.charge  );

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
