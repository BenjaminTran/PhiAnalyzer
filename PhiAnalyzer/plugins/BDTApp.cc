#include "PhiAnalyzer/PhiAnalyzer/interface/BDTApp.h"


using namespace std;

BDTApp::BDTApp(const edm::ParameterSet& iConfig)
{
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexCollName"));
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
}

BDTApp::~BDTApp()
{
}


void
BDTApp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    h_nEvt->Fill(1);
    std::vector<double> pts = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6};


    std::vector<PhiMeson> Phis;
    std::vector<kaon> bkgPKp;
    std::vector<kaon> bkgPKm;

    edm::Handle<reco::TrackCollection> trkSrc;
    edm::Handle<reco::VertexCollection> vertices;
    edm::Handle<edm::ValueMap<reco::DeDxData> > DeDx_Harm;
    try{
        utility::GetCollection<reco::TrackCollection>(iEvent,_trkSrc,trkSrc);
        utility::GetCollection<reco::VertexCollection>(iEvent,_vertexCollName,vertices);
        utility::GetCollection<edm::ValueMap<reco::DeDxData> >(iEvent,_Dedx_Harmonic2,DeDx_Harm);
    }
    catch(const std::invalid_argument& e){
        std::cerr << e.what();
        return;
    }

    utility::myVertex vertex = utility::MyVertexBuild(vertices);

    for(reco::TrackCollection::const_iterator trk = trkSrc->begin();
            trk != trkSrc->end();
            ++trk)
    {
        //if(!utility::SelectionCut(trk,vertex,false,1.0,1.0,2.4,0,5))
        if(!trk->quality(reco::TrackBase::highPurity)) continue;
        reco::TrackRef track_ref = reco::TrackRef(trkSrc,trk - trkSrc->begin());
        utility::track_combo track_bundle(trk, track_ref);

        kaon::cutVariables kaonCutVariables_;
        kaonCutVariables_.ptError  = trk->ptError();
        kaonCutVariables_.dz       = trk->dz(vertex.bestvtx);
        kaonCutVariables_.dzError  = sqrt(TMath::Power(trk->dzError(),2) + TMath::Power(vertex.bestvzError,2));
        kaonCutVariables_.dxy      = trk->dxy(vertex.bestvtx);
        kaonCutVariables_.dxyError = sqrt(TMath::Power(trk->d0Error(),2) + vertex.bestvxError*vertex.bestvyError);
        kaonCutVariables_.nhits    = trk->numberOfValidHits();
        kaonCutVariables_.chi2     = trk->chi2();
        kaonCutVariables_.chi2norm = trk->normalizedChi2();
        kaonCutVariables_.vx       = trk->vx();
        kaonCutVariables_.vy       = trk->vy();
        kaonCutVariables_.vz       = trk->vz();
        kaonCutVariables_.ndof     = trk->ndof();

        kaon K(TVector3(trk->px(), trk->py(), trk->pz()), trk->eta(), trk->phi(), kaonCutVariables_, trk->charge(), utility::getDeDx(track_bundle, DeDx_Harm), false);

        if(K.getCharge() == 1)
            bkgPKp.push_back(K);
        else if(K.getCharge() == -1)
            bkgPKm.push_back(K);
    }

    //Make tree for event
    TTree* BDTPhiTree = new TTree("tmpBDTphiTree","tmpBDTphiTree");
    BDTPhiTree->Branch( "momentum_1" , &phiKaonCandidate.momentum_1 );
    BDTPhiTree->Branch( "pt_1"       , &phiKaonCandidate.pt_1       );
    BDTPhiTree->Branch( "ptError_1"  , &phiKaonCandidate.ptError_1  );
    BDTPhiTree->Branch( "dedx_1"     , &phiKaonCandidate.dedx_1     );
    BDTPhiTree->Branch( "dz_1"       , &phiKaonCandidate.dz_1       );
    BDTPhiTree->Branch( "dzError_1"  , &phiKaonCandidate.dzError_1  );
    BDTPhiTree->Branch( "dxy_1"      , &phiKaonCandidate.dxy_1      );
    BDTPhiTree->Branch( "dxyError_1" , &phiKaonCandidate.dxyError_1 );
    BDTPhiTree->Branch( "eta_1"      , &phiKaonCandidate.eta_1      );
    BDTPhiTree->Branch( "rapidity_1" , &phiKaonCandidate.rapidity_1 );
    BDTPhiTree->Branch( "nhits_1"    , &phiKaonCandidate.nhits_1    );
    BDTPhiTree->Branch( "momentum_2" , &phiKaonCandidate.momentum_2 );
    BDTPhiTree->Branch( "pt_2"       , &phiKaonCandidate.pt_2       );
    BDTPhiTree->Branch( "ptError_2"  , &phiKaonCandidate.ptError_2  );
    BDTPhiTree->Branch( "dedx_2"     , &phiKaonCandidate.dedx_2     );
    BDTPhiTree->Branch( "dz_2"       , &phiKaonCandidate.dz_2       );
    BDTPhiTree->Branch( "dzError_2"  , &phiKaonCandidate.dzError_2  );
    BDTPhiTree->Branch( "dxy_2"      , &phiKaonCandidate.dxy_2      );
    BDTPhiTree->Branch( "dxyError_2" , &phiKaonCandidate.dxyError_2 );
    BDTPhiTree->Branch( "eta_2"      , &phiKaonCandidate.eta_2      );
    BDTPhiTree->Branch( "rapidity_2" , &phiKaonCandidate.rapidity_2 );
    BDTPhiTree->Branch( "nhits_2"    , &phiKaonCandidate.nhits_2    );

    //Setup reader
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    float local_relpterr_1    , local_relpterr_2;
    float local_dca_z1     , local_dca_z2;
    float local_dca_xy1    , local_dca_xy2;
    float local_rapidity_1 , local_rapidity_2;
    float local_nhits_1    , local_nhits_2;
    float local_dedx_1     , local_dedx_2;
    float local_eta_1      , local_eta_2;
    float local_momentum_1 , local_momentum_2;

    float local_pt_1 , local_pt_2;
    float local_ptError_1, local_ptError_2;
    float local_dzError_1 , local_dzError_2;
    float local_dxyError_1 , local_dxyError_2;
    float local_dz_1, local_dz_2;
    float local_dxy_1, local_dxy_2;

    reader->AddVariable("fabs(pt_1/ptError_1)",&local_relpterr_1);
    reader->AddVariable("fabs(dz_1/dzError_1)",&local_dca_z1);
    reader->AddVariable("fabs(dxy_1/dxyError_1)",&local_dca_xy1);
    reader->AddVariable("rapidity_1",&local_rapidity_1);
    reader->AddVariable("nhits_1",&local_nhits_1);
    reader->AddVariable("dedx_1",&local_dedx_1);
    reader->AddVariable("eta_1",&local_eta_1);
    reader->AddVariable("momentum_1",&local_momentum_1);

    reader->AddVariable("fabs(pt_2/ptError_2)",&local_relpterr_2);
    reader->AddVariable("fabs(dz_2/dzError_2)",&local_dca_z2);
    reader->AddVariable("fabs(dxy_2/dxyError_2)",&local_dca_xy2);
    reader->AddVariable("rapidity_2",&local_rapidity_2);
    reader->AddVariable("nhits_2",&local_nhits_2);
    reader->AddVariable("dedx_2",&local_dedx_2);
    reader->AddVariable("eta_2",&local_eta_2);
    reader->AddVariable("momentum_2",&local_momentum_2);

    std::ostringstream dir;
    //dir << "dataset"
    dir << "PhiAnalyzer/PhiAnalyzer/dataset/dataset_pt1.5-2/weights/Phi_BDT_BDT.weights.xml";

    reader->BookMVA("BDT method",weightFileName.c_str());


    //Build Phis
    Phis = PhiMeson::EventCombinatorialPhi(bkgPKp,bkgPKm);

    for(PhiMeson phi : Phis)
    {
        if(phi.getMass() < 1.0 || phi.getMass() > 1.04) continue;
        utility::FillTreeStruct(phiKaonCandidate, &phi);
        BDTPhiTree->Fill();
        h_masspt->Fill(phi.getMass(),phi.getPt());
    }

    //Setup branch addresses
    BDTPhiTree->SetBranchAddress("momentum_1"  , &local_momentum_1);
    BDTPhiTree->SetBranchAddress( "pt_1"       , &local_pt_1       );
    BDTPhiTree->SetBranchAddress( "ptError_1"  , &local_ptError_1  );
    BDTPhiTree->SetBranchAddress( "dedx_1"     , &local_dedx_1     );
    BDTPhiTree->SetBranchAddress( "dz_1"       , &local_dz_1       );
    BDTPhiTree->SetBranchAddress( "dzError_1"  , &local_dzError_1  );
    BDTPhiTree->SetBranchAddress( "dxy_1"      , &local_dxy_1      );
    BDTPhiTree->SetBranchAddress( "dxyError_1" , &local_dxyError_1 );
    BDTPhiTree->SetBranchAddress( "eta_1"      , &local_eta_1      );
    BDTPhiTree->SetBranchAddress( "rapidity_1" , &local_rapidity_1 );
    BDTPhiTree->SetBranchAddress( "nhits_1"    , &local_nhits_1    );
    BDTPhiTree->SetBranchAddress( "momentum_2" , &local_momentum_2 );
    BDTPhiTree->SetBranchAddress( "pt_2"       , &local_pt_2       );
    BDTPhiTree->SetBranchAddress( "ptError_2"  , &local_ptError_2  );
    BDTPhiTree->SetBranchAddress( "dedx_2"     , &local_dedx_2     );
    BDTPhiTree->SetBranchAddress( "dz_2"       , &local_dz_2       );
    BDTPhiTree->SetBranchAddress( "dzError_2"  , &local_dzError_2  );
    BDTPhiTree->SetBranchAddress( "dxy_2"      , &local_dxy_2      );
    BDTPhiTree->SetBranchAddress( "dxyError_2" , &local_dxyError_2 );
    BDTPhiTree->SetBranchAddress( "eta_2"      , &local_eta_2      );
    BDTPhiTree->SetBranchAddress( "rapidity_2" , &local_rapidity_2 );
    BDTPhiTree->SetBranchAddress( "nhits_2"    , &local_nhits_2    );

   for (Long64_t ievt=0; ievt<BDTPhiTree->GetEntries();ievt++) {

      //if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      BDTPhiTree->GetEntry(ievt);

      local_relpterr_1 = fabs(local_pt_1/local_ptError_1);
      local_dca_z1 = fabs(local_dz_1/local_dzError_1);
      local_dca_xy1 = fabs(local_dxy_1/local_dxyError_1);

      local_relpterr_2 = fabs(local_pt_2/local_ptError_2);
      local_dca_z2 = fabs(local_dz_2/local_dzError_2);
      local_dca_xy2 = fabs(local_dxy_2/local_dxyError_2);

      h_BDTresponse->Fill(reader->EvaluateMVA("BDT method"));
   }

   delete BDTPhiTree;


}

void
BDTApp::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    //edm::FileInPath fip("PhiAnalyzer/PhiAnalyzer/dataset/dataset_pt1.5-2/weights/Phi_BDT_BDT.weights.xml");
    edm::FileInPath fip("PhiAnalyzer/PhiAnalyzer/dataset/dataset_pt5.5-6/weights/Phi_BDT_BDT.weights.xml");
    weightFileName = fip.fullPath();


    h_nEvt = fs->make<TH1D>("nEvt","",10,0,10);
    h_BDTresponse = fs->make<TH1D>("BDTresp","",40,-1,1);
    h_masspt = fs->make<TH2D>("masspt","",100,1,1.05,200,0,20);
    phiKaonTree = fs->make<TTree>("BDTAppTree","BDTAppTree");


    phiKaonTree->Branch( "mass"       , &mass       );
    phiKaonTree->Branch( "pt"       , &pt       );
    phiKaonTree->Branch( "BDTresponse"   , &BDTresponse);
    //phiKaonTree->Branch( "momentum_1" , &phiKaonCandidate.momentum_1 );
    //phiKaonTree->Branch( "pt_1"       , &phiKaonCandidate.pt_1       );
    //phiKaonTree->Branch( "ptError_1"  , &phiKaonCandidate.ptError_1  );
    ////phiKaonTree->Branch( "energy_1"   , &phiKaonCandidate.energy_1   );
    //phiKaonTree->Branch( "dedx_1"     , &phiKaonCandidate.dedx_1     );
    ////phiKaonTree->Branch( "charge_1"   , &phiKaonCandidate.charge_1   );
    //phiKaonTree->Branch( "dz_1"       , &phiKaonCandidate.dz_1       );
    //phiKaonTree->Branch( "dzError_1"  , &phiKaonCandidate.dzError_1  );
    //phiKaonTree->Branch( "dxy_1"      , &phiKaonCandidate.dxy_1      );
    //phiKaonTree->Branch( "dxyError_1" , &phiKaonCandidate.dxyError_1 );
    //phiKaonTree->Branch( "eta_1"      , &phiKaonCandidate.eta_1      );
    //phiKaonTree->Branch( "rapidity_1" , &phiKaonCandidate.rapidity_1 );
    ////phiKaonTree->Branch( "phi_1"      , &phiKaonCandidate.phi_1      );
    ////phiKaonTree->Branch( "vx_1"       , &phiKaonCandidate.vx_1       );
    ////phiKaonTree->Branch( "vy_1"       , &phiKaonCandidate.vy_1       );
    ////phiKaonTree->Branch( "vz_1"       , &phiKaonCandidate.vz_1       );
    ////phiKaonTree->Branch( "px_1"       , &phiKaonCandidate.px_1       );
    ////phiKaonTree->Branch( "py_1"       , &phiKaonCandidate.py_1       );
    ////phiKaonTree->Branch( "pz_1"       , &phiKaonCandidate.pz_1       );
    ////phiKaonTree->Branch( "vzFlip_1"   , &phiKaonCandidate.vzFlip_1   );
    ////phiKaonTree->Branch( "chi2_1"     , &phiKaonCandidate.chi2_1     );
    ////phiKaonTree->Branch( "chi2norm_1" , &phiKaonCandidate.chi2norm_1 );
    ////phiKaonTree->Branch( "ndof_1"     , &phiKaonCandidate.ndof_1     );
    //phiKaonTree->Branch( "nhits_1"    , &phiKaonCandidate.nhits_1    );
    //phiKaonTree->Branch( "momentum_2" , &phiKaonCandidate.momentum_2 );
    //phiKaonTree->Branch( "pt_2"       , &phiKaonCandidate.pt_2       );
    //phiKaonTree->Branch( "ptError_2"  , &phiKaonCandidate.ptError_2  );
    ////phiKaonTree->Branch( "energy_2"   , &phiKaonCandidate.energy_2   );
    //phiKaonTree->Branch( "dedx_2"     , &phiKaonCandidate.dedx_2     );
    ////phiKaonTree->Branch( "charge_2"   , &phiKaonCandidate.charge_2   );
    //phiKaonTree->Branch( "dz_2"       , &phiKaonCandidate.dz_2       );
    //phiKaonTree->Branch( "dzError_2"  , &phiKaonCandidate.dzError_2  );
    //phiKaonTree->Branch( "dxy_2"      , &phiKaonCandidate.dxy_2      );
    //phiKaonTree->Branch( "dxyError_2" , &phiKaonCandidate.dxyError_2 );
    //phiKaonTree->Branch( "eta_2"      , &phiKaonCandidate.eta_2      );
    //phiKaonTree->Branch( "rapidity_2" , &phiKaonCandidate.rapidity_2 );
    ////phiKaonTree->Branch( "phi_2"      , &phiKaonCandidate.phi_2      );
    ////phiKaonTree->Branch( "vx_2"       , &phiKaonCandidate.vx_2       );
    ////phiKaonTree->Branch( "vy_2"       , &phiKaonCandidate.vy_2       );
    ////phiKaonTree->Branch( "vz_2"       , &phiKaonCandidate.vz_2       );
    ////phiKaonTree->Branch( "px_2"       , &phiKaonCandidate.px_2       );
    ////phiKaonTree->Branch( "py_2"       , &phiKaonCandidate.py_2       );
    ////phiKaonTree->Branch( "pz_2"       , &phiKaonCandidate.pz_2       );
    ////phiKaonTree->Branch( "vzFlip_2"   , &phiKaonCandidate.vzFlip_2   );
    ////phiKaonTree->Branch( "chi2_2"     , &phiKaonCandidate.chi2_2     );
    ////phiKaonTree->Branch( "chi2norm_2" , &phiKaonCandidate.chi2norm_2 );
    ////phiKaonTree->Branch( "ndof_2"     , &phiKaonCandidate.ndof_2     );
    //phiKaonTree->Branch( "nhits_2"    , &phiKaonCandidate.nhits_2    );

}

void
BDTApp::endJob()
{
}

void
BDTApp::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
    desc.addUntracked<edm::InputTag>("vertexCollName",edm::InputTag("offlinePrimaryVertices"));
    descriptions.add("BDTApp",desc);
}

DEFINE_FWK_MODULE(BDTApp);
