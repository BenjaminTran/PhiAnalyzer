/*
 * Making methods that are frequently used across modules
 */

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/PhiSelector.h"

namespace utility
{
    const double kaonMass = 0.493677;
    myVertex MyVertexBuild(edm::Handle<reco::VertexCollection> vertices)
    {
        const reco::Vertex & vtx = (*vertices)[0];
        math::XYZPoint vtx_point(vtx.x(),vtx.y(),vtx.z());
        myVertex myVertex_(vtx.x(),vtx.y(),vtx.z(),vtx.xError(),vtx.yError(),vtx.zError(),vtx,vtx_point);
        return myVertex_;
    }

    //Multiplicity criteria
    bool isTrackGood(reco::TrackCollection::const_iterator &track, myVertex myVtx, bool ptCut = true)
    {
        double dzvtx = track->dz(myVtx.bestvtx);
        double dxyvtx = track->dxy(myVtx.bestvtx);
        double dzerror = sqrt(track->dzError()*track->dzError()+myVtx.bestvzError*myVtx.bestvzError);
        double dxyerror = sqrt(track->d0Error()*track->d0Error()+myVtx.bestvxError*myVtx.bestvyError);

        if(!track->quality(reco::TrackBase::highPurity)) return false;
        if(fabs(track->ptError())/track->pt()>0.10) return false;
        if(fabs(dzvtx/dzerror) > 3) return false;
        if(fabs(dxyvtx/dxyerror) > 3) return false;
        if(fabs(track->eta()) > 3) return false;
        if(ptCut)
        {
            if(track->pt() <= 0.4) return false;
        }

        return true;
    }

    //Check if track passes the given constraints
    bool SelectionCut(reco::TrackCollection::const_iterator &track, myVertex myVtx, bool ptCut, double dzdca, double dxydca, double eta, double ptCutVal, int nhits)
    {
        double dzvtx = track->dz(myVtx.bestvtx);
        double dxyvtx = track->dxy(myVtx.bestvtx);
        double dzerror = sqrt(track->dzError()*track->dzError()+myVtx.bestvzError*myVtx.bestvzError);
        double dxyerror = sqrt(track->d0Error()*track->d0Error()+myVtx.bestvxError*myVtx.bestvyError);

        if(!track->quality(reco::TrackBase::highPurity)) return false;
        if(track->numberOfValidHits() < nhits) return false;
        if(fabs(track->ptError())/track->pt()>0.10) return false;
        if(fabs(dzvtx/dzerror) > dzdca) return false;
        if(fabs(dxyvtx/dxyerror) > dxydca) return false;
        if(fabs(track->eta()) > eta) return false;
        if(ptCut)
        {
            if(track->pt() <= ptCutVal) return false;
        }

        return true;

    }


    //Return the multiplicity of the event
    int Multiplicity(edm::Handle<reco::TrackCollection> tracks,
            edm::Handle<reco::VertexCollection> vertices)
    {
        int Multiplicity = 0;

        myVertex myVtx = MyVertexBuild(vertices);

        for(reco::TrackCollection::const_iterator it = tracks->begin();
                it != tracks->end();
                ++it)
        {
            if(isTrackGood(it,myVtx,true))
                Multiplicity++;
            else
                continue;
        }

        return Multiplicity;
    }

    /*
     * For constraint use tight, default, or loose. Default will always be
     * defined as the band of choice
     */
    bool AcceptTrackDeDx(track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, std::string constraint = "default")
    {
        double functionValueTop = -999;
        double functionValueBot = -999;
        double momentum = track_combo_.track->p();
        double dedx = getDeDx(track_combo_, DeDxTrack);
        //int nhits = track_combo_.track->numberOfValidHits();
        if(constraint == "tight")
        {
            functionValueTop = 0.55*(TMath::Power(1.6/momentum,2) - 2.95*TMath::Power(0.6/momentum,1)) + 3.3;
            functionValueBot = 0.55*(TMath::Power(1.15/momentum,2) - 1.7*TMath::Power(0.6/momentum,1)) + 2.7;

            if(dedx < functionValueTop && dedx > functionValueBot)
                return true;
            else
                return false;
        }
        else if(constraint == "default")
        {
            functionValueTop = 0.55*(TMath::Power(1.62/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3.6;
            functionValueBot = 0.55*(TMath::Power(1.15/momentum,2) - 1.7*TMath::Power(0.6/momentum,1)) + 2.5;

            if(dedx < functionValueTop && dedx > functionValueBot)
                return true;
            else
                return false;
        }
        else if(constraint == "loose")
        {

            functionValueTop = 0.55*(TMath::Power(1.62/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3.6;
            functionValueBot = 0.5*(TMath::Power(0.5/momentum,4) - 1*TMath::Power(0.50/momentum,2)) + 2.7;

            if(dedx < functionValueTop && dedx > functionValueBot)
                return true;
            else
                return false;
        }
        else
        {
            throw std::invalid_argument(constraint + " is not a valid dedx function constraint!");
        }
    }

    double getDeDx(utility::track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack)
    {
        double dedx_ = -1;
        const edm::ValueMap<reco::DeDxData> dedxTrack = *DeDxTrack.product();
        dedx_ = dedxTrack[track_combo_.track_ref].dEdx();

        return dedx_;
    }

    std::vector<PhiMeson> EventCombinatorialPhi(std::vector<kaon> PKp_, std::vector<kaon> PKm_)
    {
        std::vector<PhiMeson> phiCollection;
        for(kaon Pkp : PKp_)
        {
            for(kaon Pkm : PKm_)
            {
                TVector3 dau1p(Pkp.getPx(), Pkp.getPy(), Pkp.getPz());
                TVector3 dau2p(Pkm.getPx(), Pkm.getPy(), Pkm.getPz());
                TVector3 dauPsum(dau1p + dau2p);
                double mass = sqrt(TMath::Power(Pkp.getEnergy()+Pkm.getEnergy(),2) - dauPsum.Mag2());

                double p = sqrt(dauPsum.Mag2());
                double pt = dauPsum.Perp();

                TLorentzVector phiLV(dauPsum,Pkp.getEnergy() + Pkm.getEnergy());

                double rapidity = phiLV.Rapidity();
                double eta = phiLV.Eta();
                double phi = phiLV.Phi();

                PhiMeson pgf(mass,pt,eta,phi,rapidity,p);

                phiCollection.push_back(pgf);
            }
        }

        return phiCollection;
    }
}
