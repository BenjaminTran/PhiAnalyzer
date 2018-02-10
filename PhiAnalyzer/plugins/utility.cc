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

    bool AcceptTrackDeDx(track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, bool tight = true)
        {
            double functionValueTop = -999;
            double functionValueBot = -999;
            double momentum = track_combo_.track->p();
            double dedx = PhiSelector::getDeDx(track_combo_, DeDxTrack);
            int nhits = track_combo_.track->numberOfValidHits();
            if(tight)
            {
                functionValueTop = 0.55*(TMath::Power(1.6/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3.3;
                functionValueBot = 0.55*(TMath::Power(1.15/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3;
                if(dedx < functionValueTop && dedx > functionValueBot)
                    return true;
                else
                    return false;
            }
            else
            {

                functionValueTop = 0.55*(TMath::Power(1.62/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3.6;
                functionValueBot = 0.5*(TMath::Power(0.5/momentum,4) - 1*TMath::Power(0.50/momentum,2)) + 2.7;
                if(dedx < functionValueTop && dedx > functionValueBot)
                    return true;
                else
                    return false;
            }
        }
}
