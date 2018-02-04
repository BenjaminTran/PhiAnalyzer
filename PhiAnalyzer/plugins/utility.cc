/*
 * Making methods that are frequently used across modules
 */

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

namespace utility
{
    myVertex MyVertexBuild(edm::Handle<reco::VertexCollection> vertices)
    {
        const reco::Vertex & vtx = (*vertices)[0];
        math::XYZPoint vtx_point(vtx.x(),vtx.y(),vtx.z());
        myVertex myVertex_(vtx.x(),vtx.y(),vtx.z(),vtx.xError(),vtx.yError(),vtx.zError(),vtx,vtx_point);
        return myVertex_;
    }

    bool isTrackGood(reco::TrackCollection::const_iterator &track, myVertex myVtx)
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
        if(track->pt() <= 0.4) return false;

        return true;
    }


    int TrackFilter(edm::Handle<reco::TrackCollection> tracks,
            edm::Handle<reco::VertexCollection> vertices)
    {
        int Multiplicity = 0;

        myVertex myVtx = MyVertexBuild(vertices);

        for(reco::TrackCollection::const_iterator it = tracks->begin();
                it != tracks->end();
                ++it)
        {
            if(isTrackGood(it,myVtx))
                Multiplicity++;
            else
                continue;
        }

        return Multiplicity;
    }
}
