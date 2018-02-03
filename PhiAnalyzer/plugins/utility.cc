/*
 * Making methods that are frequently used across modules
 */

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

namespace utility
{
    struct myVertex{
        double bestvx = -999;
        double bestvy = -999;
        double bestvz = -999;
        double bestvxError = -999;
        double bestvyError = -999;
        double bestvzError = -999;

        reco::Vertex vtx;
        math::XYZPoint bestvtx;

        myVertex(double vx, double vy, double vz, double vxErr, double vyErr, double vzErr, reco::Vertex &vtx_, math::XYZPoint bestvtx_) :
            bestvx(vx), bestvy(vy), bestvz(vz), bestvxError(vxErr), bestvyError(vyErr), bestvzError(vzErr), vtx(vtx_), bestvtx(bestvtx_) {}

    };

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
        //double bestvx = -999.9, bestvy = -999.9, bestvz = -999.9;
        //double bestvxError = -999.9, bestvyError = -999.9, bestvzError = -999.9;
        //const reco::Vertex & vtx = (*vertices)[0];
        //bestvx = vtx.x();
        //bestvy = vtx.y();
        //bestvz = vtx.z();
        //bestvxError = vtx.xError();
        //bestvyError = vtx.yError();
        //bestvzError = vtx.zError();

        myVertex myVtx = myVertexBuild(vertices);

        for(reco::TrackCollection::const_iterator it = tracks->begin();
                it != tracks->end();
                ++it)
        {
            //double dzvtx = it->dz(bestvtx);
            //double dxyvtx = it->dxy(bestvtx);
            //double dzerror = sqrt(it->dzError()*it->dzError()+bestvzError*bestvzError);
            //double dxyerror = sqrt(it->d0Error()*it->d0Error()+bestvxError*bestvyError);

            //if(!it->quality(reco::TrackBase::highPurity)) continue;
            //if(fabs(it->ptError())/it->pt()>0.10) continue;
            //if(fabs(dzvtx/dzerror) > 3) continue;
            //if(fabs(dxyvtx/dxyerror) > 3) continue;

            //double eta = it->eta();
            //double pt  = it->pt();

            //if(fabs(eta)>2.4) continue;
            //if(pt<=0.4) continue;
            if(isTrackGood(it,myVtx))
                Multiplicity++;
            else
                continue;
        }

        return Multiplicity;
    }
}
