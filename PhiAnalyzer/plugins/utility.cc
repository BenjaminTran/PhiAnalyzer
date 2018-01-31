/*
 * Making methods that are frequently used across modules
 */

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

namespace utility
{
    int trackFilter(edm::Handle<reco::TrackCollection> tracks,
            edm::Handle<reco::VertexCollection> vertices)
    {
        int Multiplicity = 0;
        double bestvx = -999.9, bestvy = -999.9, bestvz = -999.9;
        double bestvxError = -999.9, bestvyError = -999.9, bestvzError = -999.9;
        const reco::Vertex & vtx = (*vertices)[0];
        bestvx = vtx.x();
        bestvy = vtx.y();
        bestvz = vtx.z();
        bestvxError = vtx.xError();
        bestvyError = vtx.yError();
        bestvzError = vtx.zError();

        for(reco::TrackCollection::const_iterator it = tracks->begin();
                it != tracks->end();
                ++it)
        {
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

            double dzvtx = it->dz(bestvtx);
            double dxyvtx = it->dxy(bestvtx);
            double dzerror = sqrt(it->dzError()*it->dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(it->d0Error()*it->d0Error()+bestvxError*bestvyError);

            if(!it->quality(reco::TrackBase::highPurity)) continue;
            if(fabs(it->ptError())/it->pt()>0.10) continue;
            if(fabs(dzvtx/dzerror) > 3) continue;
            if(fabs(dxyvtx/dxyerror) > 3) continue;

            double eta = it->eta();
            double pt  = it->pt();

            if(fabs(eta)>2.4) continue;
            if(pt<=0.4) continue;
            Multiplicty++;
        }

        return Multiplicity;
    }
}
