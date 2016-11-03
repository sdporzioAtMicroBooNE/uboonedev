
#include "SpacePointAnalysisAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include <cmath>
#include <algorithm>

namespace SpacePointAnalysis
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
SpacePointAnalysisAlg::SpacePointAnalysisAlg(fhicl::ParameterSet const & pset) 
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("SpacePointAnalysisAlg") << "SpacePointAnalysisAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
SpacePointAnalysisAlg::~SpacePointAnalysisAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void SpacePointAnalysisAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fLocalDirName = pset.get<std::string>("LocalDirName", std::string("wow"));
}

//----------------------------------------------------------------------------
/// Begin job method.
void SpacePointAnalysisAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    fDeltaTByPlane[0] = dir.make<TH1D>("DeltaT_WU", ";Delta T", 100, -50., 50.);
    fDeltaTByPlane[1] = dir.make<TH1D>("DeltaT_WV", ";Delta T", 100, -50., 50.);
    fDeltaTByPlane[2] = dir.make<TH1D>("DeltaT_VU", ";Delta T", 100, -50., 50.);
    
    fOverlapSmall[0]  = dir.make<TH1D>("OverlapSmall_UV", ";Overlap Fraction", 120, -1.1, 1.1);
    fOverlapSmall[1]  = dir.make<TH1D>("OverlapSmall_WU", ";Overlap Fraction", 120, -1.1, 1.1);
    fOverlapSmall[2]  = dir.make<TH1D>("OverlapSmall_WV", ";Overlap Fraction", 120, -1.1, 1.1);
    
    fOverlapLarge[0]  = dir.make<TH1D>("OverlapLarge_UV", ";Overlap Fraction", 120, -1.1, 1.1);
    fOverlapLarge[1]  = dir.make<TH1D>("OverlapLarge_WU", ";Overlap Fraction", 120, -1.1, 1.1);
    fOverlapLarge[2]  = dir.make<TH1D>("OverlapLarge_WV", ";Overlap Fraction", 120, -1.1, 1.1);

    return;
}
    
void SpacePointAnalysisAlg::fillHistograms(const TrackSpacePointMap& trackSpacePointMap, const SpacePointHitMap& spacePointHitMap) const
{
    
    return;
}
    
void SpacePointAnalysisAlg::fillHistograms(const SpacePointPtrVec& spacePointPtrVec, const SpacePointHitMap& spacePointHitMap) const
{
    // We'll need the offsets for each plane
    fViewOffsetMap.insert(std::pair<geo::View_t,float>(geo::kU,float(fDetectorProperties->GetXTicksOffset(geo::kU, 0, 0)-fDetectorProperties->TriggerOffset())));
    fViewOffsetMap.insert(std::pair<geo::View_t,float>(geo::kV,float(fDetectorProperties->GetXTicksOffset(geo::kV, 0, 0)-fDetectorProperties->TriggerOffset())));
    fViewOffsetMap.insert(std::pair<geo::View_t,float>(geo::kW,float(fDetectorProperties->GetXTicksOffset(geo::kW, 0, 0)-fDetectorProperties->TriggerOffset())));
    
    // Loop through space points
    for(const auto& spacePoint : spacePointPtrVec)
    {
        // Recover vector of 2D hits associated to this space point
        SpacePointHitMap::const_iterator spacePointHitItr = spacePointHitMap.find(spacePoint.key());
        
        if (spacePointHitItr != spacePointHitMap.end())
        {
            const HitPtrVec& hitPtrVec = spacePointHitItr->second;
            
            // Focus only on space points with 3 2D hits for now
            if (hitPtrVec.size() != 3) continue;
            
            // Organize the associated hits by view so we can access easily
            std::map<geo::View_t,const art::Ptr<recob::Hit>> viewHitMap;
            
            for(const auto& hit : hitPtrVec) viewHitMap.insert(std::pair<geo::View_t,const art::Ptr<recob::Hit>>(hit->View(),hit));
            
            // Ok, first thing to look at is the time different plane to plane
            float peakTimeU = viewHitMap.at(geo::kU)->PeakTime() - fViewOffsetMap.at(geo::kU);
            float peakTimeV = viewHitMap.at(geo::kV)->PeakTime() - fViewOffsetMap.at(geo::kV);
            float peakTimeW = viewHitMap.at(geo::kW)->PeakTime() - fViewOffsetMap.at(geo::kW);
            
            float deltaT_WU = peakTimeW - peakTimeU;
            float deltaT_WV = peakTimeW - peakTimeV;
            float deltaT_VU = peakTimeV - peakTimeU;
            
            // Make some hists
            fDeltaTByPlane[0]->Fill(deltaT_WU, 1.);
            fDeltaTByPlane[1]->Fill(deltaT_WV, 1.);
            fDeltaTByPlane[2]->Fill(deltaT_VU, 1.);
            
            // Look at overlaps
            float hitUWidth       = viewHitMap.at(geo::kU)->RMS();
            float hitVWidth       = viewHitMap.at(geo::kV)->RMS();
            float hitWWidth       = viewHitMap.at(geo::kW)->RMS();
            
            // Get overlap fractions
            float maxUpper_VU     = std::min(peakTimeV+hitVWidth,peakTimeU+hitUWidth);
            float minLower_VU     = std::max(peakTimeV-hitVWidth,peakTimeU-hitUWidth);
            float overlap_VU      = maxUpper_VU - minLower_VU;
            float overlapSmall_VU = 0.5 * overlap_VU / std::min(hitVWidth,hitUWidth);
            float overlapLarge_VU = 0.5 * overlap_VU / std::max(hitVWidth,hitUWidth);
            
            fOverlapSmall[0]->Fill(overlapSmall_VU, 1.);
            fOverlapLarge[0]->Fill(overlapLarge_VU, 1.);
            
            float maxUpper_WU     = std::min(peakTimeW+hitWWidth,peakTimeU+hitUWidth);
            float minLower_WU     = std::max(peakTimeW-hitWWidth,peakTimeU-hitUWidth);
            float overlap_WU      = maxUpper_WU - minLower_WU;
            float overlapSmall_WU = 0.5 * overlap_WU / std::min(hitWWidth,hitUWidth);
            float overlapLarge_WU = 0.5 * overlap_WU / std::max(hitWWidth,hitUWidth);
            
            fOverlapSmall[1]->Fill(overlapSmall_WU, 1.);
            fOverlapLarge[1]->Fill(overlapLarge_WU, 1.);
            
            float maxUpper_WV     = std::min(peakTimeW+hitWWidth,peakTimeV+hitVWidth);
            float minLower_WV     = std::max(peakTimeW-hitWWidth,peakTimeV-hitVWidth);
            float overlap_WV      = maxUpper_WV - minLower_WV;
            float overlapSmall_WV = 0.5 * overlap_WV / std::min(hitWWidth,hitVWidth);
            float overlapLarge_WV = 0.5 * overlap_WV / std::max(hitWWidth,hitVWidth);
            
            fOverlapSmall[2]->Fill(overlapSmall_WV, 1.);
            fOverlapLarge[2]->Fill(overlapLarge_WV, 1.);
        }
    }
    return;
}
    
// Useful for normalizing histograms
void SpacePointAnalysisAlg::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
//    double normFactor(1./numEvents);
    
//    for(size_t idx = 0; idx < 3; idx++) fHitsByWire->Scale(normFactor);
    
    return;
}

}
