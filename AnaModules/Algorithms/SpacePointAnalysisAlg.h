#ifndef SPACEPOINTANALYSISALG_H
#define SPACEPOINTANALYSISALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       SpacePointAnalysisAlg
// Module Type: producer
// File:        SpacePointAnalysisAlg.h
//
//              The intent of this module is to provide methods for
//              "analyzing" hits associated to space points
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//
// Created by Tracy Usher (usher@slac.stanford.edu) on February 19, 2016
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace SpacePointAnalysis
{
    
// The following typedefs will, obviously, be useful
using HitPtrVec          = std::vector<art::Ptr<recob::Hit>>;
using SpacePointHitMap   = std::map<size_t,HitPtrVec>;
using SpacePointPtrVec   = std::vector<art::Ptr<recob::SpacePoint>>;
using TrackSpacePointMap = std::map<size_t,SpacePointPtrVec>;
    
class SpacePointAnalysisAlg
{
public:

    // Copnstructors, destructor.
    SpacePointAnalysisAlg(fhicl::ParameterSet const & pset);
    ~SpacePointAnalysisAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&);
    void endJob(int numEvents);
    
    void fillHistograms(const TrackSpacePointMap&, const SpacePointHitMap&) const;
    void fillHistograms(const SpacePointPtrVec&, const SpacePointHitMap&)   const;
    
private:

    // Fcl parameters.
    std::string fLocalDirName;     ///< Fraction for truncated mean
    
    // Pointers to the histograms we'll create.
    TH1D*     fDeltaTByPlane[3];
    TH1D*     fOverlapSmall[3];
    TH1D*     fOverlapLarge[3];
    
    mutable std::map<geo::View_t, float>       fViewOffsetMap;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
};

} // end of namespace caldata

#endif