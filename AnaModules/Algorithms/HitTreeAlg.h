#ifndef HITTREEALG_H
#define HITTREEALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       HitTreeAlg
// Module Type: producer
// File:        HitTreeAlg.h
//
//              The intent of this module is to provide methods for
//              "analyzing" hits on waveforms
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//
// Created by Davide Porzio (salvatore.porzio@postgrad.manchester.ac.uk) on March 02, 2017
//
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace HitTree
{

  // The following typedefs will, obviously, be useful
  using HitPtrVec       = std::vector<art::Ptr<recob::Hit>>;
  using ViewHitMap      = std::map<size_t,HitPtrVec>;
  using TrackViewHitMap = std::map<int,ViewHitMap>;
  using Track = recob::Track;

  class HitTreeAlg
  {
  public:

    // Copnstructors, destructor.
    HitTreeAlg(fhicl::ParameterSet const & pset);
    ~HitTreeAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeTree(art::ServiceHandle<art::TFileService>&, const std::string&);
    void fillTree(const int, const art::Handle<std::vector<recob::Track>>, const TrackViewHitMap&);
    void fillTree(const int, const HitPtrVec&);
    double CalcLength(const Track*);

  private:

    // Fcl parameters.
    std::string fLocalDirName;     ///< Fraction for truncated mean

    // Tree variables
    TTree* fTrackTree;
    int t_event;
    double t_length;
    int t_nHits[3];
    int t_ID;
    std::vector< std::vector<float> > t_PH_v, t_PW_v;

    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
  };

} // end of namespace caldata

#endif
