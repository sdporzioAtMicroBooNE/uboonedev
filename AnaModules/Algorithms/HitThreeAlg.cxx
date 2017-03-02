#include "HitTreeAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include <cmath>
#include <algorithm>

namespace HitTree
{
  //----------------------------------------------------------------------------
  /// Constructor.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameters.
  ///
  HitTreeAlg::HitTreeAlg(fhicl::ParameterSet const & pset)
  : t_PH_v(3, std::vector<float>())
  , t_PW_v(3, std::vector<float>())
  {
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    reconfigure(pset);

    // Report.
    mf::LogInfo("HitTreeAlg") << "HitTreeAlg configured\n";
  }

  //----------------------------------------------------------------------------
  /// Destructor.
  HitTreeAlg::~HitTreeAlg()
  {

  }

  //----------------------------------------------------------------------------
  /// Reconfigure method.
  ///
  /// Arguments:
  ///
  /// pset - Fcl parameter set.
  ///
  void HitTreeAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fLocalDirName = pset.get<std::string>("LocalDirName", std::string("wow"));
  }

  //----------------------------------------------------------------------------
  /// Begin job method.
  void HitTreeAlg::initializeTree(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
  {
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    fTrackTree = dir.make<TTree>("TrackTree",    "");
    fTrackTree->Branch("Event",&t_event,"event/I");
    fTrackTree->Branch("TrackID",&t_ID,"ID/I");
    fTrackTree->Branch("TrackLength",&t_length,"length/D");
    fTrackTree->Branch("NumHits",&t_nHits,"nHits_U/I:nHits_V/I:nHits_Y/I");
    fTrackTree->Branch("Hit_PulseHeight_U",&t_PH_v[0]);
    fTrackTree->Branch("Hit_PulseHeight_V",&t_PH_v[1]);
    fTrackTree->Branch("Hit_PulseHeight_Y",&t_PH_v[2]);
    fTrackTree->Branch("Hit_PulseWidth_U",&t_PW_v[0]);
    fTrackTree->Branch("Hit_PulseWidth_V",&t_PW_v[1]);
    fTrackTree->Branch("Hit_PulseWidth_Y",&t_PW_v[2]);

    return;
  }

  void HitTreeAlg::fillTree(const int fEvent, const art::Handle<std::vector<recob::Track>> trackHandle, const TrackViewHitMap& trackViewHitMap)
  {
    t_event = fEvent;
    // Looping for each track
    for(const auto& trackHitVecMapItr : trackViewHitMap)
    {
      t_ID = trackHitVecMapItr.first;
      art::Ptr<recob::Track> track(trackHandle,t_ID);
      t_length = CalcLength(track.get());

      t_PH_v.at(0).clear();
      t_PH_v.at(1).clear();
      t_PH_v.at(2).clear();
      t_PW_v.at(0).clear();
      t_PW_v.at(1).clear();
      t_PW_v.at(2).clear();
      t_nHits[0] = 0;
      t_nHits[1] = 0;
      t_nHits[2] = 0;
      // Looping for each plane (?)
      for(const auto& viewHitPair : trackHitVecMapItr.second)
      {
        for(const auto& hitPtr : viewHitPair.second)
        {
          t_PH_v.at(viewHitPair.first).push_back(hitPtr->PeakAmplitude());
          t_PW_v.at(viewHitPair.first).push_back(hitPtr->RMS());
          t_nHits[viewHitPair.first]++;
        }
      }
      fTrackTree->Fill();
    }
    return;
  }

  void HitTreeAlg::fillTree(const int fEvent, const HitPtrVec& hitPtrVec)
  {
    // Looping for each hit in a plane (viewHitPair.second)
    for(const auto& hitPtr : hitPtrVec)
    {
      // Extract interesting hit parameters
      const geo::WireID& wireID   = hitPtr->WireID();
      size_t view     = wireID.Plane;

      t_PH_v.at(view).push_back(hitPtr->PeakAmplitude());
      t_PW_v.at(view).push_back(hitPtr->RMS());
      t_nHits[view]++;
    }
    fTrackTree->Fill();
    return;
  }

  double HitTreeAlg::CalcLength(const recob::Track* track)
  {
    double   result(0.);
    TVector3 disp(track->LocationAtPoint(0));
    TVector3 lastPoint(track->LocationAtPoint(0));
    TVector3 lastDir(0.,0.,0.);
    int      n(track->NumberTrajectoryPoints());

    for(int i = 1; i < n; ++i)
    {
      const TVector3& pos = track->LocationAtPoint(i);

      TVector3 trajDir = pos - lastPoint;

      if (trajDir.Mag2()) trajDir.SetMag(1.);

      //        if (lastDir.Dot(trajDir) >= 0.)
      //        {
      disp   -= pos;
      result += disp.Mag();
      disp    = pos;
      //        }

      lastPoint = pos;
      lastDir   = trajDir;
    }

    return result;
  }


}
