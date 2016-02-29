#ifndef CALWIREANALYSISALG_H
#define CALWIREANALYSISALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       CalWireAnalysisAlg
// Module Type: producer
// File:        CalWireAnalysisAlg.h
//
//              The intent of this module is to provide methods for
//              "analyzing" hits on waveforms
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
#include "art/Framework/Core/FindManyP.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardata/RecoBase/Wire.h"
#include "lardata/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace calwireanalysis
{
    
// The following typedefs will, obviously, be useful
using  HitPtrVec  = std::vector<art::Ptr<recob::Hit>>;
using  WirePtrVec = std::vector<art::Ptr<recob::Wire>>;
    
class CalWireAnalysisAlg
{
public:

    // Copnstructors, destructor.
    CalWireAnalysisAlg(fhicl::ParameterSet const & pset);
    ~CalWireAnalysisAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&);
    void endJob(int numEvents);
    
    void fillHistograms(const WirePtrVec&, const art::FindManyP<recob::Hit>&) const;
    
private:

    // Fcl parameters.
    std::string    fLocalDirName;     ///< Fraction for truncated mean
    raw::TDCtick_t fFirstBin;
    raw::TDCtick_t fEndBin;
    float          fThreshold;
    
    // Pointers to the histograms we'll create.
    TH1D*     fNumWiresView[3];
    TH1D*     fWireValueHist[3];
    TH1D*     fROISize[3];
    TH1D*     fNumAboveThres[3];
    TH1D*     fNumHitsROI[3];
    TH1D*     fNumROIsWire[3];
    TH1D*     fNumHitsView[3];
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
};

} // end of namespace caldata

#endif