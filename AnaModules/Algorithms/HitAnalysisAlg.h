#ifndef HITANALYSISALG_H
#define HITANALYSISALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       HitAnalysisAlg
// Module Type: producer
// File:        HitAnalysisAlg.h
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
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardata/RecoBase/Hit.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace hitanalysis
{
    
// The following typedefs will, obviously, be useful
using  HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
    
class HitAnalysisAlg
{
public:

    // Copnstructors, destructor.
    HitAnalysisAlg(fhicl::ParameterSet const & pset);
    ~HitAnalysisAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&);
    void endJob(int numEvents);
    
    void fillHistograms(const HitPtrVec& hitPtrVec) const;
    
private:

    // Fcl parameters.
    std::string fLocalDirName;     ///< Fraction for truncated mean
    
    // Pointers to the histograms we'll create.
    TH1D*     fHitsByWire[3];
    TH1D*     fPulseHeight[3];
    TH1D*     fPulseHeightSingle[3];
    TH1D*     fPulseHeightMulti[3];
    TH1D*     fChi2DOF[3];
    TH1D*     fNumDegFree[3];
    TH1D*     fChi2DOFSingle[3];
    TH1D*     fHitMult[3];
    TH1D*     fHitCharge[3];
    TH1D*     fFitWidth[3];
    TH1D*     fHitSumADC[3];
    TH2D*     fNDFVsChi2[3];
    TH2D*     fPulseHVsWidth[3];
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
};

} // end of namespace caldata

#endif