#ifndef PANDORAANALYSISALG_H
#define PANDORAANALYSISALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       PandoraAnalysisAlg
// Module Type: producer
// File:        PandoraAnalysisAlg.h
//
//              The intent of this module is to provide methods for
//              "analyzing" the output of pandora
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
#include "art/Framework/Principal/Event.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Spacepoint.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace pandoraanalysis
{
    
// The following typedefs will, obviously, be useful
using  HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
    
class PandoraAnalysisAlg
{
public:

    // Copnstructors, destructor.
    PandoraAnalysisAlg(fhicl::ParameterSet const & pset);
    ~PandoraAnalysisAlg();

    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset);
    void initializeHists(art::ServiceHandle<art::TFileService>&, const std::string&);
    
    void pandoraAnalysis(const art::Event& event) const;
    
private:
    
    // Traverse PFParticle hierarchy
    int traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>&,
                                    size_t,
                                    const art::FindManyP<recob::Track>&,
                                    const art::FindManyP<recob::Vertex>&,
                                    int&,
                                    int&) const;
    
    // The following typedefs will, obviously, be useful
    using  HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
    
    double length(const recob::Track* track) const;
    
    double projectedLength(const recob::Track* track) const;
        
    // The parameters we'll read from the .fcl file.
    std::string fHitProducerLabel;
    std::string fPFParticleProducerLabel;
    std::string fTrackProducerLabel;
    
    // Pointers to the histograms we'll create.
    // Track histograms
    TH1D*     fNPrimPFParticles;
    TH1D*     fNPrimWDaughters;
    TH1D*     fNpdg13PFParticles;
    TH1D*     fNpdg11PFParticles;
    TH1D*     fNDaughters;
    TH1D*     fNTracks;
    TH1D*     fNVertices;
    TH1D*     fNHitsPFparticles;
    TH1D*     fNHitsWDaughters;
    TH1D*     fNHitsNDaughters;
    
    TH1D*     fNumPFPartTracks;
    TH1D*     fNumPFPartTracksL;
    TH1D*     fNumFitTracks;
    TH1D*     fNumTrksPFPart;
    TH1D*     fNumTrksPFPartL;
    TH1D*     fNumFitTrksPFPart;
    TH1D*     fNumFitTrksPFPartL;
    TH1D*     fPFPartTrackLen;
    TH1D*     fPFPartTrackLenL;
    TH1D*     fPFPartEndLenL;
    TH1D*     fFitTrackLen;
    TH1D*     fFitTrackProjLen;
    TH1D*     fFitEndLen;
    TH1D*     fTrackDeltaLen;
    TH1D*     fTrackDeltaProjLen;
    TH2D*     fFitVsPFPartLen;
    TH2D*     fFitELVsTL;
    TProfile* fFitVsPFPartEff;
    TProfile* fFitVsPFNHitsEff;
    
    TH1D*     fDeltaStartPos;
    TH1D*     fDeltaEndPos;
    TH1D*     fTrackDeltaStart;
    TH1D*     fCosTracks;
    TH2D*     fDStartVsDEnd;
    
    TH1D*     fDeltaWiresTrk[3];
    TH1D*     fNumHitsTrk[3];
    TH1D*     fHitWireRatioTrk[3];
    TProfile* fRatioVsDWires[3];
    TH1D*     fTrajDispDiff;
    TH1D*     fTrajDispAng;
    TH1D*     fTrajAng;
    TH1D*     fTrajDispAngRev;
    TH1D*     fTrajDocaAll;
    TH1D*     fTrajDoca[3];
    TH1D*     fTrajStartDiff;
    TH1D*     fTrajEndDiff;
    
    TH1D*     fNumPFPartHits;
    TH1D*     fNumPFPartViews;
    TH1D*     fNumPFPartViewsL;
    TH2D*     fViewVsHits;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    art::ServiceHandle<geo::Geometry>            fGeometry;             ///< pointer to Geometry service
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;   ///< Detector properties service
};

} // end of namespace caldata

#endif