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
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Spacepoint.h"

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
    void endJob(int numEvents);

    void compareTwoTracks(const recob::Track*, const recob::Track*) const;
    
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
    
    // Keep track of histograms
    std::vector<TH1D*>                  fTH1DVec;
    std::vector<std::vector<TH1D*>>     fTH1DVecVec;   // hists by view (for example
    std::vector<TH2D*>                  fTH2DVec;
    std::vector<std::vector<TH2D*>>     fTH2DVecVec;
    std::vector<TProfile*>              fTProfileVec;
    std::vector<std::vector<TProfile*>> fTProfileVecVec;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             // pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
};

} // end of namespace caldata

#endif