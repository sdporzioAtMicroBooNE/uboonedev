////////////////////////////////////////////////////////////////////////
// Class:       SpecialTracks
// Module Type: producer
// File:        SpecialTracks_module.cc
//              The goal of this module is to output a set of tracks,
//              hits, clusters, PFParticles, etc., which correspond
//              to some special selection criteria
//
// Generated at Wed Sep 17 19:17:00 2014 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::SpecialTracks
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/geo.h"

#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/Hit.h"

#include "lardata/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"


class SpecialTracks : public art::EDProducer
{
public:
    explicit SpecialTracks(fhicl::ParameterSet const & p);
    virtual ~SpecialTracks();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

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
    
    double length(const recob::Track* track);
    
    double projectedLength(const recob::Track* track);
    
    std::string fPFParticleModuleLabel;
    std::string fHitProducerLabel;
    std::string fTrackModuleLabel;
    
    // Other variables that will be shared between different methods.
    const geo::GeometryCore* fGeometry;       // pointer to Geometry service
};


SpecialTracks::SpecialTracks(fhicl::ParameterSet const & p)
{
    fGeometry = lar::providerFrom<geo::Geometry>();
    
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<recob::PFParticle>>();
    produces< std::vector<recob::Cluster>>();
    produces< std::vector<recob::Track>>();
    produces< std::vector<recob::Vertex>>();
    produces< std::vector<recob::Hit>>();
//    produces< std::vector<recob::SpacePoint> >();
//    produces< std::vector<recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster>>();
    produces< art::Assns<recob::PFParticle, recob::Track>>();
    produces< art::Assns<recob::PFParticle, recob::Vertex>>();
//    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
//    produces< art::Assns<recob::PFParticle, recob::Seed> >();
//    produces< art::Assns<recob::Seed,       recob::Hit> >();
    produces< art::Assns<recob::Cluster,    recob::Hit>>();
    produces< art::Assns<recob::Track,      recob::Hit>>();
//    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
}

SpecialTracks::~SpecialTracks()
{
    // Clean up dynamic memory and other resources here.
}

void SpecialTracks::produce(art::Event & event)
{
    // Instatiate the output
    std::unique_ptr< std::vector<recob::PFParticle>> pfParticleVector( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::Cluster>>    clusterVector(    new std::vector<recob::Cluster>    );
    std::unique_ptr< std::vector<recob::Track>>      trackVector(      new std::vector<recob::Track>      );
    std::unique_ptr< std::vector<recob::Vertex>>     vertexVector(     new std::vector<recob::Vertex>     );
    std::unique_ptr< std::vector<recob::Hit>>        hitVector(        new std::vector<recob::Hit>        );
    
    std::unique_ptr< art::Assns<recob::Cluster,    recob::Hit>>         clusterHitAssociations(   new art::Assns<recob::Cluster,    recob::Hit>     );
    std::unique_ptr< art::Assns<recob::Track,      recob::Hit>>         trackHitAssociations(     new art::Assns<recob::Track,      recob::Hit>     );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster>>     pfPartClusAssociations(   new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Track>>       pfPartTrackAssociations(  new art::Assns<recob::PFParticle, recob::Track>   );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex>>      pfPartVertexAssociations( new art::Assns<recob::PFParticle, recob::Vertex>  );
    
    // Start by recovering the PFParticle collection from art
    art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
    event.getByLabel(fPFParticleModuleLabel, pfParticleHandle);
    
    if (pfParticleHandle.isValid())
    {
        // We need the cluster handle to find associated hits
        art::Handle<std::vector<recob::Cluster>> clusterHandle;
        event.getByLabel(fPFParticleModuleLabel, clusterHandle);
        
        // Look up the tracks that were NOT produced by the PFParticle producer
        art::Handle<std::vector<recob::Track>> trackHandle;
        event.getByLabel(fTrackModuleLabel, trackHandle);
        
        // Recover the collection of associations between PFParticles and Tracks from the PFParticle producer
        art::FindManyP<recob::Track> pfTrackAssns(pfParticleHandle, event, fPFParticleModuleLabel);
        
        // Now get the track associations from the track producer
        art::FindManyP<recob::Track> kfTrackAssns(pfParticleHandle, event, fTrackModuleLabel);
        
        // Recover the collection of associations between PFParticles and vertices from the PFParticle producer
        art::FindManyP<recob::Vertex> vertexAssns(pfParticleHandle, event, fPFParticleModuleLabel);
        
        // Next we get the clusters associated to the PFParticle
        art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, event, fPFParticleModuleLabel);
        
        // Finally we get hit cluster associations
        art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, event, fPFParticleModuleLabel);
        
        // Recover the collection of associations between tracks and hits
        art::FindManyP<recob::Hit> trackHitAssns(trackHandle, event, fTrackModuleLabel);
        
        if (pfTrackAssns.isValid() && kfTrackAssns.isValid() > 0 && vertexAssns.isValid() && clusterAssns.isValid() && clusterHitAssns.isValid() && trackHitAssns.isValid())
        {
            // Ok, loop through the PFParticles
            for(size_t pfParticleIdx = 0; pfParticleIdx < pfParticleHandle->size(); pfParticleIdx++)
            {
                bool saveThisParticle(false);
                
                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);
                
                // Start by counting hits for this PFParticle...
                int nHitsPerView[] = {0,0,0};
                int nHitsTotal(0);
                int nClusPerView[] = {0,0,0};
                
                std::vector<art::Ptr<recob::Cluster>> clusterVec = clusterAssns.at(pfParticle.key());
                
                for(auto& cluster : clusterVec)
                {
                    std::vector<art::Ptr<recob::Hit>> hitVec = clusterHitAssns.at(cluster.key());
                    
                    if (hitVec.empty()) continue;
                    
                    int view = cluster->View();
                    
                    nClusPerView[view]++;
                    nHitsPerView[view]  = hitVec.size();
                    nHitsTotal         += hitVec.size();
                }
                
                // We don't think we care about the situation where we have a few number of hits
                if (nHitsTotal > 100)
                {
                    // Look up the associated track (if one)
                    std::vector<art::Ptr<recob::Track>> trackVec = kfTrackAssns.at(pfParticle.key());
                    
                    // If there is a track then let's see if it is somehow special
                    if (trackVec.size() == 1)
                    {
                        art::Ptr<recob::Track> track = trackVec.front();
                        
                        // Recover start/end points and directions
                        const TVector3& trackStart    = track->Vertex();
                        const TVector3& trackEnd      = track->End();
                        
                        double trackLen     = length(track.get());
                        double trackProjLen = projectedLength(track.get());
                        double trackEndLen  = (trackEnd - trackStart).Mag();
                        
                        if (fabs(trackLen - trackEndLen) > 25.)
                        {
                            std::cout << "--> Track id: " << track->ID() << " has len: " << trackLen << ", proj: " << trackProjLen << ", end: " << trackEndLen << std::endl;
                            std::vector<art::Ptr<recob::Hit>> trackHits = trackHitAssns.at(track.key());
                            std::cout << "    Track has " << trackHits.size() << " hits, " << track->NumberTrajectoryPoints() << " trajectory points" << std::endl;
                            saveThisParticle = true;
                        }
                        
                    }
                    // If there is not a track then we need to understand why
                    else if (trackVec.size() > 1) saveThisParticle = true;
                }
                else
                {
                    int nNoClus = std::count(nClusPerView,nClusPerView+3,0);
                    int nNoHits = std::count(nHitsPerView,nHitsPerView+3,0);
                    
                    if (nNoClus > 0 || nNoHits > 0) saveThisParticle = true;
                    saveThisParticle = false;
                }
            
                if (saveThisParticle)
                {
                    // Save copies of everything
                    pfParticleVector->push_back(*pfParticle);
                    
                    std::cout << "++++++>> Saving PFParticle key: " << pfParticle.key() << std::endl;
                    
                    // Make a map between the old and new hits...
                    std::map<const recob::Hit*,size_t> oldToNewMap;

                    // Keep track of starting positions
                    size_t clusterVectorOrigSize(clusterVector->size());
                    size_t trackVectorOrigSize(trackVector->size());
                    size_t vertexVectorOrigSize(vertexVector->size());
                    size_t hitVectorOrigSize(hitVector->size());
                    size_t hitIdx(hitVectorOrigSize);
                    
                    // copy clusters
                    for(const auto& cluster : clusterVec)
                    {
                        clusterVector->push_back(*cluster);
                        
                        std::vector<art::Ptr<recob::Hit>> clusterHitVec = clusterHitAssns.at(cluster.key());
                        
                        for(const auto& clusHit : clusterHitVec)
                        {
                            hitVector->push_back(*clusHit);
                            
                            oldToNewMap[clusHit.get()] = hitIdx++;
                        }
                        
                        util::CreateAssn(*this, event, *clusterVector, *hitVector, *clusterHitAssociations, hitVectorOrigSize, hitVector->size());
                        
                        hitVectorOrigSize += clusterHitVec.size();
                    }
                    
                    util::CreateAssn(*this, event, *pfParticleVector, *clusterVector, *pfPartClusAssociations, clusterVectorOrigSize, clusterVector->size());
                    
                    // fit tracks
                    std::vector<art::Ptr<recob::Track>> trackVec = kfTrackAssns.at(pfParticle.key());
                    
                    for(const auto& track : trackVec)
                    {
                        trackVector->push_back(*track);
                        
                        std::cout << "    ##>> Saving track id: " << track->ID() << ", associated with PFParticle: " << pfParticle.key() << std::endl;
                        
                        std::vector<art::Ptr<recob::Hit>> trackHitVec = trackHitAssns.at(track.key());
                        
                        // Use this to keep track of hits associated to the track
                        std::vector<size_t> hitIdxVec;
                        
                        for(const auto& hit : trackHitVec) hitIdxVec.push_back(oldToNewMap[hit.get()]);
                        
                        util::CreateAssn(*this, event, *trackVector, *hitVector, *trackHitAssociations, hitIdxVec);
                    }
                    
                    util::CreateAssn(*this, event, *pfParticleVector, *trackVector, *pfPartTrackAssociations, trackVectorOrigSize, trackVector->size());
                    
                    // vertices
                    std::vector<art::Ptr<recob::Vertex>> vertexVec = vertexAssns.at(pfParticle.key());
                    for(const auto& vertex : vertexVec) vertexVector->push_back(*vertex);
                    
                    util::CreateAssn(*this, event, *pfParticleVector, *vertexVector, *pfPartVertexAssociations, vertexVectorOrigSize, vertexVector->size());
                }
            }
            
            std::cout << "!! Outputting " << pfParticleVector->size() << " PFParticles with issues of " << pfParticleHandle->size() << " input" << std::endl;
        }
    }
    
    event.put( std::move(pfParticleVector)         );
    event.put( std::move(clusterVector)            );
    event.put( std::move(trackVector)              );
    event.put( std::move(vertexVector)             );
    event.put( std::move(hitVector)                );
    event.put( std::move(clusterHitAssociations)   );
    event.put( std::move(trackHitAssociations)     );
    event.put( std::move(pfPartClusAssociations)   );
    event.put( std::move(pfPartTrackAssociations)  );
    event.put( std::move(pfPartVertexAssociations) );
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void SpecialTracks::beginJob()
{
    double xyzStart[3],xyzEnd[3];
    fGeometry->WireEndPoints(geo::WireID(0,0,2,0),xyzStart,xyzEnd);
    
    std::cout << "==> Y wire 0 pos: " << xyzStart[1] << ", " << xyzStart[2] << std::endl;
    
    fGeometry->WireEndPoints(geo::WireID(0,0,2,1),xyzStart,xyzEnd);
    
    std::cout << "==> Y wire 1 pos: " << xyzStart[1] << ", " << xyzStart[2] << std::endl;
}

void SpecialTracks::reconfigure(fhicl::ParameterSet const & p)
{
    // Implementation of optional member function here.
    fPFParticleModuleLabel  = p.get< std::string >("PFParticleModuleLabel", "pandoraCosmic");
    fHitProducerLabel       = p.get< std::string >("HitProducerLabel",      "gaushit");
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel",      "pandoraCosmicKHit");
}

void SpecialTracks::endJob() {
  // Implementation of optional member function here.
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double SpecialTracks::length(const recob::Track* track)
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
        
        if (lastDir.Dot(trajDir) >= 0.)
        {
            disp   -= pos;
            result += disp.Mag();
            disp    = pos;
        }
        
        lastPoint = pos;
        lastDir   = trajDir;
    }
    
    return result;
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double SpecialTracks::projectedLength(const recob::Track* track)
{
    double   result(0.);
    TVector3 lastPoint(track->LocationAtPoint(0));
    TVector3 lastDir(track->DirectionAtPoint(0));
    int      n(track->NumberTrajectoryPoints());
    
    for(int i = 1; i < n; ++i)
    {
        const TVector3& newPoint = track->LocationAtPoint(i);
        
        TVector3 lastToNewPoint = newPoint - lastPoint;
        double   arcLenToDoca   = lastDir.Dot(lastToNewPoint);
        
        result    += arcLenToDoca;
        lastPoint  = lastPoint + arcLenToDoca * lastDir;
        lastDir    = track->DirectionAtPoint(i);
    }
    
    return result;
}

int SpecialTracks::traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
                                             size_t                                       pfParticleIdx,
                                             const art::FindManyP<recob::Track>&          trackAssns,
                                             const art::FindManyP<recob::Vertex>&         vertexAssns,
                                             int&                                         nTracks,
                                             int&                                         nVertices) const
{
    // So far no daughters...
    int nDaughters(0);
    
    // Get pointer to PFParticle
    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);
    
    // Recover tracks/vertices associated to this PFParticle
    std::vector<art::Ptr<recob::Track>>  pfPartTrackVec  = trackAssns.at(pfParticle.key());
    std::vector<art::Ptr<recob::Vertex>> pfPartVertexVec = vertexAssns.at(pfParticle.key());
    
    nTracks    += pfPartTrackVec.size();
    nVertices  += pfPartVertexVec.size();
    nDaughters += pfParticle->Daughters().size();
    
    for(auto& daughterIdx : pfParticle->Daughters())
    {
        nDaughters += traversePFParticleHierarchy(pfParticleHandle, daughterIdx, trackAssns, vertexAssns, nTracks, nVertices);
    }
    
    return nDaughters;
}


DEFINE_ART_MODULE(SpecialTracks)
