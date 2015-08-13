////////////////////////////////////////////////////////////////////////
// Class:       CompositeCosmicTagger
// Module Type: producer
// File:        CompositeCosmicTagger_module.cc
//              The goal of this module is to synthesize a Cosmic Tag from
//              the set of existing CosmicTag's for a given object.
//              Basically, we start with a PFParticle and get the associated
//              Tracks and with this information we can collect the CosmicTag
//              objects and then see if we can find a "better" overall tag
//
// Generated at Wed Sep 17 19:17:00 2014 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::CompositeCosmicTagger
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

#include "Geometry/Geometry.h"
#include "Geometry/geo.h"

#include "RecoBase/PFParticle.h"
#include "RecoBase/PCAxis.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"

#include "AnalysisBase/CosmicTag.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"

#include "TVector3.h"


namespace cosmic
{
    class CompositeCosmicTagger;
}

class cosmic::CompositeCosmicTagger : public art::EDProducer
{
public:
    explicit CompositeCosmicTagger(fhicl::ParameterSet const & p);
    virtual ~CompositeCosmicTagger();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

private:
    // Forward declaration of object to handle sorting of PFParticles associated to Track IDs
    class SortTrackVec;
    
    std::string fPFParticleModuleLabel;
    std::string fPCAxisModuleLabel;
    std::string fTrackModuleLabel;
    
    std::string fPFCosmicTagModuleLabel;
    std::string fTKCosmicTagModuleLabel;
    std::string fTKFlashTagModuleLabel;
};


cosmic::CompositeCosmicTagger::CompositeCosmicTagger(fhicl::ParameterSet const & p)
{
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<anab::CosmicTag> >();
    produces< art::Assns<recob::PFParticle, anab::CosmicTag> >();
}

cosmic::CompositeCosmicTagger::~CompositeCosmicTagger()
{
    // Clean up dynamic memory and other resources here.
}

// This will be used to sort PFParticles in the map below
class cosmic::CompositeCosmicTagger::SortTrackVec
{
    /**
     * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
     */
public:
    SortTrackVec() {}
    
    bool operator()(const art::Ptr<recob::Track>& left, const art::Ptr<recob::Track>& right)
    {
//        size_t numHitsLeft  = fPFPartCntMap.at(left);
//        size_t numHitsRight = fPFPartCntMap.at(right);
        
        return left->NumberTrajectoryPoints() > right->NumberTrajectoryPoints();
    }
private:
//    const PFParticleMcAna::PFParticleHitCntMap& fPFPartCntMap;
};

void cosmic::CompositeCosmicTagger::produce(art::Event & evt)
{
    // Instatiate the output
    std::unique_ptr< std::vector< anab::CosmicTag >>                  cosmicTagPFParticleVector(  new std::vector<anab::CosmicTag> );
    std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag >> assnOutCosmicTagPFParticle( new art::Assns<recob::PFParticle, anab::CosmicTag>);
    
    // Recover handle for PFParticles
    art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
    evt.getByLabel( fPFParticleModuleLabel, pfParticleHandle);
    
    if (!pfParticleHandle.isValid())
    {
        evt.put( std::move(cosmicTagPFParticleVector)  );
        evt.put( std::move(assnOutCosmicTagPFParticle) );
        return;
    }
    
    // Recover the list of associated PCA axes
    art::FindManyP<recob::PCAxis> pfPartToPCAxisAssns(pfParticleHandle, evt, fPCAxisModuleLabel);
    
    // We need a handle to the Track collection so we can recover associated CosmicTags
    art::Handle<std::vector<recob::Track>> trackHandle;
    evt.getByLabel(fTrackModuleLabel, trackHandle);
    
    // Recover the list of associated Tracks
    art::FindManyP<recob::Track> pfPartToTrackAssns(pfParticleHandle, evt, fTrackModuleLabel);
    
    // Recover relations between PFParticles and CosmicTag objects
    art::FindManyP<anab::CosmicTag> pfCosmicTagAssnVec(pfParticleHandle, evt, fPFCosmicTagModuleLabel);
    
    // Recover relations between Tracks and CosmicTag objects
    art::FindManyP<anab::CosmicTag> tkCosmicTagAssnVec(trackHandle, evt, fTKCosmicTagModuleLabel);
    art::FindManyP<anab::CosmicTag> tkFlashTagAssnVec( trackHandle, evt, fTKFlashTagModuleLabel);
    
    // The outer loop is going to be over PFParticles
    for(size_t pfPartIdx = 0; pfPartIdx != pfParticleHandle->size(); pfPartIdx++)
    {
        art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfPartIdx);
        
        // Recover the PCAxis vector
        std::vector<art::Ptr<recob::PCAxis> > pcAxisVec = pfPartToPCAxisAssns.at(pfParticle.key());
        
        // Is there an axis associated to this PFParticle?
        if (pcAxisVec.empty()) continue;
        
        // *****************************************************************************************
        // For what follows below we want the "best" PCAxis object only. However, it can be that
        // there are two PCAxes for a PFParticle (depending on source) where the "first" axis will
        // be the "better" one that we want (this statement by fiat, it is defined that way in the
        // axis producer module).
        if (pcAxisVec.size() > 1 && pcAxisVec.front()->getID() > pcAxisVec.back()->getID()) std::reverse(pcAxisVec.begin(), pcAxisVec.end());
        // We need to confirm this!!
        // *****************************************************************************************
        
        // Recover the axis
        const art::Ptr<recob::PCAxis>& pcAxis = pcAxisVec.front();

        // The purpose of this section is to find end points to assign a CosmicTag object in the event
        // we don't have one associated to this PFParticle. This is a highly unlikely situation but...
        double eigenVal0 = sqrt(pcAxis->getEigenValues()[0]);
        double maxArcLen = 3. * eigenVal0;
        
        // Recover PCA end points
        TVector3 vertexPosition(pcAxis->getAvePosition()[0],pcAxis->getAvePosition()[1],pcAxis->getAvePosition()[2]);
        TVector3 vertexDirection(pcAxis->getEigenVectors()[0][0],pcAxis->getEigenVectors()[0][1],pcAxis->getEigenVectors()[0][2]);
        
        TVector3 pcAxisStart = vertexPosition - maxArcLen * vertexDirection;
        TVector3 pcAxisEnd   = vertexPosition + maxArcLen * vertexDirection;
        
        // "Track" end points in easily readable form
        float trackEndPt1_X = pcAxisStart[0];
        float trackEndPt1_Y = pcAxisStart[1];
        float trackEndPt1_Z = pcAxisStart[2];
        float trackEndPt2_X = pcAxisEnd[0];
        float trackEndPt2_Y = pcAxisEnd[1];
        float trackEndPt2_Z = pcAxisEnd[2];
        
        // Now on to the business of forming a composite tag!
        // Recover the CosmicTag associated to this PFParticle
        std::vector<art::Ptr<anab::CosmicTag>> pfCosmicTagVec = pfCosmicTagAssnVec.at(pfParticle.key());
        
        // What we think this is...
        float               cosmicScore(0.);
        anab::CosmicTagID_t tagID(anab::kNotTagged);
        
        if (!pfCosmicTagVec.empty())
        {
            // If there exists a CosmicTag for this PFParticle the we automatically
            // assign the score and tag values from that object
            cosmicScore = pfCosmicTagVec.front()->CosmicScore();
            tagID       = pfCosmicTagVec.front()->CosmicType();
            
            // If this particle has hits out of time it is pointless to do anything else
            if (tagID < 100)
            {
                // Recover tracks associated to this PFParticle
                std::vector<art::Ptr<recob::Track>> trackVec = pfPartToTrackAssns.at(pfParticle.key());
            
                // If there are no tracks then nothing to do
                if (!trackVec.empty())
                {
                    // Keep track of lowest score
                    float lowCosmicTagScore(100.);
                    float lowFlashTagScore(100.);
                    
                    // Let's sort the tracks so we can be sure the longest is the first
                    std::sort(trackVec.begin(), trackVec.end(), SortTrackVec());
                    
                    // Recover pointing of first track
                    TVector3 firstTrackStartPos = trackVec.front()->Vertex();
                    TVector3 firstTrackEndPos   = trackVec.front()->End();
                    TVector3 firstTrackStartDir = trackVec.front()->VertexDirection();
                    
                    int firstTrackHitCnt(trackVec.front()->NumberTrajectoryPoints());
                    
                    // poor man's vertexing... just keep track of two closest tracks
                    double trackMatchDist(1000.);
                    double cosThetaMatch(1.);
                    int    trackCnt(0);
                    
                    int minHitCnt = std::max(10, std::min(75,int(firstTrackHitCnt/5)));
                    
                    // Loop over tracks
                    for(auto& track : trackVec)
                    {
                        // Tracks are ordered so if we get a run we're done
                        int curTrackHitCnt(track->NumberTrajectoryPoints());
                        
                        if (curTrackHitCnt < minHitCnt) break;
                        
                        if (trackCnt++)
                        {
                            TVector3 endPointDiff = track->Vertex() - firstTrackStartPos;
                            
                            if (endPointDiff.Mag() < trackMatchDist)
                            {
                                trackMatchDist = endPointDiff.Mag();
                                cosThetaMatch  = fabs(firstTrackStartDir.Dot(track->VertexDirection()));
                            }
                            
                            endPointDiff = track->Vertex() - firstTrackEndPos;
                            
                            if (endPointDiff.Mag() < trackMatchDist)
                            {
                                trackMatchDist = endPointDiff.Mag();
                                cosThetaMatch  = fabs(firstTrackStartDir.Dot(track->VertexDirection()));
                            }
                            
                            endPointDiff = track->End() - firstTrackStartPos;
                            
                            if (endPointDiff.Mag() < trackMatchDist)
                            {
                                trackMatchDist = endPointDiff.Mag();
                                cosThetaMatch  = fabs(firstTrackStartDir.Dot(track->VertexDirection()));
                            }
                            
                            endPointDiff = track->End() - firstTrackEndPos;
                            
                            if (endPointDiff.Mag() < trackMatchDist)
                            {
                                trackMatchDist = endPointDiff.Mag();
                                cosThetaMatch  = fabs(firstTrackStartDir.Dot(track->VertexDirection()));
                            }
                        }
                        
                        // Recover CosmicTag and process
                        std::vector<art::Ptr<anab::CosmicTag>> tkCosmicTagVec = tkCosmicTagAssnVec.at(track.key());
                
                        if (!tkCosmicTagVec.empty())
                        {
                            float tkCosmicScore = tkCosmicTagVec.front()->CosmicScore();
                            
                            if (tkCosmicScore < 1.) lowCosmicTagScore = std::min(lowCosmicTagScore, tkCosmicScore);
                            else
                            {
                                lowCosmicTagScore = tkCosmicScore;
                                break;
                            }
                        }
                
                        // Recover CosmicTag and process
                        std::vector<art::Ptr<anab::CosmicTag>> tkFlashTagVec = tkFlashTagAssnVec.at(track.key());
                
                        if (!tkFlashTagVec.empty())
                        {
                            lowFlashTagScore = std::min(lowFlashTagScore, tkFlashTagVec.front()->CosmicScore());
                        }
                    }
                    
                    // At this point we can try to evaluate things...
                    // First consider when PFParticle tag is 0.5
                    //if (cosmicScore > 0.4 && cosmicScore < 1.0)
                    //{
                    //    if (lowFlashTagScore < 1.)
                    //    {
                    //        std::cout << "Track: " << pfParticle.key() << " has cosmicScore: " << cosmicScore << ", flash: " << lowFlashTagScore << std::endl;
                    //        cosmicScore = 0.1;
                    //    }
                    //}
                    if (cosmicScore > 0.4)
                    {
                        std::cout << "Track: " << pfParticle.key() << " has cosmicScore: " << cosmicScore;
                        
                        if (lowCosmicTagScore < 0.5) cosmicScore = 0.2;
                        else if (trackCnt > 1 && trackMatchDist < 10. && cosThetaMatch < 0.97 && lowCosmicTagScore < 1.)
                            cosmicScore = 0.3;
                        
                        std::cout << ", # tracks: " << trackVec.size() << ", lowCosmicTagScore: " << lowCosmicTagScore << ", new: " << cosmicScore << std::endl;
                    }
                }
            }
        }
        
        std::vector<float> endPt1;
        std::vector<float> endPt2;
        endPt1.push_back( trackEndPt1_X );
        endPt1.push_back( trackEndPt1_Y );
        endPt1.push_back( trackEndPt1_Z );
        endPt2.push_back( trackEndPt2_X );
        endPt2.push_back( trackEndPt2_Y );
        endPt2.push_back( trackEndPt2_Z );
        
        // Create the tag object for this PFParticle and make the corresponding association
        cosmicTagPFParticleVector->emplace_back( endPt1, endPt2, cosmicScore, tagID);
        
        util::CreateAssn(*this, evt, *cosmicTagPFParticleVector, pfParticle, *assnOutCosmicTagPFParticle );
    }
    
    evt.put( std::move(cosmicTagPFParticleVector)  );
    evt.put( std::move(assnOutCosmicTagPFParticle) );
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void cosmic::CompositeCosmicTagger::beginJob()
{
}

void cosmic::CompositeCosmicTagger::reconfigure(fhicl::ParameterSet const & p)
{
    // Implementation of optional member function here.
    fPFParticleModuleLabel  = p.get< std::string >("PFParticleModuleLabel");
    fPCAxisModuleLabel      = p.get< std::string >("PCAxisModuleLabel");
    fTrackModuleLabel       = p.get< std::string >("TrackModuleLabel");
    
    fPFCosmicTagModuleLabel = p.get< std::string >("PFCosmicTagModuleLabel");
    fTKCosmicTagModuleLabel = p.get< std::string >("TKCosmicTagModuleLabel");
    fTKFlashTagModuleLabel  = p.get< std::string >("TKFlashTagModuleLabel");
}

void cosmic::CompositeCosmicTagger::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(cosmic::CompositeCosmicTagger)
