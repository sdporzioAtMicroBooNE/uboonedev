////////////////////////////////////////////////////////////////////////
//
// Class:       NoiseFilterHits
// Module Type: producer
// File:        NoiseFilterHits_module.cc
//
//              The aim is to try to remove noise hits and produce a clean list of hits
//
// Configuration parameters:
//
// HitProducerLabel        - the producer of the recob::Hit objects
//
// Created by Tracy Usher (usher@slac.stanford.edu) on September 18, 2014
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/TimeService.h"
#include "Utilities/SimpleTimeService.h"

class Propagator;

class NoiseFilterHits : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit NoiseFilterHits(fhicl::ParameterSet const & pset);
    virtual ~NoiseFilterHits();

    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:

    // Fcl parameters.
    std::string              fHitProducerLabel;        ///< The full collection of hits

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
};

DEFINE_ART_MODULE(NoiseFilterHits)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
NoiseFilterHits::NoiseFilterHits(fhicl::ParameterSet const & pset) :
  fNumEvent(0)
{
    reconfigure(pset);
    produces<std::vector<recob::Hit> >();

    // Report.
    mf::LogInfo("NoiseFilterHits") << "NoiseFilterHits configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
NoiseFilterHits::~NoiseFilterHits()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void NoiseFilterHits::reconfigure(fhicl::ParameterSet const & pset)
{
    fHitProducerLabel = pset.get<std::string>("HitProducerLabel");
}

//----------------------------------------------------------------------------
/// Begin job method.
void NoiseFilterHits::beginJob()
{
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method. The goal is to produce a list of recob::Hit
/// objects which are a "clean" subset of all hits and which are believed to
/// be due to a neutrino interaction. It does this by considering input CosmicTag
/// objects, relating them to PFParticles/Tracks and removing the hits
/// associated to those objects which are believed to be Cosmic Rays.
///
void NoiseFilterHits::produce(art::Event & evt)
{
    ++fNumEvent;
    
    // Start by looking up the original hits
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitProducerLabel, hitHandle);
    
    // If there are hits then we are going to output something so get a new
    // output hit vector
    std::unique_ptr<std::vector<recob::Hit> > outputHits(new std::vector<recob::Hit>);
    
    // If there are no hits then there should be no output
    if (hitHandle.isValid())
    {
        // Loop through the hits
        for(const auto& hit : *hitHandle)
        {
            bool   keepHit(false);
            double pulseHeight(hit.PeakAmplitude());
            double RMS(hit.RMS());
            double roiLength = std::min(double(199.), double(hit.EndTick() - hit.StartTick()));
            
            if (hit.View() == geo::kU)
            {
                double cutVal = pulseHeight / 3.;
                
                if (RMS <= cutVal) keepHit = true;
            }
            else if (hit.View() == geo::kV)
            {
                double cutVal = pulseHeight / 3.;
                
                if (RMS <= cutVal) keepHit = true;
            }
            else if (hit.View() == geo::kW)
            {
                double cutVal = 2.00;
                //double cutVal = 1.75;
                
                //if (pulseHeight < 15) cutVal = -pulseHeight / 3. + 5.;
                
                if (RMS < cutVal) keepHit = true;
            }
            
            keepHit = roiLength <= 12. && pulseHeight <= 15.;
            keepHit = RMS > pulseHeight;
            
            if (keepHit) outputHits->push_back(hit);
        }
    }
    
    std::cout << "NoiseFilterHits - size of input hit collection: " << hitHandle->size() << ", size of output collection: " << outputHits->size() << std::endl;
    
    // Add tracks and associations to event.
    evt.put(std::move(outputHits));
}

//----------------------------------------------------------------------------
/// End job method.
void NoiseFilterHits::endJob()
{
}
