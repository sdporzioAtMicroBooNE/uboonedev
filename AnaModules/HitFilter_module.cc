////////////////////////////////////////////////////////////////////////
// Class:       HitFilter
// Module Type: producer
// File:        HitFilter_module.cc
//              The goal of this module is to try to filter out noise hits
//
// Generated at Tues Feb 23 19:17:00 2016 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::HitFilter
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/geo.h"

#include "lardata/RecoBase/Hit.h"

class HitFilter : public art::EDProducer
{
public:
    explicit HitFilter(fhicl::ParameterSet const & p);
    virtual ~HitFilter();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

private:
    
    std::string        fHitProducerLabel;
    std::vector<float> fMinPulseHeight;
    std::vector<float> fMinPulseSigma;
    
    // Other variables that will be shared between different methods.
    const geo::GeometryCore* fGeometry;       // pointer to Geometry service
};


HitFilter::HitFilter(fhicl::ParameterSet const & p)
{
    fGeometry = lar::providerFrom<geo::Geometry>();
    
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<recob::Hit>>();
}

HitFilter::~HitFilter()
{
    // Clean up dynamic memory and other resources here.
}

void HitFilter::produce(art::Event & event)
{
    // Instatiate the output
    std::unique_ptr< std::vector<recob::Hit>> hitVector( new std::vector<recob::Hit> );
    
    // Start by recovering the PFParticle collection from art
    art::Handle<std::vector<recob::Hit>> inputHitHandle;
    event.getByLabel(fHitProducerLabel, inputHitHandle);
    
    if (inputHitHandle.isValid())
    {
        // Loop through hits
        for(const auto& hit : *inputHitHandle)
        {
            // Extract the hit parameters we think we need
            const geo::WireID& wireID   = hit.WireID();
            float              hitPH    = std::min(hit.PeakAmplitude(),float(249.8));
            float              hitSigma = hit.RMS();
            size_t             view     = wireID.Plane;
            
            if (!(hitPH < fMinPulseHeight[view] && hitSigma < fMinPulseSigma[view]))
                hitVector->emplace_back(recob::Hit(hit));
        }
    }
    
    event.put( std::move(hitVector) );
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void HitFilter::beginJob()
{
}

void HitFilter::reconfigure(fhicl::ParameterSet const & p)
{
    // Implementation of optional member function here.
    fHitProducerLabel = p.get< std::string        >("HitProducerLabel", "gaushit");
    fMinPulseHeight   = p.get< std::vector<float> >("MinPulseHeight",   std::vector<float>() = {14., 15., 30.});
    fMinPulseSigma    = p.get< std::vector<float> >("MinPulseSigma",    std::vector<float>() = {1.8, 2.0, 1.4});
}

void HitFilter::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(HitFilter)
