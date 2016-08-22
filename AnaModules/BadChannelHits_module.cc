////////////////////////////////////////////////////////////////////////
// Class:       BadChannelHits
// Module Type: producer
// File:        BadChannelHits_module.cc
//              The goal of this module is to try to filter out noise hits
//
// Generated at Tues Feb 23 19:17:00 2016 by Tracy Usher by cloning CosmicTrackTagger
// from art v1_02_02.
// artmod -e beginJob -e reconfigure -e endJob producer trkf::BadChannelHits
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
//#include "larcore/Geometry/geo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"

class BadChannelHits : public art::EDProducer
{
public:
    explicit BadChannelHits(fhicl::ParameterSet const & p);
    virtual ~BadChannelHits();

    void produce(art::Event & e) override;

    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void endJob() override;

private:
    
    int   fMinAllowedChanStatus; ///< Don't consider channels with lower status
    
    // Other variables that will be shared between different methods.
    const geo::GeometryCore* fGeometry;       // pointer to Geometry service
};


BadChannelHits::BadChannelHits(fhicl::ParameterSet const & p)
{
    fGeometry = lar::providerFrom<geo::Geometry>();
    
    this->reconfigure(p);

    // Call appropriate Produces<>() functions here.
    produces< std::vector<recob::Hit>>();
}

BadChannelHits::~BadChannelHits()
{
    // Clean up dynamic memory and other resources here.
}

void BadChannelHits::produce(art::Event & event)
{
    // Instatiate the output
    std::unique_ptr< std::vector<recob::Hit>> hitVector( new std::vector<recob::Hit> );
    
    // Recover the channel status service
    const lariov::ChannelStatusProvider& chanFilt           = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
    const detinfo::DetectorProperties&   detectorProperties = *lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    raw::TDCtick_t start_tick(0);
    raw::TDCtick_t end_tick(detectorProperties.NumberTimeSamples());
    
    float     deltaTime            = end_tick - start_tick;
    float     peak_time            = 0.5 * deltaTime;
    float     sigma_peak_time      = deltaTime / std::sqrt(12.);
    float     rms                  = sigma_peak_time;
    float     peak_amplitude       = 100.;
    float     sigma_peak_amplitude = 10.;
    float     summedADC            = 10000.;
    float     hit_integral         = 10000.;
    float     hit_sigma_integral   = 100.;
    short int multiplicity         = 1;
    short int local_index          = 0;
    float     goodness_of_fit      = -10000.;
    int       dof                  = 1;
    
    // We'll loop over views and wires per view to build a channel ID which we can input to the channel status service
    for(size_t viewIdx = 0; viewIdx < fGeometry->Nviews(); viewIdx++)
    {
        for(size_t wireIdx = 0; wireIdx < fGeometry->Nwires(viewIdx); wireIdx++)
        {
            raw::ChannelID_t channel = fGeometry->PlaneWireToChannel(viewIdx, wireIdx);
            
            // make sure a valid channel
            if (!chanFilt.IsPresent(channel)) continue;
            
            // Is this a "bad" channel?
            if (chanFilt.Status(channel) < fMinAllowedChanStatus)
            {
                // get the WireID for this hit
                std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
                // for now, just take the first option returned from ChannelToWire
                geo::WireID wid  = wids[0];
                geo::View_t view = fGeometry->View(wid.Plane);
                
                // ok, now create the hit
                hitVector->emplace_back(recob::Hit(channel,
                                                   start_tick,
                                                   end_tick,
                                                   peak_time,
                                                   sigma_peak_time,
                                                   rms,
                                                   peak_amplitude,
                                                   sigma_peak_amplitude,
                                                   summedADC,
                                                   hit_integral,
                                                   hit_sigma_integral,
                                                   multiplicity,
                                                   local_index,
                                                   goodness_of_fit,
                                                   dof,
                                                   view,
                                                   geo::kMysteryType,
                                                   wid));
            }
        }
    }
    
    event.put( std::move(hitVector) );
    
    return;

} // end of produce
//////////////////////////////////////////////////////////////////////////////////////////////////////

void BadChannelHits::beginJob()
{
}

void BadChannelHits::reconfigure(fhicl::ParameterSet const & p)
{
    // Implementation of optional member function here.
    fMinAllowedChanStatus = p.get< int >("MinAllowedChannelStatus", 3);
}

void BadChannelHits::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(BadChannelHits)
