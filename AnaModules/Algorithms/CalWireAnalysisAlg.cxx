
#include "CalWireAnalysisAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include <cmath>
#include <algorithm>

namespace calwireanalysis
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CalWireAnalysisAlg::CalWireAnalysisAlg(fhicl::ParameterSet const & pset) 
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("CalWireAnalysisAlg") << "CalWireAnalysisAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
CalWireAnalysisAlg::~CalWireAnalysisAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void CalWireAnalysisAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fLocalDirName = pset.get<std::string   >("LocalDirName", std::string("wow"));
    fFirstBin     = pset.get<raw::TDCtick_t>("FirstBin",                      0);
    fEndBin       = pset.get<raw::TDCtick_t>("EndBin",                     6400);
    fThreshold    = pset.get<float         >("Threshold",                    8.);
}

//----------------------------------------------------------------------------
/// Begin job method.
void CalWireAnalysisAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    fNumWiresView[0]    = dir.make<TH1D>("NumWiresView0",  ";#",     240,    0.,  2400.);
    fNumWiresView[1]    = dir.make<TH1D>("NumWiresView1",  ";#",     240,    0.,  2400.);
    fNumWiresView[2]    = dir.make<TH1D>("NumWiresView2",  ";#",     340,    0.,  3400.);
    fWireValueHist[0]   = dir.make<TH1D>("WireValueHist0", ";ADC",   200,  -25.,    25.);
    fWireValueHist[1]   = dir.make<TH1D>("WireValueHist1", ";ADC",   200,  -25.,    25.);
    fWireValueHist[2]   = dir.make<TH1D>("WireValueHist2", ";ADC",   200,  -25.,    25.);
    fROISize[0]         = dir.make<TH1D>("ROISize0",       ";size",  600,    0.,   600.);
    fROISize[1]         = dir.make<TH1D>("ROISize1",       ";size",  400,    0.,   400.);
    fROISize[2]         = dir.make<TH1D>("ROISize2",       ";size",  400,    0.,   400.);
    fNumAboveThres[0]   = dir.make<TH1D>("NumAboveThres0", ";#",     100,    0.,   100.);
    fNumAboveThres[1]   = dir.make<TH1D>("NumAboveThres1", ";#",     100,    0.,   100.);
    fNumAboveThres[2]   = dir.make<TH1D>("NumAboveThres2", ";#",     100,    0.,   100.);
    fNumHitsROI[0]      = dir.make<TH1D>("NumHitsROI0",    ";#",      40,    0.,    40.);
    fNumHitsROI[1]      = dir.make<TH1D>("NumHitsROI1",    ";#",      40,    0.,    40.);
    fNumHitsROI[2]      = dir.make<TH1D>("NumHitsROI2",    ";#",      40,    0.,    40.);
    fNumROIsWire[0]     = dir.make<TH1D>("NumROIsWire0",   ";#",      50,    0.,    50.);
    fNumROIsWire[1]     = dir.make<TH1D>("NumROIsWire1",   ";#",      50,    0.,    50.);
    fNumROIsWire[2]     = dir.make<TH1D>("NumROIsWire2",   ";#",      50,    0.,    50.);
    fNumHitsView[0]     = dir.make<TH1D>("NumHitsView0",   ";#",     500,    0., 50000.);
    fNumHitsView[1]     = dir.make<TH1D>("NumHitsView1",   ";#",     500,    0., 50000.);
    fNumHitsView[2]     = dir.make<TH1D>("NumHitsView2",   ";#",     500,    0., 50000.);

    return;
}
    
void CalWireAnalysisAlg::fillHistograms(const WirePtrVec& wirePtrVec, const art::FindManyP<recob::Hit>& hitWireAssns) const
{
    // Keep track of number of hits per view
    size_t nHitsPerView[]  = {0,0,0};
    size_t nWiresPerView[] = {0,0,0};
    double maxROISize[]    = {599.9,399.9,399.9};
    
    // Loop the hits and make some plots
    for(const auto& wirePtr : wirePtrVec)
    {
        // --- Setting Channel Number and Signal type ---
        size_t channel = wirePtr->Channel();
        
        // get the WireID for this hit
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
        
        // for now, just take the first option returned from ChannelToWire
        size_t view = wids[0].Plane;

        const recob::Wire::RegionsOfInterest_t& signalROI = wirePtr->SignalROI();
        
        int numROIs(0);
        
        for(const auto& range : signalROI.get_ranges())
        {
            const std::vector<float>& signal = range.data();
            
            // ROI start time
            raw::TDCtick_t roiFirstBinTick = range.begin_index();
            
            if (roiFirstBinTick <= fFirstBin) continue;
            if (roiFirstBinTick >  fEndBin  ) continue;
            
            numROIs++;
            
            int numAboveThreshold(0);
            
            for(const auto& binValue : signal)
            {
                fWireValueHist[view]->Fill(std::max(float(-24.9),std::min(float(24.9),binValue)),1.);
                if (binValue > fThreshold) numAboveThreshold++;
                if (roiFirstBinTick >= fEndBin) break;
                roiFirstBinTick++;
            }
            
            fROISize[view]->Fill(std::min(maxROISize[view],double(signal.size())), 1.);
            fNumAboveThres[view]->Fill(std::min(99.9,double(numAboveThreshold)), 1.);
        }
        
        if (numROIs > 0 && hitWireAssns.isValid())
        {
            // How many hits associated to this ROI?
            std::vector<art::Ptr<recob::Hit>> hitRoiVec = hitWireAssns.at(wirePtr.key());
        
            nHitsPerView[view] += hitRoiVec.size();
        
            fNumHitsROI[view]->Fill(std::min(39.9,double(hitRoiVec.size())/numROIs), 1.);
            fNumROIsWire[view]->Fill(std::min(49.9,double(numROIs)), 1.);
            
            nWiresPerView[view]++;
        }
    }
    
    for (size_t view = 0; view < 3; view++)
    {
        fNumHitsView[view]->Fill(nHitsPerView[view], 1.);
        fNumWiresView[view]->Fill(nWiresPerView[view], 1.);
    }
    
    return;
}
    
// Useful for normalizing histograms
void CalWireAnalysisAlg::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
//    double normFactor(1./numEvents);
    
//    for(size_t idx = 0; idx < 3; idx++) fHitsByWire[idx]->Scale(normFactor);
    
    return;
}

}
