// PMTPickupAna_module.cc
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them. 

// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.

// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.

// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.

#ifndef PMTPickupAna_Module
#define PMTPickupAna_Module

// LArSoft includes
#include "Geometry/Geometry.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Track.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"

//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "Utilities/TimeService.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TVirtualFFT.h"

// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>

namespace PMTPickupAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class PMTPickupAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit PMTPickupAna(fhicl::ParameterSet const& pset);
    virtual ~PMTPickupAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();
    void endJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

private:
    // The following typedefs will, obviously, be useful
    typedef std::vector<const recob::Hit*>                                    HitVec;
    typedef std::vector<HitVec>                                               ViewHitVec;
    typedef std::pair<size_t,size_t>                                          WirePair;
    typedef std::vector<std::vector<std::tuple<float,float,WirePair,HitVec>>> PeakTupleVec;
    typedef std::tuple<double,art::Ptr<recob::Cluster>,const recob::Hit*>     ClusterHitMatchTuple;
    
    void buildPeakDataStructure(ViewHitVec&   viewHitVec,
                                PeakTupleVec& peakTupleVec) const;

    // The parameters we'll read from the .fcl file.
    std::string fDigitModuleLabel;           // The name of the producer that created raw digits
    std::string fHitProducerLabel;
    std::string fClusterProducerLabel;       //
    std::string fTrackProducerLabel;
    bool        fWritePedestals;             // Output new file of pedestals
    float       fTruncMeanFraction;          // Fraction for truncated mean

    // Pointers to the histograms we'll create.
    TH1D*     fHitsByWire[3];
    TProfile* fHitsByWireProf[3];
    TH1D*     fNumSpikesEvent[3];
    TH1D*     fNumHitsSpike[3];
    TH1D*     fDeltaTSpike[3];
    TH1D*     fSpreadTSpike[3];
    TH1D*     fRMSSpike[3];
    TH1D*     fPulseHeight[3];
    TH1D*     fRMS[3];
    TH1D*     fSigmaPeakTime[3];
    TH1D*     fHitIndex[3];
    TH1D*     fROIlength[3];
    TH1D*     fChi2[3];
    TH2D*     fPHvsRMS[3];
    TH2D*     fPHvsSigma[3];
    TH2D*     fROIvsPH[3];
    TH2D*     fROIvsChi[3];
    
    TH1D*     fTimeHitCount[3];
    
    TH2D*     fWireVsTime[3];
    TH2D*     fWireVsPH[3];
    TH2D*     fWireVsPHRat[3];

    TH1D*     fNumMatchedClusters[3];
    TH1D*     fNumMatchedTracks[3];
    TH1D*     fDeltaClusterT[3];
    TH1D*     fDeltaTrackT[3];
    TH1D*     fTrackComplete[3];
    TH2D*     fTrackPos[3];
    TH2D*     fTrackCompPos[3];

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fNumEvents;
    
    std::vector<std::vector<double>> fChannelPedVec;

    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry>            fGeometry;       // pointer to Geometry service
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;
    const lariov::IDetPedestalProvider&          fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg

    double                                       fElectronsToGeV; // conversion factor

}; // class PMTPickupAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
PMTPickupAna::PMTPickupAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet),
      fPedestalRetrievalAlg(art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider())

{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
PMTPickupAna::~PMTPickupAna()
{}
   
//-----------------------------------------------------------------------
void PMTPickupAna::beginJob()
{
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    //double detectorLength = fGeometry->DetLength();

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes. 

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fHitsByWire[0]         = tfs->make<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]         = tfs->make<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]         = tfs->make<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    fHitsByWireProf[0]     = tfs->make<TProfile>("HitsByWireProf0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0), 0., 250.);
    fHitsByWireProf[1]     = tfs->make<TProfile>("HitsByWireProf1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1), 0., 250.);
    fHitsByWireProf[2]     = tfs->make<TProfile>("HitsByWireProf2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2), 0., 250.);
    fNumSpikesEvent[0]     = tfs->make<TH1D>("NSpikesEvent0",  "; # spikes",      25,    0.,  25.);
    fNumSpikesEvent[1]     = tfs->make<TH1D>("NSpikesEvent1",  "; # spikes",      25,    0.,  25.);
    fNumSpikesEvent[2]     = tfs->make<TH1D>("NSpikesEvent2",  "; # spikes",      25,    0.,  25.);
    fNumHitsSpike[0]       = tfs->make<TH1D>("NHitsSpike0",    "; # hits",       250,    0., 500.);
    fNumHitsSpike[1]       = tfs->make<TH1D>("NHitsSpike1",    "; # hits",       250,    0., 500.);
    fNumHitsSpike[2]       = tfs->make<TH1D>("NHitsSpike2",    "; # hits",       250,    0., 500.);
    fDeltaTSpike[0]        = tfs->make<TH1D>("DeltaTSpike0",   "; dTicks",       100,  -10.,  10.);
    fDeltaTSpike[1]        = tfs->make<TH1D>("DeltaTSpike1",   "; dTicks",       100,  -10.,  10.);
    fDeltaTSpike[2]        = tfs->make<TH1D>("DeltaTSpike2",   "; dTicks",       100,  -10.,  10.);
    fSpreadTSpike[0]       = tfs->make<TH1D>("SpreadTSpike0",  "; dTicks",       100,    0.,  25.);
    fSpreadTSpike[1]       = tfs->make<TH1D>("SpreadTSpike1",  "; dTicks",       100,    0.,  25.);
    fSpreadTSpike[2]       = tfs->make<TH1D>("SpreadTSpike2",  "; dTicks",       100,    0.,  25.);
    fRMSSpike[0]           = tfs->make<TH1D>("RMSSpike0",      "; rms",           50,    0.,   5.);
    fRMSSpike[1]           = tfs->make<TH1D>("RMSSpike1",      "; rms",           50,    0.,   5.);
    fRMSSpike[2]           = tfs->make<TH1D>("RMSSpike2",      "; rms",           50,    0.,   5.);
    fPulseHeight[0]        = tfs->make<TH1D>("PulseHeight0",   "; PH(ADC)",      500,    0., 200.);
    fPulseHeight[1]        = tfs->make<TH1D>("PulseHeight1",   "; PH(ADC)",      500,    0., 200.);
    fPulseHeight[2]        = tfs->make<TH1D>("PulseHeight2",   "; PH(ADC)",      500,    0., 200.);
    fRMS[0]                = tfs->make<TH1D>("RMS0",           "; RMS(ticks)",    50,    0.,  10.);
    fRMS[1]                = tfs->make<TH1D>("RMS1",           "; RMS(ticks)",    50,    0.,  10.);
    fRMS[2]                = tfs->make<TH1D>("RMS2",           "; RMS(ticks)",    50,    0.,  10.);
    fSigmaPeakTime[0]      = tfs->make<TH1D>("SigPeak0",       "; sigma(ticks)", 100,    0.,  20.);
    fSigmaPeakTime[1]      = tfs->make<TH1D>("SigPeak1",       "; sigma(ticks)", 100,    0.,  20.);
    fSigmaPeakTime[2]      = tfs->make<TH1D>("SigPeak2",       "; sigma(ticks)", 100,    0.,  20.);
    fHitIndex[0]           = tfs->make<TH1D>("HitIndex0",      "; index",         20,    0.,  20.);
    fHitIndex[1]           = tfs->make<TH1D>("HitIndex1",      "; index",         20,    0.,  20.);
    fHitIndex[2]           = tfs->make<TH1D>("HitIndex2",      "; index",         20,    0.,  20.);
    fROIlength[0]          = tfs->make<TH1D>("ROIlen0",        "; # bins",       200,    0., 200.);
    fROIlength[1]          = tfs->make<TH1D>("ROIlen1",        "; # bins",       200,    0., 200.);
    fROIlength[2]          = tfs->make<TH1D>("ROIlen2",        "; # bins",       200,    0., 200.);
    fChi2[0]               = tfs->make<TH1D>("Chi20",          "; chi2",         200,    0., 100.);
    fChi2[1]               = tfs->make<TH1D>("Chi21",          "; chi2",         200,    0., 100.);
    fChi2[2]               = tfs->make<TH1D>("Chi22",          "; chi2",         200,    0., 100.);
    fPHvsRMS[0]            = tfs->make<TH2D>("PHvsRMS0",       "RMS;PH(ADC)",    200,    0., 100., 100, 0.,  20.);
    fPHvsRMS[1]            = tfs->make<TH2D>("PHvsRMS1",       "RMS;PH(ADC)",    200,    0., 100., 100, 0.,  20.);
    fPHvsRMS[2]            = tfs->make<TH2D>("PHvsRMS2",       "RMS;PH(ADC)",    200,    0., 100., 100, 0.,  20.);
    fPHvsSigma[0]          = tfs->make<TH2D>("PHvsSigma0",     "Sigma;PH(ADC)",  200,    0., 100., 100, 0.,  20.);
    fPHvsSigma[1]          = tfs->make<TH2D>("PHvsSigma1",     "Sigma;PH(ADC)",  200,    0., 100., 100, 0.,  20.);
    fPHvsSigma[2]          = tfs->make<TH2D>("PHvsSigma2",     "Sigma;PH(ADC)",  200,    0., 100., 100, 0.,  20.);
    fROIvsPH[0]            = tfs->make<TH2D>("ROIvsPH0",       "PH (ADC); ROI",  200,    0., 200., 200, 0., 100.);
    fROIvsPH[1]            = tfs->make<TH2D>("ROIvsPH1",       "PH (ADC); ROI",  200,    0., 200., 200, 0., 100.);
    fROIvsPH[2]            = tfs->make<TH2D>("ROIvsPH2",       "PH (ADC); ROI",  200,    0., 200., 200, 0., 100.);
    fROIvsChi[0]           = tfs->make<TH2D>("ROIvsChi0",      "Chi; ROI",       200,    0., 100., 200, 0., 100.);
    fROIvsChi[1]           = tfs->make<TH2D>("ROIvsChi1",      "Chi; ROI",       200,    0., 100., 200, 0., 100.);
    fROIvsChi[2]           = tfs->make<TH2D>("ROIvsChi2",      "Chi; ROI",       200,    0., 100., 200, 0., 100.);
    
    fTimeHitCount[0]       = tfs->make<TH1D>("TimeHitCnt0",    "Count;Ticks",    9600,   0., 9600.);
    fTimeHitCount[1]       = tfs->make<TH1D>("TimeHitCnt1",    "Count;Ticks",    9600,   0., 9600.);
    fTimeHitCount[2]       = tfs->make<TH1D>("TimeHitCnt2",    "Count;Ticks",    9600,   0., 9600.);
    
    fWireVsTime[0]         = tfs->make<TH2D>("WireVsTime0",    ";Wire;Ticks",    fGeometry->Nwires(0)/2, 0., fGeometry->Nwires(0), 2400, 0., 9600.);
    fWireVsTime[1]         = tfs->make<TH2D>("WireVsTime1",    ";Wire;Ticks",    fGeometry->Nwires(1)/2, 0., fGeometry->Nwires(1), 2400, 0., 9600.);
    fWireVsTime[2]         = tfs->make<TH2D>("WireVsTime2",    ";Wire;Ticks",    fGeometry->Nwires(2)/2, 0., fGeometry->Nwires(2), 2400, 0., 9600.);
    fWireVsPH[0]           = tfs->make<TH2D>("WireVsPH0",      ";Wire;PH",       fGeometry->Nwires(0)/2, 0., fGeometry->Nwires(0),  100, 0.,  100.);
    fWireVsPH[1]           = tfs->make<TH2D>("WireVsPH1",      ";Wire;PH",       fGeometry->Nwires(1)/2, 0., fGeometry->Nwires(1),  100, 0.,  100.);
    fWireVsPH[2]           = tfs->make<TH2D>("WireVsPH2",      ";Wire;PH",       fGeometry->Nwires(2)/2, 0., fGeometry->Nwires(2),  100, 0.,  100.);
    fWireVsPHRat[0]        = tfs->make<TH2D>("WireVsPHRat0",   ";Wire;PH Ratio", fGeometry->Nwires(0)/2, 0., fGeometry->Nwires(0),  100, 0.,    2.);
    fWireVsPHRat[1]        = tfs->make<TH2D>("WireVsPHRat1",   ";Wire;PH Ratio", fGeometry->Nwires(1)/2, 0., fGeometry->Nwires(1),  100, 0.,    2.);
    fWireVsPHRat[2]        = tfs->make<TH2D>("WireVsPHRat2",   ";Wire;PH Ratio", fGeometry->Nwires(2)/2, 0., fGeometry->Nwires(2),  100, 0.,    2.);
    
    fNumMatchedClusters[0] = tfs->make<TH1D>("NMatClus0",      ";# clusters",     25,    0.,  25.);
    fNumMatchedClusters[1] = tfs->make<TH1D>("NMatClus1",      ";# clusters",     25,    0.,  25.);
    fNumMatchedClusters[2] = tfs->make<TH1D>("NMatClus2",      ";# clusters",     25,    0.,  25.);
    
    fNumMatchedTracks[0]   = tfs->make<TH1D>("NMatTrack0",     ";# tracks",       25,    0.,  25.);
    fNumMatchedTracks[1]   = tfs->make<TH1D>("NMatTrack1",     ";# tracks",       25,    0.,  25.);
    fNumMatchedTracks[2]   = tfs->make<TH1D>("NMatTrack2",     ";# tracks",       25,    0.,  25.);
    
    fDeltaTrackT[0]        = tfs->make<TH1D>("DTrackT0",       ";Delta T",       200, -100., 100.);
    fDeltaTrackT[1]        = tfs->make<TH1D>("DTrackT1",       ";Delta T",       200, -100., 100.);
    fDeltaTrackT[2]        = tfs->make<TH1D>("DTrackT2",       ";Delta T",       200, -100., 100.);
    
    fDeltaClusterT[0]      = tfs->make<TH1D>("DClusterT0",     ";Delta T",       200, -100., 100.);
    fDeltaClusterT[1]      = tfs->make<TH1D>("DClusterT1",     ";Delta T",       200, -100., 100.);
    fDeltaClusterT[2]      = tfs->make<TH1D>("DClusterT2",     ";Delta T",       200, -100., 100.);
    
    fTrackComplete[0]      = tfs->make<TH1D>("TrackComplete0", ";fraction",      100,    0.,   1.);
    fTrackComplete[1]      = tfs->make<TH1D>("TrackComplete1", ";fraction",      100,    0.,   1.);
    fTrackComplete[2]      = tfs->make<TH1D>("TrackComplete2", ";fraction",      100,    0.,   1.);
    
    fTrackPos[0]           = tfs->make<TH2D>("TrackPos0",      ";z(cm);y(cm)",   200,    0.,   fGeometry->DetLength(), 100, -fGeometry->DetHalfHeight(), fGeometry->DetHalfHeight());
    fTrackPos[1]           = tfs->make<TH2D>("TrackPos1",      ";z(cm);y(cm)",   200,    0.,   fGeometry->DetLength(), 100, -fGeometry->DetHalfHeight(), fGeometry->DetHalfHeight());
    fTrackPos[2]           = tfs->make<TH2D>("TrackPos2",      ";z(cm);y(cm)",   200,    0.,   fGeometry->DetLength(), 100, -fGeometry->DetHalfHeight(), fGeometry->DetHalfHeight());
    
    fTrackCompPos[0]       = tfs->make<TH2D>("TrackCompPos0",  ";z(cm);y(cm)",   200,    0.,   fGeometry->DetLength(), 100, -fGeometry->DetHalfHeight(), fGeometry->DetHalfHeight());
    fTrackCompPos[1]       = tfs->make<TH2D>("TrackCompPos1",  ";z(cm);y(cm)",   200,    0.,   fGeometry->DetLength(), 100, -fGeometry->DetHalfHeight(), fGeometry->DetHalfHeight());
    fTrackCompPos[2]       = tfs->make<TH2D>("TrackCompPos2",  ";z(cm);y(cm)",   200,    0.,   fGeometry->DetLength(), 100, -fGeometry->DetHalfHeight(), fGeometry->DetHalfHeight());
    
    // zero out the event counter
    fNumEvents = 0;
    
    // Set up our channel pedestal vec (for output) if needed
    fChannelPedVec.clear();
    
    if (fWritePedestals) fChannelPedVec.resize(fGeometry->Nchannels());
}
   
//-----------------------------------------------------------------------
void PMTPickupAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void PMTPickupAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fDigitModuleLabel     = p.get< std::string >("DigitModuleLabel",     "daq"  );
    fHitProducerLabel     = p.get< std::string >("HitModuleLabel",       "gauss");
    fClusterProducerLabel = p.get< std::string >("ClusterProducerLabel", "cluster3d");
    fTrackProducerLabel   = p.get< std::string >("TrackProducerLabel",   "trackkalmanhit");

    return;
}

//-----------------------------------------------------------------------
void PMTPickupAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    fNumEvents++;
    
    // We'll use this to control output
    static bool dumpEm(true);
    
    // Do a quick check of the reco hit finding
    art::Handle< std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);

    if (hitHandle.isValid())
    {
        std::vector<std::vector<size_t>> channelHitVec;
        
        channelHitVec.resize(fGeometry->Nchannels());
        
        std::vector<std::vector<const recob::Hit*>> viewHitVec;
        
        viewHitVec.resize(fGeometry->Nplanes());

        // Make pass through the complete list of hits and fill our data structures
        for(size_t rdIter = 0; rdIter < hitHandle->size(); ++rdIter)
        {
            // get the reference to the current raw::RawDigit
            art::Ptr<recob::Hit> hitPtr(hitHandle, rdIter);
            
            // Keep track if the index in our vector of vectors
            channelHitVec[hitPtr->Channel()].push_back(rdIter);
            viewHitVec[hitPtr->View()].push_back(hitPtr.get());
        }

        // It is helpful if the hits are in time order (by view)
        for(auto& hitVec : viewHitVec)
            std::sort(hitVec.begin(),hitVec.end(),[](const recob::Hit* left, const recob::Hit* right) {return left->PeakTime() < right->PeakTime();});
        
        // Now run through to find spikes
        PeakTupleVec peakTupleVec;
        
        peakTupleVec.resize(fGeometry->Nplanes());
        
        // Build the data structure containing the hits associated to spikes
        buildPeakDataStructure(viewHitVec, peakTupleVec);

        // Go through and do some basic histogramming
        for(size_t viewIdx = 0; viewIdx < fGeometry->Nplanes(); viewIdx++)
        {
            // Guard against nothing to do
            if (viewHitVec[viewIdx].empty()) continue;
            
//            std::vector<const recob::Hit*>& hitVec = viewHitVec[viewIdx];
            
            fNumSpikesEvent[viewIdx]->Fill(peakTupleVec[viewIdx].size(), 1.);
            
            // Output to screen if asked for
            if (dumpEm)
            {
                std::cout << "***************************************************************************************" << std::endl;
                std::cout << "----> Run: " << fRun << ", subrun: " << fSubRun << ", event: " << fEvent << std::endl;
                for(size_t peakIdx = 0; peakIdx < peakTupleVec[viewIdx].size(); peakIdx++)
                {
                    const auto&   peakTupleVal = peakTupleVec[viewIdx][peakIdx];
                    double        avePeakTime  = std::get<0>(peakTupleVal);
                    const HitVec& spikeHitVec  = std::get<3>(peakTupleVal);
                    
                    std::cout << "****>> Spike # " << peakIdx << " has average hit time: " << avePeakTime << std::endl;
                    std::cout << "Collection plane wires with hits: ";
                    
                    for(const auto& hit : spikeHitVec)
                        std::cout << hit->WireID().Wire << " ";
                    
                    std::cout << std::endl;
                }
                
                std::cout << "***************************************************************************************" << std::endl;
                dumpEm = false;
            }
        }
        
        // ***************************************************************************************************************
        // In the beginnig we will match to clusters in the collection plane
        art::Handle<std::vector<recob::Cluster> > clusterHandle;
        event.getByLabel(fClusterProducerLabel, clusterHandle);
        
        if (clusterHandle.isValid())
        {
            // Keep track of results in the following scheme:
            std::vector<std::vector<ClusterHitMatchTuple>> clusterPeakMatchVec;
            
            clusterPeakMatchVec.resize(peakTupleVec[2].size());
            
            // We'll use this to keep track of the collection plane hits associated to a cluster
            std::map<const recob::Cluster*, HitVec> clusterHitVecMap;
            
            // max difference
            double maxDeltaTicks(200.);
            
            // Recover the collection of associations between tracks and hits
            art::FindManyP<recob::Hit> clusterAssns(clusterHandle, event, fClusterProducerLabel);
            
            // Are there PFParticles with these clusters?
            art::FindManyP<recob::PFParticle> pfParticleAssns(clusterHandle, event, fClusterProducerLabel);
        
            // Go through the clusters to recover hits for the collection plane
            for(size_t clusterIdx = 0; clusterIdx < clusterHandle->size(); clusterIdx++)
            {
                art::Ptr<recob::Cluster> cluster(clusterHandle,clusterIdx);
            
                // if not in collection plane then skip
                if (cluster->View() != geo::kW) continue;
            
                // Recover associated hits
                std::vector<art::Ptr<recob::Hit>> clusterHitPtrVec = clusterAssns.at(cluster.key());
                
                // We actually just want bare pointers...
                HitVec& clusterHitVec = clusterHitVecMap[cluster.get()];
                
                for(const auto& hit : clusterHitPtrVec)
                    clusterHitVec.push_back(hit.get());
            
                // It is helpful if the hits are in time order (by view)
                std::sort(clusterHitVec.begin(),clusterHitVec.end(),[](const recob::Hit* left, const recob::Hit* right) {return left->PeakTime() < right->PeakTime();});
            
                // Try to match the end hits to one of our spikes
                const recob::Hit* frontHit = clusterHitVec.front();
                const recob::Hit* backHit  = clusterHitVec.back();
            
                const auto& peakTupleVals = peakTupleVec[2];
            
                for(size_t peakIdx = 0; peakIdx < peakTupleVals.size(); peakIdx++)
                {
                    const auto& peakTupleVal = peakTupleVals[peakIdx];
                    
                    // require that the hit be contained in the longest range of the hits in the spike
                    const WirePair& spikeWirePair = std::get<2>(peakTupleVal);
                    double          avePeakTime = std::get<0>(peakTupleVal);
                    double          deltaFront  = frontHit->PeakTime() - avePeakTime;
                    double          deltaBack   = backHit->PeakTime()  - avePeakTime;
                    
                    if (fabs(deltaFront) >= maxDeltaTicks && fabs(deltaBack) >= maxDeltaTicks) continue;
                    
                    if (fabs(deltaFront) <  maxDeltaTicks && (frontHit->WireID().Wire >= spikeWirePair.first || frontHit->WireID().Wire <= spikeWirePair.second))
                    {
                        clusterPeakMatchVec[peakIdx].push_back(std::make_tuple(deltaFront,cluster,frontHit));
                    }
                    else if (fabs(deltaBack)  <  maxDeltaTicks && (backHit->WireID().Wire >= spikeWirePair.first || backHit->WireID().Wire <= spikeWirePair.second))
                    {
                        clusterPeakMatchVec[peakIdx].push_back(std::make_tuple(deltaBack, cluster,backHit));
                    }
                }
            }
            
            // Keep track of how many "good" matches we make
            size_t numClusterMatches(0);
            
            // Go through the spikes and sort the matches so we can see who wins
            for(size_t peakMatchIdx = 0; peakMatchIdx < clusterPeakMatchVec.size(); peakMatchIdx++)
            {
                auto& clusterMatchVec = clusterPeakMatchVec[peakMatchIdx];
                
                if (clusterMatchVec.empty())
                {
                    continue;
                }
                
                std::sort(clusterMatchVec.begin(), clusterMatchVec.end(), [](const ClusterHitMatchTuple& left, const ClusterHitMatchTuple& right) {return fabs(std::get<0>(left)) < fabs(std::get<0>(right));});
                
                // The winner is the first one
                double bestMatch = std::get<0>(clusterMatchVec.front());
            
                fDeltaClusterT[2]->Fill(std::max(-99.9,std::min(99.9,bestMatch)), 1.);
                
//                std::cout << "  Spike #" << peakMatchIdx << " has best match: " << bestMatch << " to cluster: " << std::get<1>(clusterMatchVec.front())->ID()
//                          << ", time: " << std::get<2>(clusterMatchVec.front())->PeakTime() << ", " << std::get<0>(peakTupleVec[2][peakMatchIdx]) << std::endl;
//                std::cout << "        cluster has " << clusterHitVecMap[std::get<1>(clusterMatchVec.front()).get()].size() << " hits " << std::endl;
                
//                if (pfParticleAssns.isValid())
//                {
//                    std::vector<art::Ptr<recob::PFParticle>> pfParticlePtrVec = pfParticleAssns.at(std::get<1>(clusterMatchVec.front()).key());
                    
//                    if (!pfParticlePtrVec.empty())
//                        std::cout << "        PFParticle ID: " << pfParticlePtrVec.front()->Self() << std::endl;
//                }
                
                // Let's only consider good matches
                if (fabs(bestMatch) > 75.) continue;
                
                numClusterMatches++;
                
                // Ok, grab the associated spike and plot stuff
                size_t      viewIdx      = 2;
                const auto& peakTupleVal = peakTupleVec[viewIdx][peakMatchIdx];
                
                double        avePeakTime = std::get<0>(peakTupleVal);
                double        rmsPeakTime = std::get<1>(peakTupleVal);
                const HitVec& peakHitVec  = std::get<3>(peakTupleVal);
                
                size_t numBins = peakHitVec.size();
                
                fNumHitsSpike[viewIdx]->Fill(numBins, 1.);
                fRMSSpike[viewIdx]->Fill(std::min(9.99,double(rmsPeakTime)), 1.);
                fTimeHitCount[viewIdx]->Fill(avePeakTime, 1.);
                
                double maxPulseHeight = 0.;
                
                for(const auto& hit : peakHitVec)
                    maxPulseHeight = std::max(maxPulseHeight,double(hit->PeakAmplitude()));
                
                double peakTimeHigh(0.);
                double peakTimeLow(9600.);
                
                // Loop through the hits keeping track of important quantities
                // And remember that the hits are now in wire order
                for(const auto& hit : peakHitVec)
                {
                    double peakTime      = hit->PeakTime();
                    double pulseHeight   = std::min(float(249.),  hit->PeakAmplitude());
                    double pulseHghtRat  = pulseHeight / maxPulseHeight;
                    double RMS           = std::min(float(19.9),  hit->RMS());
                    double sigmaPeakTime = std::min(float(19.9),  hit->SigmaPeakTime());
                    double roiLength     = std::min(double(199.), double(hit->EndTick() - hit->StartTick()));
                    double chi2          = std::min(float(99.),   hit->GoodnessOfFit());
                    double deltaTime     = std::max(-24.9,        std::min(24.9, peakTime-avePeakTime));
                    int    hitIndex      = std::min(int(19),      int(hit->LocalIndex()));
                    
                    peakTimeHigh = std::max(peakTimeHigh, peakTime);
                    peakTimeLow  = std::min(peakTimeLow,  peakTime);
                    
                    fHitsByWire[viewIdx]->Fill(hit->WireID().Wire);
                    fHitsByWireProf[viewIdx]->Fill(hit->WireID().Wire,1.);
                    fWireVsTime[viewIdx]->Fill(hit->WireID().Wire,peakTime,1.);
                    fWireVsPH[viewIdx]->Fill(hit->WireID().Wire,pulseHeight);
                    fWireVsPHRat[viewIdx]->Fill(hit->WireID().Wire,pulseHghtRat);
                    fDeltaTSpike[viewIdx]->Fill(deltaTime, 1.);
                    fPulseHeight[viewIdx]->Fill(pulseHeight);
                    fRMS[viewIdx]->Fill(RMS);
                    fSigmaPeakTime[viewIdx]->Fill(sigmaPeakTime);
                    fHitIndex[viewIdx]->Fill(hitIndex);
                    fPHvsRMS[viewIdx]->Fill(pulseHeight, RMS);
                    fPHvsSigma[viewIdx]->Fill(pulseHeight, sigmaPeakTime);
                    fTimeHitCount[viewIdx]->Fill(hit->PeakTime(), 1.);
                    
                    if (hitIndex == 0)
                    {
                        fROIlength[viewIdx]->Fill(roiLength);
                    }
                    
                    fChi2[viewIdx]->Fill(chi2);
                    fROIvsPH[viewIdx]->Fill(roiLength, pulseHeight);
                    fROIvsChi[viewIdx]->Fill(roiLength, chi2);
                }
                
                double deltaPeakTimeMax = peakTimeHigh - peakTimeLow;
                
                fSpreadTSpike[viewIdx]->Fill(deltaPeakTimeMax, 1.);
            }
            
            fNumMatchedClusters[2]->Fill(numClusterMatches, 1.);
         
            // *****************************************************************************************************************
            // Now do some track matching
            art::Handle<std::vector<recob::Track> > trackHandle;
            event.getByLabel(fTrackProducerLabel, trackHandle);

            if (trackHandle.isValid())
            {
                // Recover the collection of associations between tracks and hits
                art::FindManyP<recob::Hit> trackAssns(trackHandle, event, fTrackProducerLabel);
                
                // Keep track of matched tracks in a vector
                std::vector<std::pair<const recob::Track*,size_t>> matchedTrackVec;
                
                matchedTrackVec.resize(clusterPeakMatchVec.size(), std::pair<const recob::Track*,size_t>(0,0));
                
                size_t bestMatchHits(5);
        
                // The goal here is to try to associate tracks to the clusters we have already associated
                // to spikes. The program is to match clusters to tracks by comparing hits.
                // Start by making a map betweeen tracks and their (sorted) collection plane hits
                std::map<const recob::Track*,HitVec> trackHitVecMap;
                
                for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
                {
                    art::Ptr<recob::Track> track(trackHandle,trackIdx);

                    HitVec& trackHitVec = trackHitVecMap[track.get()];
                    
                    // Recover the associated hits
                    std::vector<art::Ptr<recob::Hit>> trackHitPtrVec = trackAssns.at(track.key());
            
                    for(const auto& hit : trackHitPtrVec)
                        if (hit->View() == 2) trackHitVec.push_back(hit.get());
            
                    // It is helpful if the hits are in time order (by view)
                    std::sort(trackHitVec.begin(),trackHitVec.end(),[](const recob::Hit* left, const recob::Hit* right) {return left->PeakTime() < right->PeakTime();});
                }
                
                // Now pass through this map and do some matching
                for(auto& trackHitVecMapItr : trackHitVecMap)
                {
                    const recob::Track* track       = trackHitVecMapItr.first;
                    HitVec&             trackHitVec = trackHitVecMapItr.second;
                    
                    for(size_t peakMatchIdx = 0; peakMatchIdx < clusterPeakMatchVec.size(); peakMatchIdx++)
                    {
                        auto& clusterMatchVec = clusterPeakMatchVec[peakMatchIdx];
                        
                        // If no match nothing to do
                        if (clusterMatchVec.empty()) continue;
                        
                        // Recover the cluster hits through a roundabout procedure
                        art::Ptr<recob::Cluster> cluster = std::get<1>(clusterMatchVec.front());
                        
                        // container for cluster hits (recob hit pointers)
                        HitVec clusterHitVec;
                        
                        // container for the matched hits
                        HitVec matchedHitVec;
                        
                        // Get the associated hits and fill the container
                        std::vector<art::Ptr<recob::Hit>> clusterHitPtrVec = clusterAssns.at(cluster.key());
                        
                        // Get this cluster's hits for collection plane
                        for(const auto& hit : clusterHitPtrVec)
                            if (hit->View() == 2) clusterHitVec.push_back(hit.get());
                        
                        // Sort them in same way as for track hits
                        std::sort(clusterHitVec.begin(),clusterHitVec.end(),[](const recob::Hit* left, const recob::Hit* right){return left->PeakTime() < right->PeakTime();});
                        
                        // Look for an intersection
                        std::set_intersection(clusterHitVec.begin(), clusterHitVec.end(),
                                              trackHitVec.begin(),   trackHitVec.end(),
                                              std::back_inserter(matchedHitVec),
                                              [](const recob::Hit* left, const recob::Hit* right){return left->PeakTime() < right->PeakTime();});
                        
                        // If the best match then keep it
                        if (matchedHitVec.size() > bestMatchHits)
                        {
                            bestMatchHits = matchedHitVec.size();
                            matchedTrackVec[peakMatchIdx] = std::pair<const recob::Track*,size_t>(track,bestMatchHits);
                        }
                    }
                }
                
                // Smell some roses...
                size_t numMatchedTracks(0);
                
                for(size_t peakIdx = 0; peakIdx < matchedTrackVec.size(); peakIdx++)
                {
                    if (matchedTrackVec[peakIdx].second > 0)
                    {
                        const recob::Track*   track         = matchedTrackVec[peakIdx].first;
                        const recob::Cluster* cluster       = std::get<1>(clusterPeakMatchVec[peakIdx].front()).get();

                        numMatchedTracks++;
                        
                        double trackComplete = double(trackHitVecMap[track].size()) / double(clusterHitVecMap[cluster].size());
                        fTrackComplete[2]->Fill(std::min(99.5,trackComplete), 1.);
                        
                        // How to figure out which end of the track we want...
                        // Convert the hit in question to x-z coordinates, find arc length along track
                        // to the hit's distance of closest approach to track. If less than half length
                        // of track then start of track is what we want, else end of track
                        // Recover track parameters at start of track
                        TVector3 trackPos    = track->Vertex();
                        TVector3 trackEnd    = track->End();
                        TVector3 trackDir    = track->VertexDirection();
                        TVector3 trackEndDir = track->EndDirection();
                        
                        // Consider an X-Z plane defined by the spike time... find the track end closest to this
                        // plane and then project the track to it
                        double avePeakTime    = std::get<0>(peakTupleVec[2][peakIdx]);
                        double avePeakTimePos = fDetectorProperties->ConvertTicksToX(avePeakTime,2,0,0);
                        double deltaXPos      = trackPos[0] - avePeakTimePos;
                        double deltaXPosEnd   = trackEnd[0] - avePeakTimePos;
                        
                        // Swap ends if end is closer
                        if (fabs(deltaXPosEnd) < fabs(deltaXPos))
                        {
                            trackPos  =  trackEnd;
                            trackDir  =  trackEndDir;
                            deltaXPos =  deltaXPosEnd;
                        }
                        
                        // Ok, now now want the arc length to project the track to intersect with this plane.
                        double trackToXCosTheta = trackDir[0];
                        double arcLenToPlane    = deltaXPos / trackToXCosTheta;
                        
                        // Now project track to this plane (remember track vector points in direction of track)
                        trackPos -= arcLenToPlane * trackDir;
                        
                        double zPos = std::max(0., std::min(0.995*fGeometry->DetLength(), trackPos[2]));
                        double yPos = std::max(-0.995*fGeometry->DetHalfHeight(), std::min(0.995*fGeometry->DetHalfHeight(), trackPos[1]));
                        
                        fTrackPos[2]->Fill(zPos,yPos,1.);
                        
                        // Now get parameters for the hit
                        double trackDeltaT = fDetectorProperties->ConvertXToTicks(trackPos[0],2,0,0) - avePeakTime;
                        
                        fDeltaTrackT[2]->Fill(std::max(-99.9,std::min(99.9,trackDeltaT)), 1.);
                        
                        if (trackComplete > 0.85) fTrackCompPos[2]->Fill(zPos,yPos,1.);
                    }
                }
                
                fNumMatchedTracks[2]->Fill(numMatchedTracks, 1.);
            }
        }
    }

    return;
}

void PMTPickupAna::buildPeakDataStructure(ViewHitVec&   viewHitVec,
                                          PeakTupleVec& peakTupleVec) const
{
    size_t maxGapSize(2);
    size_t peakThreshold(10);
    size_t maxWireGap(99);    // we're a bit generous here, looking for electronics problems
    
    // The purpose of this loop is to find hits which are common in time (so, a "spike")
    // Once we find them we'll pause to make a few useful histograms
    for(size_t viewIdx = 0; viewIdx < fGeometry->Nplanes(); viewIdx++)
    {
        // Guard against nothing to do
        if (viewHitVec[viewIdx].empty()) continue;
        
        HitVec& inputHitVec = viewHitVec[viewIdx];
        
        size_t lastPeakTime(inputHitVec.front()->PeakTime());
        double avePeakTime(0.);
        double avePeakTime2(0.);
        double avePeakRms(0.);
        
        HitVec peakHitVec;
        
        for(size_t hitIdx = 0; hitIdx < inputHitVec.size(); hitIdx++)
        {
            const recob::Hit* hit      = inputHitVec[hitIdx];
            size_t            curTicks = hit->PeakTime();
            
            // If the gap is more than we tolerate then eject the last peak
            // and start a new one.
            if (curTicks > lastPeakTime + maxGapSize)
            {
                // Get average peak RMS
                double peakCount   = peakHitVec.size();
                
                avePeakRms /= peakCount;
                
                // Save if over threshold
                if (peakHitVec.size() > peakThreshold && avePeakRms < 4.)
                {
                    double aveTime   = avePeakTime / peakCount;
                    double rmsTime   = avePeakTime2 / peakCount - aveTime * aveTime;
                    
                    rmsTime = std::sqrt(std::max(0.,rmsTime));
                    
                    if (rmsTime < 2.)
                    {
                        auto peakTupleVals = std::make_tuple(aveTime,rmsTime,WirePair(0,0),peakHitVec);
                        peakTupleVec[viewIdx].push_back(peakTupleVals);
                    }
                }
                
                // Set to new values
                avePeakTime  = 0.;
                avePeakTime2 = 0.;
                avePeakRms   = 0.;
                peakHitVec.clear();
            }
            
            avePeakTime  += curTicks;
            avePeakTime2 += curTicks * curTicks;
            avePeakRms   += hit->RMS();
            lastPeakTime  = curTicks;
            peakHitVec.push_back(hit);
        }
        
        // end condition
        if (peakHitVec.size() > peakThreshold)
        {
            double peakCount = peakHitVec.size();
            double aveTime   = avePeakTime / peakCount;
            double rmsTime   = avePeakTime2 / peakCount - aveTime * aveTime;
            
            rmsTime = std::sqrt(std::max(0.,rmsTime));
            
            auto peakTupleVals = std::make_tuple(aveTime,rmsTime,WirePair(0,0),peakHitVec);
            peakTupleVec[viewIdx].push_back(peakTupleVals);
        }
        
        // We want to reorder the hits associated to each spike in wire order
        // and also now we can populate our WirePair structure
        for(auto& peakTupleVal : peakTupleVec[viewIdx])
        {
            HitVec& tupleHitVec  = std::get<3>(peakTupleVal);
            std::sort(tupleHitVec.begin(),tupleHitVec.end(),[](const recob::Hit* left, const recob::Hit* right){return left->WireID().Wire < right->WireID().Wire;});
            
            // Interestingly... the first wires often see the pulses so what we want to do next is find the longest contiguous block
            std::vector<WirePair> wireStartStopVec;
            
            size_t lastWire  = tupleHitVec.front()->WireID().Wire;
            size_t startWire = lastWire;
            
            for(const auto& hit : tupleHitVec)
            {
                size_t curWire = hit->WireID().Wire;
                
                if (curWire - lastWire > maxWireGap)
                {
                    wireStartStopVec.push_back(WirePair(startWire,lastWire));
                    startWire = curWire;
                }
                
                lastWire = curWire;
            }
            
            if (lastWire > startWire) wireStartStopVec.push_back(WirePair(startWire,lastWire));
            
            std::sort(wireStartStopVec.begin(),wireStartStopVec.end(),[](const WirePair& left, const WirePair& right){return left.second-left.first > right.second-right.first;});
            
            WirePair& wirePair = std::get<2>(peakTupleVal);
            
            wirePair = wireStartStopVec.front();
        }
    }

    return;
}
    
void PMTPickupAna::endJob()
{
    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see PMTPickupAna.fcl for more information.
DEFINE_ART_MODULE(PMTPickupAna)

} // namespace PMTPickupAna

#endif // PMTPickupAna_Module
