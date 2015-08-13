// LArRawDigitAna_module.cc
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

#ifndef LArRawDigitAna_Module
#define LArRawDigitAna_Module

// LArSoft includes
#include "Geometry/Geometry.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Hit.h"

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
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>

namespace LArRawDigitAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class LArRawDigitAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit LArRawDigitAna(fhicl::ParameterSet const& pset);
    virtual ~LArRawDigitAna();

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

    // The parameters we'll read from the .fcl file.
    std::string fDigitModuleLabel;           // The name of the producer that created raw digits
    std::string fHitProducerLabel;
    bool        fUseChannelPedestals;        // Use channel-by-channel pedestals
    bool        fWritePedestals;             // Output new file of pedestals

    // Pointers to the histograms we'll create.
    TH1D*     fAdcCntHist[3];
    TH1D*     fAveValHist[3];
    TH1D*     fRmsValHist[3];
    TH1D*     fPedValHist[3];
    TProfile* fRmsValProf[3];
    
    TProfile* fFFTHist;
    
    TH1D*     fHitsByWire[3];
    TProfile* fHitsByWireProf[3];
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
    
    TH1D*     fSelPulseHeight[3];
    
    TH1D*     fTimeHitCount[3];

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fNumEvents;
    
    std::vector<short>  fPedestalVec;
    std::vector<double> fChannelPedVec;

    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry>            fGeometry;       // pointer to Geometry service
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;

    double                                       fElectronsToGeV; // conversion factor

}; // class LArRawDigitAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
LArRawDigitAna::LArRawDigitAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
LArRawDigitAna::~LArRawDigitAna()
{}
   
//-----------------------------------------------------------------------
void LArRawDigitAna::beginJob()
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
    fAdcCntHist[0] = tfs->make<TH1D>("CntUPlane", ";#adc",  200, 9000., 10000.);
    fAdcCntHist[1] = tfs->make<TH1D>("CntVPlane", ";#adc",  200, 9000., 10000.);
    fAdcCntHist[2] = tfs->make<TH1D>("CntWPlane", ";#adc",  200, 9000., 10000.);
    fAveValHist[0] = tfs->make<TH1D>("AveUPlane", ";Ave",   120,  -30.,    30.);
    fAveValHist[1] = tfs->make<TH1D>("AveVPlane", ";Ave",   120,  -30.,    30.);
    fAveValHist[2] = tfs->make<TH1D>("AveWPlane", ";Ave",   120,  -30.,    30.);
    fRmsValHist[0] = tfs->make<TH1D>("RmsUPlane", ";RMS",  1000,    0.,   100.);
    fRmsValHist[1] = tfs->make<TH1D>("RmsVPlane", ";RMS",  1000,    0.,   100.);
    fRmsValHist[2] = tfs->make<TH1D>("RmsWPlane", ";RMS",  1000,    0.,   100.);
    fPedValHist[0] = tfs->make<TH1D>("PedUPlane", ";Ped",   200,  1950,  2150.);
    fPedValHist[1] = tfs->make<TH1D>("PedVPlane", ";Ped",   200,  1950,  2150.);
    fPedValHist[2] = tfs->make<TH1D>("PedWPlane", ";Ped",   200,   350,   550.);
    
    fRmsValProf[0] = tfs->make<TProfile>("RmsUPlaneProf",  ";Wire #",  2400, 0., 2400., 0., 100.);
    fRmsValProf[1] = tfs->make<TProfile>("RmsVPlaneProf",  ";Wire #",  2400, 0., 2400., 0., 100.);
    fRmsValProf[2] = tfs->make<TProfile>("RmsWPlaneProf",  ";Wire #",  3456, 0., 3456., 0., 100.);
    
    fFFTHist       = tfs->make<TProfile>("FFT1700", "FFT1700;Tick", 9600, 0., 9600., 0., 10000.);
    
    fHitsByWire[0]     = tfs->make<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]     = tfs->make<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]     = tfs->make<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    fHitsByWireProf[0] = tfs->make<TProfile>("HitsByWireProf0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0), 0., 250.);
    fHitsByWireProf[1] = tfs->make<TProfile>("HitsByWireProf1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1), 0., 250.);
    fHitsByWireProf[2] = tfs->make<TProfile>("HitsByWireProf2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2), 0., 250.);
    fPulseHeight[0]    = tfs->make<TH1D>("PulseHeight0", "; PH(ADC)",      500, 0., 250.);
    fPulseHeight[1]    = tfs->make<TH1D>("PulseHeight1", "; PH(ADC)",      500, 0., 250.);
    fPulseHeight[2]    = tfs->make<TH1D>("PulseHeight2", "; PH(ADC)",      500, 0., 250.);
    fRMS[0]            = tfs->make<TH1D>("RMS0",         "; RMS(ticks)",   100, 0.,  20.);
    fRMS[1]            = tfs->make<TH1D>("RMS1",         "; RMS(ticks)",   100, 0.,  20.);
    fRMS[2]            = tfs->make<TH1D>("RMS2",         "; RMS(ticks)",   100, 0.,  20.);
    fSigmaPeakTime[0]  = tfs->make<TH1D>("SigPeak0",     "; sigma(ticks)", 100, 0.,  20.);
    fSigmaPeakTime[1]  = tfs->make<TH1D>("SigPeak1",     "; sigma(ticks)", 100, 0.,  20.);
    fSigmaPeakTime[2]  = tfs->make<TH1D>("SigPeak2",     "; sigma(ticks)", 100, 0.,  20.);
    fHitIndex[0]       = tfs->make<TH1D>("HitIndex0",    "; index",         10, 0.,  10.);
    fHitIndex[1]       = tfs->make<TH1D>("HitIndex1",    "; index",         10, 0.,  10.);
    fHitIndex[2]       = tfs->make<TH1D>("HitIndex2",    "; index",         10, 0.,  10.);
    fROIlength[0]      = tfs->make<TH1D>("ROIlen0",      "; # bins",       200, 0., 200.);
    fROIlength[1]      = tfs->make<TH1D>("ROIlen1",      "; # bins",       200, 0., 200.);
    fROIlength[2]      = tfs->make<TH1D>("ROIlen2",      "; # bins",       200, 0., 200.);
    fChi2[0]           = tfs->make<TH1D>("Chi20",        "; chi2",         200, 0., 100.);
    fChi2[1]           = tfs->make<TH1D>("Chi21",        "; chi2",         200, 0., 100.);
    fChi2[2]           = tfs->make<TH1D>("Chi22",        "; chi2",         200, 0., 100.);
    fPHvsRMS[0]        = tfs->make<TH2D>("PHvsRMS0",     "RMS;PH(ADC)",    200, 0., 100., 100, 0.,  20.);
    fPHvsRMS[1]        = tfs->make<TH2D>("PHvsRMS1",     "RMS;PH(ADC)",    200, 0., 100., 100, 0.,  20.);
    fPHvsRMS[2]        = tfs->make<TH2D>("PHvsRMS2",     "RMS;PH(ADC)",    200, 0., 100., 100, 0.,  20.);
    fPHvsSigma[0]      = tfs->make<TH2D>("PHvsSigma0",   "Sigma;PH(ADC)",  200, 0., 100., 100, 0.,  20.);
    fPHvsSigma[1]      = tfs->make<TH2D>("PHvsSigma1",   "Sigma;PH(ADC)",  200, 0., 100., 100, 0.,  20.);
    fPHvsSigma[2]      = tfs->make<TH2D>("PHvsSigma2",   "Sigma;PH(ADC)",  200, 0., 100., 100, 0.,  20.);
    fROIvsPH[0]        = tfs->make<TH2D>("ROIvsPH0",     "PH (ADC); ROI",  200, 0., 200., 200, 0., 100.);
    fROIvsPH[1]        = tfs->make<TH2D>("ROIvsPH1",     "PH (ADC); ROI",  200, 0., 200., 200, 0., 100.);
    fROIvsPH[2]        = tfs->make<TH2D>("ROIvsPH2",     "PH (ADC); ROI",  200, 0., 200., 200, 0., 100.);
    fROIvsChi[0]       = tfs->make<TH2D>("ROIvsChi0",    "Chi; ROI",       200, 0., 100., 200, 0., 100.);
    fROIvsChi[1]       = tfs->make<TH2D>("ROIvsChi1",    "Chi; ROI",       200, 0., 100., 200, 0., 100.);
    fROIvsChi[2]       = tfs->make<TH2D>("ROIvsChi2",    "Chi; ROI",       200, 0., 100., 200, 0., 100.);
    
    fSelPulseHeight[0] = tfs->make<TH1D>("SelPulseHgt0", "; PH(ADC)",      500, 0., 250.);
    fSelPulseHeight[1] = tfs->make<TH1D>("SelPulseHgt1", "; PH(ADC)",      500, 0., 250.);
    fSelPulseHeight[2] = tfs->make<TH1D>("SelPulseHgt2", "; PH(ADC)",      500, 0., 250.);
    
    fTimeHitCount[0]   = tfs->make<TH1D>("TimeHitCnt0",  "Count;Ticks",    9600, 0., 9600.);
    fTimeHitCount[1]   = tfs->make<TH1D>("TimeHitCnt1",  "Count;Ticks",    9600, 0., 9600.);
    fTimeHitCount[2]   = tfs->make<TH1D>("TimeHitCnt2",  "Count;Ticks",    9600, 0., 9600.);
    
    // zero out the event counter
    fNumEvents = 0;
    
    // Set up our channel pedestal vec (for output) if needed
    fChannelPedVec.clear();
    
    if (fWritePedestals) fChannelPedVec.resize(fGeometry->Nchannels(), 0.);
    
    // Set up our channel by channel pedestals if requested
    fPedestalVec.clear();
    
    if (fUseChannelPedestals)
    {
        fPedestalVec.resize(fGeometry->Nchannels(), 0.);

        std::ifstream pedestalFile("pedestalRaw.txt");
        
        size_t totalCount(0);
        int    channel;
        float  pedestal;
        
        while(pedestalFile >> channel >> pedestal)
        {
            fPedestalVec[channel] = short(pedestal+0.5);
            totalCount++;
            
            if (totalCount < 200) std::cout << "Read channel, pedestal: " << channel << ", " << pedestal << ", storing: " << fPedestalVec[channel] << std::endl;
        }
        
        if (totalCount != fGeometry->Nchannels())
        {
            std::cout << "Did not read back correct number of channels..." << totalCount << std::endl;
        }
        
        pedestalFile.close();
    }
    else
    {
        fPedestalVec.resize(4800, 2045);                   // Set pedestal value for the U/V planes
        fPedestalVec.resize(fGeometry->Nchannels(), 473);  // Set the pedestal value for the W plane
    }
}
   
//-----------------------------------------------------------------------
void LArRawDigitAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void LArRawDigitAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fDigitModuleLabel    = p.get< std::string >("DigitModuleLabel",    "daq");
    fHitProducerLabel    = p.get< std::string >("HitModuleLabel",      "gauss");
    fUseChannelPedestals = p.get<bool>         ("UseChannelPedestals", false);
    fWritePedestals      = p.get<bool>         ("WritePedestals",      false);
    
    return;
}

//-----------------------------------------------------------------------
void LArRawDigitAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    fNumEvents++;
    
    std::cout << "LArRawdigitAna - run/subrun/event: " << fRun << "/" << fSubRun << "/" << fEvent << std::endl;
    std::cout << "Geomtry expects " << fGeometry->Nwires(0) << " U wires, " << fGeometry->Nwires(1) << " V wires, " << fGeometry->Nwires(2) << std::endl;
    
    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    event.getByLabel(fDigitModuleLabel, digitVecHandle);
    
    if (!digitVecHandle->size())  return;
    
    //we're gonna probably need the time service to convert hit times to TDCs
    art::ServiceHandle<util::TimeService> timeService;
    
    raw::ChannelID_t channel     = raw::InvalidChannelID; // channel number
    unsigned int     maxChannels = fGeometry->Nchannels();
    
    // Define vectors to hold info
    std::vector<int>    channelHitVec;
    std::vector<double> channelAveVec;
    std::vector<double> channelRmsVec;
    
    channelHitVec.resize(maxChannels, 0);
    channelAveVec.resize(maxChannels, 0.);
    channelRmsVec.resize(maxChannels, 0.);
    
    // Define eyeball pedestals for raw adcs
    //short int pedestals[] = {2045, 2045, 473};
    
    std::cout << "***>> In LArRawDigitAna with RawDigit vec of size: " << digitVecHandle->size() << std::endl;
    
    // Commence looping over raw digits
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        
        channel = digitVec->Channel();
        
        bool goodChan(true);
       
        // Decode the channel and make sure we have a valid one
        std::vector<geo::WireID> wids;
        try {
            wids = fGeometry->ChannelToWire(channel);
        }
        catch(...)
        {
            //std::cout << "===>> Found illegal channel with id: " << channel << std::endl;
            goodChan = false;
        }
        
        if (channel >= maxChannels || !goodChan)
        {
            std::cout << "***>> Found channel at edge of allowed: " << channel << ", maxChannel: " << maxChannels << std::endl;
            
            unsigned int dataSize = digitVec->Samples();
            
            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
            
            // And now uncompress
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
            
            int    adcCnt   = 0;
            double aveVal   = 0.;
            double rmsVal   = 0.;
            
            // Loop over the data to get the mean using the eyeball pedestals
            for(size_t bin = 0; bin < dataSize; ++bin)
            {
                float rawAdc = rawadc[bin];
                
                adcCnt++;
                aveVal += rawAdc;
            }
            
            aveVal /= double(adcCnt);
            
            std::vector<float> adcDiffVec;
            
            adcDiffVec.resize(dataSize, 0.);
            
            // Go through again to get the rms
            for(size_t bin = 0; bin < dataSize; ++bin)
            {
                float rawAdc  = rawadc[bin];
                float adcDiff = rawAdc - aveVal;
                
                rmsVal += adcDiff * adcDiff;
                
                adcDiffVec[bin] = adcDiff;
            }
            
            rmsVal = std::sqrt(rmsVal/double(adcCnt));
            
            std::cout << "*** Channel out of range: " << channel << ", ave adc: " << aveVal << ", rms: " << rmsVal << std::endl;
            
            continue;
        }
        
        // Recover plane and wire in the plane
        unsigned int view = wids[0].Plane;
        unsigned int wire = wids[0].Wire;
        
        unsigned int dataSize = digitVec->Samples();
        
        // vector holding uncompressed adc values
        std::vector<short> rawadc(dataSize);
        
        // And now uncompress
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        
        int    adcCnt   = 0;
        double aveVal   = 0.;
        double rmsVal   = 0.;
        short  pedestal = fPedestalVec[channel];
        
        // Loop over the data to get the mean using the eyeball pedestals
        for(size_t bin = 0; bin < dataSize; ++bin)
        {
            float rawAdcLessPed = rawadc[bin] - pedestal;
            
            adcCnt++;
            aveVal += rawAdcLessPed;
        }
        
        aveVal /= double(adcCnt);
        
        std::vector<float> adcDiffVec;
        
        adcDiffVec.resize(dataSize, 0.);
        
        // Go through again to get the rms
        for(size_t bin = 0; bin < dataSize; ++bin)
        {
            float rawAdcLessPed = rawadc[bin]   - pedestal;
            float adcDiff       = rawAdcLessPed - aveVal;
            
            rmsVal += adcDiff * adcDiff;
            
            adcDiffVec[bin] = adcDiff;
        }
        
        rmsVal = std::sqrt(rmsVal/double(adcCnt));
        
        channelHitVec[channel] = adcCnt;
        channelAveVec[channel] = aveVal;
        channelRmsVec[channel] = rmsVal;
        
        rmsVal = std::min(rmsVal, 99.9);

        // Fill some histograms here
        fAdcCntHist[view]->Fill(adcCnt,    1.);
        fAveValHist[view]->Fill(aveVal,    1.);
        fRmsValHist[view]->Fill(rmsVal,    1.);
        fRmsValProf[view]->Fill(wire, rmsVal, 1.);
        
        // If we are outputting a pedestal file then do the work for that here
        if (fWritePedestals)
        {
            std::sort(adcDiffVec.begin(), adcDiffVec.end(), std::less<float>());
            
            size_t maxBin = 0.90 * dataSize;
            float  aveRawAdc(0.);
            
            for(size_t bin = 0; bin < maxBin; bin++)
            {
                aveRawAdc += adcDiffVec[bin] + aveVal + pedestal;
            }
            
            aveRawAdc /= float(maxBin);
            
            fPedValHist[view]->Fill(aveRawAdc, 1.);
            
            fChannelPedVec[channel] += aveRawAdc;
        }
/*
        // Short aside to look at FFT
        if (view == 0 && wire == 408)
        {
            int fftDataSize = dataSize;
            
            TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
            
            double fftInputArray[dataSize];
            
            for(size_t idx = 0; idx < dataSize; idx++) fftInputArray[idx] = adcDiffVec[idx];
            
            fftr2c->SetPoints(fftInputArray);
            fftr2c->Transform();
            
            // Recover the power spectrum...
            double realPart(0.);
            double imaginaryPart(0.);
            
            for(size_t idx = 0; idx < dataSize; idx++)
            {
                fftr2c->GetPointComplex(idx+1, realPart, imaginaryPart);
                
                double bin = double(idx) + 0.5;
                double power = realPart*realPart + imaginaryPart*imaginaryPart;
                
                if (power > 0.) power = std::sqrt(power);
                
                fFFTHist->Fill(bin, power, 1.);
            }

        }
 */
    }
/*
    // Spin through the channelHitVec to look for missing channels
    for(size_t idx = 0; idx < maxChannels; idx++)
    {
        std::vector<geo::WireID> wids = fGeometry->ChannelToWire(idx);
        
        if (channelRmsVec[idx] < 1.)
        {
            std::cout << "###>> channel " << idx << " (view,wire: " << wids[0].Plane << "," << wids[0].Wire << ") has rms: " << channelRmsVec[idx] << std::endl;
        }
        else if (channelRmsVec[idx] < 2.)
        {
            std::cout << "--->> channel " << idx << " (view,wire: " << wids[0].Plane << "," << wids[0].Wire << ") has rms: " << channelRmsVec[idx] << std::endl;
        }
        else if (channelRmsVec[idx] > 40.)
        {
            std::cout << ">>>>> channel " << idx << " (view,wire: " << wids[0].Plane << "," << wids[0].Wire << ") has rms = " << channelRmsVec[idx] << std::endl;
        }
        
        if (channelRmsVec[idx] > 4. && channelRmsVec[idx] < 4.5)
        {
            std::cout << "***>> Chirp: " << idx << ", (view/wire: " << wids[0].Plane << "/" << wids[0].Wire << "), rmsVal: " << channelRmsVec[idx] << std::endl;
        }
    }
*/
    
    // Do a quick check of the reco hit finding
    art::Handle< std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);

    if (hitHandle.isValid())
    {
        std::vector<std::vector<size_t>> channelHitVec;
        
        std::vector<std::vector<std::vector<const recob::Hit*>>> planeTimeHitVec;
        
        planeTimeHitVec.resize(3);
        
        planeTimeHitVec[0].resize(9600);
        planeTimeHitVec[1].resize(9600);
        planeTimeHitVec[2].resize(9600);
        
        channelHitVec.resize(maxChannels);
        
        // Commence looping over raw digits
        for(size_t rdIter = 0; rdIter < hitHandle->size(); ++rdIter)
        {
            // get the reference to the current raw::RawDigit
            art::Ptr<recob::Hit> hitPtr(hitHandle, rdIter);
            
            // Keep track if the index in our vector of vectors
            channelHitVec[hitPtr->Channel()].push_back(rdIter);
            
            // Keep track of the hit itself
            size_t view = hitPtr->View();
            size_t time = hitPtr->PeakTime();
            
            planeTimeHitVec[view][time].push_back(hitPtr.get());
        }
        
        // Now pass through again to do some histogramming
        for (size_t channel = 0; channel < maxChannels; channel++)
        {
            std::vector<geo::WireID> wids = fGeometry->ChannelToWire(channel);
            
            // Recover plane and wire in the plane
            unsigned int view = wids[0].Plane;
            unsigned int wire = wids[0].Wire;
            
            // Fill the outer level hists
            fHitsByWireProf[view]->Fill(wire, channelHitVec[channel].size());
            
            // loop through hits
            for(const auto& rdIter : channelHitVec[channel])
            {
                art::Ptr<recob::Hit> hitPtr(hitHandle, rdIter);
                
                double pulseHeight   = std::min(float(250.),  hitPtr->PeakAmplitude());
                double RMS           = std::min(float(20.),   hitPtr->RMS());
                double sigmaPeakTime = std::min(float(20.),   hitPtr->SigmaPeakTime());
                double roiLength     = std::min(double(200.), double(hitPtr->EndTick() - hitPtr->StartTick()));
                double chi2          = std::min(float(100.),  hitPtr->GoodnessOfFit());
                
                fHitsByWire[view]->Fill(wire);
                fPulseHeight[view]->Fill(pulseHeight);
                fRMS[view]->Fill(RMS);
                fSigmaPeakTime[view]->Fill(sigmaPeakTime);
                fHitIndex[view]->Fill(hitPtr->LocalIndex());
                fPHvsRMS[view]->Fill(pulseHeight, RMS);
                fPHvsSigma[view]->Fill(pulseHeight, sigmaPeakTime);
                
                if (hitPtr->LocalIndex() == 0)
                {
                    fROIlength[view]->Fill(roiLength);
                }
                
                fChi2[view]->Fill(chi2);
                fROIvsPH[view]->Fill(roiLength, pulseHeight);
                fROIvsChi[view]->Fill(roiLength, chi2);
                
                if (view == 0)
                {
                    double cutVal = 2. * pulseHeight / 5.;
                    
                    if (RMS <= cutVal) fSelPulseHeight[view]->Fill(pulseHeight);
                }
                else if (view == 2)
                {
                    double cutVal = 1.75;
                    
                    if (pulseHeight < 100.) cutVal = -pulseHeight / 50. + 2.;
                    
                    if (RMS > cutVal) fSelPulseHeight[view]->Fill(pulseHeight);
                }
            }
        }
        
        for(size_t view = 0; view < 3; view++)
        {
            for(size_t time = 0; time < 9600; time++)
            {
                fTimeHitCount[view]->Fill(time, planeTimeHitVec[view][time].size());
            }
        }
    }

    return;
}
    
void LArRawDigitAna::endJob()
{
    if (fWritePedestals)
    {
        std::ofstream pedestalFile;
        pedestalFile.open("pedestalRaw.txt");
        
        for(size_t idx = 0; idx < fGeometry->Nchannels(); idx++)
        {
            //std::vector<geo::WireID> wids = fGeometry->ChannelToWire(idx);
            
            //pedestalFile << "pedestal[" << idx << "] = " << channelPedVec[idx] << ";\n"; //     view: " << wids[0].Plane << ", wire: " << wids[0].Wire << std::endl;
            pedestalFile <<  idx << " " << fChannelPedVec[idx]/double(fNumEvents) << "\n"; //     view: " << wids[0].Plane << ", wire: " << wids[0].Wire << std::endl;
        }
        
        pedestalFile.close();
    }
    
    std::cout << "***************** Wires with no hits **********************" << std::endl;
    
    // Look for channels which did not see any hits
    for(size_t view = 0; view < fGeometry->Views().size(); view++)
    {
        for(size_t wire = 0; wire < fGeometry->Nwires(view); wire++)
        {
            double nHitsByWire = fHitsByWire[view]->GetBinContent(wire+1);
            
            if (nHitsByWire < 1)
            {
                std::cout << "  --> View: " << view << ", wire: " << wire << " has no hits" << std::endl;
            }
        }
    }
    
    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see LArRawDigitAna.fcl for more information.
DEFINE_ART_MODULE(LArRawDigitAna)

} // namespace LArRawDigitAna

#endif // LArRawDigitAna_Module
