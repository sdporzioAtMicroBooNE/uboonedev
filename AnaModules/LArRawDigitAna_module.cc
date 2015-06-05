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
    bool        fUseChannelPedestals;        // Use channel-by-channel pedestals
    bool        fWritePedestals;             // Output new file of pedestals

    // Pointers to the histograms we'll create.
    TH1D*     fAdcCntHist[3];
    TH1D*     fAveValHist[3];
    TH1D*     fRmsValHist[3];
    TH1D*     fPedValHist[3];
    TProfile* fRmsValProf[3];
    
    TProfile* fFFTHist;

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
    
    // Commence looping over raw digits
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        
        channel = digitVec->Channel();
       
        // Decode the channel and make sure we have a valid one
        std::vector<geo::WireID> wids;
        try {
            wids = fGeometry->ChannelToWire(channel);
        }
        catch(...)
        {
            //std::cout << "===>> Found illegal channel with id: " << channel << std::endl;
            continue;
        }
        
        if (channel >= maxChannels)
        {
            //std::cout << "***>> Found channel at edge of allowed: " << channel << ", maxChannel: " << maxChannels << std::endl;
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
    
    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see LArRawDigitAna.fcl for more information.
DEFINE_ART_MODULE(LArRawDigitAna)

} // namespace LArRawDigitAna

#endif // LArRawDigitAna_Module