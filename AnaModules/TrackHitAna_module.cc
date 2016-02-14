// TrackHitAna_module.cc
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

#ifndef TrackHitAna_module
#define TrackHitAna_module

// LArSoft includes
#include "Geometry/Geometry.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Track.h"
#include "RecoBase/Spacepoint.h"
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

namespace TrackHitAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class TrackHitAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit TrackHitAna(fhicl::ParameterSet const& pset);
    virtual ~TrackHitAna();

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
    using  HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
    
    double length(const recob::Track* track);

    // The parameters we'll read from the .fcl file.
    std::string fHitProducerLabel;
    std::string fPFParticleProducerLabel;
    std::string fTrackProducerLabel;

    // Pointers to the histograms we'll create.
    // Hit histograms attached to tracks
    TH1D*     fHitsByWire[3];
    
    TH1D*     fPulseHeight[3];
    TH1D*     fPulseHeightSingle[3];
    TH1D*     fPulseHeightMulti[3];
    TH1D*     fChi2DOF[3];
    TH1D*     fChi2DOFSingle[3];
    TH1D*     fHitMult[3];
    TH1D*     fHitCharge[3];
    TH1D*     fFitWidth[3];
    TH1D*     fHitSumADC[3];
    
    TH1D*     fPulseHeightP[3];
    TH1D*     fPulseHeightSingleP[3];
    TH1D*     fPulseHeightMultiP[3];
    TH1D*     fChi2DOFP[3];
    TH1D*     fChi2DOFSingleP[3];
    TH1D*     fHitMultP[3];
    TH1D*     fHitChargeP[3];
    TH1D*     fHitSumADCP[3];
    
    // Hit histograms unattached
    TH1D*     fPulseHeightAll[3];
    TH1D*     fPulseHeightSingleAll[3];
    TH1D*     fChi2DOFAll[3];
    TH1D*     fChi2DOFSingleAll[3];
    TH1D*     fHitMultAll[3];
    TH1D*     fHitChargeAll[3];
    TH1D*     fHitSumADCAll[3];
    
    // Track histograms
    TH1D*     fNumPFPartTracks;
    TH1D*     fNumPFPartTracksL;
    TH1D*     fNumFitTracks;
    TH1D*     fNumTrksPFPart;
    TH1D*     fNumTrksPFPartL;
    TH1D*     fNumFitTrksPFPart;
    TH1D*     fNumFitTrksPFPartL;
    TH1D*     fPFPartTrackLen;
    TH1D*     fPFPartTrackLenL;
    TH1D*     fPFPartEndLenL;
    TH1D*     fFitTrackLen;
    TH1D*     fFitEndLen;
    TH1D*     fTrackDeltaLen;
    TH2D*     fFitVsPFPartLen;
    TH2D*     fFitELVsTL;
    TProfile* fFitVsPFPartEff;
    
    TH1D*     fDeltaStartPos;
    TH1D*     fDeltaEndPos;
    TH1D*     fTrackDeltaStart;
    TH1D*     fCosTracks;
    TH2D*     fDStartVsDEnd;
    
    TH1D*     fDeltaWiresTrk[3];
    TH1D*     fNumHitsTrk[3];
    TH1D*     fHitWireRatioTrk[3];
    TProfile* fRatioVsDWires[3];
    TH1D*     fTrajDispDiff;
    TH1D*     fTrajDispAng;
    TH1D*     fTrajDispAngRev;
//    TH1D*     fProjLenTrk[3];
    
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
}; // class TrackHitAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
TrackHitAna::TrackHitAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet),
      fPedestalRetrievalAlg(art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider())

{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
TrackHitAna::~TrackHitAna()
{}
   
//-----------------------------------------------------------------------
void TrackHitAna::beginJob()
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
    fHitsByWire[0]            = tfs->make<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]            = tfs->make<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]            = tfs->make<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fPulseHeight[0]           = tfs->make<TH1D>("PulseHeight0",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeight[1]           = tfs->make<TH1D>("PulseHeight1",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeight[2]           = tfs->make<TH1D>("PulseHeight2",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingle[0]     = tfs->make<TH1D>("PulseHeightS0", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingle[1]     = tfs->make<TH1D>("PulseHeightS1", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingle[2]     = tfs->make<TH1D>("PulseHeightS2", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightMulti[0]      = tfs->make<TH1D>("PulseHeightM0", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightMulti[1]      = tfs->make<TH1D>("PulseHeightM1", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightMulti[2]      = tfs->make<TH1D>("PulseHeightM2", "PH (ADC)",  150,  0.,  150.);
    fChi2DOF[0]               = tfs->make<TH1D>("Chi2DOF0",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[1]               = tfs->make<TH1D>("Chi2DOF1",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[2]               = tfs->make<TH1D>("Chi2DOF2",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[0]         = tfs->make<TH1D>("Chi2DOFS0",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[1]         = tfs->make<TH1D>("Chi2DOFS1",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[2]         = tfs->make<TH1D>("Chi2DOFS2",     "Chi2DOF",   502, -1.,  250.);
    fHitMult[0]               = tfs->make<TH1D>("HitMult0",      "# hits",     15,  0.,   15.);
    fHitMult[1]               = tfs->make<TH1D>("HitMult1",      "# hits",     15,  0.,   15.);
    fHitMult[2]               = tfs->make<TH1D>("HitMult2",      "# hits",     15,  0.,   15.);
    fHitCharge[0]             = tfs->make<TH1D>("HitCharge0",    "Charge",   1000,  0., 2000.);
    fHitCharge[1]             = tfs->make<TH1D>("HitCharge1",    "Charge",   1000,  0., 2000.);
    fHitCharge[2]             = tfs->make<TH1D>("HitCharge2",    "Charge",   1000,  0., 2000.);
    fFitWidth[0]              = tfs->make<TH1D>("FitWidth0",     "Charge",    100,  0.,   10.);
    fFitWidth[1]              = tfs->make<TH1D>("FitWidth1",     "Charge",    100,  0.,   10.);
    fFitWidth[2]              = tfs->make<TH1D>("FitWidth2",     "Charge",    100,  0.,   10.);
    fHitSumADC[0]             = tfs->make<TH1D>("SumADC0",       "Sum ADC",  1000,  0.,   50.);
    fHitSumADC[1]             = tfs->make<TH1D>("SumADC1",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[2]             = tfs->make<TH1D>("SumADC2",       "Sum ADC",  1000,  0., 2000.);
    
    fPulseHeightP[0]          = tfs->make<TH1D>("PulseHeightP0",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightP[1]          = tfs->make<TH1D>("PulseHeightP1",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightP[2]          = tfs->make<TH1D>("PulseHeightP2",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingleP[0]    = tfs->make<TH1D>("PulseHeightSP0", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingleP[1]    = tfs->make<TH1D>("PulseHeightSP1", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingleP[2]    = tfs->make<TH1D>("PulseHeightSP2", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightMultiP[0]     = tfs->make<TH1D>("PulseHeightMP0", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightMultiP[1]     = tfs->make<TH1D>("PulseHeightMP1", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightMultiP[2]     = tfs->make<TH1D>("PulseHeightMP2", "PH (ADC)",  150,  0.,  150.);
    fChi2DOFP[0]              = tfs->make<TH1D>("Chi2DOFP0",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFP[1]              = tfs->make<TH1D>("Chi2DOFP1",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFP[2]              = tfs->make<TH1D>("Chi2DOFP2",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingleP[0]        = tfs->make<TH1D>("Chi2DOFSP0",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingleP[1]        = tfs->make<TH1D>("Chi2DOFSP1",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingleP[2]        = tfs->make<TH1D>("Chi2DOFSP2",     "Chi2DOF",   502, -1.,  250.);
    fHitMultP[0]              = tfs->make<TH1D>("HitMultP0",      "# hits",     15,  0.,   15.);
    fHitMultP[1]              = tfs->make<TH1D>("HitMultP1",      "# hits",     15,  0.,   15.);
    fHitMultP[2]              = tfs->make<TH1D>("HitMultP2",      "# hits",     15,  0.,   15.);
    fHitChargeP[0]            = tfs->make<TH1D>("HitChargeP0",    "Charge",   1000,  0., 2000.);
    fHitChargeP[1]            = tfs->make<TH1D>("HitChargeP1",    "Charge",   1000,  0., 2000.);
    fHitChargeP[2]            = tfs->make<TH1D>("HitChargeP2",    "Charge",   1000,  0., 2000.);
    fHitSumADCP[0]            = tfs->make<TH1D>("SumADCP0",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADCP[1]            = tfs->make<TH1D>("SumADCP1",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADCP[2]            = tfs->make<TH1D>("SumADCP2",       "Sum ADC",  1000,  0., 2000.);
    
    fPulseHeightAll[0]        = tfs->make<TH1D>("PulseHeightAll0",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightAll[1]        = tfs->make<TH1D>("PulseHeightAll1",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightAll[2]        = tfs->make<TH1D>("PulseHeightAll2",  "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingleAll[0]  = tfs->make<TH1D>("PulseHeightSAll0", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingleAll[1]  = tfs->make<TH1D>("PulseHeightSAll1", "PH (ADC)",  150,  0.,  150.);
    fPulseHeightSingleAll[2]  = tfs->make<TH1D>("PulseHeightSAll2", "PH (ADC)",  150,  0.,  150.);
    fChi2DOFAll[0]            = tfs->make<TH1D>("Chi2DOFAll0",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFAll[1]            = tfs->make<TH1D>("Chi2DOFAll1",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFAll[2]            = tfs->make<TH1D>("Chi2DOFAll2",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingleAll[0]      = tfs->make<TH1D>("Chi2DOFSAll0",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingleAll[1]      = tfs->make<TH1D>("Chi2DOFSAll1",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingleAll[2]      = tfs->make<TH1D>("Chi2DOFSAll2",     "Chi2DOF",   502, -1.,  250.);
    fHitMultAll[0]            = tfs->make<TH1D>("HitMultAll0",      "# hits",     15,  0.,   15.);
    fHitMultAll[1]            = tfs->make<TH1D>("HitMultAll1",      "# hits",     15,  0.,   15.);
    fHitMultAll[2]            = tfs->make<TH1D>("HitMultAll2",      "# hits",     15,  0.,   15.);
    fHitChargeAll[0]          = tfs->make<TH1D>("HitChargeAll0",    "Charge",   1000,  0., 2000.);
    fHitChargeAll[1]          = tfs->make<TH1D>("HitChargeAll1",    "Charge",   1000,  0., 2000.);
    fHitChargeAll[2]          = tfs->make<TH1D>("HitChargeAll2",    "Charge",   1000,  0., 2000.);
    fHitSumADCAll[0]          = tfs->make<TH1D>("SumADCAll0",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADCAll[1]          = tfs->make<TH1D>("SumADCAll1",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADCAll[2]          = tfs->make<TH1D>("SumADCAll2",       "Sum ADC",  1000,  0., 2000.);

    fNumPFPartTracks          = tfs->make<TH1D>("NPFPartTracks",    ";# Tracks",         100,    0., 100.);
    fNumPFPartTracksL         = tfs->make<TH1D>("NPFPartTracksL",   ";# Tracks",         100,    0., 100.);
    fNumFitTracks             = tfs->make<TH1D>("NFitTracks",       ";# Tracks",         100,    0., 100.);
    fNumTrksPFPart            = tfs->make<TH1D>("NTrksPFPart",      " # Tracks",          10,    0.,  10.);
    fNumTrksPFPartL           = tfs->make<TH1D>("NTrksPFPartL",     " # Tracks",          10,    0.,  10.);
    fNumFitTrksPFPart         = tfs->make<TH1D>("NFitTrksPFPart",   " # Tracks",          10,    0.,  10.);
    fNumFitTrksPFPartL        = tfs->make<TH1D>("NFitTrksPFPartL",  " # Tracks",          10,    0.,  10.);
    fPFPartTrackLen           = tfs->make<TH1D>("PFPartTrackLen",   ";track len",        250,    0., 500.);
    fPFPartTrackLenL          = tfs->make<TH1D>("PFPartTrackLenL",  ";track len",        250,    0., 500.);
    fPFPartEndLenL            = tfs->make<TH1D>("PFPartEndLenL",    ";End Point len",    250,    0., 500.);
    fFitTrackLen              = tfs->make<TH1D>("FitTrackLen",      ";track len",        250,    0., 500.);
    fFitEndLen                = tfs->make<TH1D>("FitEndLen",        ";End Point len",    250,    0., 500.);
    fTrackDeltaLen            = tfs->make<TH1D>("TrackDeltaLen",    ";Delta Len",        500, -250., 250.);
    fFitVsPFPartLen           = tfs->make<TH2D>("FitVsPFPartLen",   ";length;length",    125,    0., 250., 125, 0., 250.);
    fFitELVsTL                = tfs->make<TH2D>("FitELVsTL",        ";length;length",    125,    0., 250., 125, 0., 250.);
    fFitVsPFPartEff           = tfs->make<TProfile>("FitVsPFPart",  ";length(cm)",       100.,   0., 500.,      0., 1.1);
    
    fDeltaStartPos            = tfs->make<TH1D>("DeltaStartPos",    ";delta(cm)",        100,    0.,  50.);
    fDeltaEndPos              = tfs->make<TH1D>("DeltaEndPos",      ";delta(cm)",        100,    0.,  50.);
    fTrackDeltaStart          = tfs->make<TH1D>("TrackDeltaStart",  ";delta(cm)",        100,    0.,  50.);
    fCosTracks                = tfs->make<TH1D>("CosTracks",        ";cos(theta)",       101,    0.,   1.01);
    
    fDStartVsDEnd             = tfs->make<TH2D>("DStartVsDEnd",     ";delta;delta",       50,    0.,  50., 50, 0., 50.);
    
    fDeltaWiresTrk[0]         = tfs->make<TH1D>("DltaWiresTrk0",    ";deltaWires",       250,    0., 1000.);
    fDeltaWiresTrk[1]         = tfs->make<TH1D>("DltaWiresTrk1",    ";deltaWires",       250,    0., 1000.);
    fDeltaWiresTrk[2]         = tfs->make<TH1D>("DltaWiresTrk2",    ";deltaWires",       250,    0., 1000.);
    fNumHitsTrk[0]            = tfs->make<TH1D>("NumHitsTrk0",      ";# hits",           250,    0., 1000.);
    fNumHitsTrk[1]            = tfs->make<TH1D>("NumHitsTrk1",      ";# hits",           250,    0., 1000.);
    fNumHitsTrk[2]            = tfs->make<TH1D>("NumHitsTrk2",      ";# hits",           250,    0., 1000.);
    fHitWireRatioTrk[0]       = tfs->make<TH1D>("HitWireRat0",      ";Ratio",            100,    0.,    2.);
    fHitWireRatioTrk[1]       = tfs->make<TH1D>("HitWireRat1",      ";Ratio",            100,    0.,    2.);
    fHitWireRatioTrk[2]       = tfs->make<TH1D>("HitWireRat2",      ";Ratio",            100,    0.,    2.);
    fRatioVsDWires[0]         = tfs->make<TProfile>("RatVsDWire0",  ";DeltaWires;Ratio", 200,    0.,  200., 0., 2.);
    fRatioVsDWires[1]         = tfs->make<TProfile>("RatVsDWire1",  ";DeltaWires;Ratio", 200,    0.,  200., 0., 2.);
    fRatioVsDWires[2]         = tfs->make<TProfile>("RatVsDWire2",  ";DeltaWires;Ratio", 200,    0.,  200., 0., 2.);
    
    fTrajDispDiff             = tfs->make<TH1D>("TrajPointDisp",    ";disp",             200,  -10.,   10.);
    fTrajDispAng              = tfs->make<TH1D>("TrajPointAng",     ";cos(ang)",         100,   -1.,    1.);
    
    // zero out the event counter
    fNumEvents = 0;
}
   
//-----------------------------------------------------------------------
void TrackHitAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void TrackHitAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fHitProducerLabel        = p.get< std::string >("HitModuleLabel",          "gauss");
    fPFParticleProducerLabel = p.get< std::string >("PFParticleProducerLabel", "cluster3d");
    fTrackProducerLabel      = p.get< std::string >("TrackProducerLabel",      "trackkalmanhit");

    return;
}

//-----------------------------------------------------------------------
void TrackHitAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    fNumEvents++;
    
    // The game plan for this module is to look at hits associated to tracks
    // To do this we need a valid track collection for those we are hoping to look at
    art::Handle<std::vector<recob::Track> > trackHandle;
    event.getByLabel(fTrackProducerLabel, trackHandle);
    
    // Get a local mapping between tracks and hits
    std::map<int,std::map<size_t,HitPtrVec>> trackHitVecMap;
    
    if (trackHandle.isValid())
    {
        // Recover the collection of associations between tracks and hits
        art::FindManyP<recob::Hit> trackAssns(trackHandle, event, fTrackProducerLabel);
        
        for(size_t trackIdx = 0; trackIdx < trackHandle->size(); trackIdx++)
        {
            art::Ptr<recob::Track> track(trackHandle,trackIdx);
            
            std::map<size_t,HitPtrVec>& viewHitMap = trackHitVecMap[track.key()];
            
            // Recover the associated hits
            HitPtrVec trackHitVec = trackAssns.at(track.key());
            
            for(int viewIdx = 0; viewIdx < 3; viewIdx++)
            {
                int numHits = std::accumulate(trackHitVec.begin(),trackHitVec.end(),int(0),[viewIdx](int sum, const auto& hit){return sum += hit->View() == viewIdx ? 1 : 0;});
                
                viewHitMap[viewIdx].resize(numHits);
                
                std::copy_if(trackHitVec.begin(),trackHitVec.end(),viewHitMap[viewIdx].begin(),[viewIdx](const auto& hit){return hit->View() == viewIdx;});
            }
            
            // It is helpful if the hits are in time order (by view)
//            std::sort(trackHitVec.begin(),trackHitVec.end(),[](const recob::Hit* left, const recob::Hit* right) {return left->PeakTime() < right->PeakTime();});
        }
    }
    
    // The game plan for this module is to look at hits associated to tracks
    // To do this we need a valid track collection for those we are hoping to look at
    art::Handle<std::vector<recob::Track> > pfTrackHandle;
    event.getByLabel(fPFParticleProducerLabel, pfTrackHandle);
    
    // Get a local mapping between tracks and hits
    std::map<int,std::map<size_t,HitPtrVec>> pfTrackHitVecMap;
    
    if (pfTrackHandle.isValid())
    {
        // Recover the collection of associations between tracks and hits
        art::FindManyP<recob::Hit> trackAssns(pfTrackHandle, event, fPFParticleProducerLabel);
        
        for(size_t trackIdx = 0; trackIdx < pfTrackHandle->size(); trackIdx++)
        {
            art::Ptr<recob::Track> track(pfTrackHandle,trackIdx);
            
            std::map<size_t,HitPtrVec>& viewHitMap = pfTrackHitVecMap[track.key()];
            
            // Recover the associated hits
            HitPtrVec trackHitVec = trackAssns.at(track.key());
            
            for(int viewIdx = 0; viewIdx < 3; viewIdx++)
            {
                int numHits = std::accumulate(trackHitVec.begin(),trackHitVec.end(),int(0),[viewIdx](int sum, const auto& hit){return sum += hit->View() == viewIdx ? 1 : 0;});
                
                viewHitMap[viewIdx].resize(numHits);
                
                std::copy_if(trackHitVec.begin(),trackHitVec.end(),viewHitMap[viewIdx].begin(),[viewIdx](const auto& hit){return hit->View() == viewIdx;});
            }
            
            // It is helpful if the hits are in time order (by view)
            //            std::sort(trackHitVec.begin(),trackHitVec.end(),[](const recob::Hit* left, const recob::Hit* right) {return left->PeakTime() < right->PeakTime();});
        }
    }
    
    // Now pass through this map and do some matching
    for(const auto& trackHitVecMapItr : trackHitVecMap)
    {
        for(const auto& viewHitPair : trackHitVecMapItr.second)
        {
            const HitPtrVec& trackHitVec = viewHitPair.second;
        
            // Loop the hits and make some plots
            for(const auto& hitPtr : trackHitVec)
            {
                // Extract interesting hit parameters
                const geo::WireID& wireID   = hitPtr->WireID();
                float              chi2DOF  = std::min(hitPtr->GoodnessOfFit(),float(249.8));
//              int                numDOF   = hitPtr->DegreesOfFreedom();
                int                hitMult  = hitPtr->Multiplicity();
//              int                hitIdx   = hitPtr->LocalIndex();
                float              charge   = hitPtr->Integral();
                float              sumADC   = hitPtr->SummedADC();
                float              hitPH    = std::min(hitPtr->PeakAmplitude(),float(249.8));
                float              hitSigma = hitPtr->RMS();
//              raw::TDCtick_t     hitStart = hitPtr->StartTick();
//              raw::TDCtick_t     hitEnd   = hitPtr->EndTick();
//              float              hitTime  = hitPtr->PeakTime();
            
                size_t             view     = wireID.Plane;
                size_t             wire     = wireID.Wire;
                
                if (view != viewHitPair.first) std::cout << "*******>> View mismatch: " << view << ", " << viewHitPair.first << std::endl;
            
                fHitsByWire[view]->Fill(wire,1.);
                fPulseHeight[view]->Fill(hitPH, 1.);
                fChi2DOF[view]->Fill(chi2DOF, 1.);
                fHitMult[view]->Fill(hitMult, 1.);
                fHitCharge[view]->Fill(charge, 1.);
                fFitWidth[view]->Fill(std::min(float(9.99),hitSigma), 1.);
                fHitSumADC[view]->Fill(sumADC, 1.);
            
                if (hitMult == 1)
                {
                    fPulseHeightSingle[view]->Fill(hitPH, 1.);
                    fChi2DOFSingle[view]->Fill(chi2DOF, 1.);
                }
                else
                    fPulseHeightMulti[view]->Fill(hitPH, 1.);
            }
        }
    }
    
    // Now pass through this map and do some matching
    for(const auto& trackHitVecMapItr : pfTrackHitVecMap)
    {
        int trackIdx = trackHitVecMapItr.first;
        
        for(const auto& viewHitPair : trackHitVecMapItr.second)
        {
            const HitPtrVec& trackHitVec = viewHitPair.second;
        
            art::Ptr<recob::Track> track(pfTrackHandle, trackIdx);
        
            if (length(track.get()) < 5.) continue;
        
            // Loop the hits and make some plots
            for(const auto& hitPtr : trackHitVec)
            {
                // Extract interesting hit parameters
                const geo::WireID& wireID   = hitPtr->WireID();
                float              chi2DOF  = std::min(hitPtr->GoodnessOfFit(),float(249.8));
//              int                numDOF   = hitPtr->DegreesOfFreedom();
                int                hitMult  = hitPtr->Multiplicity();
//              int                hitIdx   = hitPtr->LocalIndex();
                float              charge   = hitPtr->Integral();
                float              sumADC   = hitPtr->SummedADC();
                float              hitPH    = std::min(hitPtr->PeakAmplitude(),float(249.8));
//              float              hitSigma = hitPtr->RMS();
//              raw::TDCtick_t     hitStart = hitPtr->StartTick();
//              raw::TDCtick_t     hitEnd   = hitPtr->EndTick();
//              float              hitTime  = hitPtr->PeakTime();
                size_t             view     = wireID.Plane;
            
                fPulseHeightP[view]->Fill(hitPH, 1.);
                fChi2DOFP[view]->Fill(chi2DOF, 1.);
                fHitMultP[view]->Fill(hitMult, 1.);
                fHitChargeP[view]->Fill(charge, 1.);
                fHitSumADCP[view]->Fill(sumADC, 1.);
            
            
                if (hitMult == 1)
                {
                    fPulseHeightSingleP[view]->Fill(hitPH, 1.);
                    fChi2DOFSingleP[view]->Fill(chi2DOF, 1.);
                }
                else
                    fPulseHeightMultiP[view]->Fill(hitPH, 1.);
            }
        }
    }
    
    // Make a pass through all hits to make contrasting plots
    art::Handle< std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);

    if (hitHandle.isValid())
    {
        for(const auto& hit : *hitHandle)
        {
            // Extract interesting hit parameters
            const geo::WireID& wireID   = hit.WireID();
            float              chi2DOF  = std::min(hit.GoodnessOfFit(),float(249.8));
            int                numDOF   = hit.DegreesOfFreedom();
            int                hitMult  = hit.Multiplicity();
            int                hitIdx   = hit.LocalIndex();
            float              charge   = hit.Integral();
            float              sumADC   = hit.SummedADC();
            float              hitPH    = std::min(hit.PeakAmplitude(),float(249.8));
            float              hitSigma = hit.RMS();
            raw::TDCtick_t     hitStart = hit.StartTick();
            raw::TDCtick_t     hitEnd   = hit.EndTick();
            float              hitTime  = hit.PeakTime();
            
            size_t             view     = wireID.Plane;
            size_t             wire     = wireID.Wire;
            
            if (chi2DOF < -10000.)
            {
                std::cout << "======> Hit view/wire: " << view << "/" << wire << ", Peak time, width: " << hitTime << ", " << hitSigma << std::endl;
                std::cout << "        chi2DOF: " << chi2DOF << ", dof: " << numDOF << ", hitMult: " << hitMult << ", hitIdx: " << hitIdx << std::endl;
                std::cout << "        charge: " << charge << ", sumADC: " << sumADC << ", hitPH: " << hitPH << ", start/end: " << hitStart << "/" << hitEnd << std::endl;
            }
            
            fPulseHeightAll[view]->Fill(hitPH, 1.);
            fChi2DOFAll[view]->Fill(chi2DOF, 1.);
            fHitMultAll[view]->Fill(hitMult, 1.);
            fHitChargeAll[view]->Fill(charge, 1.);
            fHitSumADCAll[view]->Fill(sumADC, 1.);
            
            if (hitMult == 1)
            {
                fPulseHeightSingleAll[view]->Fill(hitPH, 1.);
                fChi2DOFSingleAll[view]->Fill(chi2DOF, 1.);
            }
                
        }
    }
    
    // Now we want to try to look at tracks associated to PFParticles
    art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
    event.getByLabel(fPFParticleProducerLabel, pfParticleHandle);
    
    if (pfParticleHandle.isValid())
    {
        // Recover the collection of associations between PFParticles and Tracks from the PFParticle producer
        art::FindManyP<recob::Track> firstTrackAssns(pfParticleHandle, event, fPFParticleProducerLabel);
        
        // Now get the track associations from the track producer
        art::FindManyP<recob::Track> secondTrackAssns(pfParticleHandle, event, fTrackProducerLabel);
        
        // Also grab a track handle from the PFParticle producer
        art::Handle<std::vector<recob::Track>> pfTrackHandle;
        event.getByLabel(fPFParticleProducerLabel, pfTrackHandle);
        
        if (firstTrackAssns.size() > 0 && secondTrackAssns.size() > 0 && pfTrackHandle.isValid())
        {
            int nPandoraTracks(0);
            int nPandoraTracksL(0);
            int nFitTracks(0);
            
            for(size_t pfParticleIdx = 0; pfParticleIdx < pfParticleHandle->size(); pfParticleIdx++)
            {
                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);

                // Get tracks made with the PFParticle
                std::vector<art::Ptr<recob::Track>> pfPartTrackVec = firstTrackAssns.at(pfParticle.key());
                
                // Get tracks from fitter
                std::vector<art::Ptr<recob::Track>> fitTrackVec = secondTrackAssns.at(pfParticle.key());
                
                if (pfParticle->PdgCode() == 13)
                {
                    int nTrksPFPart(0);
                    int nTrksPFPartL(0);
                    
                    // Se are meant to be able to recover track <--> hit associations from pandora
                    art::FindManyP<recob::Hit> pfTrackHitAssns(pfTrackHandle, event, fPFParticleProducerLabel);
                    
                    for(size_t trkIdx = 0; trkIdx < pfPartTrackVec.size(); trkIdx++)
                    {
                        art::Ptr<recob::Track> track = pfPartTrackVec.at(trkIdx);
                        
                        // Recover start/end points and directions
                        const TVector3& pfPartStart    = track->Vertex();
                        const TVector3& pfPartEnd      = track->End();
                        const TVector3& pfPartStartDir = track->VertexDirection();
//                        const TVector3& pfPartEndDir   = track->EndDirection();
                        
                        double pfPartLen    = length(track.get());
                        double pfPartEndLen = (pfPartEnd - pfPartStart).Mag();
                        
                        fPFPartTrackLen->Fill(std::min(499.5,pfPartLen),1.);
                        fPFPartEndLenL->Fill(std::min(499.5,pfPartEndLen),1.);
                        nTrksPFPart++;
                        
                        if (pfPartLen > 5.)
                        {
                            nPandoraTracksL++;
                            fPFPartTrackLenL->Fill(std::min(499.5,pfPartLen));
                            nTrksPFPartL++;
                        }
                        
                        bool   foundMatch(false);
                        double deltaLen(99999.);
                        double matchLen(0.);
                        
                        // Have to make local copies for the next step...
                        TVector3  fitTrkStart;
                        TVector3  fitTrkEnd;
                        TVector3  fitTrkStartDir;
                        TVector3  fitTrkEndDir;
                        
                        std::map<size_t,HitPtrVec> viewHitMap;
//                        const recob::Track*        matchTrack = 0;
                        
                        for(const auto& fitTrk : fitTrackVec)
                        {
                            double fitTrkLen = length(fitTrk.get());
                            
                            if (fabs(pfPartLen - fitTrkLen) < fabs(deltaLen))
                            {
                                foundMatch     = true;
                                deltaLen       = pfPartLen - fitTrkLen;
                                matchLen       = fitTrkLen;
                                fitTrkStart    = fitTrk->Vertex();
                                fitTrkStartDir = fitTrk->VertexDirection();
                                fitTrkEnd      = fitTrk->End();
                                fitTrkEndDir   = fitTrk->EndDirection();
                                viewHitMap     = trackHitVecMap[fitTrk.key()];
//                                matchTrack     = fitTrk.get();
                            }
                        }
                        
                        double trkMatch(0.);
                        
                        if (foundMatch)
                        {
                            deltaLen = std::max(-249.5,std::min(249.5,deltaLen));
                            trkMatch = 1.;
                            
                            fTrackDeltaLen->Fill(deltaLen, 1.);
                            fFitVsPFPartLen->Fill(std::min(249.5,pfPartLen), std::min(249.5,matchLen), 1.);
                            
                            TVector3 forwardDiff = pfPartStart - fitTrkStart;
                            TVector3 reverseDiff = pfPartStart - fitTrkEnd;
                            
                            bool flipped(false);
                            
                            // If end of Fit track closer then reverse
                            if (reverseDiff.Mag2() < forwardDiff.Mag2())
                            {
                                TVector3 tempPos = fitTrkStart;
                                
                                fitTrkStart = fitTrkEnd;
                                fitTrkEnd   = tempPos;
                                
                                TVector3 tempDir = fitTrkStartDir;
                                
                                fitTrkStartDir = -fitTrkEndDir;
                                fitTrkEndDir   = -tempDir;
                                
                                flipped = true;
                            }
                            
                            double startDist = std::min(49.9,(pfPartStart - fitTrkStart).Mag());
                            double endDist   = std::min(49.9,(pfPartEnd   - fitTrkEnd  ).Mag());
                            double cosTrkAng = std::max(0.01,pfPartStartDir.Dot(fitTrkStartDir));
//                            double cosTrkEnd = pfPartEndDir.Dot(fitTrkEndDir);
                            
                            fDeltaStartPos->Fill(startDist, 1.);
                            fDeltaEndPos->Fill(endDist, 1.);
                            fCosTracks->Fill(cosTrkAng, 1.);

                            if (flipped) fTrackDeltaStart->Fill(endDist, 1.);
                            else         fTrackDeltaStart->Fill(startDist, 1.);
                            
                            fDStartVsDEnd->Fill(startDist, endDist, 1.);
                            
                            // Look at number of wires crossed in each plane
                            for(size_t planeIdx = 0; planeIdx < 3; planeIdx++)
                            {
                                try
                                {
                                    geo::WireID startWire = fGeometry->NearestWireID(fitTrkStart, planeIdx);
                                    geo::WireID endWire   = fGeometry->NearestWireID(fitTrkEnd,   planeIdx);
                                
                                    int    deltaWires = fabs(int(endWire.Wire - startWire.Wire));
                                    int    numHits    = viewHitMap[planeIdx].size();
                                    double hitRatio   = deltaWires > 0 ? double(numHits)/double(deltaWires) : 0.;
                                
                                    fDeltaWiresTrk[planeIdx]->Fill(deltaWires, 1.);
                                    fNumHitsTrk[planeIdx]->Fill(numHits, 1.);
                                    if (deltaWires > 7) fHitWireRatioTrk[planeIdx]->Fill(std::min(1.99,hitRatio), 1.);
                                    fRatioVsDWires[planeIdx]->Fill(std::min(deltaWires,199), std::min(1.99,hitRatio));
                                }
                                catch(...) {}
                            }
/*
                            // Quick check to look for tracks which have trajectories gone awry
                            TVector3 lastPos(matchTrack->LocationAtPoint(0));
                            TVector3 lastDir(0.,0.,0.);

                            for(size_t trajIdx = 0; trajIdx < matchTrack->NumberTrajectoryPoints(); trajIdx++)
                            {
                                const TVector3& pos = matchTrack->LocationAtPoint(trajIdx);
                                TVector3        dir = pos - lastPos;
                                
                                if (dir.Mag2() > 0.) dir.SetMag(1.);
                                
                                if (dir.Dot(lastDir) < 0.)
                                {
                                    std::cout << "****>> Trajectory flips for track: " << matchTrack->ID() << ", angle: " << dir.Dot(lastDir) << std::endl;
                                    std::cout << "       idx: " << trajIdx << ", dir: " << dir.X() << "," << dir.Y() << "," << dir.Z() << ", lastDir: " << lastDir.X() << "," << lastDir.Y() << "," << lastDir.Z() << std::endl;
                                }
                                
                                lastDir = dir;
                                lastPos = pos;
                            }
 */
                        }

                        fFitVsPFPartEff->Fill(pfPartLen, trkMatch);
                    }
                    
                    fNumTrksPFPart->Fill(nTrksPFPart, 1.);
                    fNumTrksPFPartL->Fill(nTrksPFPartL, 1.);
                    
                    int nFitTrksPFPart(fitTrackVec.size());
                    
                    for(size_t trkIdx = 0; trkIdx < fitTrackVec.size(); trkIdx++)
                    {
                        art::Ptr<recob::Track> fitTrack = fitTrackVec.at(trkIdx);

                        double fitLen = std::min(499.5,length(fitTrack.get()));
                        double fitEndLen = std::min(499.5,(fitTrack->Vertex() - fitTrack->End()).Mag());
                        
                        fFitTrackLen->Fill(fitLen,1.);
                        fFitEndLen->Fill(fitEndLen,1.);
                        fFitELVsTL->Fill(fitLen,fitEndLen,1.);
                        
                        if (fitLen > 30. && fitLen > 1.05 * fitEndLen)
                        {
                            std::cout << "*******>> Fit length " << fitLen << " > 1.05xEndLen " << fitEndLen << std::endl;
                            std::cout << "          Run: " << fRun << ", event: " << fEvent << ", track: " << fitTrack->ID() << ", # tracks: " << fitTrackVec.size() << std::endl;
                        }
                        // Quick check to look for tracks which have trajectories gone awry
                        TVector3 lastTrajPos(fitTrack->LocationAtPoint(0));
                        TVector3 lastTrajDir(fitTrack->DirectionAtPoint(0));
                        TVector3 lastPntDir(0.,0.,0.);
                        std::vector<std::tuple<size_t,double,double>> flipPointVec;
                        
                        for(size_t trajIdx = 0; trajIdx < fitTrack->NumberTrajectoryPoints(); trajIdx++)
                        {
                            const TVector3& pos  = fitTrack->LocationAtPoint(trajIdx);
                            TVector3        dir  = pos - lastTrajPos;
                            double          disp = dir.Mag();
                            
                            if (disp > 0.) dir.SetMag(1.);
                            
                            if (dir.Dot(lastPntDir) < 0.)
                            {
                                disp = -disp;
                                flipPointVec.push_back(std::make_tuple(trajIdx,dir.Dot(lastPntDir),disp));
                            }
                            
                            fTrajDispDiff->Fill(std::max(-9.99,std::min(9.99,disp)), 1.);
                            fTrajDispAng->Fill(std::max(-0.99,std::min(0.99,dir.Dot(lastPntDir))), 1.);
                            
                            lastPntDir  = dir;
                            lastTrajPos = pos;
                            lastTrajDir = fitTrack->DirectionAtPoint(trajIdx);
                        }
                        
                        std::map<size_t,HitPtrVec> viewHitMap = trackHitVecMap[fitTrack.key()];
                        int                        numTrackHits(0);
                        
                        for(const auto& itr : viewHitMap) numTrackHits += itr.second.size();
                    
                        std::cout << ">>>> Trajectory flips for track: " << fitTrack->ID() << ", with " << flipPointVec.size() << " reversals" << std::endl;
                        std::cout << "     Number trajectory points: " << fitTrack->NumberTrajectoryPoints() << ", # hits: " << numTrackHits << std::endl;
//                        std::cout << "     Points/angle/trajAng: ";
//                        for(const auto& point : flipPointVec) std::cout << " " << std::get<0>(point) << "/" << std::get<1>(point) << "/" << std::get<2>(point);
//                        std::cout << std::endl;
                    }
                    
                    fNumFitTrksPFPart->Fill(nFitTrksPFPart, 1.);
                    if (nTrksPFPartL > 0) fNumFitTrksPFPartL->Fill(nFitTrksPFPart, 1.);
                    
                    nPandoraTracks += pfPartTrackVec.size();
                    nFitTracks     += fitTrackVec.size();
                }
/*
                else if (fitTrackVec.size() > 0)
                {
                    std::cout << "++ PFParticle: " << pfParticleIdx << ", pdg: " << pfParticle->PdgCode() << ", # tracks: " << pfPartTrackVec.size() << ", fit tracks: " << fitTrackVec.size() << std::endl;
                    std::cout << "   Tracks: ";
                    for(size_t trkIdx = 0; trkIdx < pfPartTrackVec.size(); trkIdx++) std::cout << " " << length(pfPartTrackVec.at(trkIdx).get());
                    for(size_t trkIdx = 0; trkIdx < fitTrackVec.size(); trkIdx++)    std::cout << " " << length(fitTrackVec.at(trkIdx).get());
                    std::cout << std::endl;
                }
 */
            }
            
            std::cout << "~~ # pandora tracks: " << nPandoraTracks << ", long: " << nPandoraTracksL << ", # fit tracks: " << nFitTracks << std::endl;
            fNumPFPartTracks->Fill(nPandoraTracks, 1.);
            fNumPFPartTracksL->Fill(nPandoraTracksL, 1.);
            fNumFitTracks->Fill(nFitTracks, 1.);
        }
        
    }

    return;
}
    
void TrackHitAna::endJob()
{
    return;
}
    
// Length of reconstructed track.
//----------------------------------------------------------------------------
double TrackHitAna::length(const recob::Track* track)
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

// This macro has to be defined for this module to be invoked from a
// .fcl file; see TrackHitAna.fcl for more information.
DEFINE_ART_MODULE(TrackHitAna)

} // namespace TrackHitAna

#endif // TrackHitAna_module
