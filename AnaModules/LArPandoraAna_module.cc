// LArPandoraAna_module.cc
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

#ifndef LArPandoraAna_Module
#define LArPandoraAna_Module

// LArSoft includes
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/AnalysisBase/CosmicTag.h"
#include "larcore/Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "lardata/MCBase/MCHitCollection.h"

//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

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

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

namespace LArPandoraAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class LArPandoraAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit LArPandoraAna(fhicl::ParameterSet const& pset);
    virtual ~LArPandoraAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

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

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.
    typedef std::map< const recob::Hit*, std::set<int> > HitToParticleMap;
    typedef std::map< int, std::set<const recob::Hit*> > ParticleToHitMap;
    typedef std::vector< art::Ptr<recob::Hit> >          HitVector;

    // The parameters we'll read from the .fcl file.
    std::string fSimulationProducerLabel;    // The name of the producer that tracked simulated particles through the detector
    std::string fMcHitCollectionModuleLabel; //> The name of the producer of MCHits
    std::string fPFParticleProducerLabel;    // The name of the produder of the PFParticle hierarchy
    std::string fHitProducerLabel;           // The name of the producer that created hits
    std::string fClusterProducerLabel;       // The name of the producer that created clusters
    std::string fTrackProducerLabel;         // The name of the producer that created the tracks
    std::string fCosmicProducerLabel;        // The name of the producer that created cosmic tags
    std::string fFlashProducerLabel;         // The name of the producer that created flash tags
    int fSelectedPDG;                        // PDG code of particle we'll focus on

    // Pointers to the histograms we'll create. 
    TH1D* fPDGCodeHist;
    TH1D* fMomentumHist;
    TH1D* fTrackLengthHist;
      
    TH1D* fNMcParticles;
    TH1D* fEneMcPart;
    TH1D* fEneSelected;
    TH1D* fEneWithHits;
    TH1D* fMcProcess;
    
    TH1D* fMcTracksWithHits;
    TH1D* fMcMuonsWithHits;
    TH1D* fMcMuonHitsPerTrack;
    
    TH1D* fMCTrackNumHits;
    TH1D* fMCTrackPhi;
    TH1D* fMCTrackTheta;
    TH1D* fMCTrackThetaXZ;
    
    TH1D* fNumPFParticles;
    TH1D* fNumPFParticlesTagged;
    TH1D* fNumTracksRecon;
    TH1D* fNumMuonsRecon;
    TH1D* fNumPartPerTrk;
    
    TH1D* fNumMcPerPFPart;
    TH1D* fNumClusPerPFPart;
    TH1D* fNumHitsTotal;
    TH1D* fNumHitsPerPFPart;
    TH1D* fNumHitsGtOneMC;
    TH1D* fNumHitsNotMC;
    TH1D* fNumHitsNoise;
    TH1D* fNumHitsBestMC;
    TH1D* fNumNoisePlane[3];
    
    TProfile* fHitEfficVsHits;
    TProfile* fHitEfficVsTheta;
    TProfile* fHitEfficVsPhi;
    
    TProfile* fClusEfficVsHits;
    TProfile* fClusEfficVsTheta;
    TProfile* fClusEfficVsPhi;
    
    TH1D*     fNumHitsPerTrack;
    TH1D*     fNumMCHitsPerTrack;
    TH1D*     fNumClusHitsPerTrk;
    TH1D*     fNumTracksPerClus;
    TProfile* fTrackEfficVsHits;
    TProfile* fTrackEfficVsTheta;
    TProfile* fTrackEfficVsPhi;
    
    TH1D* fHitEfficiency;
    TH1D* fWrongEfficiency;
    TH1D* fHitPurity;
    TH1D* fHitEfficPurity;
    
    TH1D* fHitEfficiency0;
    TH1D* fWrongEfficiency0;
    TH1D* fHitPurity0;
    TH1D* fHitEfficPurity0;
    
    TH1D* fTotHitsPerCluster;
    TH1D* fCorrectPerCluster;
    TH1D* fWrongPerCluster;
    TH1D* fActualHits;
    
    TH1D* fCosmicScoreAll;
    TH1D* fCosmicScore;
    TH1D* fCosmicTag;
    TH1D* fFlashScoreAll;
    TH1D* fFlashScore;
    TH1D* fFlashTag;
    TH1D* fOrScore;
    TH1D* fOrTag;
    
    TH1D* fCosmicScoreNuAll;
    TH1D* fCosmicScoreNu;
    TH1D* fCosmicTagNu;
    TH1D* fFlashScoreNuAll;
    TH1D* fFlashScoreNu;
    TH1D* fFlashTagNu;
    TH1D* fOrScoreNu;
    TH1D* fOrTagNu;
    
    TH1D* fHitTotalCharge[6];
    TH1D* fHitPulseHeight[6];
    TH1D* fHitPHClose[6];
    TH1D* fHitMult[6];
    
    TH1D* fDeltaPeakTime[3];
    TH1D* fPullPeakTime[3];

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;

    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    // Note: old-style C++ arrays are considered obsolete. However,
    // to create simple n-tuples, we still need to use them. 
    double fStartXYZT[4];
    double fEndXYZT[4];
    double fStartPE[4];
    double fEndPE[4];

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*           fGeometry;       // pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;

    double                                       fElectronsToGeV; // conversion factor

}; // class LArPandoraAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
LArPandoraAna::LArPandoraAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
LArPandoraAna::~LArPandoraAna()
{}
   
//-----------------------------------------------------------------------
void LArPandoraAna::beginJob()
{
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    double detectorLength = fGeometry->DetLength(); 

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
  
    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes. 

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fPDGCodeHist        = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist       = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist    = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detectorLength);
      
    fNMcParticles       = tfs->make<TH1D>("NumMcParts", ";log(Num);", 100, -1., 9.);
    fEneMcPart          = tfs->make<TH1D>("EneMcPart",  ";log(Ene);", 100, -7., 3.);
    fEneSelected        = tfs->make<TH1D>("EneSelect",  ";log(Ene);", 100, -7., 3.);
    fEneWithHits        = tfs->make<TH1D>("EneWithHit", ";log(Ene);", 100, -7., 3.);
    fMcProcess          = tfs->make<TH1D>("McProcess",  ";Process;",  100,  0, 100);
    
    fMcTracksWithHits   = tfs->make<TH1D>("McTrksWHits",      ";# hits",      50,    0.0,    50.);
    fMcMuonsWithHits    = tfs->make<TH1D>("McMuonsWHits",     ";# hits",      50,    0.0,    50.);
    fMcMuonHitsPerTrack = tfs->make<TH1D>("McHitsPerMu",      ";# hits",      50,    0.0,     4.);
    
    fMCTrackNumHits     = tfs->make<TH1D>("MCTrackNumHits",   ";log(# hits)", 50,    0.0,    4.0);
    fMCTrackPhi         = tfs->make<TH1D>("MCTrackPhi",       ";degrees",     40, -180.0,  180.0);
    fMCTrackTheta       = tfs->make<TH1D>("MCTrackTheta",     ";theta",       20,    0.0,  180.0);
    fMCTrackThetaXZ     = tfs->make<TH1D>("MCTrackThetaXZ",   ";theta",       40, -180.0,  180.0);
    
    fNumPFParticles     = tfs->make<TH1D>("NumPFParticles",   ";# PFParticles",  100,    0.0,   100.);
    fNumTracksRecon     = tfs->make<TH1D>("NumTracksRecon",   ";# PFParticles",  100,    0.0,   100.);
    fNumMuonsRecon      = tfs->make<TH1D>("NumMuonsRecon",    ";# PFParticles",  50,    0.0,   50.);
    fNumPartPerTrk      = tfs->make<TH1D>("NumPartPerTrack",  ";# PFParticles",  10,    0.0,   10.);
    
    fNumMcPerPFPart     = tfs->make<TH1D>("NumMcPerPFPart",   ";# MCParticles",  100,   0.0,   200.);
    fNumClusPerPFPart   = tfs->make<TH1D>("NumClusPerPFPart", ";# Clusters",       5,   0.0,     5.);
    fNumHitsPerPFPart   = tfs->make<TH1D>("NumHitsPerPFPart", ";log(# Hits)",     50,   0.0,     4.);
    fNumHitsTotal       = tfs->make<TH1D>("NumHitsTotal",     ";log(# Hits)",    100,   0.0,     5.);
    fNumHitsGtOneMC     = tfs->make<TH1D>("NumHitsGtOneMC",   ";log(# Hits)",    100,   0.0,     5.);
    fNumHitsNotMC       = tfs->make<TH1D>("NumHitsNotMC",     ";log(# Hits)",    100,   0.0,     5.);
    fNumHitsNoise       = tfs->make<TH1D>("NumHitsNoise",     ";log(# Hits)",    100,   0.0,     5.);
    fNumHitsBestMC      = tfs->make<TH1D>("NumHitsBestMC",    ";log(# Hits)",    100,   0.0,     5.);
    fNumNoisePlane[0]   = tfs->make<TH1D>("NumHitsNoiseU",    ";log(# Hits)",    100,   0.0,     5.);
    fNumNoisePlane[1]   = tfs->make<TH1D>("NumHitsNoiseV",    ";log(# Hits)",    100,   0.0,     5.);
    fNumNoisePlane[2]   = tfs->make<TH1D>("NumHitsNoiseW",    ";log(# Hits)",    100,   0.0,     5.);
    
    fHitEfficVsHits     = tfs->make<TProfile>("HitEfficVsHits",  ";log(# hits)",     20,    0.5,   4.5, 0., 1.1);
    fHitEfficVsTheta    = tfs->make<TProfile>("HitEfficVsTheta", ";Theta(degrees)",  20,    0.0, 180.0, 0., 1.1);
    fHitEfficVsPhi      = tfs->make<TProfile>("HitEfficVsPhi",   ";Phi(degrees)",    40, -180.0, 180.0, 0., 1.1);
    
    fClusEfficVsHits    = tfs->make<TProfile>("ClusEfficVsHits",  ";log(# hits)",     20,    0.5,   4.5, 0., 1.1);
    fClusEfficVsTheta   = tfs->make<TProfile>("ClusEfficVsTheta", ";Theta(degrees)",  20,    0.0, 180.0, 0., 1.1);
    fClusEfficVsPhi     = tfs->make<TProfile>("ClusEfficVsPhi",   ";Phi(degrees)",    40, -180.0, 180.0, 0., 1.1);

    fNumHitsPerTrack    = tfs->make<TH1D>("NumHitsPerTrack",       ";log(# hits)",     50,    0.0,   4.0);
    fNumMCHitsPerTrack  = tfs->make<TH1D>("NumMCHitsPerTrack",     ";log(# hits)",     50,    0.0,   4.0);
    fNumClusHitsPerTrk  = tfs->make<TH1D>("NumClusHitsPerTrk",     ";log(# hits)",     50,    0.0,   4.0);
    fNumTracksPerClus   = tfs->make<TH1D>("NumTracksPerClus",      ";# Tracks",        20,    0.0,  20.0);
    fTrackEfficVsHits   = tfs->make<TProfile>("TrackEfficVsHits",  ";log(# hits)",     20,    0.5,   4.5, 0., 1.1);
    fTrackEfficVsTheta  = tfs->make<TProfile>("TrackEfficVsTheta", ";Theta(degrees)",  20,    0.0, 180.0, 0., 1.1);
    fTrackEfficVsPhi    = tfs->make<TProfile>("TrackEfficVsPhi",   ";Phi(degrees)",    40, -180.0, 180.0, 0., 1.1);
    
    fHitEfficiency      = tfs->make<TH1D>("HitEffic",         ";efficiency",  51,    0.0,    1.02);
    fWrongEfficiency    = tfs->make<TH1D>("WrongEffic",       ";efficiency",  50,    0.0,    1.0);
    fHitPurity          = tfs->make<TH1D>("HitPurity",        ";purity",      51,    0.0,    1.02);
    fHitEfficPurity     = tfs->make<TH1D>("HitEfficPurity",   ";e*p",         51,    0.0,    1.02);
    
    fHitEfficiency0     = tfs->make<TH1D>("HitEffic0",        ";efficiency",  51,    0.0,    1.02);
    fWrongEfficiency0   = tfs->make<TH1D>("WrongEffic0",      ";efficiency",  50,    0.0,    1.0);
    fHitPurity0         = tfs->make<TH1D>("HitPurity0",       ";purity",      51,    0.0,    1.02);
    fHitEfficPurity0    = tfs->make<TH1D>("HitEfficPurity0",  ";e*p",         51,    0.0,    1.02);
    
    fTotHitsPerCluster  = tfs->make<TH1D>("HitsPerCluster",   ";# hits",      50,    0.0,     4.);
    fCorrectPerCluster  = tfs->make<TH1D>("CorrectPerClus",   ";# hits",      50,    0.0,     4.);
    fWrongPerCluster    = tfs->make<TH1D>("WrongPerClus",     ";# hits",      50,    0.0,     4.);
    fActualHits         = tfs->make<TH1D>("ActualPerClus",    ";# hits",      50,    0.0,     4.);
    
    fCosmicScoreAll     = tfs->make<TH1D>("CosmicScoreAll",   "score",       101,    0.0,   1.01);
    fCosmicScore        = tfs->make<TH1D>("CosmicScore",      "score",       101,    0.0,   1.01);
    fCosmicTag          = tfs->make<TH1D>("CosmicTag",        "tag value",   400,    0.0,   400.);
    fFlashScoreAll      = tfs->make<TH1D>("FlashScoreAll",    "score",       101,    0.0,   1.01);
    fFlashScore         = tfs->make<TH1D>("FlashScore",       "score",       101,    0.0,   1.01);
    fFlashTag           = tfs->make<TH1D>("FlashTag",         "tag value",   400,    0.0,   400.);
    fOrScore            = tfs->make<TH1D>("OrScore",          "score",       101,    0.0,   1.01);
    fOrTag              = tfs->make<TH1D>("OrTag",            "tag value",   400,    0.0,   400.);
    
    fCosmicScoreNuAll   = tfs->make<TH1D>("CosmicScoreNuAll", "score",       101,    0.0,   1.01);
    fCosmicScoreNu      = tfs->make<TH1D>("CosmicScoreNu",    "score",       101,    0.0,   1.01);
    fCosmicTagNu        = tfs->make<TH1D>("CosmicTagNu",      "tag value",   400,    0.0,   400.);
    fFlashScoreNuAll    = tfs->make<TH1D>("FlashScoreNuAll",  "score",       101,    0.0,   1.01);
    fFlashScoreNu       = tfs->make<TH1D>("FlashScoreNu",     "score",       101,    0.0,   1.01);
    fFlashTagNu         = tfs->make<TH1D>("FlashTagNu",       "tag value",   400,    0.0,   400.);
    fOrScoreNu          = tfs->make<TH1D>("OrScoreNu",        "score",       101,    0.0,   1.01);
    fOrTagNu            = tfs->make<TH1D>("OrTagNu",          "tag value",   400,    0.0,   400.);
    
    fHitTotalCharge[0]  = tfs->make<TH1D>("HitTotalChargeU",  "Tot Q U",     501,    0.0,  1002.);
    fHitPulseHeight[0]  = tfs->make<TH1D>("HitPulseHeightU",  "PH U",        501,    0.0,   501.);
    fHitPHClose[0]      = tfs->make<TH1D>("HitPHCloseU",      "PH U",        501,    0.0,   50.1);
    fHitMult[0]         = tfs->make<TH1D>("HitMultU",         "Mult U",       10,    0.0,   10.0);
    fHitTotalCharge[1]  = tfs->make<TH1D>("HitTotalChargeV",  "Tot Q V",     501,    0.0,  1002.);
    fHitPulseHeight[1]  = tfs->make<TH1D>("HitPulseHeightV",  "PH V",        501,    0.0,   501.);
    fHitPHClose[1]      = tfs->make<TH1D>("HitPHCloseV",      "PH V",        501,    0.0,   50.1);
    fHitMult[1]         = tfs->make<TH1D>("HitMultV",         "Mult V",       10,    0.0,   10.0);
    fHitTotalCharge[2]  = tfs->make<TH1D>("HitTotalChargeW",  "Tot Q W",     501,    0.0,  1002.);
    fHitPulseHeight[2]  = tfs->make<TH1D>("HitPulseHeightW",  "PH W",        501,    0.0,   501.);
    fHitPHClose[2]      = tfs->make<TH1D>("HitPHCloseW",      "PH W",        501,    0.0,   50.1);
    fHitMult[2]         = tfs->make<TH1D>("HitMultW",         "Mult W",       10,    0.0,   10.0);
    
    fHitTotalCharge[3]  = tfs->make<TH1D>("NseTotalChargeU",  "Tot Q U",     501,    0.0,  1002.);
    fHitPulseHeight[3]  = tfs->make<TH1D>("NsePulseHeightU",  "PH U",        501,    0.0,   501.);
    fHitPHClose[3]      = tfs->make<TH1D>("NsePHCloseU",      "PH U",        501,    0.0,   50.1);
    fHitMult[3]         = tfs->make<TH1D>("NseMultU",         "Mult U",       10,    0.0,   10.0);
    fHitTotalCharge[4]  = tfs->make<TH1D>("NseTotalChargeV",  "Tot Q V",     501,    0.0,  1002.);
    fHitPulseHeight[4]  = tfs->make<TH1D>("NsePulseHeightV",  "PH V",        501,    0.0,   501.);
    fHitPHClose[4]      = tfs->make<TH1D>("NsePHCloseV",      "PH V",        501,    0.0,   50.1);
    fHitMult[4]         = tfs->make<TH1D>("NseMultV",         "Mult V",       10,    0.0,   10.0);
    fHitTotalCharge[5]  = tfs->make<TH1D>("NseTotalChargeW",  "Tot Q W",     501,    0.0,  1002.);
    fHitPulseHeight[5]  = tfs->make<TH1D>("NsePulseHeightW",  "PH W",        501,    0.0,   501.);
    fHitPHClose[5]      = tfs->make<TH1D>("NsePHCloseW",      "PH W",        501,    0.0,   50.1);
    fHitMult[5]         = tfs->make<TH1D>("NseMultW",         "Mult W",       10,    0.0,   10.0);
    
    fDeltaPeakTime[0]   = tfs->make<TH1D>("DeltaPeakTimeU",   "Delta T",     120,  -15.0,   15.0);
    fDeltaPeakTime[1]   = tfs->make<TH1D>("DeltaPeakTimeV",   "Delta T",     120,  -15.0,   15.0);
    fDeltaPeakTime[2]   = tfs->make<TH1D>("DeltaPeakTimeW",   "Delta T",     120,  -15.0,   15.0);
    
    fPullPeakTime[0]    = tfs->make<TH1D>("PullPeakTimeU",    "Delta T",     100,  -10.0,   10.0);
    fPullPeakTime[1]    = tfs->make<TH1D>("PullPeakTimeV",    "Delta T",     100,  -10.0,   10.0);
    fPullPeakTime[2]    = tfs->make<TH1D>("PullPeakTimeW",    "Delta T",     100,  -10.0,   10.0);
}
   
//-----------------------------------------------------------------------
void LArPandoraAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void LArPandoraAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fSimulationProducerLabel    = p.get< std::string >("SimulationLabel");
    fMcHitCollectionModuleLabel = p.get< std::string >("MCHitLabel");
    fPFParticleProducerLabel    = p.get< std::string >("PFParticleLabel");
    fHitProducerLabel           = p.get< std::string >("HitLabel");
    fClusterProducerLabel       = p.get< std::string >("ClusterProducerLabel");
    fTrackProducerLabel         = p.get< std::string >("TrackProducerLabel");
    fCosmicProducerLabel        = p.get< std::string >("CosmicProducerLabel");
    fFlashProducerLabel         = p.get< std::string >("FlashProducerLabel");
    fSelectedPDG                = p.get< int         >("PDGcode");
    
    fMcHitCollectionModuleLabel = "mchitfinder";
    
    return;
}

//-----------------------------------------------------------------------
void LArPandoraAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    
    // The first step is to attempt to recover the collection of MCHits that
    // we will need for doing our PFParticle to MC matching
    art::Handle< std::vector<sim::MCHitCollection> > mcHitCollectionHandle;
    event.getByLabel(fMcHitCollectionModuleLabel, mcHitCollectionHandle);
    
    if (!mcHitCollectionHandle.isValid()) return;
    
    // Recover this into a local stl version
    const std::vector<sim::MCHitCollection> &mcHitCollectionVec = *mcHitCollectionHandle;

    // This is the standard method of reading multiple objects
    // associated with the same event; see
    // <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
    // for more information. Define a "handle" to point to a vector of
    // the objects.
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    // Then tell the event to fill the vector with all the objects of
    // that type produced by a particular producer.
    event.getByLabel(fSimulationProducerLabel, particleHandle);
    
    // Let's recover the MCTruth objects
    art::FindOneP<simb::MCTruth> mcTruthAssns(particleHandle, event, fSimulationProducerLabel);

    // Get all the simulated channels for the event. These channels
    // include the energy deposited for each track.
    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel(fSimulationProducerLabel, simChannelHandle);

    // The MCParticle objects are not necessarily in any particular
    // order. Since we may have to search the list of particles, let's
    // put them into a sorted map that will make searching fast and
    // easy. To save both space and time, the map will not contain a
    // copy of the MCParticle, but a pointer to it.
    std::map< int, const simb::MCParticle* > particleMap;
    
    int    nMcParticles = particleHandle->size();
    double logNumPart   = nMcParticles > 0 ? std::log10(nMcParticles) : -1.;
      
    fNMcParticles->Fill(logNumPart, 1.);
    
    // Before starting to loop through the particles, we are going to want to
    // build a mapping between hits and track id's
    // Start by recovering info from the event store
    art::Handle< std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);

    //we're gonna probably need the time service to convert hit times to TDCs
    const auto* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
    
    // Let us try to build up a mapping of tracks to "hits" on those tracks
    // Those "hits" will be the MCHit objects contained on a single channel per track
    // Let's see how this might work.
    // First, we'll have a map from channel to a set of MCHit points (to make them unique)
    typedef std::map<const unsigned int, std::set<const sim::MCHit*> > ChannelToMCHitMap;
    
    // Then we make a map to our tracks
    typedef std::map<int, ChannelToMCHitMap > TrackToChannelHitMap;
    
    TrackToChannelHitMap trackToChannelHitMap;
    
    // Loop through the MCHitCollection vector to fill out this map
    for(const auto& mcHitCollection : mcHitCollectionVec)
    {
        unsigned int channel = mcHitCollection.Channel();
        
        for(const auto& mcHit : mcHitCollection)
        {
            int trackID = mcHit.PartTrackId();
            
            trackToChannelHitMap[trackID][channel].insert(&mcHit);
        }
    }
    
    // We'll need the offsets for each plane
    std::map<geo::View_t, double> viewOffsetMap;
    
    viewOffsetMap[geo::kU] = 0.0; //fDetectorProperties->GetXTicksOffset(geo::kU, 0, 0)-fDetectorProperties->TriggerOffset();
    viewOffsetMap[geo::kV] = 0.0; //fDetectorProperties->GetXTicksOffset(geo::kV, 0, 0)-fDetectorProperties->TriggerOffset();
    viewOffsetMap[geo::kW] = 0.0; //fDetectorProperties->GetXTicksOffset(geo::kW, 0, 0)-fDetectorProperties->TriggerOffset();
    
    // our ultimate goal here
    HitToParticleMap hitToParticleMap;
    ParticleToHitMap particleToHitMap;
    
    int numHitsTotal(int(hitHandle->size()));
    int numHitsNoise[3] = {0,0,0};
    
    typedef std::map<int, std::vector<const sim::MCHit*> >  TrackIdToMcHitMap;
    typedef std::map<const recob::Hit*, TrackIdToMcHitMap > HitToTrackIdToMcHitMap;
    
    HitToTrackIdToMcHitMap hitToTrackIdToMcHitMap;
    
    // Ok, so this loop obviously takes the MC information and builds two maps
    // 1) a map from a Hit2D object to the track ID's that made it
    // 2) the reverse map, going from track ID to Hit2D object
    for (unsigned int iHit = 0, iHitEnd = hitHandle->size(); iHit < iHitEnd; ++iHit)
    {
        art::Ptr<recob::Hit> hit(hitHandle, iHit);
        
        const geo::WireID& wireId = hit->WireID();
        
        unsigned int channel = fGeometry->PlaneWireToChannel(wireId.Plane, wireId.Wire, wireId.TPC, wireId.Cryostat);
        
        const std::vector<sim::MCHit>& mcHitVec = mcHitCollectionVec.at(channel);
        
        int start_tdc = timeService->TPCTick2TDC( hit->StartTick() - viewOffsetMap[hit->View()]);
        int end_tdc   = timeService->TPCTick2TDC( hit->EndTick()   - viewOffsetMap[hit->View()]);
        
        sim::MCHit startTime;
        sim::MCHit endTime;
        
        startTime.SetTime(start_tdc, 0);
        endTime.SetTime(end_tdc, 0);
        
        std::vector<sim::MCHit>::const_iterator startItr = std::lower_bound(mcHitVec.begin(), mcHitVec.end(), startTime);
        std::vector<sim::MCHit>::const_iterator endItr   = std::upper_bound(startItr,         mcHitVec.end(), endTime);
        
        // Can it be that a match is not made?
        if (startItr != endItr)
        {
            while(startItr != endItr)
            {
                int trackID = (*startItr++).PartTrackId();
            
                hitToParticleMap[hit.get()].insert(trackID);
                particleToHitMap[trackID].insert(hit.get());
                
                hitToTrackIdToMcHitMap[hit.get()][trackID].push_back(&*startItr);
            }
        
            // Keep track of hit charge/pulse height
            double totalCharge = std::min(double(hit->Integral()),      1001.);
            double pulseHeight = std::min(double(hit->PeakAmplitude()),  500.);
            int    hitMult     = hit->Multiplicity();
            int    hitIndex    = hit->LocalIndex();
            
            fHitTotalCharge[hit->View()]->Fill(totalCharge, 1.);
            fHitPulseHeight[hit->View()]->Fill(pulseHeight, 1.);
            fHitMult[hit->View()]->Fill(hitMult, 1.);
            if (hitIndex == 0) fHitPHClose[hit->View()]->Fill(pulseHeight, 1.);
            
//            if (pulseHeight < 5.)
//            {
//                std::cout << "***>> Found good hit with small pulse height: " << pulseHeight << ", total Q: " << totalCharge << std::endl;
//                std::cout << *hit << std::endl;
//                std::cout << "----- MCHitVec size: " << mcHitVec.size() << ", start_tdc: " << start_tdc << ", end_tdc: " << end_tdc << std::endl;
                
//                for(const auto& mcHitTmp : mcHitVec)
//                {
//                    std::cout << "----> hit peak time: " << mcHitTmp.PeakTime() << ", Charge: " << mcHitTmp.Charge(false) << ", PH: " << mcHitTmp.Charge(true) << std::endl;
//                }
//            }
        }
        else
        {
            // Keep track of hit charge/pulse height for noise hits
            double totalCharge = std::min(double(hit->Integral()),      1001.);
            double pulseHeight = std::min(double(hit->PeakAmplitude()),  500.);
            int    hitMult     = hit->Multiplicity();
            int    hitIndex    = hit->LocalIndex();
            int    planeIdx    = hit->View() + 3;
            
            fHitTotalCharge[planeIdx]->Fill(totalCharge, 1.);
            fHitPulseHeight[planeIdx]->Fill(pulseHeight, 1.);
            fHitMult[planeIdx]->Fill(hitMult, 1.);
            if (hitIndex == 0) fHitPHClose[planeIdx]->Fill(pulseHeight, 1.);
            numHitsNoise[hit->View()]++;
            
            if (pulseHeight > 12.)
            {
                std::cout << "***>> Found noise hit with large pulse height: " << pulseHeight << ", total Q: " << totalCharge << std::endl;
                std::cout << *hit << std::endl;
                std::cout << "----- MCHitVec size: " << mcHitVec.size() << ", start_tdc: " << start_tdc << ", end_tdc: " << end_tdc << std::endl;
                
//                for(const auto& mcHitTmp : mcHitVec)
//                {
//                    std::cout << "----> hit peak time: " << mcHitTmp.PeakTime() << ", Charge: " << mcHitTmp.Charge(false) << ", PH: " << mcHitTmp.Charge(true) << std::endl;
//                }
            }
        }
    }
    
    for(const auto& hitItr : hitToTrackIdToMcHitMap)
    {
        if (hitItr.second.size() == 1)
        {
            const TrackIdToMcHitMap& trackIdToMcHitMap = hitItr.second;
            int                      trackId           = trackIdToMcHitMap.begin()->first;
            const recob::Hit*        hit               = hitItr.first;
            const geo::WireID&       wireId            = hit->WireID();
            int                      view              = wireId.Plane;
            
            unsigned int channel = fGeometry->PlaneWireToChannel(wireId.Plane, wireId.Wire, wireId.TPC, wireId.Cryostat);
            
//            std::cout << "--> Hit with track id: " << trackId << ", with mchit count: " << trackIdToMcHitMap.begin()->second.size() << std::endl;
//            std::cout << "    mchit full count: " << trackToChannelHitMap[trackId][channel].size() << std::endl;
//            std::cout << "    Hit Peak time: " << timeService->TPCTick2TDC(hitItr.first->PeakTime()) << ", PH: " << hitItr.first->PeakAmplitude() << ", TotCharge: " << hitItr.first->SummedADC() << std::endl;
//            std::cout << "    View: " << hitItr.first->View() << ", Wire: " << hitItr.first->WireID() << std::endl;
            
            double totCharge(0.);
            double maxPH(0.);
            double peakTime(0.);
            
//            for (const auto& mcHit : trackIdToMcHitMap.begin()->second)
            for (const auto& mcHit : trackToChannelHitMap[trackId][channel])
            {
//                std::cout << "    hit time: " << mcHit->PeakTime() << ", Charge: " << mcHit->Charge(false) << ", PH: " << mcHit->Charge(true) << std::endl;
                totCharge += mcHit->Charge(false);
                
                if (mcHit->Charge(true) > maxPH)
                {
                    maxPH    = mcHit->Charge(true);
                    peakTime = mcHit->PeakTime();
                }
            }

            double deltaPeakTime = timeService->TPCTick2TDC(hit->PeakTime() - viewOffsetMap[hit->View()]) - peakTime;
            double hitSigma      = 4.*hit->SigmaPeakTime();
            double pullPeakTime  = deltaPeakTime / hitSigma;
            
            fDeltaPeakTime[view]->Fill(deltaPeakTime, 1.);
            fPullPeakTime[view]->Fill(pullPeakTime, 1.);
            
//            std::cout << "    total charge: " << totCharge << ", max PH: " << maxPH << ", peak Time: " << peakTime << ", delta: " << deltaPeakTime << std::endl;
        }
    }
                     
    // Keep track of hits before proceeding
    fNumHitsTotal->Fill(std::log10(double(numHitsTotal)), 1.);
    fNumNoisePlane[0]->Fill(std::log10(double(numHitsNoise[0])), 1.);
    fNumNoisePlane[1]->Fill(std::log10(double(numHitsNoise[1])), 1.);
    fNumNoisePlane[2]->Fill(std::log10(double(numHitsNoise[2])), 1.);
    fNumHitsNoise->Fill(std::log10(double(numHitsNoise[0]+numHitsNoise[1]+numHitsNoise[2])), 1.);
    
    std::cout << "Total hits: " << numHitsTotal << ", noise U: " << numHitsNoise[0] << ", V: " << numHitsNoise[1] << ", W: " << numHitsNoise[2] << std::endl;
    
    // Some counters for the loop below
    int nMuonsWithHits(0);

    // There are two purposes to looping through the MCParticles now
    // 1) we want to build a map of track ID to MCParticle
    // 2) while we do that we can count the number of CR muons that make hits (hence are tracks)
    for ( auto const& particle : (*particleHandle) )
    {
        // For the methods you can call to get particle information,
        // see ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
        int trackID = particle.TrackId();
        
        // Temporary
//        std::cout << particle << std::endl;

        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[trackID] = &particle;

        // Histogram the PDG code of every particle in the event.
        int pdgCode = particle.PdgCode();
        fPDGCodeHist->Fill( pdgCode );
        
        // Keep track of energy
        double logEne = std::log10(particle.E());
        fEneMcPart->Fill(logEne, 1.);
        
        std::string partProcess = particle.Process();

        // Just keeping tabs on decay processes here
        if (partProcess != "primary")
        {
            if (pdgCode == 11)
            {
                fMcProcess->Fill(partProcess.c_str(), 1.);
                //std::cout << "Found electron, process = " << partProcess << std::endl;
            }
        }
        
        // No point on proceeding from here if no hits produced
        if (particleToHitMap.find(trackID) == particleToHitMap.end()) continue;
        
        // Keep track of the number of Hit2D objects created by this track
        int nParticleHit2Ds(particleToHitMap[trackID].size());
        
        // Tag the cosmic ray muons here
        if ((pdgCode == 13 || pdgCode == -13) && particle.Process() == "primary")
        {
            nMuonsWithHits++;
            fMcMuonHitsPerTrack->Fill(log10(double(nParticleHit2Ds)), 1.);
        }

        // For this example, we want to fill the n-tuples and histograms
        // only with information from the primary particles in the
        // event, whose PDG codes match a value supplied in the .fcl file.
        if ( particle.Process() == "primary"  &&  pdgCode == fSelectedPDG )
        {
            // A particle has a trajectory, consisting of a set of
            // 4-positions and 4-mommenta.
            size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

            // For trajectories, as for vectors and arrays, the
            // first point is #0, not #1.
            int last = numberTrajectoryPoints - 1;
            const TLorentzVector& positionStart = particle.Position(0);
            const TLorentzVector& positionEnd   = particle.Position(last);
            const TLorentzVector& momentumStart = particle.Momentum(0);
            const TLorentzVector& momentumEnd   = particle.Momentum(last);
            
            fEneSelected->Fill(logEne, 1.);
            if (!particleToHitMap[trackID].empty()) fEneWithHits->Fill(logEne, 1.);

            // Make a histogram of the starting momentum.
            fMomentumHist->Fill( momentumStart.P() );

            // Fill arrays with the 4-values. (Don't be fooled by
            // the name of the method; it just puts the numbers from
            // the 4-vector into the array.)
            positionStart.GetXYZT( fStartXYZT );
            positionEnd.GetXYZT( fEndXYZT );
            momentumStart.GetXYZT( fStartPE );
            momentumEnd.GetXYZT( fEndPE );
        } // if primary and PDG selected by user
    } // loop over all particles in the event.
    
    fMcMuonsWithHits->Fill(nMuonsWithHits, 1.);
    
    // Ok, at this point we're ready to recover PFParticles and the assorted necessities for the next step
    // For sure we need a boatload of stuff here...
    // Start with some useful typdefs
    typedef std::set<art::Ptr<recob::Hit> >                        Hit2DSet;              // Need to count only unique hits
    typedef std::map< int, Hit2DSet >                              TrackIDToHit2DMap;
    typedef std::vector<const recob::PFParticle*>                  PFParticleVec;
    typedef std::map< int, PFParticleVec >                         TrackIDToPFParticleVecMap;
    typedef std::map<const recob::PFParticle*, TrackIDToHit2DMap > PFParticleToTrackHit2DMap;
    
    // Now define the maps relating pfparticles to tracks
    PFParticleToTrackHit2DMap pfParticleToTrackHitMap;
    TrackIDToPFParticleVecMap trackIDToPFParticleMap;
    TrackIDToPFParticleVecMap muonIDToClusterMap;
    
    // Something to keep track of number of hits associated to a cluster
    std::map<const recob::PFParticle*, int> pfParticleHitCntMap;
    
    // Recover the PFParticles, the main products for our next major loop
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    event.getByLabel(fPFParticleProducerLabel, pfParticleHandle);
    
    // Without a valid collection of PFParticles there is nothing to do here
    if (!pfParticleHandle.isValid()) return;
    
    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    event.getByLabel(fPFParticleProducerLabel, clusterHandle);
    
    // If there are no clusters then something is really wrong
    if (!clusterHandle.isValid()) return;
    
    // Now retrieve a handle to the fit tracks
    art::Handle<std::vector<recob::Track> > trackHandle;
    event.getByLabel(fTrackProducerLabel, trackHandle);
    
    if (!trackHandle.isValid()) return;
    
    // Recover the collection of associations between PFParticles and clusters, this will
    // be the mechanism by which we actually deal with clusters
    art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, event, fPFParticleProducerLabel);
    
    // Recover the collection of associations between PFParticles and tracks, this will
    // be the mechanism by which we actually deal with tracks
    art::FindManyP<recob::Track> trackAssns(pfParticleHandle, event, fTrackProducerLabel); //fPFParticleProducerLabel);
    
    // Recover two sets of associations for linking cosmic tags to the tracks
    art::FindManyP<anab::CosmicTag> cosmicAssns(trackHandle, event, fCosmicProducerLabel);
    art::FindManyP<anab::CosmicTag> flashAssns(trackHandle,  event, fFlashProducerLabel); 
    
    // Likewise, recover the collection of associations to hits
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, event, fPFParticleProducerLabel);
    
    // Keep track of number of pfparticles
    fNumPFParticles->Fill(pfParticleHandle->size(), 1.);
    
    // The goal of this loop is to recover the list of 2D hits associated with a given PFParticle
    // (through associations with 2D clusters) and develop some maps which allow us to associate
    // the PFParticle to a MC particle.
    for (const auto& pfParticle : (*pfParticleHandle))
    {
        // First a storage container for all the hits
        art::PtrVector<recob::Hit> hits;
        
        // Now we get the 2D clusters associated to this PFParticle
        std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(pfParticle.Self());
        
        // Nominally, one believes there will be 3 2D clusters so we need to loop through the list
        for(const auto& cluster : clusterVec)
        {
            // Recover the 2D hits associated to a given cluster
            std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster->ID());
            
            // Quick mismatch test
            geo::View_t viewToMatch = cluster->View();
            
            for(const auto& hit : hitVec)
            {
                geo::View_t view = hit->View();
                
                if (viewToMatch == geo::kUnknown) viewToMatch = view;
                
                if (view != viewToMatch)
                {
                    std::cout << "*****>>> Found a view mismatch: " << viewToMatch << " wanted, found " << view << std::endl;
                }
            }
            
            // And insert them into our list
            hits.insert(hits.end(), hitVec.begin(), hitVec.end());
        }
        
        fNumClusPerPFPart->Fill(clusterVec.size(), 1.);
        
        // Keep track of this to hopefully save some time later
        pfParticleHitCntMap[&pfParticle] = hits.size();

        // To categorize the PFParticle (associate to MCParticle), we will
        // create and fill an instance of a TrackIDToHitMap.
        // So, for a given PFParticle we will have a map of the track ID's (MCParticles)
        // and the hits that they contributed energy too
        pfParticleToTrackHitMap[&pfParticle] = TrackIDToHit2DMap();
        TrackIDToHit2DMap& trackIdToHit2DMap = pfParticleToTrackHitMap[&pfParticle];
        
        // Something to count MCParticles contributing here
        std::set<int> trackIdCntSet;
        int           nMultiParticleHits(0);
        int           nNoParticleHits(0);
        
        // To fill out this map we now loop over the recovered Hit2D's and stuff into the map
        for (const auto& hit : hits)
        {
            // Given the Hit2D, recover the list of asscociated track ID's (MCParticles)
            HitToParticleMap::iterator hitToParticleMapItr = hitToParticleMap.find(hit.get());
            
            if (hitToParticleMapItr != hitToParticleMap.end())
            {
                // More than one MCParticle can contribute energy to make a given hit
                // So loop over the track ID's for this hit
                for(const auto& trackId : hitToParticleMapItr->second)
                {
                    trackIdToHit2DMap[trackId].insert(hit);
                    trackIdCntSet.insert(trackId);
                }
                
                if (hitToParticleMapItr->second.size() > 1) nMultiParticleHits++;
            }
            else
            {
                //std::cout << "Hit does not have entry in HitToParticleMap, peak time: " << hit2D->getHit()->PeakTime() << std::endl;
                nNoParticleHits++;
            }
        }
        
        // Make sure something happened here...
        if (!trackIdToHit2DMap.empty())
        {
            int    bestTrackId(-99999);
            size_t bestTrackCnt(0);
        
            // Make a quick loop through the map to do some majority logic matching
            // We only care about find the best match so no worries about sorting
            // Take advantage of this loop to count the number of unique 2D hits
            for(const auto& trackItr : trackIdToHit2DMap)
            {
                // Majority logic matching
                if (trackItr.second.size() > bestTrackCnt)
                {
                    bestTrackId  = trackItr.first;
                    bestTrackCnt = trackItr.second.size();
                }
            }
        
            if (bestTrackId != -99999)
            {
                trackIDToPFParticleMap[bestTrackId].push_back(&pfParticle);
            
                std::map< int, const simb::MCParticle* >::const_iterator partMapItr = particleMap.find(bestTrackId);
            
                if (partMapItr != particleMap.end())
                {
                    const simb::MCParticle* mcPart = partMapItr->second;
                    
                    if (mcPart && (mcPart->PdgCode() == 13 || mcPart->PdgCode() == -13) && mcPart->Process() == "primary")
                        muonIDToClusterMap[bestTrackId].push_back(&pfParticle);
                }
                
                // Do a little histogramming
                fNumMcPerPFPart->Fill(trackIdCntSet.size(), 1.);
                fNumHitsPerPFPart->Fill(std::log10(double(hits.size())), 1.);
                fNumHitsBestMC->Fill(std::log10(double(particleToHitMap[bestTrackId].size())), 1.);
                if (nMultiParticleHits > 0.) fNumHitsGtOneMC->Fill(std::log10(double(nMultiParticleHits)), 1.);
                if (nNoParticleHits    > 0.) fNumHitsNoise->Fill(std::log10(double(nNoParticleHits)), 1.);
            }
        }
    } // end of loop over the PFParticle collection
    
    // More counting...
    fNumTracksRecon->Fill(trackIDToPFParticleMap.size(), 1.);
    fNumMuonsRecon->Fill(muonIDToClusterMap.size(), 1.);
    
    // Always a handy thing to have hidden in your code:
    const double radToDegrees = 180. / 3.14159265;
    
    // Ok, at this point we have all the tools we need to proceed with trying to get some track/hit efficiencies
    // Our strategy will be to loop over the MCParticle collection again, considering only those particles we're interested
    // in (CR muons for now) and, of those, only those which leave hits in the TPC (since the CR Muon overlay will have
    // particles which make no hits in the TPC).
    for (const auto& particle : (*particleHandle))
    {
        // For the methods you can call to get particle information,
        // see ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
        int bestTrackID = particle.TrackId();
        
        // Recover the particle's identity
        int trackPDGCode = particle.PdgCode();
        
        // Only want to keep primary cosmic rays for now
        if (particle.Process() != "primary"  || !(trackPDGCode == 13 || trackPDGCode == -13)) continue;
        
        // Did this mc particle leave hits in the TPC?
        ParticleToHitMap::iterator particleToHitItr = particleToHitMap.find(bestTrackID);
        
        // No hits no work
        if (particleToHitItr == particleToHitMap.end())
        {
            if (0) std::cout << "** bestTrackID = " << bestTrackID << std::endl;
            continue;
        }
        
        // Let's get the total number of "true" hits that are created by this MCParticle
        int nTrueMcHits = particleToHitItr->second.size();
        
        // We need to make sure there was some activity in the TPC
        if (nTrueMcHits < 1) continue;
        
        // Set baseline for MC muons
        TVector3 mcPartDir(particle.Momentum().Vect());
        
        mcPartDir.SetMag(1.);
        
        double mcLogNumHits = log10(double(nTrueMcHits));
        double mcPhi        = radToDegrees * mcPartDir.Phi();
        double mcTheta      = radToDegrees * mcPartDir.Theta();
        double mcMomInXZ    = sqrt(mcPartDir.X()*mcPartDir.X() + mcPartDir.Z()*mcPartDir.Z());
        double cosMcThetaXZ = mcMomInXZ > 0. ? mcPartDir.X() / mcMomInXZ : 0.;
        double mcThetaXZ    = radToDegrees * acos(cosMcThetaXZ);
        
        if (particle.Momentum().Z() < 0.) mcThetaXZ = -mcThetaXZ;

        // Fill hists for all MCParticles that leave hits in TPC
        fMCTrackNumHits->Fill(mcLogNumHits, 1.);
        fMCTrackPhi->Fill(mcPhi, 1.);
        fMCTrackTheta->Fill(mcTheta, 1.);
        fMCTrackThetaXZ->Fill(mcThetaXZ, 1.);
        
        // Now check to see that a PFParticle is associated with this MCParticle (trackID)
        TrackIDToPFParticleVecMap::iterator trackToPFParticleItr = trackIDToPFParticleMap.find(bestTrackID);
        
        // If no PFParticle then skip the rest
        if (trackToPFParticleItr == trackIDToPFParticleMap.end())
        {
            // Effiencies are zero when hits in TPC but not PFParticle!
            fHitEfficVsHits->Fill(mcLogNumHits, 0.);
            fHitEfficVsTheta->Fill(mcTheta, 0.);
            fHitEfficVsPhi->Fill(mcPhi, 0.);
            fClusEfficVsHits->Fill(mcLogNumHits, 0.);
            fClusEfficVsTheta->Fill(mcTheta, 0.);
            fClusEfficVsPhi->Fill(mcPhi, 0.);
            fTrackEfficVsHits->Fill(mcLogNumHits, 0.);
            fTrackEfficVsTheta->Fill(mcTheta, 0.);
            fTrackEfficVsPhi->Fill(mcPhi, 0.);
            fHitEfficiency->Fill(0., 1.);
            fHitPurity->Fill(0., 1.);
            fHitEfficPurity->Fill(0., 1.);
            
            std::cout << "***>> Ana no PFParticle for CR mu with " << nTrueMcHits << " hits, theta/phi= " << mcTheta << "/" << mcPhi << std::endl;
            
            continue;
        }
        
        // Here our goal is to count the number of hits associated to a PFParticle by the "best" MCParticle
        // (which is the one which produced the most number of hits the PFParticle has grouped together)
        const recob::PFParticle* pfParticle(0);
        int                      bestCnt(0);
        int                      nPFParticles(0);

        // Loop over the PFParticles associated to this track ID (MCParticle
        for(const auto& tmpPart : trackToPFParticleItr->second)
        {
            // This is to get the list of hits this PFParticle has broken out by track id
            TrackIDToHit2DMap::iterator tmpPartTrackItr = pfParticleToTrackHitMap[tmpPart].find(bestTrackID);
            
            if (tmpPartTrackItr != pfParticleToTrackHitMap[tmpPart].end())
            {
                int trackHitCnt = tmpPartTrackItr->second.size();
                
                if (trackHitCnt > bestCnt)
                {
                    bestCnt   = trackHitCnt;
                    pfParticle = tmpPart;
                }
                
                if (trackHitCnt > int(nTrueMcHits / 10)) nPFParticles++;
            }
        }
        
        fNumPartPerTrk->Fill(nPFParticles, 1.);
        
        // Now we can determine:
        // 1) The total number of hits in this cluster
        // 2) The number of hits belonging to the matched mc track
        // 3) The number of hits on that track
        // 4) The number of hits not belonging to this track which are on the cluster (the impurities)
        size_t nTotalClusterHits = pfParticleHitCntMap[pfParticle];
        size_t nMatchedHits      = bestCnt;
        size_t nWrongClusterHits = nTotalClusterHits - nMatchedHits;
        
        // MicroBooNE definitions of efficiency and purity:
        // efficiency E = number of true hits in the cluster / number of true hits from MC particle
        // purity P     = number of true hits in the cluster / all hits in the cluster
        double hitEfficiency = double(nMatchedHits)      / double(nTrueMcHits);
        double wrongEff      = double(nWrongClusterHits) / double(nTrueMcHits);
        double hitPurity     = double(nMatchedHits)      / double(nTotalClusterHits);
        
        fHitEfficiency->Fill(hitEfficiency, 1.);
        fWrongEfficiency->Fill(wrongEff, 1.);
        fHitPurity->Fill(hitPurity, 1.);
        fHitEfficPurity->Fill(hitEfficiency*hitPurity, 1.);
        
        fHitEfficiency0->Fill(hitEfficiency, 1.);
        fWrongEfficiency0->Fill(wrongEff, 1.);
        fHitPurity0->Fill(hitPurity, 1.);
        fHitEfficPurity0->Fill(hitEfficiency*hitPurity, 1.);
        
        fTotHitsPerCluster->Fill(log10(double(nTotalClusterHits)), 1.);
        fCorrectPerCluster->Fill(log10(double(nMatchedHits)), 1.);
        fWrongPerCluster->Fill(log10(double(nWrongClusterHits)), 1.);
        fActualHits->Fill(log10(double(nTrueMcHits)), 1.);

        // We can also now associate to tracks and see how many we have fit for this particle
        
        // Now we get the 2D clusters associated to this PFParticle
        if (trackAssns.isValid() && pfParticle->IsPrimary())
        {
            std::vector<art::Ptr<recob::Track> > trackVec = trackAssns.at(pfParticle->Self());
        
            fNumTracksPerClus->Fill(trackVec.size(), 1.);
        
            if (trackVec.size() > 0)
            {
                art::Ptr<recob::Track>& track(trackVec.front());
                
                fNumHitsPerTrack->Fill(std::log10(track->NumberTrajectoryPoints()), 1.);
                fNumMCHitsPerTrack->Fill(mcLogNumHits, 1.);
                fNumClusHitsPerTrk->Fill(std::log10(nTotalClusterHits), 1.);
                
                if (hitEfficiency > 0.6)
                {
                    fTrackEfficVsHits->Fill(mcLogNumHits, 1.);
                    fTrackEfficVsTheta->Fill(mcTheta, 1.);
                    fTrackEfficVsPhi->Fill(mcPhi, 1.);
                }
                
                // Try recovering associations to cosmic tags
                if (cosmicAssns.isValid() && cosmicAssns.size() > 0)
                {
                    std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track->ID());
                    
                    if (!cosmicVec.empty())
                    {
                        art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                        
                        float               cosmicTagScore(cosmicTag->CosmicScore());
                        anab::CosmicTagID_t cosmicTagType(cosmicTag->CosmicType());
                        
                        fCosmicScoreAll->Fill(cosmicTagScore, 1.);
                        
                        if (cosmicTagType != 100)
                        {
                            fCosmicScore->Fill(cosmicTagScore, 1.);
                            fCosmicTag->Fill(cosmicTagType, 1.);
                        
                            // This is for the -or- case
                            if (cosmicTagScore == 0.5)
                            {
                                if (flashAssns.isValid() && flashAssns.size() > 0)
                                {
                                    std::vector<art::Ptr<anab::CosmicTag> > flashVec = flashAssns.at(track->ID());
                                
                                    if (!flashVec.empty())
                                    {
                                        art::Ptr<anab::CosmicTag>& flashTag(flashVec.front());
                                    
                                        float               flashTagScore(flashTag->CosmicScore());
                                        anab::CosmicTagID_t flashTagType(flashTag->CosmicType());
                                    
                                        fOrScore->Fill(flashTagScore, 1.);
                                        fOrTag->Fill(flashTagType, 1.);
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Try recovering associations to flash tags
                if (flashAssns.isValid() && flashAssns.size() > 0)
                {
                    std::vector<art::Ptr<anab::CosmicTag> > flashVec = flashAssns.at(track->ID());
                
                    if (!flashVec.empty())
                    {
                        art::Ptr<anab::CosmicTag>& flashTag(flashVec.front());
                    
                        float               flashTagScore(flashTag->CosmicScore());
                        anab::CosmicTagID_t flashTagType(flashTag->CosmicType());
                    
                        fFlashScoreAll->Fill(flashTagScore, 1.);
                        
                        if (flashTagType != 100)
                        {
                            fFlashScore->Fill(flashTagScore, 1.);
                            fFlashTag->Fill(flashTagType, 1.);
                        }
                    }
                }
            }
        }
        
        // What is the right way to handle this? It seems not right to say you "found" a track
        // If you have fewer than 50% or so of the hits...
        if (hitEfficiency > 0.6)
        {
            fClusEfficVsHits->Fill(mcLogNumHits, 1.);
            fClusEfficVsTheta->Fill(mcTheta, 1.);
            fClusEfficVsPhi->Fill(mcPhi, 1.);
        }
        else
        {
            fClusEfficVsHits->Fill(mcLogNumHits, 0.);
            fClusEfficVsTheta->Fill(mcTheta, 0.);
            fClusEfficVsPhi->Fill(mcPhi, 0.);
            
            std::cout << "***>> Ana hitEffic < 0.6 for CR mu with " << nTrueMcHits << " hits, theta/phi= " << mcTheta << "/" << mcPhi << ", hitEffic = " << hitEfficiency << std::endl;
        }
        
        fHitEfficVsHits->Fill(mcLogNumHits, hitEfficiency);
        fHitEfficVsTheta->Fill(mcTheta, hitEfficiency);
        fHitEfficVsPhi->Fill(mcPhi, hitEfficiency);
    }
    
    // One last task worth pursuing is to see if we can pick out the neutrino interaction and do some categorization of it.
    bool foundNeutrino(false);
    bool foundPFParticle(false);
    
    const simb::MCTruth* theAnswerIsOutThere(0);
    std::vector<const simb::MCParticle*> mcPartVec;
    
    for (size_t particleIdx = 0; particleIdx < particleHandle->size(); particleIdx++)
    {
        art::Ptr<simb::MCParticle> particle(particleHandle, particleIdx);
        
        // Focus on primaries for the moment
//        if (particle->Process() != "primary") continue;
        if (particle->NumberTrajectoryPoints() < 3) continue;
        
        bool isNeutrino(false);
        
        try
        {
            art::Ptr<simb::MCTruth> mcTruth = mcTruthAssns.at(particleIdx);
        
            if (mcTruth->Origin() != simb::kBeamNeutrino) continue;
            
            isNeutrino = true;
            
            theAnswerIsOutThere = mcTruth.get();
        }
        catch(...)
        {
            isNeutrino = false;
        }
        
        if (isNeutrino)
        {
            foundNeutrino = true;
            
            // For the methods you can call to get particle information,
            // see ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
            int bestTrackID = particle->TrackId();
            
            // Recover the particle's identity
//          int trackPDGCode = particle.PdgCode();
            
            // Did this mc particle leave hits in the TPC?
            ParticleToHitMap::iterator particleToHitItr = particleToHitMap.find(bestTrackID);
            
            // No hits no work
            if (particleToHitItr == particleToHitMap.end()) continue;
            
            // Let's get the total number of "true" hits that are created by this MCParticle
            int nTrueMcHits = particleToHitItr->second.size();
            
            // We need to make sure there was some activity in the TPC
            if (nTrueMcHits < 1) continue;
            
            mcPartVec.push_back(particle.get());

            // Now check to see that a PFParticle is associated with this MCParticle (trackID)
            TrackIDToPFParticleVecMap::iterator trackToPFParticleItr = trackIDToPFParticleMap.find(bestTrackID);
            
            // If no PFParticle then skip the rest
            if (trackToPFParticleItr == trackIDToPFParticleMap.end())
            {
                continue;
            }
            
            // Here our goal is to count the number of hits associated to a PFParticle by the "best" MCParticle
            // (which is the one which produced the most number of hits the PFParticle has grouped together)
            const recob::PFParticle* pfParticle(0);
            int                      bestCnt(0);
            int                      nPFParticles(0);
            
            // Loop over the PFParticles associated to this track ID (MCParticle
            for(const auto& tmpPart : trackToPFParticleItr->second)
            {
                // This is to get the list of hits this PFParticle has broken out by track id
                TrackIDToHit2DMap::iterator tmpPartTrackItr = pfParticleToTrackHitMap[tmpPart].find(bestTrackID);
            
                if (tmpPartTrackItr != pfParticleToTrackHitMap[tmpPart].end())
                {
                    int trackHitCnt = tmpPartTrackItr->second.size();
                
                    if (trackHitCnt > bestCnt)
                    {
                        bestCnt   = trackHitCnt;
                        pfParticle = tmpPart;
                    }
                
                    if (trackHitCnt > int(nTrueMcHits / 10)) nPFParticles++;
                }
            }
            
            // I don't think this can happen...
            if (!pfParticle) continue;
            
            // Only consider primaries for now
            if (!pfParticle->IsPrimary())
            {
                std::cout << "PFParticle is not primary, trackHitCnt = " << bestCnt << std::endl;
                continue;
            }
            
            foundPFParticle = true;
            
            // Our immediate goal here is to find out if we tagged this PFParticle as a cosmic ray
            if (trackAssns.isValid())
            {
                std::vector<art::Ptr<recob::Track> > trackVec = trackAssns.at(pfParticle->Self());
            
                if (trackVec.size() > 0)
                {
                    // For our purposes here we only need the first track since all will be associated
                    // to the same CosmicTag
                    art::Ptr<recob::Track>& track(trackVec.front());
                
                    // Try recovering associations to cosmic tags
                    if (cosmicAssns.isValid() && cosmicAssns.size() > 0)
                    {
                        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track->ID());
                    
                        if (!cosmicVec.empty())
                        {
                            art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                        
                            float               cosmicTagScore(cosmicTag->CosmicScore());
                            anab::CosmicTagID_t cosmicTagType(cosmicTag->CosmicType());
                        
                            fCosmicScoreNuAll->Fill(cosmicTagScore, 1.);
                            
                            if (cosmicTagType != 100)
                            {
                                fCosmicScoreNu->Fill(cosmicTagScore, 1.);
                                fCosmicTagNu->Fill(cosmicTagType, 1.);
                            
                                // This is for the -or- case
                                if (cosmicTagScore == 0.5)
                                {
                                    if (flashAssns.isValid() && flashAssns.size() > 0)
                                    {
                                        std::vector<art::Ptr<anab::CosmicTag> > flashVec = flashAssns.at(track->ID());
                                    
                                        if (!flashVec.empty())
                                        {
                                            art::Ptr<anab::CosmicTag>& flashTag(flashVec.front());
                                        
                                            float               flashTagScore(flashTag->CosmicScore());
                                            anab::CosmicTagID_t flashTagType(flashTag->CosmicType());
                                        
                                            fOrScoreNu->Fill(flashTagScore, 1.);
                                            fOrTagNu->Fill(flashTagType, 1.);
                                        }
                                    }
                                }
                            }
                            else if (cosmicTagScore == 1.)
                            {
                                std::cout << "***===> Geometry tagger mistag neutrino: tag id " << cosmicTagType << ", run " << fRun << ", event " << fEvent << std::endl;
                            }
                        }
                    }
                
                    // Try recovering associations to flash tags
                    if (flashAssns.isValid() && flashAssns.size() > 0)
                    {
                        std::vector<art::Ptr<anab::CosmicTag> > flashVec = flashAssns.at(track->ID());
                    
                        if (!flashVec.empty())
                        {
                            art::Ptr<anab::CosmicTag>& flashTag(flashVec.front());
                        
                            float               flashTagScore(flashTag->CosmicScore());
                            anab::CosmicTagID_t flashTagType(flashTag->CosmicType());
                            
                            fFlashScoreNuAll->Fill(flashTagScore, 1.);
                            
                            if (flashTagType != 100)
                            {
                                fFlashScoreNu->Fill(flashTagScore, 1.);
                                fFlashTagNu->Fill(flashTagType, 1.);
                            
                                if (flashTagScore == 1.)
                                {
                                    std::cout << "***===> Flash tagger mistag neutrino: tag id " << flashTagType << ", run " << fRun << ", event " << fEvent << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (!foundPFParticle)
    {
        std::cout << "@@@===> Did not find PFParticle, foundNeutrino = " << foundNeutrino << " for event " << fEvent << std::endl;
        
        if (theAnswerIsOutThere) std::cout << "@@@===> MCTruth: " << *theAnswerIsOutThere << std::endl;
        
        for(const auto& particle : mcPartVec)
        {
            std::cout << "@@@===> " << *particle << std::endl;
        }
    }
    

    return;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see LArPandoraAna.fcl for more information.
DEFINE_ART_MODULE(LArPandoraAna)

} // namespace LArPandoraAna

#endif // LArPandoraAna_Module
