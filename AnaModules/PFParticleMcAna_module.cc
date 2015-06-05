/**
 *  @file   PFParticleMcAna_module.cc
 *
 *  @brief  Ana analysis module aimed at comparing PFParticles to MC Truth
 *          with a specific emphasis on comparisons to single particles
 *
 *          Original version created by modifying AnalysisExample
 */

#ifndef PFParticleMcAna_module
#define PFParticleMcAna_module

// LArSoft includes
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/PCAxis.h"
#include "AnalysisBase/CosmicTag.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "MCBase/MCHitCollection.h"

//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "Utilities/TimeService.h"
#include "Utilities/AssociationUtil.h"

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

namespace PFParticleMcAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class PFParticleMcAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit PFParticleMcAna(fhicl::ParameterSet const& pset);
    virtual ~PFParticleMcAna();

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
    
    // Calculate lengths
    double length(const recob::Track* track);
    double length(const simb::MCParticle& part, double dx,
                  TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                  unsigned int tpc = 0, unsigned int cstat = 0);
    
    // This method is meant to be called at the start of the event
    void PrepareEvent(const art::Event& evt);

    // Useful typedefs to keep from getting lost in stl'ese
    typedef std::map< const recob::Hit*, std::set<int> > HitToParticleMap;    // Maps recob::Hits to track ID's (hence MCParticles)
    typedef std::map< int, std::set<const recob::Hit*> > ParticleToHitMap;    // Maps track IDs to recob Hits
    typedef std::vector< art::Ptr<recob::Hit> >          HitVector;
    
    // Define a function to fill the above data structures
    void MakeHitParticleMaps(const std::vector<sim::MCHitCollection>& mcHitCollectionVec,
                             const std::vector<recob::Hit>&           recoHitVec,
                             HitToParticleMap&                        hitToParticleMap,
                             ParticleToHitMap&                        particleToHitMap);
    
    // More useful typedefs, this time for relating PFParticles to Track IDs and hits
    typedef std::set<art::Ptr<recob::Hit> >                        Hit2DSet;                   // Need to count only unique hits
    typedef std::map< int, Hit2DSet >                              TrackIDToHit2DMap;          // Maps Track ID's to lists of hits
    typedef std::vector<const recob::PFParticle*>                  PFParticleVec;              // Typedef a vector of PFParticles
    typedef std::map< int, PFParticleVec >                         TrackIDToPFParticleVecMap;  // Maps Track ID's to the PFParticles associated to them
    typedef std::map<const recob::PFParticle*, TrackIDToHit2DMap > PFParticleToTrackHit2DMap;  // Maps PFParticles to a map of associated Track ID's to their hits
    typedef std::map<const recob::PFParticle*, int>                PFParticleHitCntMap;        // Allows us to count reco hits per PFParticle directly

    // Forward declaration of object to handle sorting of PFParticles associated to Track IDs
    class SortPFParticleVec;
    
    // Define a function to fill all of the above
    void MakePFParticleMaps(const art::Event&                                   event,
                            const art::Handle<std::vector<recob::PFParticle> >& pfParticleHandle,
                            const HitToParticleMap&                             hitToParticleMap,
                            PFParticleToTrackHit2DMap&                          pfParticleToTrackHit2DMap,
                            TrackIDToPFParticleVecMap&                          trackIDtoPFParticleVecMap,
                            PFParticleHitCntMap&                                pfParticleHitCntMap);

    // The parameters we'll read from the .fcl file.
    std::string fSimulationProducerLabel;    //> The name of the producer that tracked simulated particles through the detector
    std::string fMcHitCollectionModuleLabel; //> The name of the producer of MCHits
    std::string fPFParticleProducerLabel;    //> The name of the produder of the PFParticle hierarchy
    std::string fHitProducerLabel;           //> The name of the producer that created hits
    std::string fClusterProducerLabel;       //> The name of the producer that created clusters
    std::string fTrackProducerLabel;         //> The name of the producer that created the tracks
    std::string fCosmicProducerLabel;        //> The name of the producer that created cosmic tags
    std::string fFlashProducerLabel;         //> The name of the producer that created flash tags
    int fSelectedPDG;                        //> PDG code of particle we'll focus on

    // Pointers to the histograms we'll create. 

    // Maximum number of entries per row in our tuple
    int      fMaxEntries;
    
    // The variables that will go into the n-tuple.
    // Start with basic run information
    Int_t    fEvent;
    Int_t    fRun;
    Int_t    fSubRun;

    // Some general information on the number of hits
    Int_t    fNumHits;
    Int_t    fNoiseHits;
    Int_t    fNegTrackIds;
    
    // Monte Carlo information
    Int_t    fNumMcParticles;
    std::vector<Int_t>   fPDGCode;
    std::vector<Float_t> fMcPartStartX;
    std::vector<Float_t> fMcPartStartY;
    std::vector<Float_t> fMcPartStartZ;
    std::vector<Float_t> fMcPartStartDirX;
    std::vector<Float_t> fMcPartStartDirY;
    std::vector<Float_t> fMcPartStartDirZ;
    std::vector<Float_t> fMcPartEndX;
    std::vector<Float_t> fMcPartEndY;
    std::vector<Float_t> fMcPartEndZ;
    std::vector<Float_t> fMcPartEndDirX;
    std::vector<Float_t> fMcPartEndDirY;
    std::vector<Float_t> fMcPartEndDirZ;
    std::vector<Float_t> fMcPartEne;
    std::vector<Float_t> fMcPartMom;
    std::vector<Float_t> fMcPartMass;
    std::vector<Float_t> fMcPartTrackLen;
    std::vector<Int_t>   fMcPartTrackID;
    std::vector<Int_t>   fMcPartNumRecoHits;
    std::vector<Int_t>   fMcPartNumPFParts;
    std::vector<Int_t>   fMcPartBestPFPart;
    std::vector<Char_t>  fMcPartPrimary;
    
    // PFParticle information (meant to match the above)
    Int_t    fNumPFParticles;
    std::vector<Int_t>   fBestMcTrackID;
    std::vector<Int_t>   fNumMcTrackIDs;
    std::vector<Int_t>   fNumRecoHitsTotal;
    std::vector<Int_t>   fNumRecoHitsMatched;
    std::vector<Int_t>   fPFPartNumTracks;
    std::vector<Float_t> fPCAStartX;
    std::vector<Float_t> fPCAStartY;
    std::vector<Float_t> fPCAStartZ;
    std::vector<Float_t> fPCAEigenVal1;
    std::vector<Float_t> fPCAEigenVal2;
    std::vector<Float_t> fPCAEigenVal3;
    std::vector<Float_t> fPCAAxis1DirX;
    std::vector<Float_t> fPCAAxis1DirY;
    std::vector<Float_t> fPCAAxis1DirZ;
    std::vector<Float_t> fPCAAxis2DirX;
    std::vector<Float_t> fPCAAxis2DirY;
    std::vector<Float_t> fPCAAxis2DirZ;
    std::vector<Float_t> fPCAAxis3DirX;
    std::vector<Float_t> fPCAAxis3DirY;
    std::vector<Float_t> fPCAAxis3DirZ;
    
    std::vector<std::string> fProcessNameVec;
    
    // Associated reco track information
    std::vector<Int_t>   fNumTrackHits;
    std::vector<Float_t> fFitTrackLen;
    std::vector<Float_t> fTrackStartX;
    std::vector<Float_t> fTrackStartY;
    std::vector<Float_t> fTrackStartZ;
    std::vector<Float_t> fTrackStartDirX;
    std::vector<Float_t> fTrackStartDirY;
    std::vector<Float_t> fTrackStartDirZ;
    std::vector<Float_t> fTrackEndX;
    std::vector<Float_t> fTrackEndY;
    std::vector<Float_t> fTrackEndZ;
    std::vector<Float_t> fTrackEndDirX;
    std::vector<Float_t> fTrackEndDirY;
    std::vector<Float_t> fTrackEndDirZ;
    
    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service
    double                            fElectronsToGeV; // conversion factor
    
    // Define our output TTree here
    TTree*    fAnaTree;

}; // class PFParticleMcAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
PFParticleMcAna::PFParticleMcAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fAnaTree(0)
{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
PFParticleMcAna::~PFParticleMcAna()
{}
   
//-----------------------------------------------------------------------
void PFParticleMcAna::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Create our TTree and populate with basic information
    fAnaTree = tfs->make<TTree>("mcanalysis", "LAr Reco");
    fAnaTree->Branch("run",                  &fRun,                   "run/I");
    fAnaTree->Branch("event",                &fEvent,                 "event/I");
    fAnaTree->Branch("hits",                 &fNumHits,               "hits/I");
    fAnaTree->Branch("NoiseHits",            &fNoiseHits,             "NoiseHits/I");
    fAnaTree->Branch("NegativeTrackIDs",     &fNegTrackIds,           "NegativeTrackIDs/I");
    
    fMaxEntries = 300;
    
    // Start with the Monte Carlo information
    fAnaTree->Branch("NumMcParticles",       &fNumMcParticles,        "NumMcParticles/I");
    
    fPDGCode.resize(fMaxEntries, 0);
    fMcPartStartX.resize(fMaxEntries, 0.);
    fMcPartStartY.resize(fMaxEntries, 0.);
    fMcPartStartZ.resize(fMaxEntries, 0.);
    fMcPartStartDirX.resize(fMaxEntries, 0.);
    fMcPartStartDirY.resize(fMaxEntries, 0.);
    fMcPartStartDirZ.resize(fMaxEntries, 0.);
    fMcPartEndX.resize(fMaxEntries, 0.);
    fMcPartEndY.resize(fMaxEntries, 0.);
    fMcPartEndZ.resize(fMaxEntries, 0.);
    fMcPartEndDirX.resize(fMaxEntries, 0.);
    fMcPartEndDirY.resize(fMaxEntries, 0.);
    fMcPartEndDirZ.resize(fMaxEntries, 0.);
    fMcPartEne.resize(fMaxEntries, 0.);
    fMcPartMom.resize(fMaxEntries, 0.);
    fMcPartMass.resize(fMaxEntries, 0.);
    fMcPartTrackLen.resize(fMaxEntries, 0.);
    fMcPartTrackID.resize(fMaxEntries, 0);
    fMcPartNumRecoHits.resize(fMaxEntries, 0);
    fMcPartNumPFParts.resize(fMaxEntries, 0);
    fMcPartBestPFPart.resize(fMaxEntries, 0);
    fMcPartPrimary.resize(fMaxEntries, 0);

    fAnaTree->Branch("McPartPDGCode",            fPDGCode.data(),            "McPartPDGCode[NumMcParticles]/I");
    fAnaTree->Branch("McPartStartX",             fMcPartStartX.data(),       "McPartStartX[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartY",             fMcPartStartY.data(),       "McPartStartY[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartZ",             fMcPartStartZ.data(),       "McPartStartZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartDirX",          fMcPartStartDirX.data(),    "McPartStartDirX[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartDirY",          fMcPartStartDirY.data(),    "McPartStartDirY[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartDirZ",          fMcPartStartDirZ.data(),    "McPartStartDirZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndX",               fMcPartEndX.data(),         "McPartEndX[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndY",               fMcPartEndY.data(),         "McPartEndY[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndZ",               fMcPartEndZ.data(),         "McPartEndZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndDirX",            fMcPartEndDirX.data(),      "McPartEndDirX[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndDirY",            fMcPartEndDirY.data(),      "McPartEndDirY[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndDirZ",            fMcPartEndDirZ.data(),      "McPartEndDirZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartEne",                fMcPartEne.data(),          "McPartEne[NumMcParticles]/F");
    fAnaTree->Branch("McPartMom",                fMcPartMom.data(),          "McPartMom[NumMcParticles]/F");
    fAnaTree->Branch("McPartMass",               fMcPartMass.data(),         "McPartMass[NumMcParticles]/F");
    fAnaTree->Branch("McPartTrackLen",           fMcPartTrackLen.data(),     "McPartTrackLen[NumMcParticles]/F");
    fAnaTree->Branch("McPartTrackID",            fMcPartTrackID.data(),      "McPartTrackID[NumMcParticles]/I");
    fAnaTree->Branch("McPartNumRecoHits",        fMcPartNumRecoHits.data(),  "McPartNumRecoHits[NumMcParticles]/I");
    fAnaTree->Branch("McPartNumPFParts",         fMcPartNumPFParts.data(),   "McPartNumPFParts[NumMcParticles]/I");
    fAnaTree->Branch("McPartBestPFPart",         fMcPartBestPFPart.data(),   "McPartBestPFPart[NumMcParticles]/I");
    fAnaTree->Branch("McPartPrimary",            fMcPartPrimary.data(),      "McPartPrimary[NumMcParticles]/B");
    
    
    fProcessNameVec.resize(fMaxEntries, "processname  ");
    
    fAnaTree->Branch("McPartProcess",           &fProcessNameVec);
    
    // Now create PFParticle values to match the above
    fAnaTree->Branch("NumPFParticles",          &fNumPFParticles,            "NumPFParticles/I");

    fBestMcTrackID.resize(fMaxEntries, 0.);
    fNumMcTrackIDs.resize(fMaxEntries, 0.);
    fNumRecoHitsTotal.resize(fMaxEntries, 0.);
    fNumRecoHitsMatched.resize(fMaxEntries, 0.);
    fPFPartNumTracks.resize(fMaxEntries, 0.);
    fPCAStartX.resize(fMaxEntries, 0.);
    fPCAStartY.resize(fMaxEntries, 0.);
    fPCAStartZ.resize(fMaxEntries, 0.);
    fPCAEigenVal1.resize(fMaxEntries, 0.);
    fPCAAxis1DirX.resize(fMaxEntries, 0.);
    fPCAAxis1DirY.resize(fMaxEntries, 0.);
    fPCAAxis1DirZ.resize(fMaxEntries, 0.);
    fPCAAxis2DirX.resize(fMaxEntries, 0.);
    fPCAAxis2DirY.resize(fMaxEntries, 0.);
    fPCAAxis2DirZ.resize(fMaxEntries, 0.);
    fPCAAxis3DirX.resize(fMaxEntries, 0.);
    fPCAAxis3DirY.resize(fMaxEntries, 0.);
    fPCAAxis3DirZ.resize(fMaxEntries, 0.);
    
    fAnaTree->Branch("PFPartBestMcTrackID",      fBestMcTrackID.data(),      "PFPartBestMcTrackID[NumPFParticles]/I");
    fAnaTree->Branch("PFPartNumMcTrackIDs",      fNumMcTrackIDs.data(),      "PFPartNumMcTrackIDs[NumPFParticles]/I");
    fAnaTree->Branch("PFPartNumRecoHitsTotal",   fNumRecoHitsTotal.data(),   "PFPartNumRecoHitsTotal[NumPFParticles]/I");
    fAnaTree->Branch("PFPartNumRecoHitsMatched", fNumRecoHitsMatched.data(), "PFPartNumRecoHitsMatched[NumPFParticles]/I");
    fAnaTree->Branch("PFPartNumTracks",          fPFPartNumTracks.data(),    "PFPartNumTracks[NumPFParticles]/I");
    fAnaTree->Branch("PFPartPCAStartX",          fPCAStartX.data(),          "PFPartPCAStartX[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAStartY",          fPCAStartY.data(),          "PFPartPCAStartY[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAStartZ",          fPCAStartZ.data(),          "PFPartPCAStartZ[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAEigenVal1",       fPCAEigenVal1.data(),       "PFPartPCAEigenVal1[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAEigenVal2",       fPCAEigenVal2.data(),       "PFPartPCAEigenVal2[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAEigenVal3",       fPCAEigenVal3.data(),       "PFPartPCAEigenVal3[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis1DirX",       fPCAAxis1DirX.data(),       "PFPartPCAAxis1DirX[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis1DirY",       fPCAAxis1DirY.data(),       "PFPartPCAAxis1DirY[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis1DirZ",       fPCAAxis1DirZ.data(),       "PFPartPCAAxis1DirZ[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis2DirX",       fPCAAxis2DirX.data(),       "PFPartPCAAxis2DirX[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis2DirY",       fPCAAxis2DirY.data(),       "PFPartPCAAxis2DirY[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis2DirZ",       fPCAAxis2DirZ.data(),       "PFPartPCAAxis2DirZ[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis3DirX",       fPCAAxis3DirX.data(),       "PFPartPCAAxis3DirX[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis3DirY",       fPCAAxis3DirY.data(),       "PFPartPCAAxis3DirY[NumPFParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis3DirZ",       fPCAAxis3DirZ.data(),       "PFPartPCAAxis3DirZ[NumPFParticles]/F");
    
    fNumTrackHits.resize(fMaxEntries, 0.);
    fFitTrackLen.resize(fMaxEntries, 0.);
    fTrackStartX.resize(fMaxEntries, 0.);
    fTrackStartY.resize(fMaxEntries, 0.);
    fTrackStartZ.resize(fMaxEntries, 0.);
    fTrackStartDirX.resize(fMaxEntries, 0.);
    fTrackStartDirY.resize(fMaxEntries, 0.);
    fTrackStartDirZ.resize(fMaxEntries, 0.);
    fTrackEndX.resize(fMaxEntries, 0.);
    fTrackEndY.resize(fMaxEntries, 0.);
    fTrackEndZ.resize(fMaxEntries, 0.);
    fTrackEndDirX.resize(fMaxEntries, 0.);
    fTrackEndDirY.resize(fMaxEntries, 0.);
    fTrackEndDirZ.resize(fMaxEntries, 0.);
    
    fAnaTree->Branch("TrackNumRecoHits",         fNumTrackHits.data(),       "TrackNumRecoHits[NumPFParticles]/I");
    fAnaTree->Branch("TrackLength",              fFitTrackLen.data(),        "TrackLength[NumPFParticles]/F");
    fAnaTree->Branch("TrackStartX",              fTrackStartX.data(),        "TrackStartX[NumPFParticles]/F");
    fAnaTree->Branch("TrackStartY",              fTrackStartY.data(),        "TrackStartY[NumPFParticles]/F");
    fAnaTree->Branch("TrackStartZ",              fTrackStartZ.data(),        "TrackStartZ[NumPFParticles]/F");
    fAnaTree->Branch("TrackStartDirX",           fTrackStartDirX.data(),     "TrackStartDirX[NumPFParticles]/F");
    fAnaTree->Branch("TrackStartDirY",           fTrackStartDirY.data(),     "TrackStartDirY[NumPFParticles]/F");
    fAnaTree->Branch("TrackStartDirZ",           fTrackStartDirZ.data(),     "TrackStartDirZ[NumPFParticles]/F");
    fAnaTree->Branch("TrackEndX",                fTrackEndX.data(),          "TrackEndX[NumPFParticles]/F");
    fAnaTree->Branch("TrackEndY",                fTrackEndY.data(),          "TrackEndY[NumPFParticles]/F");
    fAnaTree->Branch("TrackEndZ",                fTrackEndZ.data(),          "TrackEndZ[NumPFParticles]/F");
    fAnaTree->Branch("TrackEndDirX",             fTrackEndDirX.data(),       "TrackEndDirX[NumPFParticles]/F");
    fAnaTree->Branch("TrackEndDirY",             fTrackEndDirY.data(),       "TrackEndDirY[NumPFParticles]/F");
    fAnaTree->Branch("TrackEndDirZ",             fTrackEndDirZ.data(),       "TrackEndDirZ[NumPFParticles]/F");
}
   
//-----------------------------------------------------------------------
void PFParticleMcAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void PFParticleMcAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fSimulationProducerLabel    = p.get< std::string >("SimulationLabel");
    fMcHitCollectionModuleLabel = p.get< std::string >("McHitFinderLabel");
    fPFParticleProducerLabel    = p.get< std::string >("PFParticleLabel");
    fHitProducerLabel           = p.get< std::string >("HitLabel");
    fClusterProducerLabel       = p.get< std::string >("ClusterProducerLabel");
    fTrackProducerLabel         = p.get< std::string >("TrackProducerLabel");
    fCosmicProducerLabel        = p.get< std::string >("CosmicProducerLabel");
    fFlashProducerLabel         = p.get< std::string >("FlashProducerLabel");
    fSelectedPDG                = p.get< int         >("PDGcode");
    
    return;
}
    
//-----------------------------------------------------------------------
void PFParticleMcAna::PrepareEvent(const art::Event &evt)
{
    fRun                   = evt.run();
    fEvent                 = evt.id().event();
    fSubRun                = evt.subRun();
    fNumHits               = 0;
    fNoiseHits             = 0;
    fNegTrackIds           = 0;
    fNumMcParticles        = 0;
    fNumPFParticles        = 0;
    
    // Now clear the tuple variables too
    fPDGCode.assign(fMaxEntries, 0);
    fMcPartEndX.assign(fMaxEntries, 0.);
    fMcPartStartY.assign(fMaxEntries, 0.);
    fMcPartStartZ.assign(fMaxEntries, 0.);
    fMcPartStartDirX.assign(fMaxEntries, 0.);
    fMcPartStartDirY.assign(fMaxEntries, 0.);
    fMcPartStartDirZ.assign(fMaxEntries, 0.);
    fMcPartEne.assign(fMaxEntries, 0.);
    fMcPartMom.assign(fMaxEntries, 0.);
    fMcPartMass.assign(fMaxEntries, 0.);
    fMcPartEndDirX.assign(fMaxEntries, 0.);
    fMcPartEndDirY.assign(fMaxEntries, 0.);
    fMcPartEndDirZ.assign(fMaxEntries, 0.);
    fMcPartTrackLen.assign(fMaxEntries, 0.);
    fMcPartTrackID.assign(fMaxEntries, 0);
    fMcPartNumRecoHits.assign(fMaxEntries, 0);
    fMcPartNumPFParts.assign(fMaxEntries, 0);
    fMcPartBestPFPart.assign(fMaxEntries, 0);
    fMcPartPrimary.assign(fMaxEntries, 0);
    
    fMcPartEndX.assign(fMaxEntries, 0.);
    fMcPartEndY.assign(fMaxEntries, 0.);
    fMcPartEndZ.assign(fMaxEntries, 0.);
    
    // PFParticle information (meant to match the above)
    fBestMcTrackID.assign(fMaxEntries, 0);
    fNumMcTrackIDs.assign(fMaxEntries, 0);
    fNumRecoHitsTotal.assign(fMaxEntries, 0);
    fNumRecoHitsMatched.assign(fMaxEntries, 0);
    fPFPartNumTracks.assign(fMaxEntries, 0);
    fPCAStartX.assign(fMaxEntries, 0.);
    fPCAStartY.assign(fMaxEntries, 0.);
    fPCAStartZ.assign(fMaxEntries, 0.);
    fPCAEigenVal1.assign(fMaxEntries, 0.);
    fPCAEigenVal2.assign(fMaxEntries, 0.);
    fPCAEigenVal3.assign(fMaxEntries, 0.);
    fPCAAxis1DirX.assign(fMaxEntries, 0.);
    fPCAAxis1DirY.assign(fMaxEntries, 0.);
    fPCAAxis1DirZ.assign(fMaxEntries, 0.);
    fPCAAxis2DirX.assign(fMaxEntries, 0.);
    fPCAAxis2DirY.assign(fMaxEntries, 0.);
    fPCAAxis2DirZ.assign(fMaxEntries, 0.);
    fPCAAxis3DirX.assign(fMaxEntries, 0.);
    fPCAAxis3DirY.assign(fMaxEntries, 0.);
    fPCAAxis3DirZ.assign(fMaxEntries, 0.);
    
    // Associated reco track information
    fNumTrackHits.assign(fMaxEntries, 0);
    fFitTrackLen.assign(fMaxEntries, 0.);
    fTrackStartX.assign(fMaxEntries, 0.);
    fTrackStartY.assign(fMaxEntries, 0.);
    fTrackStartZ.assign(fMaxEntries, 0.);
    fTrackStartDirX.assign(fMaxEntries, 0.);
    fTrackStartDirY.assign(fMaxEntries, 0.);
    fTrackStartDirZ.assign(fMaxEntries, 0.);
    fTrackEndX.assign(fMaxEntries, 0.);
    fTrackEndY.assign(fMaxEntries, 0.);
    fTrackEndZ.assign(fMaxEntries, 0.);
    fTrackEndDirX.assign(fMaxEntries, 0.);
    fTrackEndDirY.assign(fMaxEntries, 0.);
    fTrackEndDirZ.assign(fMaxEntries, 0.);
}

//-----------------------------------------------------------------------
void PFParticleMcAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    this->PrepareEvent(event);
    
    // The first step is to attempt to recover the collection of MCHits that
    // we will need for doing our PFParticle to MC matching
    art::Handle< std::vector<sim::MCHitCollection> > mcHitCollectionHandle;
    event.getByLabel(fMcHitCollectionModuleLabel, mcHitCollectionHandle);
    
    if (!mcHitCollectionHandle.isValid())
    {
        mf::LogDebug("PFParticleMcAna") << "===>> NO McHitColllection found for run: " << fRun << ", event: " << fEvent << std::endl;
        fAnaTree->Fill();
        return;
    }
    
    // Recover this into a local stl version
    const std::vector<sim::MCHitCollection> &mcHitCollectionVec = *mcHitCollectionHandle;
    
    // Recover the PFParticles, the main products for our next major loop
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    event.getByLabel(fPFParticleProducerLabel, pfParticleHandle);
    
    // no point continuing if no PFParticles
    if (!pfParticleHandle.isValid())
    {
        mf::LogDebug("PFParticleMcAna") << "===>> NO PFParticle collection found for run: " << fRun << ", event: " << fEvent << std::endl;
        fAnaTree->Fill();
        return;
    }

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

    // The MCParticle objects are not necessarily in any particular
    // order. Since we may have to search the list of particles, let's
    // put them into a sorted map that will make searching fast and
    // easy. To save both space and time, the map will not contain a
    // copy of the MCParticle, but a pointer to it.
    std::map< int, const simb::MCParticle* > particleMap;
    
    // Before starting to loop through the particles, we are going to want to
    // build a mapping between hits and track id's
    // Start by recovering info from the event store
    art::Handle< std::vector<recob::Hit> > hitHandle;
    event.getByLabel(fHitProducerLabel, hitHandle);
    
    // our ultimate goal here are the following two maps:
    // 1) a map between a recob::Hit and the track ID (meaning MCParticle)
    // 2) a map between the trackID (MCParticle) and the hits it created
    HitToParticleMap hitToParticleMap;
    ParticleToHitMap particleToHitMap;
    
    // Recover the list of hits in stl vector form
    const std::vector<recob::Hit>& recoHitVec = *hitHandle;
    
    // Build out our MCParticle <---> reco Hit maps
    MakeHitParticleMaps(mcHitCollectionVec, recoHitVec, hitToParticleMap, particleToHitMap);

    // This loop will build the map between track ID and the MCParticle related to it
    for ( auto const& particle : (*particleHandle) )
    {
        int trackID = particle.TrackId();

        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[trackID] = &particle;
    } // loop over all particles in the event.
    
    // Now define the maps relating pfparticles to tracks
    TrackIDToPFParticleVecMap trackIDToPFParticleMap;
    PFParticleToTrackHit2DMap pfParticleToTrackHitMap;
    
    // Something to keep track of number of hits associated to a cluster
    std::map<const recob::PFParticle*, int> pfParticleHitCntMap;

    // Call the method for building out the PFParticle data structures
    MakePFParticleMaps(event, pfParticleHandle, hitToParticleMap, pfParticleToTrackHitMap, trackIDToPFParticleMap, pfParticleHitCntMap);
    
    // Always a handy thing to have hidden in your code:
//    const double radToDegrees = 180. / 3.14159265;
    art::ServiceHandle<util::LArProperties> larProperties;
    
    // Recover the collection of associations between PFParticles and tracks, this will
    // be the mechanism by which we actually deal with tracks
    art::FindManyP<recob::Track> trackAssns(pfParticleHandle, event, fTrackProducerLabel);
    
    // Recover the collection of associations between PFParticles and PCAxis objects
    art::FindManyP<recob::PCAxis> pcaxisAssns(pfParticleHandle, event, fPFParticleProducerLabel);
    
    // Now retrieve a handle to the fit tracks, we will only use this when we are sure it is valid
    //art::Handle<std::vector<recob::Track> > trackHandle;
    //event.getByLabel(fTrackProducerLabel, trackHandle);
    
    // Recover two sets of associations for linking cosmic tags to the tracks
    //art::FindManyP<anab::CosmicTag> cosmicAssns(trackHandle, event, fCosmicProducerLabel);
    //art::FindManyP<anab::CosmicTag> flashAssns(trackHandle,  event, fFlashProducerLabel);

    // We will loop through all of the MCParticles and look for those particles that produce hits
    for (size_t particleIdx = 0; particleIdx < particleHandle->size(); particleIdx++)
    {
        art::Ptr<simb::MCParticle> particle(particleHandle, particleIdx);

        // Recover the track ID, for historical reasons we call it "best"
        int bestTrackID = particle->TrackId();
        
        // Did this mc particle leave hits in the TPC?
        ParticleToHitMap::iterator particleToHitItr = particleToHitMap.find(bestTrackID);
        
        // Let's get the total number of "true" hits that are created by this MCParticle
        // Count number of hits in each view
        int nTrueMcHits(0);
        
        if (particleToHitItr != particleToHitMap.end())  nTrueMcHits = particleToHitItr->second.size();
        
        // It is pointless to go through the rest of the loop if there are no hits to analyze
        if (nTrueMcHits < 1) continue;
        
        // Recover the particle's identity
        int trackPDGCode = particle->PdgCode();
        
        // Can we check the range of this particle?
        if(particle->E() < 0.001*particle->Mass() + 0.05)
        {
            if (fabs(trackPDGCode) == 13 && particle->Process() == "primary")
            {
                mf::LogDebug("PFParticleMcAna") << "~~~~> Primary muon below energy threshold, rejecting" << std::endl;
            }
        }
        
        // Calculate the x offset due to nonzero mc particle time.
        double mctime = particle->T();                                    // nsec
        double mcdx   = mctime * 1.e-3 * larProperties->DriftVelocity();  // cm
            
        // Calculate the length of this mc particle inside the fiducial volume.
        TVector3 mcstart;
        TVector3 mcend;
        TVector3 mcstartmom;
        TVector3 mcendmom;
        
        double mcTrackLen = length(*particle, mcdx, mcstart, mcend, mcstartmom, mcendmom);
        
        // Set the momentum vector to a unit direction vector
        if (mcstartmom.Mag() > 0.) mcstartmom.SetMag(1.);
        if (mcendmom.Mag()   > 0.) mcendmom.SetMag(1.);
        
        std::string process   = particle->Process();
        bool        isPrimary = process == "primary";

        // Define the various quantities we want to fill if we find a PFParticle associated to this MCParticle
        const recob::PFParticle* pfParticle(0);
        const recob::Track*      tTrack(0);
        int                      bestCnt(0);
        double                   trackLen(0.);
        int                      nTotalClusterHits(0);
        int                      numPFParticles(0);
        int                      numPFPartTracks(0);
        int                      bestPFParticleID(-1);
        int                      numTrackHits(0);
        TVector3                 pcaxisPos(0.,0.,-10.);
        TVector3                 pcaEigenVals(0.,0.,0.);
        TVector3                 pcaxis1Dir(0.,0.,1.);
        TVector3                 pcaxis2Dir(0.,0.,1.);
        TVector3                 pcaxis3Dir(0.,0.,1.);
        
        // Now check to see that a track is associated with this MCParticle (trackID)
        TrackIDToPFParticleVecMap::iterator trackIDToPFParticleItr = trackIDToPFParticleMap.find(bestTrackID);
        
        // Proceed here if we have at least a chance to match to PFParticle
        if (trackIDToPFParticleItr != trackIDToPFParticleMap.end())
        {
            // We want to keep track of this...
            numPFParticles = trackIDToPFParticleItr->second.size();
            
            // Ok, it is possible for one MCParticle to have more than one PFParticle (think of a broken cluster) so we
            // need to loop through the possible PFParticles (and fully expect this to normally be only one PFParticle)
            for(const auto& tmpPart : trackIDToPFParticleItr->second)
            {
                // For the given PFParticle, recover the element corresponding to the current MCParticle Track ID
                TrackIDToHit2DMap::iterator tmpPartTrackItr = pfParticleToTrackHitMap[tmpPart].find(bestTrackID);

                // In theory this condition is satisfied at least once per loop
                if (tmpPartTrackItr != pfParticleToTrackHitMap[tmpPart].end())
                {
                    // Recover the hit count
                    int trackHitCnt = tmpPartTrackItr->second.size();
            
                    // For the PFParticles that might be associated to the single MCParticle, select the one with the most hits
                    if (trackHitCnt > bestCnt)
                    {
                        bestCnt    = trackHitCnt;
                        pfParticle = tmpPart;
                    }
                }
            }

            // If we found a matching PFParticle then look up the associated information
            if (pfParticle)
            {
                bestPFParticleID = pfParticle->Self();
        
                // Start by looking for an associated fit track(s)
                if (trackAssns.isValid())
                {
                    std::vector<art::Ptr<recob::Track> > tTrackVec = trackAssns.at(pfParticle->Self());
                    
                    if (!tTrackVec.empty())
                    {
                        art::Ptr<recob::Track> trackPtr = tTrackVec[0];
                        size_t                 nPoints  = trackPtr->NumberTrajectoryPoints();
                    
                        for(size_t idx = 1; idx < tTrackVec.size(); idx++)
                        {
                            // In the event more than one track was associated to the PFParticle then
                            // accept the one with the most trajectory points (which is NOT the same
                            // thing as the most recob::Hits!)
                            if (tTrackVec[idx]->NumberTrajectoryPoints() > nPoints)
                            {
                                trackPtr = tTrackVec[idx];
                                nPoints  = trackPtr->NumberTrajectoryPoints();
                            }
                        }
                        

                        art::PtrVector<recob::Track> trackPtrVec;
                         
                        trackPtrVec.push_back(trackPtr);

                        // Ok, now we want to get the number of actual recob::Hits in the track
                        art::FindManyP<recob::Hit> trackHitAssns(trackPtrVec, event, fTrackProducerLabel);
                         
                        if (trackHitAssns.isValid())
                        {
                            numTrackHits = trackHitAssns.at(0).size();
                        }
                        
                        tTrack        = trackPtr.get();
                        trackLen      = length(tTrack);
                    }
                    
                    numPFPartTracks = tTrackVec.size();
                }
                
                // Look for an associated PCAxis
                if (pcaxisAssns.isValid())
                {
                    std::vector<art::Ptr<recob::PCAxis> > pcaxisVec = pcaxisAssns.at(pfParticle->Self());
                    
                    if (!pcaxisVec.empty())
                    {
                        // Take the first element, in theory this should correspond to the skeleton axis
                        const art::Ptr<recob::PCAxis> pcaxis = pcaxisVec[0];
                        
                        pcaxisPos    = TVector3(pcaxis->getAvePosition()[0],pcaxis->getAvePosition()[1],pcaxis->getAvePosition()[2]);
                        pcaEigenVals = TVector3(pcaxis->getEigenValues()[0],pcaxis->getEigenValues()[1],pcaxis->getEigenValues()[2]);
                        pcaxis1Dir   = TVector3(pcaxis->getEigenVectors()[0][0],pcaxis->getEigenVectors()[0][1],pcaxis->getEigenVectors()[0][2]);
                        pcaxis2Dir   = TVector3(pcaxis->getEigenVectors()[1][0],pcaxis->getEigenVectors()[1][1],pcaxis->getEigenVectors()[1][2]);
                        pcaxis3Dir   = TVector3(pcaxis->getEigenVectors()[2][0],pcaxis->getEigenVectors()[2][1],pcaxis->getEigenVectors()[2][2]);
                    }
                }
        
                // MicroBooNE definitions of efficiency and purity:
                // efficiency E = number of true hits in the cluster / number of true hits from MC particle  ("completeness")
                // purity P     = number of true hits in the cluster / all hits in the cluster
                nTotalClusterHits = pfParticleHitCntMap[pfParticle];
//                hitEfficiency     = double(nMatchedHits) / double(nTrueMcHits);
//                hitPurity         = double(nMatchedHits) / double(nTotalClusterHits);
            }
        }
        
        // Ok, decoding done, now transfer to our output ntuple
        // Start with MC...
        if (fNumMcParticles < fMaxEntries)
        {
            fPDGCode[fNumMcParticles]           = trackPDGCode;
            fMcPartStartX[fNumMcParticles]      = mcstart.X();
            fMcPartStartY[fNumMcParticles]      = mcstart.Y();
            fMcPartStartZ[fNumMcParticles]      = mcstart.Z();
            fMcPartStartDirX[fNumMcParticles]   = mcstartmom.X();
            fMcPartStartDirY[fNumMcParticles]   = mcstartmom.Y();
            fMcPartStartDirZ[fNumMcParticles]   = mcstartmom.Z();
            fMcPartEndX[fNumMcParticles]        = mcend.X();
            fMcPartEndY[fNumMcParticles]        = mcend.Y();
            fMcPartEndZ[fNumMcParticles]        = mcend.Z();
            fMcPartEndDirX[fNumMcParticles]     = mcendmom.X();
            fMcPartEndDirY[fNumMcParticles]     = mcendmom.Y();
            fMcPartEndDirZ[fNumMcParticles]     = mcendmom.Z();
            fMcPartEne[fNumMcParticles]         = particle->E();
            fMcPartMass[fNumMcParticles]        = particle->Mass();
            fMcPartTrackLen[fNumMcParticles]    = mcTrackLen;
            fMcPartTrackID[fNumMcParticles]     = bestTrackID;
            fMcPartNumRecoHits[fNumMcParticles] = nTrueMcHits;
            fMcPartNumPFParts[fNumMcParticles]  = numPFParticles;
            fMcPartBestPFPart[fNumMcParticles]  = bestPFParticleID;
            fProcessNameVec[fNumMcParticles]    = process;
            fMcPartPrimary[fNumMcParticles]     = isPrimary;
            
            fNumMcParticles++;
        }
        
        // Now do PFParticle corresponding to the above
        if (fNumPFParticles < fMaxEntries)
        {
            fBestMcTrackID[fNumPFParticles]      = bestTrackID;
            fNumMcTrackIDs[fNumPFParticles]      = numPFParticles;
            fNumRecoHitsTotal[fNumPFParticles]   = nTotalClusterHits;
            fNumRecoHitsMatched[fNumPFParticles] = bestCnt;
            fPFPartNumTracks[fNumPFParticles]    = numPFPartTracks;
            fPCAStartX[fNumPFParticles]          = pcaxisPos[0];
            fPCAStartY[fNumPFParticles]          = pcaxisPos[1];
            fPCAStartZ[fNumPFParticles]          = pcaxisPos[2];
            fPCAEigenVal1[fNumPFParticles]       = pcaEigenVals[0];
            fPCAEigenVal2[fNumPFParticles]       = pcaEigenVals[1];
            fPCAEigenVal3[fNumPFParticles]       = pcaEigenVals[2];
            fPCAAxis1DirX[fNumPFParticles]       = pcaxis1Dir[0];
            fPCAAxis1DirY[fNumPFParticles]       = pcaxis1Dir[1];
            fPCAAxis1DirZ[fNumPFParticles]       = pcaxis1Dir[2];
            fPCAAxis2DirX[fNumPFParticles]       = pcaxis2Dir[0];
            fPCAAxis2DirY[fNumPFParticles]       = pcaxis2Dir[1];
            fPCAAxis2DirZ[fNumPFParticles]       = pcaxis2Dir[2];
            fPCAAxis3DirX[fNumPFParticles]       = pcaxis3Dir[0];
            fPCAAxis3DirY[fNumPFParticles]       = pcaxis3Dir[1];
            fPCAAxis3DirZ[fNumPFParticles]       = pcaxis3Dir[2];
            
            fNumTrackHits[fNumPFParticles]       = numTrackHits;
            fFitTrackLen[fNumPFParticles]        = trackLen;
            
            if (tTrack)
            {
                fTrackStartX[fNumPFParticles]    = tTrack->Vertex().X();
                fTrackStartY[fNumPFParticles]    = tTrack->Vertex().Y();
                fTrackStartZ[fNumPFParticles]    = tTrack->Vertex().Z();
                fTrackStartDirX[fNumPFParticles] = tTrack->VertexDirection().X();
                fTrackStartDirY[fNumPFParticles] = tTrack->VertexDirection().Y();
                fTrackStartDirZ[fNumPFParticles] = tTrack->VertexDirection().Z();
                fTrackEndX[fNumPFParticles]      = tTrack->End().X();
                fTrackEndY[fNumPFParticles]      = tTrack->End().Y();
                fTrackEndZ[fNumPFParticles]      = tTrack->End().Z();
                fTrackEndDirX[fNumPFParticles]   = tTrack->EndDirection().X();
                fTrackEndDirY[fNumPFParticles]   = tTrack->EndDirection().Y();
                fTrackEndDirZ[fNumPFParticles]   = tTrack->EndDirection().Z();
            }
            
            fNumPFParticles++;
        }
        
        // We can't proceed if we have run off the end of our allocated space
        if (fNumMcParticles >= fMaxEntries) break;
    }

    fAnaTree->Fill();
    
    return;
}

// This builds out the hit to particle maps
//----------------------------------------------------------------------------
void PFParticleMcAna::MakeHitParticleMaps(const std::vector<sim::MCHitCollection>& mcHitCollectionVec,
                                          const std::vector<recob::Hit>&           recoHitVec,
                                          HitToParticleMap&                        hitToParticleMap,
                                          ParticleToHitMap&                        particleToHitMap)
{
    //we're gonna probably need the time service to convert hit times to TDCs
    art::ServiceHandle<util::TimeService> timeService;
    
    // Ok, so this loop obviously takes the MC information and builds two maps
    // 1) a map from a Hit2D object to the track ID's that made it
    // 2) the reverse map, going from track ID to Hit2D object
    for (const recob::Hit& hit : recoHitVec)
    {
        const geo::WireID& wireId = hit.WireID();
        
        unsigned int channel = fGeometry->PlaneWireToChannel(wireId.Plane, wireId.Wire, wireId.TPC, wireId.Cryostat);
        
        if (channel < mcHitCollectionVec.size())
        {
            const std::vector<sim::MCHit>& mcHitVec = mcHitCollectionVec[channel];
            
            int start_tdc = timeService->TPCTick2TDC( hit.StartTick() );
            int end_tdc   = timeService->TPCTick2TDC( hit.EndTick()   );
            
            sim::MCHit startTime;
            sim::MCHit endTime;
            
            startTime.SetTime(start_tdc, 0);
            endTime.SetTime(end_tdc, 0);
            
            std::vector<sim::MCHit>::const_iterator startItr = std::lower_bound(mcHitVec.begin(), mcHitVec.end(), startTime);
            std::vector<sim::MCHit>::const_iterator endItr   = std::upper_bound(startItr,         mcHitVec.end(), endTime);
            
            if (startItr != mcHitVec.end())
            {
                while(startItr != endItr)
                {
                    int trackID = (*startItr++).PartTrackId();
                    
                    if (trackID < 0)
                    {
                        trackID = -trackID;
                        fNegTrackIds++;
                    }
                    
                    hitToParticleMap[&hit].insert(trackID);
                    particleToHitMap[trackID].insert(&hit);
                }
            }
            else fNoiseHits++;
        }
        else fNoiseHits++;
    }

    return;
}
    
// This will be used to sort PFParticles in the map below
class PFParticleMcAna::SortPFParticleVec
{
    /**
     * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
     */
public:
    SortPFParticleVec(const PFParticleMcAna::PFParticleHitCntMap& pfPartCntMap) : fPFPartCntMap(pfPartCntMap) {}
    
    bool operator()(const recob::PFParticle* left, const recob::PFParticle* right)
    {
        size_t numHitsLeft  = fPFPartCntMap.at(left);
        size_t numHitsRight = fPFPartCntMap.at(right);
        
        return numHitsLeft > numHitsRight;
    }
private:
    const PFParticleMcAna::PFParticleHitCntMap& fPFPartCntMap;
};

// Length of reconstructed track.
//----------------------------------------------------------------------------
void PFParticleMcAna::MakePFParticleMaps(const art::Event&                                   event,
                                         const art::Handle<std::vector<recob::PFParticle> >& pfParticleHandle,
                                         const HitToParticleMap&                             hitToParticleMap,
                                         PFParticleToTrackHit2DMap&                          pfParticleToTrackHit2DMap,
                                         TrackIDToPFParticleVecMap&                          trackIDToPFParticleVecMap,
                                         PFParticleHitCntMap&                                pfParticleHitCntMap)
{
    // Recover pfparticle to cluster associations
    art::FindManyP<recob::Cluster> pfParticleClusterAssns(pfParticleHandle, event, fPFParticleProducerLabel);

    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    event.getByLabel(fClusterProducerLabel, clusterHandle);
    
    // Now recover the collection of associations to hits
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, event, fClusterProducerLabel);
    
    // Commence looping over PFParticles
    for(size_t pfPartIdx = 0; pfPartIdx < pfParticleHandle->size(); pfPartIdx++)
    {
        // Get our PFParticle
        art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfPartIdx);
        
        // Get the list of clusters associated to this PFParticle
        std::vector<art::Ptr<recob::Cluster> > clusterVec = pfParticleClusterAssns.at(pfPartIdx);
        
        // Number 2D hits this PFParticle
        int numPfParticleHits(0);
        
        // To categorize the PFParticle (associate to MCParticle), we will
        // create and fill an instance of a TrackIDToHitMap.
        // So, for a given PFParticle we will have a map of the track ID's (MCParticles)
        // and the hits that they contributed energy too
        pfParticleToTrackHit2DMap[pfParticle.get()] = TrackIDToHit2DMap();
        TrackIDToHit2DMap& trackIdToHit2DMap        = pfParticleToTrackHit2DMap[pfParticle.get()];
        
        // Loop over the associated clusters
        for(const auto& cluster : clusterVec)
        {
            std::vector<art::Ptr<recob::Hit> > hitVec;
            
            try {
                hitVec = clusterHitAssns.at(cluster->ID());
            }
            catch(...)
            {
                continue;
            }
            
            // Keep track of this to hopefully save some time later
            pfParticleHitCntMap[pfParticle.get()] += hitVec.size();
            
            // Something to count MCParticles contributing here
            std::set<int> trackIdCntSet;
            int           nMultiParticleHits(0);
            int           nNoParticleHits(0);
            
            numPfParticleHits += hitVec.size();
            
            // To fill out this map we now loop over the recovered Hit2D's and stuff into the map
            for (const auto& hit : hitVec)
            {
                // Given the Hit2D, recover the list of asscociated track ID's (MCParticles)
                HitToParticleMap::const_iterator hitToParticleMapItr = hitToParticleMap.find(hit.get());
                
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
                    nNoParticleHits++;
                }
            }
        } // end loop over clusters
        
        // Make sure something happened here...
        if (!trackIdToHit2DMap.empty())
        {
            // Now spin through the trackIdToHit2DMap to build the reverse map of above, taking track ID's to PFParticles...
            for (const auto& trackItr : trackIdToHit2DMap)
            {
                trackIDToPFParticleVecMap[trackItr.first].push_back(pfParticle.get());
            }
        }
        else
        {
            mf::LogDebug("PFParticleMcAna") << "***>> No PFParticle to MCParticle match made for run " << fRun << ", event: " << fEvent << std::endl;
        }
    } // end of loop over the PFParticle collection
    
    // One last spin through to sort the PFParticles in the track ID to PFParticle map
    for (auto& trackItr : trackIDToPFParticleVecMap)
    {
        std::sort(trackItr.second.begin(), trackItr.second.end(), SortPFParticleVec(pfParticleHitCntMap));
    }
    
    return;
}
    
// Length of reconstructed track.
//----------------------------------------------------------------------------
double PFParticleMcAna::length(const recob::Track* track)
{
    double   result(0.);
    TVector3 disp(track->LocationAtPoint(0));
    int      n(track->NumberTrajectoryPoints());
    
    for(int i = 1; i < n; ++i)
    {
        const TVector3& pos = track->LocationAtPoint(i);
        
        disp   -= pos;
        result += disp.Mag();
        disp    = pos;
    }
    
    return result;
}

// Length of MC particle.
//----------------------------------------------------------------------------
double PFParticleMcAna::length(const simb::MCParticle& part, double dx,
                              TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                              unsigned int tpc, unsigned int cstat)
{
    // Get services.
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    // Get fiducial volume boundary.
    
    //double xmin = 0.;
    //double xmax = 2.*geom->DetHalfWidth();
    //double ymin = -geom->DetHalfHeight();
    //double ymax = geom->DetHalfHeight();
    //double zmin = 0.;
    //double zmax = geom->DetLength();
    double xmin = 0.; //-2.*geom->DetHalfWidth();
    double xmax = 2.*geom->DetHalfWidth(); //2.*2.*geom->DetHalfWidth();
    double ymin = -(1.1*geom->DetHalfHeight());
    double ymax =  (1.1*geom->DetHalfHeight());
    double zmin = 0.;
    double zmax = geom->DetLength();
    
    double readOutWindowSize = detprop->ReadOutWindowSize();
    
    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;
    
    for(int i = 0; i < n; ++i)
    {
        TVector3 pos = part.Position(i).Vect();
        
        // Make fiducial cuts.  Require the particle to be within the physical volume of
        // the tpc, and also require the apparent x position to be within the expanded
        // readout frame.
        
        if(pos.X() >= xmin &&
           pos.X() <= xmax &&
           pos.Y() >= ymin &&
           pos.Y() <= ymax &&
           pos.Z() >= zmin &&
           pos.Z() <= zmax)
        {
            pos[0] += dx;
            double ticks = detprop->ConvertXToTicks(pos[0], 0, 0, 0);
            if(ticks >= 0. && ticks < readOutWindowSize)
            {
                if(first)
                {
                    start = pos;
                    startmom = part.Momentum(i).Vect();
                }
                else
                {
                    disp -= pos;
                    result += disp.Mag();
                }
                first = false;
                disp = pos;
                end = pos;
                endmom = part.Momentum(i).Vect();
            }
        }
    }
    
    return result;
}

// This macro has to be defined for this module to be invoked from a
// .fcl file; see PFParticleMcAna.fcl for more information.
DEFINE_ART_MODULE(PFParticleMcAna)

} // namespace PFParticleMcAna

#endif // PFParticleMcAna_module
