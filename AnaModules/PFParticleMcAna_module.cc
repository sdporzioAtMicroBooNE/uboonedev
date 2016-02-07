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
#include "Utilities/TimeService.h"
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
    
    // Find distance to nearest TPC boundary
    double distanceToTPCEdge(const TVector3& position) const;
    
    // This method is meant to be called at the start of the event
    void PrepareEvent(const art::Event& evt, int numColumns);

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
    typedef std::set<art::Ptr<recob::Hit> >                           Hit2DSet;                   // Need to count only unique hits
    typedef std::map< int, Hit2DSet >                                 TrackIDToHit2DMap;          // Maps Track ID's to lists of hits
    typedef std::vector<art::Ptr<recob::PFParticle>>                  PFParticleVec;              // Typedef a vector of PFParticles
    typedef std::map< int, PFParticleVec >                            TrackIDToPFParticleVecMap;  // Maps Track ID's to the PFParticles associated to them
    typedef std::map<art::Ptr<recob::PFParticle>, TrackIDToHit2DMap > PFParticleToTrackHit2DMap;  // Maps PFParticles to a map of associated Track ID's to their hits
    typedef std::map<art::Ptr<recob::PFParticle>, int>                PFParticleHitCntMap;        // Allows us to count reco hits per PFParticle directly

    // Forward declaration of object to handle sorting of PFParticles associated to Track IDs
    class SortPFParticleVec;
    
    // Define a function to fill all of the above
    void MakePFParticleMaps(const art::Event&                                   event,
                            const art::Handle<std::vector<recob::PFParticle> >& pfParticleHandle,
                            const HitToParticleMap&                             hitToParticleMap,
                            PFParticleToTrackHit2DMap&                          pfParticleToTrackHit2DMap,
                            TrackIDToPFParticleVecMap&                          trackIDtoPFParticleVecMap,
                            PFParticleHitCntMap&                                pfParticleHitCntMap);
    
    // Forward declaration of object to handle sorting of PFParticles associated to Track IDs
    class SortKTrackVec;
    
    // More useful typedefs, this time for fit tracks
    typedef std::vector<art::Ptr<recob::Track>>                  KTrackVec;              // Typedef a vector of PFParticles
    typedef std::map< int, KTrackVec >                           TrackIDToKTrackVecMap;  // Maps Track ID's to the PFParticles associated to them
    typedef std::map<art::Ptr<recob::Track>, TrackIDToHit2DMap > KTrackToTrackHit2DMap;  // Maps PFParticles to a map of associated Track ID's to their hits
    typedef std::map<art::Ptr<recob::Track>, int>                KTrackHitCntMap;        // Allows us to count reco hits per PFParticle directly
    
    // Define a function to fill all of the above
    void MakeKTrackMaps(const art::Event&                              event,
                        const art::Handle<std::vector<recob::Track> >& trackHandle,
                        const HitToParticleMap&                        hitToParticleMap,
                        KTrackToTrackHit2DMap&                         kTrackToTrackHit2DMap,
                        TrackIDToKTrackVecMap&                         trackIDtoKTrackVecMap,
                        KTrackHitCntMap&                               kTrackHitCntMap);

    // The parameters we'll read from the .fcl file.
    std::string fSimulationProducerLabel;    //> The name of the producer that tracked simulated particles through the detector
    std::string fMcHitCollectionModuleLabel; //> The name of the producer of MCHits
    std::string fPFParticleProducerLabel;    //> The name of the produder of the PFParticle hierarchy
    std::string fHitProducerLabel;           //> The name of the producer that created hits
    std::string fClusterProducerLabel;       //> The name of the producer that created clusters
    std::string fTrackProducerLabel;         //> The name of the producer that created the tracks
    std::string fPFCosmicProducerLabel;      //> The name of the producer that created the PFParticle cosmic tags
    std::string fCosmicProducerLabel;        //> The name of the producer that created cosmic tags
    std::string fFlashProducerLabel;         //> The name of the producer that created flash tags

    // Maximum number of entries per row in our tuple
    int         fMaxEntries;
    
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
    std::vector<Int_t>       fPDGCode;
    std::vector<Float_t>     fMcPartStartX;
    std::vector<Float_t>     fMcPartStartY;
    std::vector<Float_t>     fMcPartStartZ;
    std::vector<Float_t>     fMcPartStartDirX;
    std::vector<Float_t>     fMcPartStartDirY;
    std::vector<Float_t>     fMcPartStartDirZ;
    std::vector<Float_t>     fMcPartEndX;
    std::vector<Float_t>     fMcPartEndY;
    std::vector<Float_t>     fMcPartEndZ;
    std::vector<Float_t>     fMcPartEndDirX;
    std::vector<Float_t>     fMcPartEndDirY;
    std::vector<Float_t>     fMcPartEndDirZ;
    std::vector<Float_t>     fMcPartEne;
    std::vector<Float_t>     fMcPartMom;
    std::vector<Float_t>     fMcPartMass;
    std::vector<Float_t>     fMcPartTrackLen;
    std::vector<Int_t>       fMcPartTrackID;
    std::vector<Int_t>       fMcParentTrackID;
    std::vector<Int_t>       fMcPartNumRecoHits;
    std::vector<Int_t>       fMcPartNumUniqueHits;
    std::vector<Int_t>       fMcPartNumPFParts;
    std::vector<Int_t>       fMcPartBestPFPart;
    std::vector<UShort_t>    fMcPartPrimary;
    std::vector<UShort_t>    fMcNeutrinoDaughter;
    
    std::vector<std::string> fProcessNameVec;
    std::vector<std::string> fParProcNameVec;
    
    // PFParticle information (meant to match the above)
    Int_t    fNumPFParticles;
    std::vector<Int_t>       fPFParticleID;
    std::vector<Int_t>       fBestMcTrackID;
    std::vector<Int_t>       fNumMcTrackIDs;
    std::vector<Int_t>       fNumRecoHitsTotal;
    std::vector<Int_t>       fNumRecoHitsMatched;
    std::vector<Int_t>       fPFPartNumTracks;
    std::vector<Float_t>     fPCAStartX;
    std::vector<Float_t>     fPCAStartY;
    std::vector<Float_t>     fPCAStartZ;
    std::vector<Float_t>     fPCAEigenVal1;
    std::vector<Float_t>     fPCAEigenVal2;
    std::vector<Float_t>     fPCAEigenVal3;
    std::vector<Float_t>     fPCAAxis1DirX;
    std::vector<Float_t>     fPCAAxis1DirY;
    std::vector<Float_t>     fPCAAxis1DirZ;
    std::vector<Float_t>     fPCAAxis2DirX;
    std::vector<Float_t>     fPCAAxis2DirY;
    std::vector<Float_t>     fPCAAxis2DirZ;
    std::vector<Float_t>     fPCAAxis3DirX;
    std::vector<Float_t>     fPCAAxis3DirY;
    std::vector<Float_t>     fPCAAxis3DirZ;
    
    std::vector<Float_t>     fPFCompositeScore;
    std::vector<Float_t>     fPFCompositeTag;
    
    std::vector<Float_t>     fPFCosmicGeoScore;
    std::vector<Float_t>     fPFCosmicGeoTag;
    
    // Associated reco track information
    std::vector<Int_t>       fTrackID;
    std::vector<Int_t>       fNumKTracks;
    std::vector<Int_t>       fTrackBestMCTrkId;
    std::vector<Int_t>       fNumTrackHitsTotal;
    std::vector<Int_t>       fNumTrackHitsMatched;
    std::vector<Float_t>     fFitTrackLen;
    std::vector<Float_t>     fTrackStartX;
    std::vector<Float_t>     fTrackStartY;
    std::vector<Float_t>     fTrackStartZ;
    std::vector<Float_t>     fTrackStartDirX;
    std::vector<Float_t>     fTrackStartDirY;
    std::vector<Float_t>     fTrackStartDirZ;
    std::vector<Float_t>     fTrackEndX;
    std::vector<Float_t>     fTrackEndY;
    std::vector<Float_t>     fTrackEndZ;
    std::vector<Float_t>     fTrackEndDirX;
    std::vector<Float_t>     fTrackEndDirY;
    std::vector<Float_t>     fTrackEndDirZ;
    std::vector<Float_t>     fDistToEdge;
    
    std::vector<Float_t>     fCosmicGeoScore;
    std::vector<Float_t>     fCosmicGeoTag;
    std::vector<Float_t>     fFlashScore;
    
    std::vector<Float_t>     fCosmicStartX;
    std::vector<Float_t>     fCosmicStartY;
    std::vector<Float_t>     fCosmicStartZ;
    std::vector<Float_t>     fCosmicEndX;
    std::vector<Float_t>     fCosmicEndY;
    std::vector<Float_t>     fCosmicEndZ;
    
    // Other variables that will be shared between different methods.
    art::ServiceHandle<geo::Geometry>            fGeometry;       // pointer to Geometry service
    art::ServiceHandle<util::TimeService>        fTimeService;
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;
    double                                       fElectronsToGeV; // conversion factor
    
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
    fMcParentTrackID.resize(fMaxEntries, -1.);
    fMcPartNumRecoHits.resize(fMaxEntries, 0);
    fMcPartNumUniqueHits.resize(fMaxEntries, 0);
    fMcPartNumPFParts.resize(fMaxEntries, 0);
    fMcPartBestPFPart.resize(fMaxEntries, 0);
    fMcPartPrimary.resize(fMaxEntries, 0);
    fMcNeutrinoDaughter.resize(fMaxEntries, 0);

    fAnaTree->Branch("McPartPDGCode",            fPDGCode.data(),             "McPartPDGCode[NumMcParticles]/I");
    fAnaTree->Branch("McPartStartX",             fMcPartStartX.data(),        "McPartStartX[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartY",             fMcPartStartY.data(),        "McPartStartY[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartZ",             fMcPartStartZ.data(),        "McPartStartZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartDirX",          fMcPartStartDirX.data(),     "McPartStartDirX[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartDirY",          fMcPartStartDirY.data(),     "McPartStartDirY[NumMcParticles]/F");
    fAnaTree->Branch("McPartStartDirZ",          fMcPartStartDirZ.data(),     "McPartStartDirZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndX",               fMcPartEndX.data(),          "McPartEndX[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndY",               fMcPartEndY.data(),          "McPartEndY[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndZ",               fMcPartEndZ.data(),          "McPartEndZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndDirX",            fMcPartEndDirX.data(),       "McPartEndDirX[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndDirY",            fMcPartEndDirY.data(),       "McPartEndDirY[NumMcParticles]/F");
    fAnaTree->Branch("McPartEndDirZ",            fMcPartEndDirZ.data(),       "McPartEndDirZ[NumMcParticles]/F");
    fAnaTree->Branch("McPartEne",                fMcPartEne.data(),           "McPartEne[NumMcParticles]/F");
    fAnaTree->Branch("McPartMom",                fMcPartMom.data(),           "McPartMom[NumMcParticles]/F");
    fAnaTree->Branch("McPartMass",               fMcPartMass.data(),          "McPartMass[NumMcParticles]/F");
    fAnaTree->Branch("McPartTrackLen",           fMcPartTrackLen.data(),      "McPartTrackLen[NumMcParticles]/F");
    fAnaTree->Branch("McPartTrackID",            fMcPartTrackID.data(),       "McPartTrackID[NumMcParticles]/I");
    fAnaTree->Branch("McParentTrackID",          fMcParentTrackID.data(),     "McParentTrackID[NumMcParticles]/I");
    fAnaTree->Branch("McPartNumRecoHits",        fMcPartNumRecoHits.data(),   "McPartNumRecoHits[NumMcParticles]/I");
    fAnaTree->Branch("McPartNumUniqueHits",      fMcPartNumUniqueHits.data(), "McPartNumUniqueHits[NumMcParticles]/I");
    fAnaTree->Branch("McPartNumPFParts",         fMcPartNumPFParts.data(),    "McPartNumPFParts[NumMcParticles]/I");
    fAnaTree->Branch("McPartBestPFPart",         fMcPartBestPFPart.data(),    "McPartBestPFPart[NumMcParticles]/I");
    fAnaTree->Branch("McPartPrimary",            fMcPartPrimary.data(),       "McPartPrimary[NumMcParticles]/s");
    fAnaTree->Branch("McNeutrinoDaughter",       fMcNeutrinoDaughter.data(),  "McNeutrinoDaughter[NumMcParticles]/s");
    
    
    fProcessNameVec.resize(fMaxEntries, "processname  ");
    fParProcNameVec.resize(fMaxEntries, "processname  ");
    
    fAnaTree->Branch("McPartProcess",           &fProcessNameVec);
    fAnaTree->Branch("McPartParProc",           &fParProcNameVec);
    
    // Now create PFParticle values to match the above
    fAnaTree->Branch("NumPFParticles",          &fNumPFParticles,            "NumPFParticles/I");

    fPFParticleID.resize(fMaxEntries, -1);
    fBestMcTrackID.resize(fMaxEntries, -1);
    fNumMcTrackIDs.resize(fMaxEntries, 0);
    fNumRecoHitsTotal.resize(fMaxEntries, 0);
    fNumRecoHitsMatched.resize(fMaxEntries, 0);
    fPFPartNumTracks.resize(fMaxEntries, 0);
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
    fPFCosmicGeoScore.resize(fMaxEntries, 0.);
    fPFCosmicGeoTag.resize(fMaxEntries, 0.);
    
    fAnaTree->Branch("PFParticleID",             fPFParticleID.data(),       "PFParticleID[NumMcParticles]/I");
    fAnaTree->Branch("PFPartBestMcTrackID",      fBestMcTrackID.data(),      "PFPartBestMcTrackID[NumMcParticles]/I");
    fAnaTree->Branch("PFPartNumMcTrackIDs",      fNumMcTrackIDs.data(),      "PFPartNumMcTrackIDs[NumMcParticles]/I");
    fAnaTree->Branch("PFPartNumRecoHitsTotal",   fNumRecoHitsTotal.data(),   "PFPartNumRecoHitsTotal[NumMcParticles]/I");
    fAnaTree->Branch("PFPartNumRecoHitsMatched", fNumRecoHitsMatched.data(), "PFPartNumRecoHitsMatched[NumMcParticles]/I");
    fAnaTree->Branch("PFPartNumTracks",          fPFPartNumTracks.data(),    "PFPartNumTracks[NumMcParticles]/I");
    fAnaTree->Branch("PFPartPCAStartX",          fPCAStartX.data(),          "PFPartPCAStartX[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAStartY",          fPCAStartY.data(),          "PFPartPCAStartY[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAStartZ",          fPCAStartZ.data(),          "PFPartPCAStartZ[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAEigenVal1",       fPCAEigenVal1.data(),       "PFPartPCAEigenVal1[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAEigenVal2",       fPCAEigenVal2.data(),       "PFPartPCAEigenVal2[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAEigenVal3",       fPCAEigenVal3.data(),       "PFPartPCAEigenVal3[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis1DirX",       fPCAAxis1DirX.data(),       "PFPartPCAAxis1DirX[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis1DirY",       fPCAAxis1DirY.data(),       "PFPartPCAAxis1DirY[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis1DirZ",       fPCAAxis1DirZ.data(),       "PFPartPCAAxis1DirZ[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis2DirX",       fPCAAxis2DirX.data(),       "PFPartPCAAxis2DirX[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis2DirY",       fPCAAxis2DirY.data(),       "PFPartPCAAxis2DirY[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis2DirZ",       fPCAAxis2DirZ.data(),       "PFPartPCAAxis2DirZ[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis3DirX",       fPCAAxis3DirX.data(),       "PFPartPCAAxis3DirX[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis3DirY",       fPCAAxis3DirY.data(),       "PFPartPCAAxis3DirY[NumMcParticles]/F");
    fAnaTree->Branch("PFPartPCAAxis3DirZ",       fPCAAxis3DirZ.data(),       "PFPartPCAAxis3DirZ[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicGeoScore",         fPFCosmicGeoScore.data(),   "PFCosmicGeoScore[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicGeoTag",           fPFCosmicGeoTag.data(),     "PFCosmicGeoTag[NumMcParticles]/F");
    
    fCosmicStartX.resize(fMaxEntries, 0.);
    fCosmicStartY.resize(fMaxEntries, 0.);
    fCosmicStartZ.resize(fMaxEntries, 0.);
    fCosmicEndX.resize(fMaxEntries, 0.);
    fCosmicEndY.resize(fMaxEntries, 0.);
    fCosmicEndZ.resize(fMaxEntries, 0.);
    
    fAnaTree->Branch("PFCosmicStartX",           fCosmicStartX.data(),       "PFCosmicStartX[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicStartY",           fCosmicStartY.data(),       "PFCosmicStartY[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicStartZ",           fCosmicStartZ.data(),       "PFCosmicStartZ[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicEndX",             fCosmicEndX.data(),         "PFCosmicEndX[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicEndY",             fCosmicEndY.data(),         "PFCosmicEndY[NumMcParticles]/F");
    fAnaTree->Branch("PFCosmicEndZ",             fCosmicEndZ.data(),         "PFCosmicEndZ[NumMcParticles]/F");
    
    fTrackID.resize(fMaxEntries, -1);
    fNumKTracks.resize(fMaxEntries, 0);
    fTrackBestMCTrkId.resize(fMaxEntries, 0);
    fNumTrackHitsMatched.resize(fMaxEntries, 0);
    fNumTrackHitsTotal.resize(fMaxEntries, 0);
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
    fDistToEdge.resize(fMaxEntries, 0.);
    fCosmicGeoScore.resize(fMaxEntries, 0.);
    fCosmicGeoTag.resize(fMaxEntries, 0.);
    fFlashScore.resize(fMaxEntries, 0.);
    
    fAnaTree->Branch("TrackID",                  fTrackID.data(),              "TrackID[NumMcParticles]/I");
    fAnaTree->Branch("TrackNumPerMCParticle",    fNumKTracks.data(),           "TrackNumPerMCParticle[NumMcParticles]/I");
    fAnaTree->Branch("TrackBestMCTrkID",         fTrackBestMCTrkId.data(),     "TrackBestMCTrkId[NumMcParticles]/I");
    fAnaTree->Branch("TrackNumRecoHitsMatched",  fNumTrackHitsMatched.data(),  "TrackNumRecoHitsMatched[NumMcParticles]/I");
    fAnaTree->Branch("TrackNumRecoHitsTotal",    fNumTrackHitsTotal.data(),    "TrackNumRecoHitsTotal[NumMcParticles]/I");
    fAnaTree->Branch("TrackLength",              fFitTrackLen.data(),          "TrackLength[NumMcParticles]/F");
    fAnaTree->Branch("TrackStartX",              fTrackStartX.data(),          "TrackStartX[NumMcParticles]/F");
    fAnaTree->Branch("TrackStartY",              fTrackStartY.data(),          "TrackStartY[NumMcParticles]/F");
    fAnaTree->Branch("TrackStartZ",              fTrackStartZ.data(),          "TrackStartZ[NumMcParticles]/F");
    fAnaTree->Branch("TrackStartDirX",           fTrackStartDirX.data(),       "TrackStartDirX[NumMcParticles]/F");
    fAnaTree->Branch("TrackStartDirY",           fTrackStartDirY.data(),       "TrackStartDirY[NumMcParticles]/F");
    fAnaTree->Branch("TrackStartDirZ",           fTrackStartDirZ.data(),       "TrackStartDirZ[NumMcParticles]/F");
    fAnaTree->Branch("TrackEndX",                fTrackEndX.data(),            "TrackEndX[NumMcParticles]/F");
    fAnaTree->Branch("TrackEndY",                fTrackEndY.data(),            "TrackEndY[NumMcParticles]/F");
    fAnaTree->Branch("TrackEndZ",                fTrackEndZ.data(),            "TrackEndZ[NumMcParticles]/F");
    fAnaTree->Branch("TrackEndDirX",             fTrackEndDirX.data(),         "TrackEndDirX[NumMcParticles]/F");
    fAnaTree->Branch("TrackEndDirY",             fTrackEndDirY.data(),         "TrackEndDirY[NumMcParticles]/F");
    fAnaTree->Branch("TrackEndDirZ",             fTrackEndDirZ.data(),         "TrackEndDirZ[NumMcParticles]/F");
    fAnaTree->Branch("TrackDistToEdge",          fDistToEdge.data(),           "TrackDistToEdge[NumMcParticles]/F");
    fAnaTree->Branch("CosmicGeoScore",           fCosmicGeoScore.data(),       "CosmicGeoScore[NumMcParticles]/F");
    fAnaTree->Branch("CosmicGeoTag",             fCosmicGeoTag.data(),         "CosmicGeoTag[NumMcParticles]/F");
    fAnaTree->Branch("FlashScore",               fFlashScore.data(),           "FlashScore[NumMcParticles]/F");
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
    fPFCosmicProducerLabel      = p.get< std::string >("PFCosmicProducerLabel");
    fCosmicProducerLabel        = p.get< std::string >("CosmicProducerLabel");
    fFlashProducerLabel         = p.get< std::string >("FlashProducerLabel");
    fMaxEntries                 = p.get< int         >("MaxEntries", 300);

    std::cout << "PFParticleMcAna: geometry has half width: " << fGeometry->DetHalfWidth() << ", half height: " << fGeometry->DetHalfHeight() << ", length: " << fGeometry->DetLength() << std::endl;
    
    return;
}
    
//-----------------------------------------------------------------------
void PFParticleMcAna::PrepareEvent(const art::Event &evt, int numColumns)
{
    fRun                   = evt.run();
    fEvent                 = evt.id().event();
    fSubRun                = evt.subRun();
    fNumHits               = 0;
    fNoiseHits             = 0;
    fNegTrackIds           = 0;
    fNumMcParticles        = 0;
    fNumPFParticles        = 0;
    
    // Set the number of entries for each vector
    // But no more than the maximum requested when we started
    size_t maxEntries = numColumns < fMaxEntries ? numColumns : fMaxEntries;
    
    // Now clear the tuple variables too
    fPDGCode.assign(maxEntries, 0);
    fMcPartEndX.assign(maxEntries, 0.);
    fMcPartStartY.assign(maxEntries, 0.);
    fMcPartStartZ.assign(maxEntries, 0.);
    fMcPartStartDirX.assign(maxEntries, 0.);
    fMcPartStartDirY.assign(maxEntries, 0.);
    fMcPartStartDirZ.assign(maxEntries, 0.);
    fMcPartEne.assign(maxEntries, 0.);
    fMcPartMom.assign(maxEntries, 0.);
    fMcPartMass.assign(maxEntries, 0.);
    fMcPartEndDirX.assign(maxEntries, 0.);
    fMcPartEndDirY.assign(maxEntries, 0.);
    fMcPartEndDirZ.assign(maxEntries, 0.);
    fMcPartTrackLen.assign(maxEntries, 0.);
    fMcPartTrackID.assign(maxEntries, 0);
    fMcParentTrackID.assign(maxEntries, -1.);
    fMcPartNumRecoHits.assign(maxEntries, 0);
    fMcPartNumUniqueHits.assign(maxEntries, 0);
    fMcPartNumPFParts.assign(maxEntries, 0);
    fMcPartBestPFPart.assign(maxEntries, 0);
    fMcPartPrimary.assign(maxEntries, 0);
    fMcNeutrinoDaughter.assign(maxEntries, 0);
    
    fMcPartEndX.assign(maxEntries, 0.);
    fMcPartEndY.assign(maxEntries, 0.);
    fMcPartEndZ.assign(maxEntries, 0.);
    
    // PFParticle information (meant to match the above)
    fPFParticleID.assign(maxEntries, -1);
    fBestMcTrackID.assign(maxEntries, -1);
    fNumMcTrackIDs.assign(maxEntries, 0);
    fNumRecoHitsTotal.assign(maxEntries, 0);
    fNumRecoHitsMatched.assign(maxEntries, 0);
    fPFPartNumTracks.assign(maxEntries, 0);
    fPCAStartX.assign(maxEntries, 0.);
    fPCAStartY.assign(maxEntries, 0.);
    fPCAStartZ.assign(maxEntries, 0.);
    fPCAEigenVal1.assign(maxEntries, 0.);
    fPCAEigenVal2.assign(maxEntries, 0.);
    fPCAEigenVal3.assign(maxEntries, 0.);
    fPCAAxis1DirX.assign(maxEntries, 0.);
    fPCAAxis1DirY.assign(maxEntries, 0.);
    fPCAAxis1DirZ.assign(maxEntries, 0.);
    fPCAAxis2DirX.assign(maxEntries, 0.);
    fPCAAxis2DirY.assign(maxEntries, 0.);
    fPCAAxis2DirZ.assign(maxEntries, 0.);
    fPCAAxis3DirX.assign(maxEntries, 0.);
    fPCAAxis3DirY.assign(maxEntries, 0.);
    fPCAAxis3DirZ.assign(maxEntries, 0.);
    fPFCosmicGeoScore.assign(maxEntries, 0.);
    fPFCosmicGeoTag.assign(maxEntries, 0.);
    
    fCosmicStartX.assign(maxEntries, 0.);
    fCosmicStartY.assign(maxEntries, 0.);
    fCosmicStartZ.assign(maxEntries, 0.);
    fCosmicEndX.assign(maxEntries, 0.);
    fCosmicEndY.assign(maxEntries, 0.);
    fCosmicEndZ.assign(maxEntries, 0.);
    
    // Associated reco track information
    fTrackID.assign(maxEntries, -1);
    fTrackBestMCTrkId.assign(maxEntries, 0);
    fNumKTracks.assign(maxEntries, 0);
    fNumTrackHitsMatched.assign(maxEntries, 0);
    fNumTrackHitsTotal.assign(maxEntries, 0);
    fFitTrackLen.assign(maxEntries, 0.);
    fTrackStartX.assign(maxEntries, 0.);
    fTrackStartY.assign(maxEntries, 0.);
    fTrackStartZ.assign(maxEntries, 0.);
    fTrackStartDirX.assign(maxEntries, 0.);
    fTrackStartDirY.assign(maxEntries, 0.);
    fTrackStartDirZ.assign(maxEntries, 0.);
    fTrackEndX.assign(maxEntries, 0.);
    fTrackEndY.assign(maxEntries, 0.);
    fTrackEndZ.assign(maxEntries, 0.);
    fTrackEndDirX.assign(maxEntries, 0.);
    fTrackEndDirY.assign(maxEntries, 0.);
    fTrackEndDirZ.assign(maxEntries, 0.);
    fDistToEdge.assign(maxEntries, 0.);
    fCosmicGeoScore.assign(maxEntries, 0.);
    fCosmicGeoTag.assign(maxEntries, 0.);
    fFlashScore.assign(maxEntries, 0.);

    art::ServiceHandle<geo::Geometry>            fGeometry;       // pointer to Geometry service
    art::ServiceHandle<util::TimeService>        fTimeService;
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;
}

//-----------------------------------------------------------------------
void PFParticleMcAna::analyze(const art::Event& event)
{
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
    std::map<art::Ptr<recob::PFParticle>, int> pfParticleHitCntMap;

    // Call the method for building out the PFParticle data structures
    MakePFParticleMaps(event, pfParticleHandle, hitToParticleMap, pfParticleToTrackHitMap, trackIDToPFParticleMap, pfParticleHitCntMap);
    
    // Now retrieve a handle to the fit tracks which we need to build the fit track MCParticle maps
    art::Handle<std::vector<recob::Track> > trackHandle;
    event.getByLabel(fTrackProducerLabel, trackHandle);
    
    // Now define the maps relating fit tracks (we'll call KTracks) to MCParticle tracks
    TrackIDToKTrackVecMap trackIDToKTrackMap;
    KTrackToTrackHit2DMap kTrackToTrackHitMap;
    
    // Something to keep track of number of hits associated to a KTrack
    std::map<art::Ptr<recob::Track>, int> kTrackHitCntMap;
    
    // Call the method for building out the KTrack data structures
    MakeKTrackMaps(event, trackHandle, hitToParticleMap, kTrackToTrackHitMap, trackIDToKTrackMap, kTrackHitCntMap);
    
    // Always a handy thing to have hidden in your code:
//    const double radToDegrees = 180. / 3.14159265;
    art::ServiceHandle<util::LArProperties> larProperties;
    
    // Recover the collection of associations between PFParticles and tracks, this will
    // be the mechanism by which we actually deal with tracks
    art::FindManyP<recob::Track> trackAssns(pfParticleHandle, event, fTrackProducerLabel);
    
    // Recover the collection of associations between PFParticles and PCAxis objects
    art::FindManyP<recob::PCAxis> pcaxisAssns(pfParticleHandle, event, fPFParticleProducerLabel);
    
    // Recover two sets of associations for linking cosmic tags to the tracks
    art::FindManyP<anab::CosmicTag> pfCosmicAssns(  pfParticleHandle, event, fPFCosmicProducerLabel);
    art::FindManyP<anab::CosmicTag> pfTkCosmicAssns(trackHandle,      event, fPFCosmicProducerLabel);
    art::FindManyP<anab::CosmicTag> cosmicAssns(    trackHandle,      event, fCosmicProducerLabel);
    art::FindManyP<anab::CosmicTag> flashAssns(     trackHandle,      event, fFlashProducerLabel);
    
    // Ok, final step is to initialize the ntuple variables, in particular now that we know how many columns
    this->PrepareEvent(event, particleToHitMap.size());
    
    // We will loop through all of the MCParticles and look for those particles that produce hits
    for (size_t particleIdx = 0; particleIdx < particleHandle->size(); particleIdx++)
    {
        art::Ptr<simb::MCParticle> particle(particleHandle, particleIdx);

        // Recover the track ID, for historical reasons we call it "best"
        int particleTrackID = particle->TrackId();
        
        // Did this mc particle leave hits in the TPC?
        ParticleToHitMap::iterator particleToHitItr = particleToHitMap.find(particleTrackID);
        
        // Let's get the total number of "true" hits that are created by this MCParticle
        // Count number of hits in each view
        size_t nTrueMcHits(0);
        size_t nUniqueTrueMcHits(0);
        
        if (particleToHitItr != particleToHitMap.end())
        {
            nTrueMcHits = particleToHitItr->second.size();
            
            for(const auto& hit : particleToHitItr->second)
            {
                if (hitToParticleMap[hit].size() < 2) nUniqueTrueMcHits++;
            }
        }
        
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
        
        // The following is meant to get the correct offset for drawing the particle trajectory
        // In particular, the cosmic rays will not be correctly placed without this
        double g4Ticks(fTimeService->TPCG4Time2Tick(particle->T())+fDetectorProperties->GetXTicksOffset(0,0,0)-fDetectorProperties->TriggerOffset());
        double xOffset(fDetectorProperties->ConvertTicksToX(g4Ticks, 0, 0, 0));
        
        // Calculate the length of this mc particle inside the fiducial volume.
        TVector3 mcstart;
        TVector3 mcend;
        TVector3 mcstartmom;
        TVector3 mcendmom;
        
        double mcTrackLen = length(*particle, xOffset, mcstart, mcend, mcstartmom, mcendmom);
        
        // Apparently it can happen that we completely miss the TPC
        if (mcTrackLen == 0.) continue;
        
//        if (mcTrackLen > 0.3 || nTrueMcHits > 3)
//        {
//            std::cout << "****>> Number true hits: " << nTrueMcHits << ", range: " << mcTrackLen << std::endl;
//            std::cout << *particle << std::endl;
//        }
        
        // Set the momentum vector to a unit direction vector
        if (mcstartmom.Mag() > 0.) mcstartmom.SetMag(1.);
        if (mcendmom.Mag()   > 0.) mcendmom.SetMag(1.);
        
        std::string process   = particle->Process();
        bool        isPrimary = process == "primary";
        
        // We also need to determine ultimate parentage so we can discern true neutrino decay products from CR
        bool isNeutrino(false);
        int  parentTrackIdx(-1);
        
        try
        {
            art::Ptr<simb::MCTruth> mcTruth = mcTruthAssns.at(particleIdx);
            
//            std::cout << *mcTruth << std::endl;
            
            simb::Origin_t particleOrigin = mcTruth->Origin();
            if (particleOrigin != simb::kCosmicRay) isNeutrino = true;
        }
        catch(...)
        {
            isNeutrino = false;
        }
        
        // See if we can find the mother to this particle in our particle map
        std::map< int, const simb::MCParticle*>::iterator particleItr = particleMap.find(particle->Mother());
        const simb::MCParticle* lastParticle = particle.get();
        
        // While we have a valid mother particle, traverse up the particle tree
        // Note: The sign we are at the top of the tree is that the mother particle index is zero and
        // there will be no entry in the index-particle map for that case.
        while(particleItr != particleMap.end())
        {
            lastParticle = particleItr->second;
            particleItr  = particleMap.find(lastParticle->Mother());
        }
        
        std::string parentProcess = lastParticle->Process();
        
        parentTrackIdx = lastParticle->TrackId();

        // Define the various quantities we want to fill if we find a PFParticle associated to this MCParticle
        const recob::Track* tTrack(0);
        int                 bestCnt(0);
        double              trackLen(0.);
        int                 nTotalClusterHits(0);
        int                 numPFParticles(0);         // PFParticle
        int                 numPFPartTracks(0);        // PFParticle
        int                 bestPFParticleID(-1);      // PFParticle
        int                 bestPFParticleMCTrkID(-1); // PFParticle
        float               trackDistToEdge(9999.);
        float               pfCosmicTagVal(0.);
        float               pfCosmicTypeVal(0.);
        float               cosmicTagVal(0.);
        float               cosmicTypeVal(0.);
        float               flashTagVal(0.);
        TVector3            pcaxisPos(0.,0.,-10.);
        TVector3            pcaEigenVals(0.,0.,0.);
        TVector3            pcaxis1Dir(0.,0.,1.);
        TVector3            pcaxis2Dir(0.,0.,1.);
        TVector3            pcaxis3Dir(0.,0.,1.);
        TVector3            pfCosmicStart(0.,0.,0.);
        TVector3            pfCosmicEnd(0.,0.,0.);
        
        int                 nTotalTrackHits(0);
        int                 bestKTrackCnt(0);
        int                 numKTracks(0);
        int                 bestKTrackID(-1);
        int                 bestKTrackMCTrkId(-1);
        
        // Now check to see that a PFParticle is associated with this MCParticle (trackID)
        TrackIDToPFParticleVecMap::iterator trackIDToPFParticleItr = trackIDToPFParticleMap.find(particleTrackID);
        
        // Proceed here if we have at least a chance to match to PFParticle
        if (trackIDToPFParticleItr != trackIDToPFParticleMap.end())
        {
            // We want to keep track of this...
            numPFParticles = trackIDToPFParticleItr->second.size();
            
            // What happens if a secondary with multiple PFParticles?
            // We want to try to avoid the case where a single hit can be confused
            // between two MCParticles, say crossed tracks or something.
            int topParticleTrackID(particleTrackID);

            if (numPFParticles > 1 && particle->Process() != "primary")
            {
                int parentId = particle->Mother();
                
                if (parentId)
                {
                    std::map< int, const simb::MCParticle*>::iterator parentItr = particleMap.find(parentId);
                    
                    if (parentItr != particleMap.end())
                    {
                        if (parentItr->second->Mother() == 0)
                        {
                            TrackIDToPFParticleVecMap::iterator parentIDToPFParticleItr = trackIDToPFParticleMap.find(parentId);
                            
                            if (parentIDToPFParticleItr != trackIDToPFParticleMap.end()) topParticleTrackID = parentId;
                        }
                    }
                    //else std::cout << *particle << std::endl;
                }
            }
            
            // Ok, it is possible for one MCParticle to have more than one PFParticle (think of a broken cluster) so we
            // need to loop through the possible PFParticles (and fully expect this to normally be only one PFParticle)
            for(const auto& tmpPart : trackIDToPFParticleItr->second)
            {
                // For the given PFParticle, recover the element corresponding to the current MCParticle Track ID
                TrackIDToHit2DMap::iterator tmpPartTrackItr = pfParticleToTrackHitMap[tmpPart].find(topParticleTrackID);

                // In theory this condition is satisfied at least once per loop
                if (tmpPartTrackItr != pfParticleToTrackHitMap[tmpPart].end())
                {
                    // Recover the hit count
                    int trackHitCnt = tmpPartTrackItr->second.size();

                    // Most hits wins!
                    if (trackHitCnt > bestCnt)
                    {
                        bestCnt          = trackHitCnt;
                        //pfBestPurity     = purity;
                        bestPFParticleID = tmpPart.key();
                    }
                }
            }

            // If we found a matching PFParticle then look up the associated information
            if (bestPFParticleID > -1)
            {
                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,bestPFParticleID);
                
                // Recover tracks associated to this PFParticle
                std::vector<art::Ptr<recob::Track>> trackVec = trackAssns.at(pfParticle.key());
                
                numKTracks = trackVec.size();
                
                for(const auto& trackPtr : trackVec)
                {
                    // For the given PFParticle, recover the element corresponding to the current MCParticle Track ID
                    TrackIDToHit2DMap::iterator tmpPartTrackItr = kTrackToTrackHitMap[trackPtr].find(topParticleTrackID);
                    
                    // In theory this condition is satisfied at least once per loop
                    if (tmpPartTrackItr != kTrackToTrackHitMap[trackPtr].end())
                    {
                        // Recover the hit count
                        int trackHitCnt = tmpPartTrackItr->second.size();
                        
                        // Most hits are the winner
                        if (trackHitCnt > bestKTrackCnt)
                        {
                            bestKTrackCnt = trackHitCnt;
                            bestKTrackID  = trackPtr.key();
                        }
                    }
                }

                // Start by looking for an associated fit track(s)
                if (bestKTrackID > -1)
                {
                    art::Ptr<recob::Track> trackPtr(trackHandle,bestKTrackID);
                    
                    tTrack          = trackPtr.get();
                    trackLen        = length(tTrack);
                    nTotalTrackHits = kTrackHitCntMap[trackPtr];
                    
                    TrackIDToHit2DMap::iterator trackIDToHit2DMapItr = kTrackToTrackHitMap[trackPtr].begin();
                    int                             bestKTrackHitCnt = trackIDToHit2DMapItr->second.size();
                    
                    bestKTrackMCTrkId = trackIDToHit2DMapItr->first;
                    
                    while(++trackIDToHit2DMapItr != kTrackToTrackHitMap[trackPtr].end())
                    {
                        if (bestKTrackHitCnt < int(trackIDToHit2DMapItr->second.size()))
                        {
                            bestKTrackHitCnt  = trackIDToHit2DMapItr->second.size();
                            bestKTrackMCTrkId = trackIDToHit2DMapItr->first;
                        }
                    }
                    
                    // While here let's check CR tagging status of this track
                    if (cosmicAssns.isValid())
                    {
                        std::vector<art::Ptr<anab::CosmicTag>> cosmicTagVec = cosmicAssns.at(trackPtr.key());
                        
                        if (!cosmicTagVec.empty())
                        {
                            art::Ptr<anab::CosmicTag> cosmicTag = cosmicTagVec.front();
                        
                            cosmicTagVal  = cosmicTag->CosmicScore();
                            cosmicTypeVal = cosmicTag->CosmicType();
                        }
                    }
                    
                    if (flashAssns.isValid())
                    {
                        std::vector<art::Ptr<anab::CosmicTag>> cosmicTagVec = flashAssns.at(trackPtr.key());
                        
                        if (!cosmicTagVec.empty())
                        {
                            art::Ptr<anab::CosmicTag> cosmicTag = cosmicTagVec.front();
                        
                            flashTagVal  = cosmicTag->CosmicScore();
                        }
                    }
                    
                    if (pfTkCosmicAssns.isValid()) // && !(pfCosmicAssns.isValid() && pfCosmicAssns.size() > 0))
                    {
                        std::vector<art::Ptr<anab::CosmicTag>> pfTkCosmicTagVec = pfTkCosmicAssns.at(trackPtr.key());
                        
                        if (!pfTkCosmicTagVec.empty())
                        {
                            art::Ptr<anab::CosmicTag> cosmicTag = pfTkCosmicTagVec.front();
                        
                            pfCosmicTagVal  = cosmicTag->CosmicScore();
                            pfCosmicTypeVal = cosmicTag->CosmicType();
                        }
                    }
                    
                    double beginDistToEdge = distanceToTPCEdge(tTrack->Vertex());
                    double endDistToEdge   = distanceToTPCEdge(tTrack->End());
                    
                    trackDistToEdge = std::min(beginDistToEdge, endDistToEdge);
                }
                
                // How many tracks associated to this PFParticle?
                if (trackAssns.isValid())
                {
                    std::vector<art::Ptr<recob::Track>> trackVec = trackAssns.at(pfParticle.key());
                    numPFPartTracks = trackVec.size();
                }
                
                // Look for an associated PCAxis
                if (pcaxisAssns.isValid())
                {
                    std::vector<art::Ptr<recob::PCAxis> > pcaxisVec = pcaxisAssns.at(pfParticle.key());
                    
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
                
                // Were there cosmic tag associations with this PFParticle?
                if (pfCosmicAssns.isValid())
                {
                    std::vector<art::Ptr<anab::CosmicTag>> pfCosmicTagVec = pfCosmicAssns.at(pfParticle.key());
                    
                    if (!pfCosmicTagVec.empty())
                    {
                        art::Ptr<anab::CosmicTag> cosmicTag = pfCosmicTagVec.front();
                    
                        pfCosmicTagVal  = cosmicTag->CosmicScore();
                        pfCosmicTypeVal = cosmicTag->CosmicType();
                        
                        pfCosmicStart   = TVector3(cosmicTag->EndPoint1()[0],cosmicTag->EndPoint1()[1],cosmicTag->EndPoint1()[2]);
                        pfCosmicEnd     = TVector3(cosmicTag->EndPoint2()[0],cosmicTag->EndPoint2()[1],cosmicTag->EndPoint2()[2]);
                    }
                }
                
                // Find the MCParticle which contributes the most to this PFParticle
                TrackIDToHit2DMap::iterator trackIDToHit2DMapItr     = pfParticleToTrackHitMap[pfParticle].begin();
                int                             bestPFParticleHitCnt = trackIDToHit2DMapItr->second.size();
                
                bestPFParticleMCTrkID = trackIDToHit2DMapItr->first;
                
                while(++trackIDToHit2DMapItr != pfParticleToTrackHitMap[pfParticle].end())
                {
                    if (bestPFParticleHitCnt < int(trackIDToHit2DMapItr->second.size()))
                    {
                        bestPFParticleHitCnt  = trackIDToHit2DMapItr->second.size();
                        bestPFParticleMCTrkID = trackIDToHit2DMapItr->first;
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
            fPDGCode[fNumMcParticles]             = trackPDGCode;
            fMcPartStartX[fNumMcParticles]        = mcstart.X();
            fMcPartStartY[fNumMcParticles]        = mcstart.Y();
            fMcPartStartZ[fNumMcParticles]        = mcstart.Z();
            fMcPartStartDirX[fNumMcParticles]     = mcstartmom.X();
            fMcPartStartDirY[fNumMcParticles]     = mcstartmom.Y();
            fMcPartStartDirZ[fNumMcParticles]     = mcstartmom.Z();
            fMcPartEndX[fNumMcParticles]          = mcend.X();
            fMcPartEndY[fNumMcParticles]          = mcend.Y();
            fMcPartEndZ[fNumMcParticles]          = mcend.Z();
            fMcPartEndDirX[fNumMcParticles]       = mcendmom.X();
            fMcPartEndDirY[fNumMcParticles]       = mcendmom.Y();
            fMcPartEndDirZ[fNumMcParticles]       = mcendmom.Z();
            fMcPartEne[fNumMcParticles]           = particle->E();
            fMcPartMass[fNumMcParticles]          = particle->Mass();
            fMcPartTrackLen[fNumMcParticles]      = mcTrackLen;
            fMcPartTrackID[fNumMcParticles]       = particleTrackID;
            fMcParentTrackID[fNumMcParticles]     = parentTrackIdx;
            fMcPartNumRecoHits[fNumMcParticles]   = nTrueMcHits;
            fMcPartNumUniqueHits[fNumMcParticles] = nUniqueTrueMcHits;
            fMcPartNumPFParts[fNumMcParticles]    = numPFParticles;
            fMcPartBestPFPart[fNumMcParticles]    = bestPFParticleID;
            fProcessNameVec[fNumMcParticles]      = process;
            fParProcNameVec[fNumMcParticles]      = parentProcess;
            fMcPartPrimary[fNumMcParticles]       = isPrimary;
            fMcNeutrinoDaughter[fNumMcParticles]  = isNeutrino;
            
            fNumMcParticles++;
        }
        
        // Now do PFParticle corresponding to the above
        if (fNumPFParticles < fMaxEntries)
        {
            fPFParticleID[fNumPFParticles]       = bestPFParticleID;
            fBestMcTrackID[fNumPFParticles]      = bestPFParticleMCTrkID;
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
            
            fPFCosmicGeoScore[fNumPFParticles]   = pfCosmicTagVal;
            fPFCosmicGeoTag[fNumPFParticles]     = pfCosmicTypeVal;
            
            fCosmicStartX[fNumPFParticles]       = pfCosmicStart.X();
            fCosmicStartY[fNumPFParticles]       = pfCosmicStart.Y();
            fCosmicStartZ[fNumPFParticles]       = pfCosmicStart.Z();
            fCosmicEndX[fNumPFParticles]         = pfCosmicEnd.X();
            fCosmicEndY[fNumPFParticles]         = pfCosmicEnd.Y();
            fCosmicEndZ[fNumPFParticles]         = pfCosmicEnd.Z();
            
            if (tTrack)
            {
                fTrackID[fNumPFParticles]             = bestKTrackID;
                fNumKTracks[fNumPFParticles]          = numKTracks;
                fNumTrackHitsMatched[fNumPFParticles] = bestKTrackCnt;
                fNumTrackHitsTotal[fNumPFParticles]   = nTotalTrackHits;
                fFitTrackLen[fNumPFParticles]         = trackLen;
                fTrackBestMCTrkId[fNumPFParticles]    = bestKTrackMCTrkId;
                fTrackStartX[fNumPFParticles]         = tTrack->Vertex().X();
                fTrackStartY[fNumPFParticles]         = tTrack->Vertex().Y();
                fTrackStartZ[fNumPFParticles]         = tTrack->Vertex().Z();
                fTrackStartDirX[fNumPFParticles]      = tTrack->VertexDirection().X();
                fTrackStartDirY[fNumPFParticles]      = tTrack->VertexDirection().Y();
                fTrackStartDirZ[fNumPFParticles]      = tTrack->VertexDirection().Z();
                fTrackEndX[fNumPFParticles]           = tTrack->End().X();
                fTrackEndY[fNumPFParticles]           = tTrack->End().Y();
                fTrackEndZ[fNumPFParticles]           = tTrack->End().Z();
                fTrackEndDirX[fNumPFParticles]        = tTrack->EndDirection().X();
                fTrackEndDirY[fNumPFParticles]        = tTrack->EndDirection().Y();
                fTrackEndDirZ[fNumPFParticles]        = tTrack->EndDirection().Z();
                fDistToEdge[fNumPFParticles]          = trackDistToEdge;
                fCosmicGeoScore[fNumPFParticles]      = cosmicTagVal;
                fCosmicGeoTag[fNumPFParticles]        = cosmicTypeVal;
                fFlashScore[fNumPFParticles]          = flashTagVal;
            }
            
            fNumPFParticles++;
        }
        
        // We can't proceed if we have run off the end of our allocated space
        if (fNumMcParticles >= fMaxEntries)
        {
            std::cout << "==> PFParticleMcAna exceeding particle with hits count: " << fNumMcParticles << std::endl;
            break;
        }
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
            
            //int start_tdc = timeService->TPCTick2TDC( hit.StartTick() );
            //int end_tdc   = timeService->TPCTick2TDC( hit.EndTick()   );
            
            int start_tdc    = timeService->TPCTick2TDC( hit.StartTick() );
            int end_tdc      = timeService->TPCTick2TDC( hit.EndTick()   );
            int hitStart_tdc = timeService->TPCTick2TDC( hit.PeakTime() - 3.*hit.SigmaPeakTime() );
            int hitEnd_tdc   = timeService->TPCTick2TDC( hit.PeakTime() + 3.*hit.SigmaPeakTime() );
            
            start_tdc = std::max(start_tdc, hitStart_tdc);
            end_tdc   = std::min(end_tdc,   hitEnd_tdc  );
            
            sim::MCHit startTime;
            sim::MCHit endTime;
            
            startTime.SetTime(start_tdc, 0);
            endTime.SetTime(end_tdc, 0);
            
            std::vector<sim::MCHit>::const_iterator startItr = std::lower_bound(mcHitVec.begin(), mcHitVec.end(), startTime);
            std::vector<sim::MCHit>::const_iterator endItr   = std::upper_bound(startItr,         mcHitVec.end(), endTime);
            
            if (startItr != mcHitVec.end())
            {
                double totCharge = 0.;
//                size_t itrDist   = std::distance(startItr, endItr);
                
                for(auto tmpItr = startItr; tmpItr != endItr; tmpItr++) totCharge += (*tmpItr).Charge();
                
                if (!(totCharge > 0.))
                {
//                    std::cout << "******** zero charge hit? *******" << std::endl;
//                    std::cout << hit << std::endl;
                    continue;
                }
                
                while(startItr != endItr)
                {
                    double fracCharge = (*startItr).Charge() / totCharge;
                    int    trackID    = (*startItr++).PartTrackId();
                    
                    // skip if no real contribution to the charge?
                    if (fracCharge < 0.01)
                    {
//                        std::cout << ">>>> skipping due to low fraction charge contribution: " << hit.SummedADC() << ", " << totCharge << ", " << fracCharge << ", " << itrDist << std::endl;
                        continue;
                    }
                    
                    if (trackID < 0)
                    {
                        trackID = -trackID;
                        fNegTrackIds++;
                    }
                    
//                    if (trackID == 34)
//                    {
//                        std::cout << "**Hit matching found trackID == 34, channel: " << channel << ", hit: " << hit << std::endl;
//                    }
                    
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
    
    bool operator()(const art::Ptr<recob::PFParticle>& left, const art::Ptr<recob::PFParticle>& right)
    {
        size_t numHitsLeft  = fPFPartCntMap.at(left);
        size_t numHitsRight = fPFPartCntMap.at(right);
        
        return numHitsLeft > numHitsRight;
    }
private:
    const PFParticleMcAna::PFParticleHitCntMap& fPFPartCntMap;
};

// Build maps for PFParticles
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
        std::vector<art::Ptr<recob::Cluster> > clusterVec = pfParticleClusterAssns.at(pfParticle.key());
        
        // Number 2D hits this PFParticle
        int numPfParticleHits(0);
        
        // To categorize the PFParticle (associate to MCParticle), we will
        // create and fill an instance of a TrackIDToHitMap.
        // So, for a given PFParticle we will have a map of the track ID's (MCParticles)
        // and the hits that they contributed energy too
        pfParticleToTrackHit2DMap[pfParticle] = TrackIDToHit2DMap();
        TrackIDToHit2DMap& trackIdToHit2DMap  = pfParticleToTrackHit2DMap[pfParticle];
        
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
            pfParticleHitCntMap[pfParticle] += hitVec.size();
            
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
                trackIDToPFParticleVecMap[trackItr.first].push_back(pfParticle);
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
    
// This will be used to sort PFParticles in the map below
class PFParticleMcAna::SortKTrackVec
{
    /**
     * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
     */
public:
    SortKTrackVec(const PFParticleMcAna::KTrackHitCntMap& kTrackCntMap) : fKTrackCntMap(kTrackCntMap) {}
    
    bool operator()(const art::Ptr<recob::Track>& left, const art::Ptr<recob::Track>& right)
    {
        size_t numHitsLeft  = fKTrackCntMap.at(left);
        size_t numHitsRight = fKTrackCntMap.at(right);
        
        return numHitsLeft > numHitsRight;
    }
private:
    const PFParticleMcAna::KTrackHitCntMap& fKTrackCntMap;
};
    
// Build maps for PFParticles
//----------------------------------------------------------------------------
void PFParticleMcAna::MakeKTrackMaps(const art::Event&                              event,
                                     const art::Handle<std::vector<recob::Track> >& trackHandle,
                                     const HitToParticleMap&                        hitToParticleMap,
                                     KTrackToTrackHit2DMap&                         kTrackToTrackHit2DMap,
                                     TrackIDToKTrackVecMap&                         trackIDtoKTrackVecMap,
                                     KTrackHitCntMap&                               kTrackHitCntMap)
{
    // Recover track to hit associations
    art::FindManyP<recob::Hit> kTrackHitAssns(trackHandle, event, fTrackProducerLabel);
    
    // Commence looping over Track
    for(size_t kTrackIdx = 0; kTrackIdx < trackHandle->size(); kTrackIdx++)
    {
        // Get our PFParticle
        art::Ptr<recob::Track> kTrack(trackHandle,kTrackIdx);
        
        // To categorize the fit track (associate to MCParticle), we will
        // create and fill an instance of a TrackIDToHitMap.
        // So, for a given PFParticle we will have a map of the track ID's (MCParticles)
        // and the hits that they contributed energy too
        kTrackToTrackHit2DMap[kTrack]         = TrackIDToHit2DMap();
        TrackIDToHit2DMap& trackIdToHit2DMap  = kTrackToTrackHit2DMap[kTrack];
        
        // Try to recover the hits associated to this track
        std::vector<art::Ptr<recob::Hit> > hitVec = kTrackHitAssns.at(kTrack.key());
        
        // Keep track of this to hopefully save some time later
        kTrackHitCntMap[kTrack] = hitVec.size();
            
        // Something to count MCParticles contributing here
        std::set<int> trackIdCntSet;
        int           nMultiParticleHits(0);
        int           nNoParticleHits(0);
        
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
        
        // Make sure something happened here...
        if (!trackIdToHit2DMap.empty())
        {
            // Now spin through the trackIdToHit2DMap to build the reverse map of above, taking track ID's to PFParticles...
            for (const auto& trackItr : trackIdToHit2DMap)
            {
                trackIDtoKTrackVecMap[trackItr.first].push_back(kTrack);
            }
        }
        else
        {
            mf::LogDebug("PFParticleMcAna") << "***>> No Track to KTrack match made for run " << fRun << ", event: " << fEvent << std::endl;
        }
    } // end of loop over the Track collection
    
    // One last spin through to sort the KTracks in the track ID to PFParticle map
    for (auto& trackItr : trackIDtoKTrackVecMap)
    {
        std::sort(trackItr.second.begin(), trackItr.second.end(), SortKTrackVec(kTrackHitCntMap));
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
    
    // Get fiducial volume boundary.
    double xmin = -0.1;
    double xmax = 2.*fGeometry->DetHalfWidth() + 0.1;
    double ymin = -fGeometry->DetHalfHeight() + 0.1;
    double ymax =  fGeometry->DetHalfHeight() + 0.1;
    double zmin = -0.1;
    double zmax = fGeometry->DetLength() + 0.1;
    
    double readOutWindowSize = fDetectorProperties->ReadOutWindowSize();
    
    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;
    
    // The following for debugging purposes
    int findTrackID(-1);
    
    if (part.TrackId() == findTrackID) std::cout << ">>> length, mcpart: " << part << std::endl;
    
    // Loop over the complete collection of trajectory points
    for(int i = 0; i < n; ++i)
    {
        TVector3 pos = part.Position(i).Vect();
        
        // Make fiducial cuts.
        // There are two sets here:
        // 1) We check the original x,y,z position of the trajectory points and require they be
        //    within the confines of the physical TPC
        // 2) We then check the timing of the presumed hit and make sure it lies within the
        //    readout window for this simulation
        
        if(pos.X() >= xmin &&
           pos.X() <= xmax &&
           pos.Y() >= ymin &&
           pos.Y() <= ymax &&
           pos.Z() >= zmin &&
           pos.Z() <= zmax)
        {
            pos[0] += dx;
            double ticks = fDetectorProperties->ConvertXToTicks(pos[0], 0, 0, 0);
            
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
            
            if (part.TrackId() == findTrackID)
            {
                std::cout << ">>> Track #" << findTrackID << ", pos: " << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ", ticks: " << ticks << ", nearest Y wire: ";
                double worldLoc[] = {pos.X(),pos.Y(),pos.Z()};
                geo::WireID wireID = fGeometry->NearestWireID(worldLoc, 2);
                std::cout << wireID << std::endl;
            }
        }
    }
    
    return result;
}
    
double PFParticleMcAna::distanceToTPCEdge(const TVector3& position) const
{
    double distToEdge(fGeometry->DetLength());
    
    // The coordinate system has x = 0 at anode plane, y = 0 in center and z = 0 at upstream boundary...
    // Get distance to X edge first
    double xDistance = std::max(fGeometry->DetHalfWidth() - fabs(position.X() - fGeometry->DetHalfWidth()), 0.);
    
    distToEdge = std::min(xDistance, distToEdge);
    
    // Get the distance to the Y edge next
    double yDistance = std::max(fGeometry->DetHalfHeight() - fabs(position.Y()), 0.);
    
    distToEdge = std::min(yDistance, distToEdge);
    
    // And now Z
    double zDistance = std::max(0.5*fGeometry->DetLength() - fabs(position.Z() - 0.5*fGeometry->DetLength()), 0.);

    return std::min(zDistance, distToEdge);
}


// This macro has to be defined for this module to be invoked from a
// .fcl file; see PFParticleMcAna.fcl for more information.
DEFINE_ART_MODULE(PFParticleMcAna)

} // namespace PFParticleMcAna

#endif // PFParticleMcAna_module
