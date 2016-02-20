// LArNeutrinoAna_module.cc
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

#ifndef LArNeutrinoAna_module
#define LArNeutrinoAna_module

// LArSoft includes
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/AnalysisBase/CosmicTag.h"
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
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

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

namespace LArNeutrinoAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class LArNeutrinoAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit LArNeutrinoAna(fhicl::ParameterSet const& pset);
    virtual ~LArNeutrinoAna();

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

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;

    // Other variables that will be shared between different methods.
    const geo::GeometryCore* fGeometry;       // pointer to Geometry service
    double                   fElectronsToGeV; // conversion factor
    
    // Define histograms here
    TH1D*     fMuonRange;
    TH1D*     fProtonRange;
    TH1D*     fOtherRange;
    
    TH1D*     fMuonRangeMC;
    TH1D*     fProtonRangeMC;
    TH1D*     fOtherRangeMC;
    
    TH1D*     fMuonDeltaRange;
    TH1D*     fProtonDeltaRange;
    TH1D*     fOtherDeltaRange;
    
    TProfile* fMuonEffVsRange;
    TProfile* fProtonEffVsRange;
    TProfile* fOtherEffVsRange;
    
    TH1D*     fMuonNumHits;
    TH1D*     fProtonNumHits;
    TH1D*     fOtherNumHits;
    
    TH1D*     fMuonNumHitsMC;
    TH1D*     fProtonNumHitsMC;
    TH1D*     fOtherNumHitsMC;
    
    TProfile* fMuonEffVsHits;
    TProfile* fProtonEffVsHits;
    TProfile* fOtherEffVsHits;
    
    TH1D*     fMuonHitEfficiency;
    TH1D*     fMuonHitPurity;
    TH1D*     fMuonHitEfficPurity;
    
    TH1D*     fProtonHitEfficiency;
    TH1D*     fProtonHitPurity;
    TH1D*     fProtonHitEfficPurity;
    
    TH1D*     fOtherHitEfficiency;
    TH1D*     fOtherHitPurity;
    TH1D*     fOtherHitEfficPurity;
    
    TH1D*     fNumMcProngs;
    TH1D*     fNumPrimaryMuons;
    TH1D*     fNumPrimaryProtons;
    TH1D*     fNumTracks;
    TH1D*     fNumMuonTracks;
    TH1D*     fNumProtonTracks;
    
    TH2D*     fMuonTLvsMCL;
    TH2D*     fProtonTLvsMCL;
    TH2D*     fOtherTLvsMCL;
    
    TH1D*     fPdgCodeOther;
    TH1D*     fPdgCodeOtherMiss;
    
    TH1D*     fNoiseHits;

}; // class LArNeutrinoAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
LArNeutrinoAna::LArNeutrinoAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    fGeometry = lar::providerFrom<geo::Geometry>();
    
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
LArNeutrinoAna::~LArNeutrinoAna()
{}
   
//-----------------------------------------------------------------------
void LArNeutrinoAna::beginJob()
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
    fMuonRangeMC          = tfs->make<TH1D>("MuonRangeMC",           "; Range",       100,    0., 0.5  * fGeometry->DetLength());
    fProtonRangeMC        = tfs->make<TH1D>("ProtonRangeMC",         "; Range",       100,    0., 0.1  * fGeometry->DetLength());
    fOtherRangeMC         = tfs->make<TH1D>("OtherRangeMC",          "; Range",       100,    0., 0.05 * fGeometry->DetLength());
    
    fMuonRange            = tfs->make<TH1D>("MuonRange",             "; Range",       100,    0., 0.5  * fGeometry->DetLength());
    fProtonRange          = tfs->make<TH1D>("ProtonRange",           "; Range",       100,    0., 0.1  * fGeometry->DetLength());
    fOtherRange           = tfs->make<TH1D>("OtherRange",            "; Range",       100,    0., 0.05 * fGeometry->DetLength());
    
    fMuonDeltaRange       = tfs->make<TH1D>("MuonDeltaRange",        "; Delta Range", 200, -100., 100.);
    fProtonDeltaRange     = tfs->make<TH1D>("ProtonDeltaRange",      "; Delta Range", 200,  -50.,  50.);
    fOtherDeltaRange      = tfs->make<TH1D>("OtherDeltaRange",       "; Delta Range", 200,  -50.,  50.);
    
    fMuonEffVsRange       = tfs->make<TProfile>("MuonEffVsRange",    ";Range",         40,    0.,  0.5  * fGeometry->DetLength(), 0., 1.1);
    fProtonEffVsRange     = tfs->make<TProfile>("ProtonEffVsRange",  ";Range",         25,    0.,  0.1  * fGeometry->DetLength(), 0., 1.1);
    fOtherEffVsRange      = tfs->make<TProfile>("OtherEffVsRange",   ";Range",         20,    0.,  0.05 * fGeometry->DetLength(), 0., 1.1);
    
    fMuonNumHits          = tfs->make<TH1D>("MuonNumHits",           ";log(# hits)",   20,    0.5,   4.5);
    fProtonNumHits        = tfs->make<TH1D>("ProtonNumHits",         ";log(# hits)",   20,    0.5,   4.5);
    fOtherNumHits         = tfs->make<TH1D>("OtherNumHits",          ";log(# hits)",   20,    0.5,   4.5);
    
    fMuonNumHitsMC        = tfs->make<TH1D>("MuonNumHitsMC",         ";log(# hits)",   20,    0.5,   4.5);
    fProtonNumHitsMC      = tfs->make<TH1D>("ProtonNumHitsMC",       ";log(# hits)",   20,    0.5,   4.5);
    fOtherNumHitsMC       = tfs->make<TH1D>("OtherNumHitsMC",        ";log(# hits)",   20,    0.5,   4.5);
    
    fMuonEffVsHits        = tfs->make<TProfile>("MuonEffVsHits",     ";log(# hits)",   20,    0.5,   4.5, 0., 1.1);
    fProtonEffVsHits      = tfs->make<TProfile>("ProtonEffVsHits",   ";log(# hits)",   20,    0.5,   4.5, 0., 1.1);
    fOtherEffVsHits       = tfs->make<TProfile>("OtherEffVsHits",    ";log(# hits)",   20,    0.5,   4.5, 0., 1.1);
    
    fMuonHitEfficiency    = tfs->make<TH1D>("MuonHitEffic",          ";efficiency",    51,    0.0,    1.02);
    fMuonHitPurity        = tfs->make<TH1D>("MuonHitPurity",         ";purity",        51,    0.0,    1.02);
    fMuonHitEfficPurity   = tfs->make<TH1D>("MuonHitEfficPurity",    ";e*p",           51,    0.0,    1.02);
    
    fProtonHitEfficiency  = tfs->make<TH1D>("ProtonHitEffic",        ";efficiency",    51,    0.0,    1.02);
    fProtonHitPurity      = tfs->make<TH1D>("ProtonHitPurity",       ";purity",        51,    0.0,    1.02);
    fProtonHitEfficPurity = tfs->make<TH1D>("ProtonHitEfficPurity",  ";e*p",           51,    0.0,    1.02);
    
    fOtherHitEfficiency   = tfs->make<TH1D>("OtherHitEffic",         ";efficiency",    51,    0.0,    1.02);
    fOtherHitPurity       = tfs->make<TH1D>("OtherHitPurity",        ";purity",        51,    0.0,    1.02);
    fOtherHitEfficPurity  = tfs->make<TH1D>("OtherHitEfficPurity",   ";e*p",           51,    0.0,    1.02);
    
    fMuonTLvsMCL          = tfs->make<TH2D>("MuonTLvsMCL",           "Reco Length vs. MC Truth Length", 100, 0., 0.5 * detectorLength, 100, 0., 0.5 * detectorLength);
    fProtonTLvsMCL        = tfs->make<TH2D>("ProtonTLvsMCL",         "Reco Length vs. MC Truth Length", 100, 0., 0.2 * detectorLength, 100, 0., 0.2 * detectorLength);
    fOtherTLvsMCL         = tfs->make<TH2D>("OtherTLvsMCL",          "Reco Length vs. MC Truth Length", 100, 0., 0.1 * detectorLength, 100, 0., 0.1 * detectorLength);
    
    fNumMcProngs          = tfs->make<TH1D>("NumMcProngs",           ";# prongs",     10, 0., 10.);
    fNumPrimaryMuons      = tfs->make<TH1D>("NumPrimMuons",          ";# muons",      10, 0., 10.);
    fNumPrimaryProtons    = tfs->make<TH1D>("NumPrimProtons",        ";# protons",    10, 0., 10.);
    fNumTracks            = tfs->make<TH1D>("NumTracks",             ";# tracks",     10, 0., 10.);
    fNumMuonTracks        = tfs->make<TH1D>("NumMuonTracks",         ";# tracks",     10, 0., 10.);
    fNumProtonTracks      = tfs->make<TH1D>("NumProtonTracks",       ";# tracks",     10, 0., 10.);
    
    fPdgCodeOther         = tfs->make<TH1D>("PdgCodeOther",          ";pdg code",    400, -200., 200.);
    fPdgCodeOtherMiss     = tfs->make<TH1D>("PdgCodeOtherM",         ";pdg code",    400, -200., 200.);
    
    fNoiseHits            = tfs->make<TH1D>("NoiseHits",             ";# hits",      100, 0., 100.);
}
   
//-----------------------------------------------------------------------
void LArNeutrinoAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void LArNeutrinoAna::reconfigure(fhicl::ParameterSet const& p)
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
    return;
}

//-----------------------------------------------------------------------
void LArNeutrinoAna::analyze(const art::Event& event)
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

    //we're gonna probably need the time service to convert hit times to TDCs
    const auto* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
      
    // our ultimate goal here
    HitToParticleMap hitToParticleMap;
    ParticleToHitMap particleToHitMap;
    
    // Keep track of the number of noise hits
    int nNoiseHits(0);
    
    // Ok, so this loop obviously takes the MC information and builds two maps
    // 1) a map from a Hit2D object to the track ID's that made it
    // 2) the reverse map, going from track ID to Hit2D object
    for (unsigned int iHit = 0, iHitEnd = hitHandle->size(); iHit < iHitEnd; ++iHit)
    {
        art::Ptr<recob::Hit> hit(hitHandle, iHit);
        
        const geo::WireID& wireId = hit->WireID();
        
        unsigned int channel = fGeometry->PlaneWireToChannel(wireId.Plane, wireId.Wire, wireId.TPC, wireId.Cryostat);
        
        const std::vector<sim::MCHit>& mcHitVec = mcHitCollectionVec.at(channel);
        
        int start_tdc = timeService->TPCTick2TDC( hit->StartTick() );
        int end_tdc   = timeService->TPCTick2TDC( hit->EndTick()   );
        
        sim::MCHit startTime;
        sim::MCHit endTime;
        
        startTime.SetTime(start_tdc, 0);
        endTime.SetTime(end_tdc, 0);
        
        std::vector<sim::MCHit>::const_iterator startItr = std::lower_bound(mcHitVec.begin(), mcHitVec.end(), startTime);
        std::vector<sim::MCHit>::const_iterator endItr   = std::upper_bound(startItr,         mcHitVec.end(), endTime);
        
        while(startItr != endItr)
        {
            int trackID = (*startItr++).PartTrackId();
            
            hitToParticleMap[hit.get()].insert(trackID);
            particleToHitMap[trackID].insert(hit.get());
        }
        
//        else nNoiseHits++;
    }
    
    // Record for posteriety
    fNoiseHits->Fill(nNoiseHits);

    // This loop will build the map between track ID and the MCParticle related to it
    for ( auto const& particle : (*particleHandle) )
    {
        // For the methods you can call to get particle information,
        // see ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
        int trackID = particle.TrackId();

        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[trackID] = &particle;
    } // loop over all particles in the event.
    
    // Ok, at this point we're ready to recover Tracks and the assorted necessities for the next step
    // For sure we need a boatload of stuff here...
    // Start with some useful typdefs
    typedef std::set<art::Ptr<recob::Hit> >                   Hit2DSet;              // Need to count only unique hits
    typedef std::map< int, Hit2DSet >                         TrackIDToHit2DMap;
    typedef std::vector<const recob::Track*>                  TrackVec;
    typedef std::map< int, TrackVec >                         TrackIDToTrackVecMap;
    typedef std::map<const recob::Track*, TrackIDToHit2DMap > TrackToTrackHit2DMap;
    
    // Now define the maps relating pfparticles to tracks
    TrackIDToTrackVecMap trackIDToTrackMap;
    TrackToTrackHit2DMap trackToTrackHitMap;
    
    // Something to keep track of number of hits associated to a cluster
    std::map<const recob::Track*, int> trackHitCntMap;
    
    // Recover the PFParticles, the main products for our next major loop
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    event.getByLabel(fPFParticleProducerLabel, pfParticleHandle);
    
    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    event.getByLabel(fPFParticleProducerLabel, clusterHandle);
    
    // Now retrieve a handle to the fit tracks
    art::Handle<std::vector<recob::Track> > trackHandle;
    event.getByLabel(fTrackProducerLabel, trackHandle);

    // If we have reconstructed tracks then we fill out the tables relating to hits
    // The valididy ot a trackHandle implies the validity of the PFParticle and Cluster handles
    if (trackHandle.isValid())
    {
        // Recover track hit associations
        art::FindManyP<recob::Hit> trackHitAssns(trackHandle, event, fTrackProducerLabel);
    
        // The goal of this loop is to recover the list of 2D hits associated with a given Track
        // (through associations with 2D clusters) and develop some maps which allow us to associate
        // the Track to a MC particle.
        for (const auto& track : (*trackHandle))
        {
            // Recover the hits associated to this track
            std::vector<art::Ptr<recob::Hit> > hits = trackHitAssns.at(track.ID());
        
            // Keep track of this to hopefully save some time later
            trackHitCntMap[&track] = hits.size();

            // To categorize the PFParticle (associate to MCParticle), we will
            // create and fill an instance of a TrackIDToHitMap.
            // So, for a given PFParticle we will have a map of the track ID's (MCParticles)
            // and the hits that they contributed energy too
            trackToTrackHitMap[&track] = TrackIDToHit2DMap();
            TrackIDToHit2DMap& trackIdToHit2DMap = trackToTrackHitMap[&track];
        
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
                    trackIDToTrackMap[bestTrackId].push_back(&track);
                }
            }
        } // end of loop over the PFParticle collection
    }
    
    // Always a handy thing to have hidden in your code:
//    const double radToDegrees = 180. / 3.14159265;
    auto const* larProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // One last task worth pursuing is to see if we can pick out the neutrino interaction and do some categorization of it.
//    const simb::MCTruth* theAnswerIsOutThere(0);
    std::vector<const simb::MCParticle*> mcPartVec;
    std::vector<const simb::MCParticle*> primaryMuonVec;
    std::vector<const simb::MCParticle*> primaryProtonVec;
    std::vector<const recob::Track*>     trackVec;
    std::vector<const recob::Track*>     muonTrackVec;
    std::vector<const recob::Track*>     protonTrackVec;
    
    for (size_t particleIdx = 0; particleIdx < particleHandle->size(); particleIdx++)
    {
        art::Ptr<simb::MCParticle> particle(particleHandle, particleIdx);
        
        // See if we can quickly eliminate junk
        if (particle->NumberTrajectoryPoints() < 2) continue;
        
        bool isNeutrino(false);
        
        try
        {
            art::Ptr<simb::MCTruth> mcTruth = mcTruthAssns.at(particleIdx);
        
            if (mcTruth->Origin() != simb::kBeamNeutrino) continue;
            
            isNeutrino = true;
            
//            theAnswerIsOutThere = mcTruth.get();
        }
        catch(...)
        {
            isNeutrino = false;
        }
        
        if (isNeutrino)
        {
            // For the methods you can call to get particle information,
            // see ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
            int bestTrackID = particle->TrackId();
            
            // Recover the particle's identity
           int trackPDGCode = particle->PdgCode();
            
            // Did this mc particle leave hits in the TPC?
            ParticleToHitMap::iterator particleToHitItr = particleToHitMap.find(bestTrackID);
            
            // No hits no work
            if (particleToHitItr == particleToHitMap.end()) continue;
            
            // Can we check the range of this particle?
            if(particle->E() < 0.001*particle->Mass() + 0.05) continue;
            
            // Let's get the total number of "true" hits that are created by this MCParticle
            int nTrueMcHits = particleToHitItr->second.size();
            
            // Count number of hits in each view
            int           nHitsPerView[3] = {0,0,0};
            std::set<int> uniqueWiresPerView[3] = {std::set<int>(),std::set<int>(),std::set<int>()};
            float         maxTimePerView[3] = {0., 0., 0.};
            float         minTimePerView[3] = {10000., 10000., 10000.};
            
            for (const auto& hit : particleToHitItr->second)
            {
                int view = hit->View();
                
                nHitsPerView[view]++;
                uniqueWiresPerView[view].insert(hit->WireID().Wire);
                maxTimePerView[view] = std::max(maxTimePerView[view], hit->PeakTime());
                minTimePerView[view] = std::min(minTimePerView[view], hit->PeakTime());
            }
            
            // We need to make sure there was some activity in the TPC
            if (uniqueWiresPerView[0].size() < 4 || uniqueWiresPerView[1].size() < 4 || uniqueWiresPerView[2].size() < 4) continue;
            
//            std::cout << "--- PDG code: " << trackPDGCode << std::endl;
            
//            for(int idx = 0; idx < 3; idx++)
//            {
//                std::cout << "Hits: " << nHitsPerView[idx] << ", wires: " << uniqueWiresPerView[idx].size() << ", max Time: " << maxTimePerView[idx]
//                << ", min Time: " << minTimePerView[idx] << std::endl;
//            }
            
            // Calculate the x offset due to nonzero mc particle time.
            double mctime = particle->T();                                  // nsec
            double mcdx   = mctime * 1.e-3 * larProperties->DriftVelocity();  // cm
                
            // Calculate the length of this mc particle inside the fiducial volume.
            TVector3 mcstart;
            TVector3 mcend;
            TVector3 mcstartmom;
            TVector3 mcendmom;
            
            double mcTrackLen = length(*particle, mcdx, mcstart, mcend, mcstartmom, mcendmom);
            
            // Require particle travel at least three wires in Y-Z plane
            if (mcTrackLen < 0.91) continue;
            
            mcPartVec.push_back(particle.get());
            
            double logNumMcHits(std::log10(nTrueMcHits));
            
            bool isPrimaryMuon   = fabs(trackPDGCode) ==   13 && particle->Process() == "primary";
            bool isPrimaryProton = fabs(trackPDGCode) == 2212 && particle->Process() == "primary";
            
            if      (isPrimaryMuon)   primaryMuonVec.push_back(particle.get());
            else if (isPrimaryProton) primaryProtonVec.push_back(particle.get());
            
            mf::LogDebug("LArNeutrinoAna") << "***>> Found Neutrino product: " << *particle << std::endl;

            // Now check to see that a track is associated with this MCParticle (trackID)
            TrackIDToTrackVecMap::iterator trackIDToTrackItr = trackIDToTrackMap.find(bestTrackID);
            
            // If no track then skip the rest
            if (trackIDToTrackItr == trackIDToTrackMap.end())
            {
                // Of course, we must record the results first, and what we do depends on the particle...
                if (isPrimaryMuon)
                {
                    mcTrackLen = std::min(mcTrackLen, 0.5*fGeometry->DetLength()-1.);

                    fMuonRangeMC->Fill(mcTrackLen, 1.);
                    fMuonRange->Fill(0., 1.);
                    fMuonEffVsRange->Fill(mcTrackLen, 0.);
                    fMuonNumHits->Fill(0.1, 1.);
                    fMuonNumHitsMC->Fill(logNumMcHits, 1.);
                    fMuonEffVsHits->Fill(logNumMcHits, 0.);
//                    fMuonHitEfficiency->Fill(0., 1.);
//                    fMuonHitPurity->Fill(0., 1.);
//                    fMuonHitEfficPurity->Fill(0., 1.);
                    fMuonTLvsMCL->Fill(mcTrackLen, 0., 1.);
                    
                    mf::LogInfo("LArNeutrinoAna") << "***>> Primary muon missing, event # " << fEvent << ", mcTrackLen: " << mcTrackLen << std::endl;
                }
                else if (isPrimaryProton)
                {
                    mcTrackLen = std::min(mcTrackLen, 0.1*fGeometry->DetLength()-1.);
                    
                    fProtonRangeMC->Fill(mcTrackLen, 1.);
                    fProtonRange->Fill(0., 1.);
                    fProtonEffVsRange->Fill(mcTrackLen,  0.);
                    fProtonNumHits->Fill(0.1, 1.);
                    fProtonNumHitsMC->Fill(logNumMcHits, 1.);
                    fProtonEffVsHits->Fill(logNumMcHits, 0.);
//                    fProtonHitEfficiency->Fill(0., 1.);
//                    fProtonHitPurity->Fill(0., 1.);
//                    fProtonHitEfficPurity->Fill(0., 1.);
                    fProtonTLvsMCL->Fill(mcTrackLen, 0., 1.);
                    mf::LogInfo("LArNeutrinoAna") << "***>> Primary proton missing, event # " << fEvent << ", mcTrackLen: " << mcTrackLen << std::endl;
                }
                else
                {
                    mcTrackLen = std::min(mcTrackLen, 0.05*fGeometry->DetLength()-1.);
                    
                    // Should add pdgcode plot to see what things are...
                    fOtherRangeMC->Fill(mcTrackLen, 1.);
                    fOtherRange->Fill(0., 1.);
                    fOtherEffVsRange->Fill(mcTrackLen,  0.);
                    fOtherNumHits->Fill(0.1, 1.);
                    fOtherNumHitsMC->Fill(logNumMcHits, 1.);
                    fOtherEffVsHits->Fill(logNumMcHits, 0.);
//                    fOtherHitEfficiency->Fill(0., 1.);
//                    fOtherHitPurity->Fill(0., 1.);
//                    fOtherHitEfficPurity->Fill(0., 1.);
                    fOtherTLvsMCL->Fill(mcTrackLen, 0., 1.);
                    fPdgCodeOtherMiss->Fill(std::max(-199.5,std::min(199.5,double(trackPDGCode))), 1.);
                }
                
                continue;
            }
            
            // Here our goal is to count the number of hits associated to a Track by the "best" MCParticle
            // (which is the one which produced the most number of hits the Track has grouped together)
            const recob::Track* track(0);
            int                 bestCnt(0);
            int                 nTracks(0);
            
            // Loop over the Tracks associated to this track ID (MCParticle)
            for(const auto& tmpPart : trackIDToTrackItr->second)
            {
                // This is to get the list of hits this PFParticle has broken out by track id
                TrackIDToHit2DMap::iterator tmpPartTrackItr = trackToTrackHitMap[tmpPart].find(bestTrackID);
            
                if (tmpPartTrackItr != trackToTrackHitMap[tmpPart].end())
                {
                    int trackHitCnt = tmpPartTrackItr->second.size();
                
                    if (trackHitCnt > bestCnt)
                    {
                        bestCnt = trackHitCnt;
                        track   = tmpPart;
                    }
                
                    if (trackHitCnt > int(nTrueMcHits / 10)) nTracks++;
                }
            }
            
            // I don't think this can happen...
            if (!track) continue;
            
            trackVec.push_back(track);
            
            // Now we can determine:
            // 1) The total number of hits in this cluster
            // 2) The number of hits belonging to the matched mc track
            // 3) The number of hits on that track
            // 4) The number of hits not belonging to this track which are on the cluster (the impurities)
            size_t nTotalClusterHits = trackHitCntMap[track];
            size_t nMatchedHits      = bestCnt;
//            size_t nWrongClusterHits = nTotalClusterHits - nMatchedHits;
            
            // MicroBooNE definitions of efficiency and purity:
            // efficiency E = number of true hits in the cluster / number of true hits from MC particle
            // purity P     = number of true hits in the cluster / all hits in the cluster
            double hitEfficiency = double(nMatchedHits)      / double(nTrueMcHits);
//            double wrongEff      = double(nWrongClusterHits) / double(nTrueMcHits);
            double hitPurity     = double(nMatchedHits)      / double(nTotalClusterHits);
            
            double trackLen   = length(track);
            double logNumHits = std::log10(nTotalClusterHits);
            
            if (isPrimaryMuon)
            {
                double deltaRange = std::min(99.9, std::max(mcTrackLen - trackLen, -99.9));
                
                mcTrackLen = std::min(mcTrackLen, 0.5*fGeometry->DetLength()-1.);
                trackLen   = std::min(trackLen,   0.5*fGeometry->DetLength()-1.);
                
                fMuonRangeMC->Fill(mcTrackLen, 1.);
                fMuonRange->Fill(trackLen, 1.);
                fMuonDeltaRange->Fill(deltaRange,  1.);
                fMuonEffVsRange->Fill(mcTrackLen,  1.);
                fMuonNumHits->Fill(logNumHits, 1.);
                fMuonNumHitsMC->Fill(logNumMcHits, 1.);
                fMuonEffVsHits->Fill(logNumMcHits, 1.);
                fMuonHitEfficiency->Fill(hitEfficiency, 1.);
                fMuonHitPurity->Fill(hitPurity, 1.);
                fMuonHitEfficPurity->Fill(hitEfficiency*hitPurity, 1.);
                fMuonTLvsMCL->Fill(mcTrackLen, trackLen, 1.);
                muonTrackVec.push_back(track);
            }
            else if (isPrimaryProton)
            {
                double deltaRange = std::min(49.9, std::max(mcTrackLen - trackLen, -49.9));
                
                mcTrackLen = std::min(mcTrackLen, 0.1*fGeometry->DetLength()-1.);
                trackLen   = std::min(trackLen,   0.1*fGeometry->DetLength()-1.);
                
                fProtonRangeMC->Fill(mcTrackLen, 1.);
                fProtonRange->Fill(trackLen, 1.);
                fProtonDeltaRange->Fill(deltaRange,  1.);
                fProtonEffVsRange->Fill(mcTrackLen,  1.);
                fProtonNumHits->Fill(logNumHits, 1.);
                fProtonNumHitsMC->Fill(logNumMcHits, 1.);
                fProtonEffVsHits->Fill(logNumMcHits, 1.);
                fProtonHitEfficiency->Fill(hitEfficiency, 1.);
                fProtonHitPurity->Fill(hitPurity, 1.);
                fProtonHitEfficPurity->Fill(hitEfficiency*hitPurity, 1.);
                fProtonTLvsMCL->Fill(mcTrackLen, trackLen, 1.);
                protonTrackVec.push_back(track);
            }
            else
            {
                double deltaRange = std::min(49.9, std::max(mcTrackLen - trackLen, -49.9));
                
                mcTrackLen = std::min(mcTrackLen, 0.05*fGeometry->DetLength()-1.);
                trackLen   = std::min(trackLen,   0.05*fGeometry->DetLength()-1.);
                
                fOtherRangeMC->Fill(mcTrackLen, 1.);
                fOtherRange->Fill(trackLen, 1.);
                fOtherDeltaRange->Fill(deltaRange,  1.);
                fOtherEffVsRange->Fill(mcTrackLen,  1.);
                fOtherNumHits->Fill(logNumHits, 1.);
                fOtherNumHitsMC->Fill(logNumMcHits, 1.);
                fOtherEffVsHits->Fill(logNumMcHits, 1.);
                fOtherHitEfficiency->Fill(hitEfficiency, 1.);
                fOtherHitPurity->Fill(hitPurity, 1.);
                fOtherHitEfficPurity->Fill(hitEfficiency*hitPurity, 1.);
                fOtherTLvsMCL->Fill(mcTrackLen, trackLen, 1.);
                fPdgCodeOther->Fill(std::max(-199.5,std::min(199.5,double(trackPDGCode))), 1.);
            }
        }
    }
    
    int nMcProngs       = mcPartVec.size();
    int nPrimaryMuons   = primaryMuonVec.size();
    int nPrimaryProtons = primaryProtonVec.size();
    int nTracks         = trackVec.size();
    int nMuonTracks     = muonTrackVec.size();
    int nProtonTracks   = protonTrackVec.size();
    
    fNumMcProngs->Fill(nMcProngs);
    fNumPrimaryMuons->Fill(nPrimaryMuons);
    fNumPrimaryProtons->Fill(nPrimaryProtons);
    fNumTracks->Fill(nTracks);
    fNumMuonTracks->Fill(nMuonTracks);
    fNumProtonTracks->Fill(nProtonTracks);

    return;
}
    
// Length of reconstructed track.
//----------------------------------------------------------------------------
double LArNeutrinoAna::length(const recob::Track* track)
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
double LArNeutrinoAna::length(const simb::MCParticle& part, double dx,
                              TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                              unsigned int tpc, unsigned int cstat)
{
    // Get services.
    
    const auto* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    const auto* geom    = lar::providerFrom<geo::Geometry>();
    
    // Get fiducial volume boundary.
    
    double xmin = 0.;
    double xmax = 2.*geom->DetHalfWidth();
    double ymin = -geom->DetHalfHeight();
    double ymax = geom->DetHalfHeight();
    double zmin = 0.;
    double zmax = geom->DetLength();
    
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
            if(ticks >= 0. && ticks < detprop->ReadOutWindowSize())
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
// .fcl file; see LArNeutrinoAna.fcl for more information.
DEFINE_ART_MODULE(LArNeutrinoAna)

} // namespace LArNeutrinoAna

#endif // LArNeutrinoAna_module
