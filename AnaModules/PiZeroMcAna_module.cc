// PiZeroMcAna_module.cc
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

#ifndef PiZeroMcAna_Module
#define PiZeroMcAna_Module

// LArSoft includes
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/PCAxis.h"
#include "lardata/AnalysisBase/CosmicTag.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "lardata/MCBase/MCHitCollection.h"

//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
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

namespace PiZeroMcAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class PiZeroMcAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit PiZeroMcAna(fhicl::ParameterSet const& pset);
    virtual ~PiZeroMcAna();

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

    // The variables that will go into the n-tuple.
    int      fEvent;
    int      fRun;
    int      fSubRun;
    
    // Some general information on the number of hits
    Int_t    fNumHits;
    Int_t    fNoiseHits;
    Int_t    fNegTrackIds;

    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    // Note: old-style C++ arrays are considered obsolete. However,
    // to create simple n-tuples, we still need to use them. 
    double fStartXYZT[4];
    double fEndXYZT[4];
    double fStartPE[4];
    double fEndPE[4];

    // Other variables that will be shared between different methods.
    const geo::GeometryCore*            fGeometry;             // pointer to Geometry service
    const detinfo::DetectorProperties*  fDetectorProperties;   ///< Detector properties service
    double                              fElectronsToGeV; // conversion factor

}; // class PiZeroMcAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
PiZeroMcAna::PiZeroMcAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
PiZeroMcAna::~PiZeroMcAna()
{}
   
//-----------------------------------------------------------------------
void PiZeroMcAna::beginJob()
{
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    double detectorLength = fGeometry->DetLength(); 

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;
    
    // Temporary
    fNegTrackIds = 0;
    fNumHits     = 0;
    fNoiseHits   = 0;
  
    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes. 

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fPDGCodeHist        = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist       = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist    = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detectorLength);

}
   
//-----------------------------------------------------------------------
void PiZeroMcAna::beginRun(const art::Run& /*run*/)
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void PiZeroMcAna::reconfigure(fhicl::ParameterSet const& p)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fSimulationProducerLabel    = p.get< std::string >("SimulationLabel");
    fMcHitCollectionModuleLabel = p.get< std::string >("MCHitFinderLabel");
    fPFParticleProducerLabel    = p.get< std::string >("PFParticleLabel");
    fHitProducerLabel           = p.get< std::string >("HitLabel");
    fClusterProducerLabel       = p.get< std::string >("ClusterProducerLabel");
    fTrackProducerLabel         = p.get< std::string >("TrackProducerLabel");
    
    fMcHitCollectionModuleLabel = "mchitfinder";
    
    return;
}

//-----------------------------------------------------------------------
void PiZeroMcAna::analyze(const art::Event& event)
{
    // The first step is to attempt to recover the collection of MCHits that
    // we will need for doing our PFParticle to MC matching
    art::Handle< std::vector<sim::MCHitCollection> > mcHitCollectionHandle;
    event.getByLabel(fMcHitCollectionModuleLabel, mcHitCollectionHandle);
    
    if (!mcHitCollectionHandle.isValid())
    {
        mf::LogDebug("PFParticleMcAna") << "===>> NO McHitColllection found for run: " << fRun << ", event: " << fEvent << std::endl;
//        fAnaTree->Fill();
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
//        fAnaTree->Fill();
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
    // While we are at it, identify the two primary gammas and get some parameters
    for ( auto const& particle : (*particleHandle) )
    {
        int trackID = particle.TrackId();
        
        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[trackID] = &particle;
        
        // Recover particle id
        int pdgCode  = particle.PdgCode();
        int motherID = particle.Mother();
        
        if (pdgCode == 22 && motherID == 1)
        {
            TVector3 startPos(particle.Vx(),particle.Vy(),particle.Vz());
            TVector3 endPos(particle.EndX(),particle.EndY(),particle.EndZ());
            double   dist = (endPos-startPos).Mag();
            
            std::cout << "***** Gamma, dist to conversion: " << dist << std::endl;
            std::cout << particle << std::endl;
        }
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
    
    // Recover the collection of associations between PFParticles and clusters, this will
    // be the mechanism by which we actually deal with clusters
    art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, event, fPFParticleProducerLabel);
    
    // Recover the collection of associations between PFParticles and tracks, this will
    // be the mechanism by which we actually deal with tracks
    art::FindManyP<recob::Track> trackAssns(pfParticleHandle, event, fTrackProducerLabel);
    
    // Recover the collection of associations between PFParticles and PCAxis objects
    art::FindManyP<recob::PCAxis> pcaxisAssns(pfParticleHandle, event, fPFParticleProducerLabel);
    
    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle<std::vector<recob::Cluster> > clusterHandle;
    event.getByLabel(fPFParticleProducerLabel, clusterHandle);
    
    // Ok, now we can get the hit associations
    art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, event, fPFParticleProducerLabel);
    
    // PiZero analysis:
    // Loop through PFParticles and match to gamma
    for(size_t pfPartIdx = 0; pfPartIdx < pfParticleHandle->size(); pfPartIdx++)
    {
        // Get our PFParticle
        art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfPartIdx);

        // We want to count the total number of hits associated to this PFParticle
        int nTotHits2D(0);
        
        // Recover the associations to clusters so we can do some hit counting
        std::vector<art::Ptr<recob::Cluster> > clusterVec = clusterAssns.at(pfParticle.key());
        
        // Nominally, one believes there will be 3 2D clusters so we need to loop through the list
        for(const auto& cluster : clusterVec)
        {
            // Recover the 2D hits associated to a given cluster
            std::vector<art::Ptr<recob::Hit> > hitVec = clusterHitAssns.at(cluster->ID());
            
            // Count em
            nTotHits2D += hitVec.size();
        }
        
        // Recover MCParticle track ID associated to this PFParticle
        const TrackIDToHit2DMap& trackIdToHitMap = pfParticleToTrackHitMap[pfParticle];
        
        if (trackIdToHitMap.empty()) continue;
        
        int    partTrackId = trackIdToHitMap.begin()->first;
        size_t hitCount    = trackIdToHitMap.begin()->second.size();
        
        if (trackIdToHitMap.size() > 1)
        {
            for(const auto& mapItr : trackIdToHitMap)
            {
                if (mapItr.second.size() > hitCount)
                {
                    partTrackId = mapItr.first;
                    hitCount    = mapItr.second.size();
                }
            }
        }
        
        // Chase this id back up to find the gamma from the pi-zero
        const simb::MCParticle* mcMatched  = particleMap[partTrackId];
        const simb::MCParticle* mcParticle = mcMatched;
        
        while(!(mcParticle->PdgCode() == 22 && mcParticle->Mother() == 1))
        {
            int trackId = mcParticle->Mother();
            
            mcParticle = particleMap[trackId];
        }
        
        // Recover the PCA axes associated to this PFParticle
        std::vector<art::Ptr<recob::PCAxis>> pcAxisVec = pcaxisAssns.at(pfParticle.key());
        
        // empty?
        if (pcAxisVec.empty()) continue;
        
        // The order of axes in the returned association vector is arbitrary... the "first" axis is
        // better and we can divine that by looking at the axis id's (the best will have been made first)
        if (pcAxisVec.size() > 1 && pcAxisVec.front()->getID() > pcAxisVec.back()->getID()) std::reverse(pcAxisVec.begin(), pcAxisVec.end());
        
        const art::Ptr<recob::PCAxis>& pca = pcAxisVec.front();
        
        // We need the mean position
//        const double*              avePosition        = pca->getAvePosition();
//        double                     primaryEigenValue  = pca->getEigenValues()[0];
        const std::vector<double>& primaryEigenVector = pca->getEigenVectors()[0];
        
        TVector3 pcAxisDir(primaryEigenVector[0],primaryEigenVector[1],primaryEigenVector[2]);
        TVector3 gammaDir(mcParticle->Px(),mcParticle->Py(),mcParticle->Pz());
        
        gammaDir.SetMag(1.);
        
        double cosTheta = std::max(-1.,std::min(1.,gammaDir.Dot(pcAxisDir)));
        
        if (cosTheta < 0.) cosTheta = -cosTheta;
        
        double angleToGamma = 180. * std::acos(cosTheta) / 3.14159;
        double energyFrac   = mcMatched->E() / mcParticle->E();
        
        std::cout << ">>> PFParticle Idx: " << pfPartIdx << ", tot hits: " << nTotHits2D << ", matched hits: " << hitCount << ", cos(theta): " << cosTheta << ", angle: " << angleToGamma << std::endl;
        std::cout << "    Matched P: " << mcMatched->P() << ", E: " << mcMatched->E() << ", E frac: " << energyFrac << ", gamma ID: " << mcParticle->TrackId() << std::endl;

    }

    return;
}
    
// This builds out the hit to particle maps
//----------------------------------------------------------------------------
void PiZeroMcAna::MakeHitParticleMaps(const std::vector<sim::MCHitCollection>& mcHitCollectionVec,
                                      const std::vector<recob::Hit>&           recoHitVec,
                                      HitToParticleMap&                        hitToParticleMap,
                                      ParticleToHitMap&                        particleToHitMap)
{
    //we're gonna probably need the time service to convert hit times to TDCs
    const auto& timeService = lar::providerFrom<detinfo::DetectorClocksService>();
    
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
                    if (fracCharge < 0.01) continue;
                    
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
class PiZeroMcAna::SortPFParticleVec
{
    /**
     * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
     */
public:
    SortPFParticleVec(const PiZeroMcAna::PFParticleHitCntMap& pfPartCntMap) : fPFPartCntMap(pfPartCntMap) {}
    
    bool operator()(const art::Ptr<recob::PFParticle>& left, const art::Ptr<recob::PFParticle>& right)
    {
        size_t numHitsLeft  = fPFPartCntMap.at(left);
        size_t numHitsRight = fPFPartCntMap.at(right);
        
        return numHitsLeft > numHitsRight;
    }
private:
    const PiZeroMcAna::PFParticleHitCntMap& fPFPartCntMap;
};

// Build maps for PFParticles
//----------------------------------------------------------------------------
void PiZeroMcAna::MakePFParticleMaps(const art::Event&                                   event,
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
            mf::LogDebug("PiZeroMcAna") << "***>> No PFParticle to MCParticle match made for run " << fRun << ", event: " << fEvent << std::endl;
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
class PiZeroMcAna::SortKTrackVec
{
    /**
     * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
     */
public:
    SortKTrackVec(const PiZeroMcAna::KTrackHitCntMap& kTrackCntMap) : fKTrackCntMap(kTrackCntMap) {}
    
    bool operator()(const art::Ptr<recob::Track>& left, const art::Ptr<recob::Track>& right)
    {
        size_t numHitsLeft  = fKTrackCntMap.at(left);
        size_t numHitsRight = fKTrackCntMap.at(right);
        
        return numHitsLeft > numHitsRight;
    }
private:
    const PiZeroMcAna::KTrackHitCntMap& fKTrackCntMap;
};

// Build maps for PFParticles
//----------------------------------------------------------------------------
void PiZeroMcAna::MakeKTrackMaps(const art::Event&                              event,
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
            mf::LogDebug("PiZeroMcAna") << "***>> No Track to KTrack match made for run " << fRun << ", event: " << fEvent << std::endl;
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
double PiZeroMcAna::length(const recob::Track* track)
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
double PiZeroMcAna::length(const simb::MCParticle& part, double dx,
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

double PiZeroMcAna::distanceToTPCEdge(const TVector3& position) const
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
// .fcl file; see PiZeroMcAna.fcl for more information.
DEFINE_ART_MODULE(PiZeroMcAna)

} // namespace PiZeroMcAna

#endif // PiZeroMcAna_Module
