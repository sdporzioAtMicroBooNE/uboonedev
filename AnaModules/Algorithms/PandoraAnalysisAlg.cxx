
#include "PandoraAnalysisAlg.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include <cmath>
#include <algorithm>

namespace pandoraanalysis
{
    
enum class TH1DLabels : size_t
{
    NPrimPFParticles,
    NPrimWDaughters,
    Npdg13PFParticles,
    Npdg11PFParticles,
    NDaughters,
    NTracks,
    NVertices,
    NHitsPFparticles,
    NHitsWDaughters,
    NHitsNDaughters,
    NumPFPartTracks,
    NumPFPartTracksL,
    NumFitTracks,
    NumTrksPFPart,
    NumTrksPFPartL,
    NumFitTrksPFPart,
    NumFitTrksPFPartL,
    PFPartTrackLen,
    PFPartTrackLenL,
    PFPartEndLenL,
    FitTrackLen,
    FitTrackProjLen,
    FitEndLen,
    TrackDeltaLen,
    TrackDeltaProjLen,
    DeltaStartPos,
    DeltaEndPos,
    TrackDeltaStart,
    CosTracks,
    DeltaWiresTrk,
    NumHitsTrk,
    HitWireRatioTrk,
    TrajDispDiff,
    TrajDispAng,
    TrajAng,
    TrajDocaAll,
    TrajDoca,
    TrajStartDiff,
    TrajEndDiff,
    NumPFPartHits,
    NumPFPartViews,
    NumPFPartViewsL,
    AngleToWire,
    AngleToAxis,
    Count
};
    
enum class TH2DLabels : size_t
{
    FitVsPFPartLen,
    FitELVsTL,
    DStartVsDEnd,
    ViewVsHits,
    HitPHVsAngle,
    Count
};

enum class TProfileLabels : size_t
{
    FitVsPFPartEff,
    FitVsPFNHitsEff,
    RatioVsDWires,
    Count
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
PandoraAnalysisAlg::PandoraAnalysisAlg(fhicl::ParameterSet const & pset) 
{
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);

    // Report.
    mf::LogInfo("PandoraAnalysisAlg") << "PandoraAnalysisAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
PandoraAnalysisAlg::~PandoraAnalysisAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void PandoraAnalysisAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fHitProducerLabel        = pset.get< std::string >("HitModuleLabel",          "gauss");
    fPFParticleProducerLabel = pset.get< std::string >("PFParticleProducerLabel", "cluster3d");
    fTrackProducerLabel      = pset.get< std::string >("TrackProducerLabel",      "trackkalmanhit");
}

//----------------------------------------------------------------------------
/// Begin job method.
void PandoraAnalysisAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());
    
    fTH1DVec.resize(size_t(TH1DLabels::Count),0);
    fTH1DVecVec.resize(size_t(TH1DLabels::Count),std::vector<TH1D*>(0));
    fTH2DVec.resize(size_t(TH2DLabels::Count),0);
    fTH2DVecVec.resize(size_t(TH2DLabels::Count),std::vector<TH2D*>(0));
    fTProfileVec.resize(size_t(TProfileLabels::Count),0);
    fTProfileVecVec.resize(size_t(TProfileLabels::Count),std::vector<TProfile*>(0));

    fTH1DVec[size_t(TH1DLabels::NPrimPFParticles)]            = dir.make<TH1D>("NPrimPFParticles", ";# particles",      100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::NPrimWDaughters)]             = dir.make<TH1D>("PrimWDaughters",   ";# particles",      100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::Npdg13PFParticles)]           = dir.make<TH1D>("Npdg13PFParts",    ";# particles",      100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::Npdg11PFParticles)]           = dir.make<TH1D>("Npdg11PFParts",    ";# particles",      100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::NDaughters)]                  = dir.make<TH1D>("NpDaughters",      ";# daughters",       20,    0.,   20.);
    fTH1DVec[size_t(TH1DLabels::NTracks)]                     = dir.make<TH1D>("NTracks",          ";# tracks",          20,    0.,   20.);
    fTH1DVec[size_t(TH1DLabels::NVertices)]                   = dir.make<TH1D>("NVertices",        ";# vertices",        20,    0.,   20.);
    fTH1DVec[size_t(TH1DLabels::NHitsPFparticles)]            = dir.make<TH1D>("NHitsPFParticles", ";# hits",           200,    0., 2000.);
    fTH1DVec[size_t(TH1DLabels::NHitsWDaughters)]             = dir.make<TH1D>("NHitsWDaughters",  ";# hits",           200,    0., 2000.);
    fTH1DVec[size_t(TH1DLabels::NHitsNDaughters)]             = dir.make<TH1D>("NHitsNDaughters",  ";# hits",           200,    0., 2000.);
    
    fTH1DVec[size_t(TH1DLabels::NumPFPartTracks)]             = dir.make<TH1D>("NPFPartTracks",    ";# Tracks",         100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::NumPFPartTracksL)]            = dir.make<TH1D>("NPFPartTracksL",   ";# Tracks",         100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::NumFitTracks)]                = dir.make<TH1D>("NFitTracks",       ";# Tracks",         100,    0.,  100.);
    fTH1DVec[size_t(TH1DLabels::NumTrksPFPart)]               = dir.make<TH1D>("NTrksPFPart",      " # Tracks",          10,    0.,   10.);
    fTH1DVec[size_t(TH1DLabels::NumTrksPFPartL)]              = dir.make<TH1D>("NTrksPFPartL",     " # Tracks",          10,    0.,   10.);
    fTH1DVec[size_t(TH1DLabels::NumFitTrksPFPart)]            = dir.make<TH1D>("NFitTrksPFPart",   " # Tracks",          10,    0.,   10.);
    fTH1DVec[size_t(TH1DLabels::NumFitTrksPFPartL)]           = dir.make<TH1D>("NFitTrksPFPartL",  " # Tracks",          10,    0.,   10.);
    fTH1DVec[size_t(TH1DLabels::PFPartTrackLen)]              = dir.make<TH1D>("PFPartTrackLen",   ";track len",        250,    0.,  500.);
    fTH1DVec[size_t(TH1DLabels::PFPartTrackLenL)]             = dir.make<TH1D>("PFPartTrackLenL",  ";track len",        250,    0.,  500.);
    fTH1DVec[size_t(TH1DLabels::PFPartEndLenL)]               = dir.make<TH1D>("PFPartEndLenL",    ";End Point len",    250,    0.,  500.);
    fTH1DVec[size_t(TH1DLabels::FitTrackLen)]                 = dir.make<TH1D>("FitTrackLen",      ";track len",        250,    0.,  500.);
    fTH1DVec[size_t(TH1DLabels::FitTrackProjLen)]             = dir.make<TH1D>("FitTrackProjLen",  ";track len",        250,    0.,  500.);
    fTH1DVec[size_t(TH1DLabels::FitEndLen)]                   = dir.make<TH1D>("FitEndLen",        ";End Point len",    250,    0.,  500.);
    fTH1DVec[size_t(TH1DLabels::TrackDeltaLen)]               = dir.make<TH1D>("TrackDeltaLen",    ";Delta Len",        500, -250.,  250.);
    fTH1DVec[size_t(TH1DLabels::TrackDeltaProjLen)]           = dir.make<TH1D>("TrackDeltaProjLen",";Delta Len",        500, -250.,  250.);
    
    fTH2DVec[size_t(TH2DLabels::FitVsPFPartLen)]              = dir.make<TH2D>("FitVsPFPartLen",   ";length;length",    100,    0.,  500., 100, 0., 500.);
    fTH2DVec[size_t(TH2DLabels::FitELVsTL)]                   = dir.make<TH2D>("FitELVsTL",        ";length;length",    100,    0.,  500., 100, 0., 500.);
    fTProfileVec[size_t(TProfileLabels::FitVsPFPartEff)]      = dir.make<TProfile>("FitVsPFPart",  ";length(cm)",       100.,   0.,  500.,      0., 1.1);
    fTProfileVec[size_t(TProfileLabels::FitVsPFNHitsEff)]     = dir.make<TProfile>("FitVsPFNHits", ";# hits",           100.,   0., 3000.,      0., 1.1);
    
    fTH1DVec[size_t(TH1DLabels::DeltaStartPos)]               = dir.make<TH1D>("DeltaStartPos",    ";delta(cm)",        100,    0.,   50.);
    fTH1DVec[size_t(TH1DLabels::DeltaEndPos)]                 = dir.make<TH1D>("DeltaEndPos",      ";delta(cm)",        100,    0.,   50.);
    fTH1DVec[size_t(TH1DLabels::TrackDeltaStart)]             = dir.make<TH1D>("TrackDeltaStart",  ";delta(cm)",        100,    0.,   50.);
    fTH1DVec[size_t(TH1DLabels::CosTracks)]                   = dir.make<TH1D>("CosTracks",        ";cos(theta)",       101,    0.,    1.01);
    
    fTH2DVec[size_t(TH2DLabels::DStartVsDEnd)]                = dir.make<TH2D>("DStartVsDEnd",     ";delta;delta",       50,    0.,   50., 50, 0., 50.);
    
    fTH1DVecVec[size_t(TH1DLabels::DeltaWiresTrk)].resize(3,0);
    fTH1DVecVec[size_t(TH1DLabels::DeltaWiresTrk)][0]         = dir.make<TH1D>("DltaWiresTrk0",    ";deltaWires",       250,    0., 1000.);
    fTH1DVecVec[size_t(TH1DLabels::DeltaWiresTrk)][1]         = dir.make<TH1D>("DltaWiresTrk1",    ";deltaWires",       250,    0., 1000.);
    fTH1DVecVec[size_t(TH1DLabels::DeltaWiresTrk)][2]         = dir.make<TH1D>("DltaWiresTrk2",    ";deltaWires",       250,    0., 1000.);
    
    fTH1DVecVec[size_t(TH1DLabels::NumHitsTrk)].resize(3,0);
    fTH1DVecVec[size_t(TH1DLabels::NumHitsTrk)][0]            = dir.make<TH1D>("NumHitsTrk0",      ";# hits",           250,    0., 1000.);
    fTH1DVecVec[size_t(TH1DLabels::NumHitsTrk)][1]            = dir.make<TH1D>("NumHitsTrk1",      ";# hits",           250,    0., 1000.);
    fTH1DVecVec[size_t(TH1DLabels::NumHitsTrk)][2]            = dir.make<TH1D>("NumHitsTrk2",      ";# hits",           250,    0., 1000.);
    
    fTH1DVecVec[size_t(TH1DLabels::HitWireRatioTrk)].resize(3,0);
    fTH1DVecVec[size_t(TH1DLabels::HitWireRatioTrk)][0]       = dir.make<TH1D>("HitWireRat0",      ";Ratio",            100,    0.,    2.);
    fTH1DVecVec[size_t(TH1DLabels::HitWireRatioTrk)][1]       = dir.make<TH1D>("HitWireRat1",      ";Ratio",            100,    0.,    2.);
    fTH1DVecVec[size_t(TH1DLabels::HitWireRatioTrk)][2]       = dir.make<TH1D>("HitWireRat2",      ";Ratio",            100,    0.,    2.);
    
    fTProfileVecVec[size_t(TProfileLabels::RatioVsDWires)].resize(3,0);
    fTProfileVecVec[size_t(TProfileLabels::RatioVsDWires)][0] = dir.make<TProfile>("RatVsDWire0",  ";DeltaWires;Ratio", 200,    0.,  200., 0., 2.);
    fTProfileVecVec[size_t(TProfileLabels::RatioVsDWires)][1] = dir.make<TProfile>("RatVsDWire1",  ";DeltaWires;Ratio", 200,    0.,  200., 0., 2.);
    fTProfileVecVec[size_t(TProfileLabels::RatioVsDWires)][2] = dir.make<TProfile>("RatVsDWire2",  ";DeltaWires;Ratio", 200,    0.,  200., 0., 2.);
    
    fTH1DVec[size_t(TH1DLabels::TrajDispDiff)]                = dir.make<TH1D>("TrajPointDisp",    ";disp",             200,   -5.,    5.);
    fTH1DVec[size_t(TH1DLabels::TrajDispAng)]                 = dir.make<TH1D>("TrajPointAng",     ";cos(ang)",         100,   -1.,    1.);
    fTH1DVec[size_t(TH1DLabels::TrajAng)]                     = dir.make<TH1D>("TrajAng",          ";cos(ang)",         100,   -1.,    1.);
    fTH1DVec[size_t(TH1DLabels::TrajDocaAll)]                 = dir.make<TH1D>("TrajDocaAll",      ";doca",             100,    0.,   10.);
    
    fTH1DVecVec[size_t(TH1DLabels::TrajDoca)].resize(3,0);
    fTH1DVecVec[size_t(TH1DLabels::TrajDoca)][0]              = dir.make<TH1D>("TrajDoca_0",       ";doca",             100,    0.,   10.);
    fTH1DVecVec[size_t(TH1DLabels::TrajDoca)][1]              = dir.make<TH1D>("TrajDoca_1",       ";doca",             100,    0.,   10.);
    fTH1DVecVec[size_t(TH1DLabels::TrajDoca)][2]              = dir.make<TH1D>("TrajDoca_2",       ";doca",             100,    0.,   10.);
    
    fTH1DVec[size_t(TH1DLabels::TrajStartDiff)]               = dir.make<TH1D>("TrajDeltaStart",   ";delta(cm)",        100,    0.,   50.);
    fTH1DVec[size_t(TH1DLabels::TrajEndDiff)]                 = dir.make<TH1D>("TrajDeltaEnd",     ";delta(cm)",        100,    0.,   50.);
    
    fTH1DVec[size_t(TH1DLabels::NumPFPartHits)]               = dir.make<TH1D>("NPFPartHits",      ";# hits",           300,    0., 3000.);
    fTH1DVec[size_t(TH1DLabels::NumPFPartViews)]              = dir.make<TH1D>("NPFPartViews",     ";# views",            5,    0.,    5.);
    fTH1DVec[size_t(TH1DLabels::NumPFPartViewsL)]             = dir.make<TH1D>("NPFPartViewsL",    ";# views",            5,    0.,    5.);
    fTH2DVec[size_t(TH2DLabels::ViewVsHits)]                  = dir.make<TH2D>("ViewVsHits",       ";# hits;# views",   100,    0., 3000., 5, 0., 5.);

    fTH1DVec[size_t(TH1DLabels::AngleToWire)]                 = dir.make<TH1D>("AngleToWire",     ";angle(deg)",        180,    0.,  180.);
    fTH1DVec[size_t(TH1DLabels::AngleToAxis)]                 = dir.make<TH1D>("AngleToX",        ";angle(deg)",        180,    0.,  180.);
    
    fTH2DVecVec[size_t(TH2DLabels::HitPHVsAngle)].resize(3,0);
    fTH2DVecVec[size_t(TH2DLabels::HitPHVsAngle)][0]          = dir.make<TH2D>("HitPHVsAngle0",   ";angle(deg);PH",      90,    0.,  180., 150, 0.,  150.);
    fTH2DVecVec[size_t(TH2DLabels::HitPHVsAngle)][1]          = dir.make<TH2D>("HitPHVsAngle1",   ";angle(deg);PH",      90,    0.,  180., 150, 0.,  150.);
    fTH2DVecVec[size_t(TH2DLabels::HitPHVsAngle)][2]          = dir.make<TH2D>("HitPHVsAngle2",   ";angle(deg);PH",      90,    0.,  180., 150, 0.,  150.);
    
    return;
}
    
void PandoraAnalysisAlg::pandoraAnalysis(const art::Event& event) const
{
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
    
    // Now we want to try to look at tracks associated to PFParticles
    art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
    event.getByLabel(fPFParticleProducerLabel, pfParticleHandle);
    art::Handle<std::vector<recob::Cluster>> clusterHandle;
    event.getByLabel(fPFParticleProducerLabel, clusterHandle);
    
    if (pfParticleHandle.isValid() && clusterHandle.isValid())
    {
        // Recover the collection of associations between PFParticles and Tracks from the PFParticle producer
        art::FindManyP<recob::Track> firstTrackAssns(pfParticleHandle, event, fPFParticleProducerLabel);
        
        // Now get the track associations from the track producer
        art::FindManyP<recob::Track> secondTrackAssns(pfParticleHandle, event, fTrackProducerLabel);
        
        // Recover the collection of associations between PFParticles and vertices from the PFParticle producer
        art::FindManyP<recob::Vertex> vertexAssns(pfParticleHandle, event, fPFParticleProducerLabel);
        
        // Next we get the clusters associated to the PFParticle
        art::FindManyP<recob::Cluster> clusterAssns(pfParticleHandle, event, fPFParticleProducerLabel);
        
        // Finally we get hit cluster associations
        art::FindManyP<recob::Hit> clusterHitAssns(clusterHandle, event, fPFParticleProducerLabel);
        
        // Recover the collection of associations between tracks and hits
        art::FindManyP<recob::Hit> trackHitAssns(trackHandle, event, fTrackProducerLabel);
        
        // Recover the collection of associations between tracks and space points
        art::FindManyP<recob::SpacePoint> spacePointAssns(trackHandle, event, fTrackProducerLabel);
        
        // Recover the collection of associations between space points and hits
        art::FindManyP<recob::Hit> hitSpacePointAssns(trackHandle, event, fTrackProducerLabel);
        
        // Also grab a track handle from the PFParticle producer
        art::Handle<std::vector<recob::Track>> pfTrackHandle;
        event.getByLabel(fPFParticleProducerLabel, pfTrackHandle);
        
        if (firstTrackAssns.size() > 0 && secondTrackAssns.size() > 0 && pfTrackHandle.isValid() && clusterAssns.isValid())
        {
            int nPandoraTracks(0);
            int nPandoraTracksL(0);
            int nFitTracks(0);
            int nPrimaryPFParticles(0);
            int nPdg13PFParticles(0);
            int nPdg11PFParticles(0);
            int nPrimWDaughters(0);
            
            for(size_t pfParticleIdx = 0; pfParticleIdx < pfParticleHandle->size(); pfParticleIdx++)
            {
                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);
                
                // We would like # hits for this PFParticle and number views
                int nHitsPerView[] = {0,0,0};
                int nHitsTotal(0);
                int nViews(0);
                int nClusPerView[] = {0,0,0};
                
                std::vector<art::Ptr<recob::Cluster>> clusterVec = clusterAssns.at(pfParticle.key());
                
                for(auto& cluster : clusterVec)
                {
                    std::vector<art::Ptr<recob::Hit>> hitVector = clusterHitAssns.at(cluster.key());
                    
                    if (hitVector.empty()) continue;
                    
                    int view = cluster->View();
                    
                    nClusPerView[view]++;
                    nHitsPerView[view]  = hitVector.size();
                    nHitsTotal         += hitVector.size();
                }
                
                for(int viewIdx = 0; viewIdx < 3; viewIdx++)
                {
                    int nClus = nClusPerView[viewIdx];
                    
                    if (nClus > 1) std::cout << "=======> # clusters: " << nClus << std::endl;
                    
                    if (nHitsPerView[viewIdx] > 0) nViews++;
                }
                
                // Get the daughter count
                if (pfParticle->Parent() == recob::PFParticle::kPFParticlePrimary)
                {
                    int nTracks(0);
                    int nVertices(0);
                    int daughterCount = traversePFParticleHierarchy(pfParticleHandle, pfParticleIdx, firstTrackAssns, vertexAssns, nTracks, nVertices);
                    
                    nPrimaryPFParticles++;
                    
                    if (pfParticle->PdgCode() == 13)
                    {
                        nPdg13PFParticles++;
                        fTH1DVec[size_t(TH1DLabels::NHitsPFparticles)]->Fill(nHitsTotal, 1.);
                    }
                    else nPdg11PFParticles++;
                    
                    if (daughterCount > 0)
                    {
                        nPrimWDaughters++;
                        fTH1DVec[size_t(TH1DLabels::NHitsWDaughters)]->Fill(nHitsTotal, 1.);
                    }
                    else fTH1DVec[size_t(TH1DLabels::NHitsNDaughters)]->Fill(nHitsTotal, 1.);
                    
                    fTH1DVec[size_t(TH1DLabels::NDaughters)]->Fill(daughterCount, 1.);
                    fTH1DVec[size_t(TH1DLabels::NTracks)]->Fill(nTracks, 1.);
                    fTH1DVec[size_t(TH1DLabels::NVertices)]->Fill(nVertices, 1.);
                }
                
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
                        
                        double pfPartLen    = projectedLength(track.get()); //length(track.get());
                        double pfPartEndLen = (pfPartEnd - pfPartStart).Mag();
                        
                        fTH1DVec[size_t(TH1DLabels::PFPartTrackLen)]->Fill(std::min(499.5,pfPartLen),1.);
                        fTH1DVec[size_t(TH1DLabels::PFPartEndLenL)]->Fill(std::min(499.5,pfPartEndLen),1.);
                        nTrksPFPart++;
                        
                        if (pfPartLen > 5.)
                        {
                            nPandoraTracksL++;
                            fTH1DVec[size_t(TH1DLabels::PFPartTrackLenL)]->Fill(std::min(499.5,pfPartLen));
                            nTrksPFPartL++;
                        }
                        
                        bool   foundMatch(false);
                        double deltaLen(99999.);
                        double deltaProjLen(99999.);
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
                            double fitTrkLen     = length(fitTrk.get());
                            double fitTrkProjLen = projectedLength(fitTrk.get());
                            double fitTrkEndLen  = (fitTrk->End() - fitTrk->Vertex()).Mag();
                            
                            if (fabs(pfPartLen - fitTrkLen) < fabs(deltaLen))
                            {
                                foundMatch     = true;
                                deltaLen       = fitTrkEndLen - fitTrkLen; //pfPartLen - fitTrkLen;
                                deltaProjLen   = pfPartLen - fitTrkProjLen;
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
                            deltaLen     = std::max(-249.5,std::min(249.5,deltaLen));
                            deltaProjLen = std::max(-249.5,std::min(249.5,deltaProjLen));
                            trkMatch     = 1.;
                            
                            fTH1DVec[size_t(TH1DLabels::TrackDeltaLen)]->Fill(deltaLen, 1.);
                            fTH1DVec[size_t(TH1DLabels::TrackDeltaProjLen)]->Fill(deltaProjLen, 1.);
                            fTH2DVec[size_t(TH2DLabels::FitVsPFPartLen)]->Fill(std::min(499.5,pfPartLen), std::min(499.5,matchLen), 1.);
                            
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
                            
                            fTH1DVec[size_t(TH1DLabels::DeltaStartPos)]->Fill(startDist, 1.);
                            fTH1DVec[size_t(TH1DLabels::DeltaEndPos)]->Fill(endDist, 1.);
                            fTH1DVec[size_t(TH1DLabels::CosTracks)]->Fill(cosTrkAng, 1.);
                            
                            if (flipped) fTH1DVec[size_t(TH1DLabels::TrackDeltaStart)]->Fill(endDist, 1.);
                            else         fTH1DVec[size_t(TH1DLabels::TrackDeltaStart)]->Fill(startDist, 1.);
                            
                            fTH2DVec[size_t(TH2DLabels::DStartVsDEnd)]->Fill(startDist, endDist, 1.);
                            
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
                                    
                                    fTH1DVecVec[size_t(TH1DLabels::DeltaWiresTrk)][planeIdx]->Fill(deltaWires, 1.);
                                    fTH1DVecVec[size_t(TH1DLabels::NumHitsTrk)][planeIdx]->Fill(numHits, 1.);
                                    if (deltaWires > 7) fTH1DVecVec[size_t(TH1DLabels::HitWireRatioTrk)][planeIdx]->Fill(std::min(1.99,hitRatio), 1.);
                                    fTProfileVecVec[size_t(TProfileLabels::RatioVsDWires)][planeIdx]->Fill(std::min(deltaWires,199), std::min(1.99,hitRatio));
                                }
                                catch(...) {}
                            }
                        }
                        
                        fTProfileVec[size_t(TProfileLabels::FitVsPFPartEff)]->Fill(pfPartLen, trkMatch);
                        fTProfileVec[size_t(TProfileLabels::FitVsPFNHitsEff)]->Fill(nHitsTotal, trkMatch);
                    }
                    
                    fTH1DVec[size_t(TH1DLabels::NumTrksPFPart)]->Fill(nTrksPFPart, 1.);
                    fTH1DVec[size_t(TH1DLabels::NumTrksPFPartL)]->Fill(nTrksPFPartL, 1.);
                    
                    int nFitTrksPFPart(fitTrackVec.size());
                    
                    for(size_t trkIdx = 0; trkIdx < fitTrackVec.size(); trkIdx++)
                    {
                        art::Ptr<recob::Track> fitTrack = fitTrackVec.at(trkIdx);
                        
                        double fitLen    = std::min(499.5,length(fitTrack.get()));
                        double projLen   = std::min(499.5,projectedLength(fitTrack.get()));
                        double fitEndLen = std::min(499.5,(fitTrack->Vertex() - fitTrack->End()).Mag());
                        
                        fTH1DVec[size_t(TH1DLabels::FitTrackLen)]->Fill(fitLen,1.);
                        fTH1DVec[size_t(TH1DLabels::FitTrackProjLen)]->Fill(projLen,1.);
                        fTH1DVec[size_t(TH1DLabels::FitEndLen)]->Fill(fitEndLen,1.);
                        fTH2DVec[size_t(TH2DLabels::FitELVsTL)]->Fill(fitEndLen,fitLen,1.);
                        
                        // Quick check to look for tracks which have trajectories gone awry
                        TVector3 lastTrajPos(fitTrack->LocationAtPoint(0));
                        TVector3 lastPoint(fitTrack->LocationAtPoint(0));
                        TVector3 lastTrajDir(fitTrack->DirectionAtPoint(0));
                        TVector3 lastPntDir(0.,0.,0.);
                        
                        TVector3 trajStartDiff = fitTrack->Vertex() - fitTrack->LocationAtPoint(1);
                        TVector3 trajEndDiff   = fitTrack->End()    - fitTrack->LocationAtPoint(fitTrack->NumberTrajectoryPoints() - 2);
                        
                        fTH1DVec[size_t(TH1DLabels::TrajStartDiff)]->Fill(std::min(trajStartDiff.Mag(),49.9), 1.);
                        fTH1DVec[size_t(TH1DLabels::TrajEndDiff)]->Fill(std::min(trajEndDiff.Mag(),49.9), 1.);
                        
                        // Make the assumption that the hit order in the association corresponds to the trajectory point order
                        std::vector<art::Ptr<recob::Hit>> trackHitVector = trackHitAssns.at(fitTrack.key());
                        
                        // What is the preferred view for this track?
                        std::map<size_t,HitPtrVec>& viewHitMap  = trackHitVecMap[fitTrack.key()];
                        int                         numTrackHits(0);
                        
                        std::vector<std::pair<size_t,size_t>> viewCountVec;
                        
                        for(const auto& viewHit : viewHitMap)
                        {
                            viewCountVec.push_back(std::make_pair(viewHit.first,viewHit.second.size()));
                            numTrackHits += viewHit.second.size();
                        }
                        
                        std::sort(viewCountVec.begin(),viewCountVec.end(),[](const auto& left,const auto& right){return left.second > right.second;});
                        
                        std::vector<size_t> bestViewMap;
                        
                        bestViewMap.resize(viewCountVec.size(), 0);
                        
                        for(size_t idx = 0; idx < viewCountVec.size(); idx++) bestViewMap[viewCountVec[idx].first] = idx;
                        
                        // Let's try to keep track of the doca's... looking for a preferred direction...
                        std::vector<std::vector<double>> docaVecByView;
                        
                        for(size_t trajIdx = 1; trajIdx < fitTrack->NumberTrajectoryPoints(); trajIdx++)
                        {
                            art::Ptr<recob::Hit> trackHit = trackHitVector.at(trajIdx);
                            size_t               view     = trackHit->View();
                            const TVector3&      trajPos  = fitTrack->LocationAtPoint(trajIdx);
                            const TVector3&      trajDir  = fitTrack->DirectionAtPoint(trajIdx);
                            TVector3             dir      = trajPos - lastTrajPos;
                            double               disp     = dir.Mag();
                            
                            const geo::WireID    wireID   = trackHit->WireID();
                            const geo::WireGeo&  wireGeo  = fGeometry->Wire(wireID);
                            
                            TVector3             wireVec  = wireGeo.Direction();
                            
                            // Want to break out y-z components from angle to x axis
                            TVector3 trajYZ(0.,trajDir.Y(),trajDir.Z());
                            
                            if (trajYZ.Mag() > 0.) trajYZ.SetMag(1.);
                            
                            double cosAngToWire = trajYZ.Dot(wireVec);
                            double angleToWire  = 180. * acos(cosAngToWire) / M_PI;
                            double angleToAxis  = 180. * acos(trajDir.X()) / M_PI;
                            
//                            if (abs(angleToAxis - 90.) < 0.25)
//                            {
//                                std::cout << "* angleToAxis: " << angleToAxis << ", trajDir: " << trajDir.X() << "," << trajDir.Y() << "," << trajDir.Z() << std::endl;
//                            }
                           
                            fTH1DVec[size_t(TH1DLabels::AngleToWire)]->Fill(angleToWire, 1.);
                            fTH1DVec[size_t(TH1DLabels::AngleToAxis)]->Fill(angleToAxis, 1.);
                            fTH2DVecVec[size_t(TH2DLabels::HitPHVsAngle)][view]->Fill(angleToWire, std::min(trackHit->PeakAmplitude(),float(249.8)), 1.);
                            
                            if (disp > 0.) dir.SetMag(1.);
                            
                            if (dir.Dot(lastPntDir) < 0.) disp = -disp;
                            
                            // Now get doca of next trajectory point
                            TVector3 lastToNewPoint = trajPos - lastPoint;
                            double   arcLenToDoca   = lastTrajDir.Dot(lastToNewPoint);
                            
                            lastPoint = lastPoint + arcLenToDoca * lastTrajDir;
                            
                            TVector3 docaVec = trajPos - lastPoint;
                            double   doca    = docaVec.Mag();
                            double   cosTraj = lastTrajDir.Dot(trajDir);
                            
                            fTH1DVec[size_t(TH1DLabels::TrajDispDiff)]->Fill(std::max(-4.99,std::min(4.99,disp)), 1.);
                            fTH1DVec[size_t(TH1DLabels::TrajDispAng)]->Fill(std::max(-0.99,std::min(0.99,dir.Dot(lastPntDir))), 1.);
                            fTH1DVec[size_t(TH1DLabels::TrajAng)]->Fill(std::max(-0.99,std::min(0.99,cosTraj)), 1.);
                            fTH1DVec[size_t(TH1DLabels::TrajDocaAll)]->Fill(std::min(doca,9.99), 1.);
                            fTH1DVecVec[size_t(TH1DLabels::TrajDoca)][bestViewMap[view]]->Fill(std::min(doca,9.99), 1.);
                            
                            lastPntDir  = dir;
                            lastTrajPos = trajPos;
                            lastTrajDir = trajDir;
                        }
                    }
                    
                    fTH1DVec[size_t(TH1DLabels::NumFitTrksPFPart)]->Fill(nFitTrksPFPart, 1.);
                    if (nTrksPFPartL > 0) fTH1DVec[size_t(TH1DLabels::NumFitTrksPFPartL)]->Fill(nFitTrksPFPart, 1.);
                    
                    if (nTrksPFPartL > 0 && nFitTrksPFPart < 1)
                    {
                        fTH1DVec[size_t(TH1DLabels::NumPFPartHits)]->Fill(nHitsTotal, 1.);
                        fTH1DVec[size_t(TH1DLabels::NumPFPartViews)]->Fill(nViews, 1.);
                        if (nHitsTotal > 500) fTH1DVec[size_t(TH1DLabels::NumPFPartViewsL)]->Fill(nViews, 1.);
                        fTH2DVec[size_t(TH2DLabels::ViewVsHits)]->Fill(nHitsTotal, nViews, 1.);
                    }
                    
                    nPandoraTracks += pfPartTrackVec.size();
                    nFitTracks     += fitTrackVec.size();
                }
            }
            
//            std::cout << "~~ # pandora tracks: " << nPandoraTracks << ", long: " << nPandoraTracksL << ", # fit tracks: " << nFitTracks << std::endl;
            fTH1DVec[size_t(TH1DLabels::NumPFPartTracks)]->Fill(nPandoraTracks, 1.);
            fTH1DVec[size_t(TH1DLabels::NumPFPartTracksL)]->Fill(nPandoraTracksL, 1.);
            fTH1DVec[size_t(TH1DLabels::NumFitTracks)]->Fill(nFitTracks, 1.);
            
            fTH1DVec[size_t(TH1DLabels::NPrimPFParticles)]->Fill(nPrimaryPFParticles, 1.);
            fTH1DVec[size_t(TH1DLabels::NPrimWDaughters)]->Fill(nPrimWDaughters, 1.);
            fTH1DVec[size_t(TH1DLabels::Npdg13PFParticles)]->Fill(nPdg13PFParticles, 1.);
            fTH1DVec[size_t(TH1DLabels::Npdg11PFParticles)]->Fill(nPdg11PFParticles, 1.);
        }
        
    }
    
    return;
}

// Compare two tracks
//----------------------------------------------------------------------------
void PandoraAnalysisAlg::compareTwoTracks(const recob::Track* track1, const recob::Track* track2) const
{
    return;
}
    
// Useful for normalizing histograms
void PandoraAnalysisAlg::endJob(int numEvents)
{
    return;
}
    
// Length of reconstructed track.
//----------------------------------------------------------------------------
double PandoraAnalysisAlg::length(const recob::Track* track) const
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
        
        //        if (lastDir.Dot(trajDir) >= 0.)
        //        {
        disp   -= pos;
        result += disp.Mag();
        disp    = pos;
        //        }
        
        lastPoint = pos;
        lastDir   = trajDir;
    }
    
    return result;
}

// Length of reconstructed track.
//----------------------------------------------------------------------------
double PandoraAnalysisAlg::projectedLength(const recob::Track* track) const
{
    double   result(0.);
    TVector3 lastPoint(track->LocationAtPoint(0));
    TVector3 lastDir(track->DirectionAtPoint(0));
    int      n(track->NumberTrajectoryPoints());
    
    for(int i = 1; i < n; ++i)
    {
        const TVector3& newPoint = track->LocationAtPoint(i);
        
        TVector3 lastToNewPoint = newPoint - lastPoint;
        double   arcLenToDoca   = lastDir.Dot(lastToNewPoint);
        
        result    += arcLenToDoca;
        lastPoint  = lastPoint + arcLenToDoca * lastDir;
        lastDir    = track->DirectionAtPoint(i);
    }
    
    return result;
}

int PandoraAnalysisAlg::traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
                                                    size_t                                       pfParticleIdx,
                                                    const art::FindManyP<recob::Track>&          trackAssns,
                                                    const art::FindManyP<recob::Vertex>&         vertexAssns,
                                                    int&                                         nTracks,
                                                    int&                                         nVertices) const
{
    // So far no daughters...
    int nDaughters(0);
    
    // Get pointer to PFParticle
    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);
    
    // Recover tracks/vertices associated to this PFParticle
    std::vector<art::Ptr<recob::Track>>  pfPartTrackVec  = trackAssns.at(pfParticle.key());
    std::vector<art::Ptr<recob::Vertex>> pfPartVertexVec = vertexAssns.at(pfParticle.key());
    
    nTracks    += pfPartTrackVec.size();
    nVertices  += pfPartVertexVec.size();
    nDaughters += pfParticle->Daughters().size();
    
    for(auto& daughterIdx : pfParticle->Daughters())
    {
        nDaughters += traversePFParticleHierarchy(pfParticleHandle, daughterIdx, trackAssns, vertexAssns, nTracks, nVertices);
    }
    
    return nDaughters;
}
}
