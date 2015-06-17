/**
 *  @file   CGALTest_module.cc
 *
 *          Class:       CGALTest
 *          Module Type: Producer
 * 
 *  @brief  Producer module to test various features available in CGAL (see cgal.org)
 *
 *  @author usher@slac.stanford.edu 
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/HalfEdge.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/PCAxis.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Vertex.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/convex_hull_2.h"

#include "CGAL/point_generators_3.h"
#include "CGAL/Delaunay_triangulation_2.h"
#include "CGAL/Polyhedron_3.h"
#include "CGAL/convex_hull_3_to_polyhedron_3.h"
#include "CGAL/algorithm.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>
#include <functional>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_CGALTest
{
    
/**
 *  @brief  Definition of the CGALTest class
 */
class CGALTest : public art::EDProducer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset - reference to the parameters used by this module and its algorithms
     */
    CGALTest(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~CGALTest();

    /**
     *  @brief declare the standard art functions that we'll implement in this producer module
     */
    void beginJob();
    void endJob();
    void produce(art::Event &evt);
    void reconfigure(fhicl::ParameterSet const &pset);

private:
    
    // Useful typdefs (I hope)
    typedef CGAL::Exact_predicates_inexact_constructions_kernel CGALKernel;
    typedef CGALKernel::Point_2                                 Point_2;
    typedef std::vector<Point_2>                                CGALPointsVec;
    //typedef std::map<const Point_2, const recob::SpacePoint*>   CGALPointsToSpacePointsMap;
    
    typedef CGALKernel::Point_3                                 Point_3;
    typedef std::list<Point_3>                                  CGALPoints3List;
    
    typedef CGAL::Delaunay_triangulation_2<CGALKernel>          Delaunay;
    //typedef CGAL::Triangulation_3<CGALKernel>                   Delaunay;
    typedef Delaunay::Point                                     Point;
    typedef Delaunay::Vertex_handle                             Vertex_handle;
    typedef Delaunay::Edge                                      Edge_handle;
    typedef Delaunay::Edge_circulator                           Edge_circulator;
    typedef Delaunay::Edge_iterator                             Edge_iterator;
//    typedef Delaunay::Cell_handle                               Cell_handle;
//    typedef CGAL::Triple<Cell_handle, int, int>                 EdgeTriple;
    
    typedef std::map<const Point, const recob::SpacePoint*>     CGALPointsToSpacePointsMap;
    typedef std::list<Point>                                    DelaunayPointList;
    
    typedef CGAL::Polyhedron_3<CGALKernel>                      Polyhedron_3;

    /**
     *  @brief  Event Preparation 
     * 
     *  @param  evt  the ART event 
     */
    void PrepareEvent(const art::Event &evt);  

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    /**
     *   Algorithm parameters
     */
    bool                      fEnableMonitoring;      ///< Turn on monitoring of this algorithm
    std::string               fHitfinderProducer;     ///< Producer of the reco hits
    std::string               fPFParticleProducer;    ///< Producer of the PFParticles

    /**
     *   Tree variables for output
     */
    TTree*                    fRecoTree;              ///<
    int                       fRun;                   ///<
    int                       fEvent;                 ///<
    int                       fHits;                  ///< Keeps track of the number of hits seen
    float                     fTotalTime;             ///< Keeps track of total execution time
    
    /** 
     *   Other useful variables
     */
    geo::Geometry*            fGeometry;              ///<  pointer to the Geometry service
    util::DetectorProperties* fDetector;              ///<  Pointer to the detector properties
};

DEFINE_ART_MODULE(CGALTest)

} // namespace lar_CGALTest

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_CGALTest {

CGALTest::CGALTest(fhicl::ParameterSet const &pset)
{
    this->reconfigure(pset);

    produces< std::vector<recob::SpacePoint>>();
    produces< std::vector<recob::HalfEdge>>();
    produces< art::Assns<recob::PFParticle, recob::SpacePoint>>();
    produces< art::Assns<recob::PFParticle, recob::HalfEdge>>();
    produces< art::Assns<recob::SpacePoint, recob::HalfEdge>>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CGALTest::~CGALTest()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CGALTest::reconfigure(fhicl::ParameterSet const &pset)
{
    fHitfinderProducer  = pset.get<std::string>("HitFinderProducer",  "gaushit");
    fPFParticleProducer = pset.get<std::string>("PFParticleProducer", "cluster3d");
    fEnableMonitoring   = pset.get<bool>       ("EnableMonitoring",   false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CGALTest::beginJob()
{
    /**
     *  @brief beginJob will be tasked with initializing monitoring, in necessary, but also to init the 
     *         geometry and detector services (and this probably needs to go in a "beginEvent" method?)
     */
    if (fEnableMonitoring)
        this->InitializeMonitoring();
    
    art::ServiceHandle<geo::Geometry>            geometry;
    art::ServiceHandle<util::DetectorProperties> detectorProperties;
    
    fGeometry = &*geometry;
    fDetector = &*detectorProperties;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CGALTest::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CGALTest::produce(art::Event &evt)
{
    /**
     *  @brief Producer method for reovering the 2D hits and driving the 3D reconstruciton
     */
    mf::LogInfo("CGALTest") << " *** CGALTest::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] Starting Now! *** " << std::endl;

    // Set up for monitoring the timing... at some point this should be removed in favor of
    // external profilers
    cet::cpu_timer theClockTotal;
    
    if (fEnableMonitoring)
    {
        theClockTotal.start();
    }
    
    // This really only does anything if we are monitoring since it clears our tree variables
    this->PrepareEvent(evt);
    
    mf::LogDebug("CGALTest") << " *** CGALTest::ProduceArtClusters() *** " << std::endl;
    
    std::unique_ptr< std::vector<recob::SpacePoint>> artSpacePointVector( new std::vector<recob::SpacePoint> );
    std::unique_ptr< std::vector<recob::HalfEdge>>   artHalfEdgeVector( new std::vector<recob::HalfEdge> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint>>  artPFPartSPAssociations(     new art::Assns<recob::PFParticle, recob::SpacePoint>   );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::HalfEdge>>    artPFPartHEAssociations(     new art::Assns<recob::PFParticle, recob::HalfEdge>     );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::HalfEdge>>    artSpacePointHEAssociations( new art::Assns<recob::SpacePoint, recob::HalfEdge>     );
    
    // Ok, recover the PFParticles from which all else derives
    art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel(fPFParticleProducer, pfParticleHandle);
    
    if (pfParticleHandle.isValid())
    {
        // Now recover the PCA axes and SpacePoints
        art::FindMany<recob::PCAxis>     pcAxisAssnVec(    pfParticleHandle, evt, fPFParticleProducer);
        art::FindMany<recob::SpacePoint> spacePointAssnVec(pfParticleHandle, evt, fPFParticleProducer);

        // We need these for building SpacePoints associated with the convex hull
        int    spacePointIdx(0);
        int    halfEdgeIdx(0);

        // If no valid space point associations then nothing to do
        if (spacePointAssnVec.isValid() && pcAxisAssnVec.isValid())
        {
            // Ok, loop over PFParticles and begin the process!
            for(size_t idx = 0; idx < pfParticleHandle->size(); idx++)
            {
                if (idx > 0) break;
                
                // Recover cluster
                const art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, idx);
                
                // Consider only "primary" PFParticles
                if (!pfParticle->IsPrimary()) continue;
                
                // Look up the associated PCA
                std::vector<const recob::PCAxis*> pcaVec(pcAxisAssnVec.at(pfParticle->Self()));
                
                if (pcaVec.empty()) continue;
                
                // And the space points
                const std::vector<const recob::SpacePoint*>& hitsVec(spacePointAssnVec.at(pfParticle->Self()));
                
                if (hitsVec.empty()) continue;
                
                // Get the primary axes
                const recob::PCAxis& pca(*pcaVec.front());
                
                // Ok, the idea here is that one wants to determine the 3D space position of the projection of the space point
                // onto the plane defined by the two primary PCA axes. The method is to compute the distance vector from the plane
                // to the space point in question... this is done by getting a vector from the plane's origin to the space point
                // and then taking the projection onto the third axis. One then simply subtracts this in the direction of the plane
                // normal from the space point position.
                // Got that? Here we go!
                TVector3 planeOrigin(pca.getAvePosition()[0],     pca.getAvePosition()[1],     pca.getAvePosition()[2]);
                TVector3 axis1Dir(   pca.getEigenVectors()[0][0], pca.getEigenVectors()[0][1], pca.getEigenVectors()[0][2]);
                TVector3 axis2Dir(   pca.getEigenVectors()[1][0], pca.getEigenVectors()[1][1], pca.getEigenVectors()[1][2]);
                TVector3 axis3Dir(   pca.getEigenVectors()[2][0], pca.getEigenVectors()[2][1], pca.getEigenVectors()[2][2]);
                
                // Define the points vector to be filled
//                CGALPointsVec   cgalPointsVec;
//                CGALPoints3List   cgalPoints3List;
                DelaunayPointList delaunayPointList;
                
                // This will be a useful map
                CGALPointsToSpacePointsMap cgalPointsToSpacePointsMap;

                for(const auto spacePoint : hitsVec)
                {
                    // Skeleton points for now
//                    if (spacePoint->Chisq() != -1. && spacePoint->Chisq() != -3. && spacePoint->Chisq() != -4.) continue;
                    
                    // Recover the hull vertex position
                    const double* xyz = spacePoint->XYZ();
                    
                    TVector3 spacePointPos(xyz[0],xyz[1],xyz[2]);
                    TVector3 pcaToSpVec = spacePointPos - planeOrigin;
                    
                    // Projections (distance)
                    double proj1Dist = pcaToSpVec.Dot(axis1Dir);
                    double proj2Dist = pcaToSpVec.Dot(axis2Dir);
//                    double proj3Dist = pcaToSpVec.Dot(axis3Dir);

//                    cgalPointsVec.push_back(Point_2(proj1Dist, proj2Dist));
//                    cgalPoints3List.push_back(Point_3(xyz[0],xyz[1],xyz[2]));
//                    cgalPointsToSpacePointsMap[cgalPointsVec.back()] = spacePoint;
//                    delaunayPointList.emplace_back(Point(xyz[0] - planeOrigin.X(), xyz[1] - planeOrigin.Y(), xyz[2] - planeOrigin.Z()));
                    delaunayPointList.emplace_back(Point(proj1Dist, proj2Dist));
                    cgalPointsToSpacePointsMap[delaunayPointList.back()] = spacePoint;
                }
/*
                for(double xPos = -50.; xPos < 51.; xPos += 20.)
                {
                    for(double yPos = -50.; yPos < 51.; yPos += 20.)
                    {
                        for(double zPos = -50.; zPos < 51.; zPos += 20.)
                        {
                            delaunayPointList.emplace_back(Point(xPos,yPos,zPos));
                        }
                    }
                }
*/
 /*
                // define a container for the results
                CGALPointsVec convexHullVec;
                
                // Calculate it!
                CGAL::convex_hull_2(cgalPointsVec.begin(), cgalPointsVec.end(), std::back_inserter(convexHullVec));
                
                // Keep track of current start for space points
                std::vector<art::Ptr<recob::SpacePoint> > spacePointPtrVec;
                
                // See if we can recover the original space points
                for(const auto& point2 : convexHullVec)
                {
                    const recob::SpacePoint* spacePoint = cgalPointsToSpacePointsMap[point2];
                    
                    if (!spacePoint)
                    {
                        std::cout << "crap!" << std::endl;
                        continue;
                    }
                    
                    // Ok, clone from the original space point
                    artSpacePointVector->push_back(recob::SpacePoint(spacePoint->XYZ(),spacePoint->ErrXYZ(),spacePoint->Chisq(),spacePointIdx));
                    
                    // We also want to make associations so we need to recover an art::Ptr to the SpacePoint
                    art::ProductID spacePointId = getProductID<std::vector<recob::SpacePoint> >(evt);
                    art::Ptr<recob::SpacePoint> spacePointPtr(spacePointId, spacePointIdx, evt.productGetter(spacePointId));
                    spacePointPtrVec.push_back(spacePointPtr);
                    
                    spacePointIdx++;
                }
                
                // Make associations to the SpacePoints
                util::CreateAssn(*this, evt, pfParticle, spacePointPtrVec, *artPFPartSPAssociations);
*/
                if (delaunayPointList.size() > 3)
                {
                    // Try the 3D stuff here...
                    // Start with the Delaunay Triangulation?
                    Delaunay T;
                    T.insert(delaunayPointList.begin(),delaunayPointList.end());
                    
                    std::list<Vertex_handle> vertexList;
//                    std::list<Cell_handle>   cellList;
                    
//                    T.incident_vertices(T.infinite_vertex(), std::back_inserter(vertexList));
//                    T.incident_cells(T.infinite_vertex(),    std::back_inserter(cellList));
                    
//                    const Vertex_handle& vertexFront = vertexList.front();
                    
//                    std::cout << "Delaunay triangulation started with " << delaunayPointList.size() << " points and returned " << vertexList.size() << " vertices, " << cellList.size() << " cells" << ", # edges: " << T.tds().number_of_edges() << std::endl;
//                    std::cout << "  ==> # vertices: " << T.number_of_vertices() << ", faces: " << T.number_of_cells() << ", finite: " << T.number_of_finite_cells()
//                              << ", edges: " << T.number_of_edges() << ", finite: " << T.number_of_finite_edges() << std::endl;
//                    std::cout << "  --> first vertex: " << vertexFront->point() << std::endl;
/*
                    // Trial loop through cells to find all possible edges...
                    for(const auto& cellHandle : cellList)
                    {
                        for(int firstVtxIdx = 0; firstVtxIdx < 4; firstVtxIdx++)
                        {
                            const Vertex_handle firstVtxHandle = cellHandle->vertex(firstVtxIdx);
                            
                            TVector3 firstPoint(firstVtxHandle->point()[0],firstVtxHandle->point()[1],firstVtxHandle->point()[2]);
                            
                            for(int secondVtxIdx = firstVtxIdx + 1; secondVtxIdx < 4; secondVtxIdx++)
                            {
                                const Vertex_handle secondVtxHandle = cellHandle->vertex(secondVtxIdx);
                                
                                TVector3 secondPoint(secondVtxHandle->point()[0],secondVtxHandle->point()[1],secondVtxHandle->point()[2]);
                                TVector3 deltaPoint = secondPoint - firstPoint;
                                
                                std::cout << "  --> cell idxs: " << firstVtxIdx << ", " << secondVtxIdx << ", dist: " << deltaPoint.Mag() << std::endl;
                            }
                        }
                    }
*/
                    // Get the vertices/space points... keep a map to go back and forth
                    std::map<Vertex_handle, const recob::SpacePoint*> vertexToSpacePointMap;
                    
                    // Keep track of current start for space points
                    std::vector<art::Ptr<recob::SpacePoint> > spacePointPtrVec;
                    
                    double spacePointErr[] = {1., 0., 1., 0., 0., 1.};
                    
                    artSpacePointVector->reserve(delaunayPointList.size()+1);
                    
                    // Traverse ALL of the vertices in the triangulation and make space points, etc.
                    int numVertices(0);
                    
                    art::ProductID spacePointId = getProductID<std::vector<recob::SpacePoint> >(evt);
                    
                    for(auto vertexItr = T.tds().vertices_begin(); vertexItr != T.tds().vertices_end(); vertexItr++)
                    {
                        // Ignore the "infinite" vertex for now
                        //if (vertexItr == T.tds().vertices_begin()) continue;
                        
                        double spacePointPos[] = {planeOrigin.X(), planeOrigin.Y(), planeOrigin.Z()};
                        double spacePointChi(0.);
                        
                        CGALPointsToSpacePointsMap::iterator mapItr = cgalPointsToSpacePointsMap.find(vertexItr->point());
                        
                        if (mapItr != cgalPointsToSpacePointsMap.end())
                        {
                            const recob::SpacePoint* spacePoint = mapItr->second;
                            spacePointPos[0] = spacePoint->XYZ()[0];
                            spacePointPos[1] = spacePoint->XYZ()[1];
                            spacePointPos[2] = spacePoint->XYZ()[2];
                            spacePointChi    = spacePoint->Chisq();
                        }
                        else
                            std::cout << "==> Creating infinite space point" << std::endl;
                        
                        //double spacePointPos[] = {vertexItr->point().x() + planeOrigin.X(), vertexItr->point().y() + planeOrigin.Y(), vertexItr->point().z() + planeOrigin.Z()};
                        artSpacePointVector->push_back(recob::SpacePoint(spacePointPos,spacePointErr,spacePointChi,spacePointIdx));
                        
                        // We also want to make associations so we need to recover an art::Ptr to the SpacePoint
                        art::Ptr<recob::SpacePoint> spacePointPtr(spacePointId, spacePointIdx, evt.productGetter(spacePointId));
                        spacePointPtrVec.push_back(spacePointPtr);
                        
                        spacePointIdx++;
                        
                        vertexToSpacePointMap[vertexItr] = &artSpacePointVector->back();
                        
                        numVertices++;
                    }
                    
                    std::cout << "--> Dealaunay point list size: " << delaunayPointList.size() << ", # vertices: " << numVertices << std::endl;
                    
                    // Set up to make the half edges which will also want associations
                    std::vector<art::Ptr<recob::HalfEdge> > halfEdgePtrVec;
                    art::ProductID halfEdgeId = getProductID<std::vector<recob::HalfEdge> >(evt);
                    
                    // Are we using all of the space points?
                    std::set<int> usedSpacePointSet;
                    
                    // repeat this loop...
                    for(auto edgeCirculator = T.tds().edges_begin(); edgeCirculator != T.tds().edges_end(); edgeCirculator++)
                    {
//                    for(auto vertexItr = T.tds().vertices_begin(); vertexItr != T.tds().vertices_end(); vertexItr++)
//                    {
                        // Ignore the "infinite" vertex for now
//                        if (vertexItr == T.tds().vertices_begin()) continue;
                        
                        // Traverse the edges incident to this vertex here?
//                        Edge_circulator edgeCirculator = T.incident_edges(vertexItr);
//                        Edge_circulator done(edgeCirculator);
//                        std::list<EdgeTriple> edgeList;
//                        T.incident_edges(vertexItr, std::back_inserter(edgeList));
                        
//                        for(const auto edgeTriplet : edgeList)
//                        do
//                        {
                            if (edgeCirculator->first->is_valid())
//                            if (edgeTriplet.first->is_valid())
                            {
//                                const Vertex_handle firstVtxHandle  = edgeTriplet.first->vertex(edgeTriplet.second);
//                                const Vertex_handle secondVtxHandle = edgeTriplet.first->vertex(edgeTriplet.third);
                                const Vertex_handle firstVtxHandle  = edgeCirculator->first->vertex(edgeCirculator->first->ccw(edgeCirculator->second));
                                const Vertex_handle secondVtxHandle = edgeCirculator->first->vertex(edgeCirculator->first->cw(edgeCirculator->second));
                                
                                if (firstVtxHandle == T.infinite_vertex() || secondVtxHandle == T.infinite_vertex())
                                {
                                    if (firstVtxHandle  == T.infinite_vertex()) std::cout << "**> Skipping first VTX is infinite" << std::endl;
                                    if (secondVtxHandle == T.infinite_vertex()) std::cout << "**> Skipping second VTX is infinite" << std::endl;
//                                    continue;
                                }
                                
                                const recob::SpacePoint* firstSpacePoint(nullptr);
                                const recob::SpacePoint* secondSpacePoint(nullptr);
                                
                                auto firstVtxToSPItr = vertexToSpacePointMap.find(firstVtxHandle);
                                if (firstVtxToSPItr != vertexToSpacePointMap.end()) firstSpacePoint = firstVtxToSPItr->second;
                                
                                auto secondVtxToSPItr  = vertexToSpacePointMap.find(secondVtxHandle);
                                if (secondVtxToSPItr != vertexToSpacePointMap.end()) secondSpacePoint = secondVtxToSPItr->second;
                                
                                if (!firstSpacePoint || !secondSpacePoint)
                                {
                                    std::cout << " **> Found null vtx pointer! " << firstSpacePoint << ", " << secondSpacePoint << std::endl;
                                    continue;
                                }
                                
                                TVector3 firstPoint( firstSpacePoint->XYZ()[0],  firstSpacePoint->XYZ()[1],  firstSpacePoint->XYZ()[2]);
                                TVector3 secondPoint(secondSpacePoint->XYZ()[0], secondSpacePoint->XYZ()[1], secondSpacePoint->XYZ()[2]);
                                TVector3 deltaPoint = secondPoint - firstPoint;
                                
                                double halfEdgeLength = deltaPoint.Mag();
                                
                                if (halfEdgeLength == 0.)
                                {
                                    std::cout << "**> Edge length: " << halfEdgeLength << ", 1st Vtx pos: " << firstVtxHandle->point() << ", 2nd Vtx pos: " << secondVtxHandle->point() << std::endl;
                                    continue;
                                }
                                
                                usedSpacePointSet.insert(firstSpacePoint->ID());
                                usedSpacePointSet.insert(secondSpacePoint->ID());
                                
                                artHalfEdgeVector->push_back(recob::HalfEdge(halfEdgeLength, firstSpacePoint->ID(), secondSpacePoint->ID(), halfEdgeIdx));
                                
                                // We also want to make associations so we need to recover an art::Ptr to the SpacePoint
                                art::Ptr<recob::HalfEdge> halfEdgePtr(halfEdgeId, halfEdgeIdx, evt.productGetter(halfEdgeId));
                                halfEdgePtrVec.push_back(halfEdgePtr);
                                
                                halfEdgeIdx++;
                            }
                            else
                                std::cout << "**> edge triplet cell was not valid?" << std::endl;
//                        } while(++edgeCirculator != done);
                        
                    }
                    // Make associations to the SpacePoints
                    util::CreateAssn(*this, evt, pfParticle, spacePointPtrVec, *artPFPartSPAssociations);
                    
                    std::cout << "**> used space point set size: " << usedSpacePointSet.size() << std::endl;
/*
                    // Set up to make the half edges which will also want associations
                    std::vector<art::Ptr<recob::HalfEdge> > halfEdgePtrVec;
                    art::ProductID halfEdgeId = getProductID<std::vector<recob::HalfEdge> >(evt);
                    
                    // Traverse all of the edges
                    for(auto edgeTripletItr = T.tds().edges_begin(); edgeTripletItr != T.tds().edges_end(); edgeTripletItr++)
                    {
                        if (edgeTripletItr->first->is_valid())
                        {
                            const Vertex_handle firstVtxHandle  = edgeTripletItr->first->vertex(edgeTripletItr->second);
                            const Vertex_handle secondVtxHandle = edgeTripletItr->first->vertex(edgeTripletItr->third);
                            
                            const recob::SpacePoint* firstSpacePoint(nullptr);
                            const recob::SpacePoint* secondSpacePoint(nullptr);
                            
                            auto firstVtxToSPItr = vertexToSpacePointMap.find(firstVtxHandle);
                            if (firstVtxToSPItr != vertexToSpacePointMap.end()) firstSpacePoint = firstVtxToSPItr->second;
                            
                            auto secondVtxToSPItr  = vertexToSpacePointMap.find(secondVtxHandle);
                            if (secondVtxToSPItr != vertexToSpacePointMap.end()) secondSpacePoint = secondVtxToSPItr->second;
                            
                            if (!firstSpacePoint || !secondSpacePoint)
                            {
                                std::cout << " **> Found null vtx pointer! " << firstSpacePoint << ", " << secondSpacePoint << std::endl;
                                continue;
                            }
                            
                            TVector3 firstPoint(firstVtxHandle->point()[0],firstVtxHandle->point()[1],firstVtxHandle->point()[2]);
                            TVector3 secondPoint(secondVtxHandle->point()[0],secondVtxHandle->point()[1],secondVtxHandle->point()[2]);
                            TVector3 deltaPoint = secondPoint - firstPoint;
                            
                            double halfEdgeLength = deltaPoint.Mag();
                            
                            if (halfEdgeLength > 10.)
                            {
                                std::cout << "**> Edge length: " << halfEdgeLength << ", 1st Vtx pos: " << firstVtxHandle->point() << ", 2nd Vtx pos: " << secondVtxHandle->point() << std::endl;
                                continue;
                            }
                            
                            artHalfEdgeVector->push_back(recob::HalfEdge(halfEdgeLength, firstSpacePoint->ID(), secondSpacePoint->ID(), halfEdgeIdx));
                            
                            // We also want to make associations so we need to recover an art::Ptr to the SpacePoint
                            art::Ptr<recob::HalfEdge> halfEdgePtr(halfEdgeId, halfEdgeIdx, evt.productGetter(halfEdgeId));
                            halfEdgePtrVec.push_back(halfEdgePtr);
                            
                            halfEdgeIdx++;
                        }
                    }
*/
                    // Make associations to the SpacePoints
                    util::CreateAssn(*this, evt, pfParticle, halfEdgePtrVec, *artPFPartHEAssociations);
/*
                    if (T.dimension() == 3 && vertexList.size() > 3)
                    {
                        // Now try the quick hull algorithm
                        Polyhedron_3 convexHull;
                
                        CGAL::convex_hull_3_to_polyhedron_3(T, convexHull);
                
                        std::cout << "   --> 3D convex hull has " << convexHull.size_of_vertices() << " vertices, "
                                  << convexHull.size_of_facets() << " facets, "
                                  << convexHull.size_of_halfedges() << " half edges, "
                                  << convexHull.size_of_border_halfedges() << " border half edges" << std::endl;
                        
                        // Get the vertices/space points... keep a map to go back and forth
                        std::map<Polyhedron_3::Vertex_handle, const recob::SpacePoint*> vertexToSpacePointMap;
                        
                        // Keep track of current start for space points
                        std::vector<art::Ptr<recob::SpacePoint> > spacePointPtrVec;
                        
                        double spacePointErr[] = {1., 0., 1., 0., 0., 1.};
                        double spacePointChi(0.);
                        
                        artSpacePointVector->reserve(convexHull.size_of_vertices());
                        
                        for(Polyhedron_3::Vertex_iterator vertexItr = convexHull.vertices_begin(); vertexItr != convexHull.vertices_end(); vertexItr++)
                        {
                            std::cout << "  --> Vertex degree: " << vertexItr->vertex_degree() << ", point: " << vertexItr->point() << std::endl;
                            
                            double spacePointPos[] = {vertexItr->point().x(), vertexItr->point().y(), vertexItr->point().z()};
                            artSpacePointVector->push_back(recob::SpacePoint(spacePointPos,spacePointErr,spacePointChi,spacePointIdx));
                            
                            // We also want to make associations so we need to recover an art::Ptr to the SpacePoint
                            art::ProductID spacePointId = getProductID<std::vector<recob::SpacePoint> >(evt);
                            art::Ptr<recob::SpacePoint> spacePointPtr(spacePointId, spacePointIdx, evt.productGetter(spacePointId));
                            spacePointPtrVec.push_back(spacePointPtr);
                            
                            spacePointIdx++;
                            
                            vertexToSpacePointMap[vertexItr] = &artSpacePointVector->back();
                        }
                        
                        // Make associations to the SpacePoints
                        util::CreateAssn(*this, evt, pfParticle, spacePointPtrVec, *artPFPartSPAssociations);
                        
                        // Set up to make the half edges which will also want associations
                        std::vector<art::Ptr<recob::HalfEdge> > halfEdgePtrVec;
                        art::ProductID halfEdgeId = getProductID<std::vector<recob::HalfEdge> >(evt);
                        
                        for(Polyhedron_3::Facet_iterator facetItr = convexHull.facets_begin(); facetItr != convexHull.facets_end(); facetItr++)
                        {
                            //const Polyhedron_3::Facet_handle& facetHandle = *facetItr;
                            
                            size_t facetDegree = facetItr->facet_degree();
                            
                            std::cout << "   == Facet degree: " << facetDegree;
                            
                            Polyhedron_3::Halfedge_handle halfEdge_handle = facetItr->halfedge();
                            
                            do
                            {
                                std::cout << ", vertex degree: " << halfEdge_handle->vertex_degree();
                                
                                Polyhedron_3::Vertex_handle firstVertex = halfEdge_handle->prev()->vertex();
                                Polyhedron_3::Vertex_handle nextVertex  = halfEdge_handle->vertex();
                                
                                const recob::SpacePoint* firstSpacePoint(nullptr);
                                const recob::SpacePoint* nextSpacePoint(nullptr);
                                
                                auto firstVtxToSPItr = vertexToSpacePointMap.find(firstVertex);
                                if (firstVtxToSPItr != vertexToSpacePointMap.end()) firstSpacePoint = firstVtxToSPItr->second;
                                
                                auto nextVtxToSPItr  = vertexToSpacePointMap.find(nextVertex);
                                if (nextVtxToSPItr != vertexToSpacePointMap.end()) nextSpacePoint = nextVtxToSPItr->second;
                                
                                if (firstSpacePoint && nextSpacePoint)
                                {
                                    std::cout << ", vtx ids: " << firstSpacePoint->ID() << ", " << nextSpacePoint->ID();
                                }
                                else
                                {
                                    std::cout << ", vtx ptrs: " << firstSpacePoint << ", " << nextSpacePoint;
                                }
                                
                                TVector3 firstPoint(firstSpacePoint->XYZ()[0],firstSpacePoint->XYZ()[1],firstSpacePoint->XYZ()[2]);
                                TVector3 nextPoint( nextSpacePoint->XYZ()[0], nextSpacePoint->XYZ()[1], nextSpacePoint->XYZ()[2]);
                                TVector3 deltaPoint = nextPoint - firstPoint;
                                
                                double halfEdgeLength = deltaPoint.Mag();
                                
                                artHalfEdgeVector->push_back(recob::HalfEdge(halfEdgeLength, firstSpacePoint->ID(), nextSpacePoint->ID(), halfEdgeIdx));
                                
                                // We also want to make associations so we need to recover an art::Ptr to the SpacePoint
                                art::Ptr<recob::HalfEdge> halfEdgePtr(halfEdgeId, halfEdgeIdx, evt.productGetter(halfEdgeId));
                                halfEdgePtrVec.push_back(halfEdgePtr);
                                
                                halfEdgeIdx++;
                                
                                halfEdge_handle = halfEdge_handle->next();
                            } while(halfEdge_handle != facetItr->halfedge());
                            
                            std::cout << std::endl;
                        }
                        
                        // Make associations to the SpacePoints
                        util::CreateAssn(*this, evt, pfParticle, halfEdgePtrVec, *artPFPartHEAssociations);
                    }
 */
                }
            }
            
        }
    }
    
    // Finaly done, now output everything to art
    evt.put(std::move(artSpacePointVector));
    evt.put(std::move(artHalfEdgeVector));
    evt.put(std::move(artPFPartSPAssociations));
    evt.put(std::move(artPFPartHEAssociations));
    evt.put(std::move(artSpacePointHEAssociations));
    
    // If monitoring then deal with the fallout
    if (fEnableMonitoring)
    {
        theClockTotal.stop();
        
        fRun                   = evt.run();
        fEvent                 = evt.id().event();
        fTotalTime             = theClockTotal.accumulated_real_time();
//        m_hits                  = static_cast<int>(clusterHit2DMasterVec.size());
        fRecoTree->Fill();
        
        mf::LogDebug("CGALTest") << "*** CGALTest total time: " << fTotalTime << std::endl;
    }
    
    // Will we ever get here? ;-)
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CGALTest::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    fRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    fRecoTree->Branch("run",                  &fRun,                   "run/I");
    fRecoTree->Branch("event",                &fEvent,                 "event/I");
    fRecoTree->Branch("hits",                 &fHits,                  "hits/I");
    fRecoTree->Branch("totalTime",            &fTotalTime,             "time/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CGALTest::PrepareEvent(const art::Event &evt)
{
    fRun                   = evt.run();
    fEvent                 = evt.id().event();
    fHits                  = 0;
    fTotalTime             = 0.f;
}

} // namespace lar_CGALTest
