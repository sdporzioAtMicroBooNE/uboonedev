
#include "HitAnalysisAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include <cmath>
#include <algorithm>

namespace hitanalysis
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
HitAnalysisAlg::HitAnalysisAlg(fhicl::ParameterSet const & pset) 
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("HitAnalysisAlg") << "HitAnalysisAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
HitAnalysisAlg::~HitAnalysisAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void HitAnalysisAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fLocalDirName = pset.get<std::string>("LocalDirName", std::string("wow"));
}

//----------------------------------------------------------------------------
/// Begin job method.
void HitAnalysisAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs, const std::string& dirName)
{
    // Make a directory for these histograms
    art::TFileDirectory dir = tfs->mkdir(dirName.c_str());

    fHitsByWire[0]            = dir.make<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]            = dir.make<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]            = dir.make<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fPulseHeight[0]           = dir.make<TH1D>("PulseHeight0",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeight[1]           = dir.make<TH1D>("PulseHeight1",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeight[2]           = dir.make<TH1D>("PulseHeight2",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[0]     = dir.make<TH1D>("PulseHeightS0", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[1]     = dir.make<TH1D>("PulseHeightS1", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[2]     = dir.make<TH1D>("PulseHeightS2", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[0]      = dir.make<TH1D>("PulseHeightM0", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[1]      = dir.make<TH1D>("PulseHeightM1", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[2]      = dir.make<TH1D>("PulseHeightM2", "PH (ADC)",  300,  0.,  150.);
    fChi2DOF[0]               = dir.make<TH1D>("Chi2DOF0",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[1]               = dir.make<TH1D>("Chi2DOF1",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[2]               = dir.make<TH1D>("Chi2DOF2",      "Chi2DOF",   502, -1.,  250.);
    fNumDegFree[0]            = dir.make<TH1D>("NumDegFree0",   "NDF",       100,  0.,  100.);
    fNumDegFree[1]            = dir.make<TH1D>("NumDegFree1",   "NDF",       100,  0.,  100.);
    fNumDegFree[2]            = dir.make<TH1D>("NumDegFree2",   "NDF",       100,  0.,  100.);
    fChi2DOFSingle[0]         = dir.make<TH1D>("Chi2DOFS0",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[1]         = dir.make<TH1D>("Chi2DOFS1",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[2]         = dir.make<TH1D>("Chi2DOFS2",     "Chi2DOF",   502, -1.,  250.);
    fHitMult[0]               = dir.make<TH1D>("HitMult0",      "# hits",     15,  0.,   15.);
    fHitMult[1]               = dir.make<TH1D>("HitMult1",      "# hits",     15,  0.,   15.);
    fHitMult[2]               = dir.make<TH1D>("HitMult2",      "# hits",     15,  0.,   15.);
    fHitCharge[0]             = dir.make<TH1D>("HitCharge0",    "Charge",   1000,  0., 2000.);
    fHitCharge[1]             = dir.make<TH1D>("HitCharge1",    "Charge",   1000,  0., 2000.);
    fHitCharge[2]             = dir.make<TH1D>("HitCharge2",    "Charge",   1000,  0., 2000.);
    fFitWidth[0]              = dir.make<TH1D>("FitWidth0",     "Width",     100,  0.,   10.);
    fFitWidth[1]              = dir.make<TH1D>("FitWidth1",     "Width",     100,  0.,   10.);
    fFitWidth[2]              = dir.make<TH1D>("FitWidth2",     "Width",     100,  0.,   10.);
    fHitSumADC[0]             = dir.make<TH1D>("SumADC0",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[1]             = dir.make<TH1D>("SumADC1",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[2]             = dir.make<TH1D>("SumADC2",       "Sum ADC",  1000,  0., 2000.);
    
    fNDFVsChi2[0]             = dir.make<TH2D>("NDFVsChi20",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    fNDFVsChi2[1]             = dir.make<TH2D>("NDFVsChi21",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    fNDFVsChi2[2]             = dir.make<TH2D>("NDFVsChi22",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    
    fPulseHVsWidth[0]         = dir.make<TH2D>("PHVsWidth0",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fPulseHVsWidth[1]         = dir.make<TH2D>("PHVsWidth1",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fPulseHVsWidth[2]         = dir.make<TH2D>("PHVsWidth2",    ";PH;Width", 100,  0.,  100., 100,  0., 10.);

    return;
}
    
void HitAnalysisAlg::fillHistograms(const HitPtrVec& hitPtrVec) const
{
    // Keep track of number of hits per view
    size_t nHitsPerView[] = {0,0,0};
    
    // Loop the hits and make some plots
    for(const auto& hitPtr : hitPtrVec)
    {
        // Extract interesting hit parameters
        const geo::WireID& wireID   = hitPtr->WireID();
        float              chi2DOF  = std::min(hitPtr->GoodnessOfFit(),float(249.8));
        int                numDOF   = hitPtr->DegreesOfFreedom();
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
        
        nHitsPerView[view]++;
        
        fHitsByWire[view]->Fill(wire,1.);
        fPulseHeight[view]->Fill(hitPH, 1.);
        fChi2DOF[view]->Fill(chi2DOF, 1.);
        fNumDegFree[view]->Fill(numDOF, 1.);
        fHitMult[view]->Fill(hitMult, 1.);
        fHitCharge[view]->Fill(charge, 1.);
        fFitWidth[view]->Fill(std::min(float(9.99),hitSigma), 1.);
        fHitSumADC[view]->Fill(sumADC, 1.);
        fNDFVsChi2[view]->Fill(numDOF, chi2DOF, 1.);
//        fPulseHVsWidth[view]->Fill(std::min(float(99.9),hitPH),std::min(float(9.99),hitSigma), 1.);
        
        if (hitMult == 1)
        {
            fPulseHeightSingle[view]->Fill(hitPH, 1.);
            fChi2DOFSingle[view]->Fill(chi2DOF, 1.);
            fPulseHVsWidth[view]->Fill(std::min(float(99.9),hitPH),std::min(float(9.99),hitSigma), 1.);
        }
        else
            fPulseHeightMulti[view]->Fill(hitPH, 1.);
        
//        if (hitPH < 15. && hitSigma > 9.) std::cout << "==> Hit view: " << view << ", wire: " << wire << ", time: " << hitPtr->PeakTime() << " has PW/PH: " << hitSigma << "/" << hitPH << ", Mult: " << hitMult << ", chi: " << chi2DOF << std::endl;
    }
    
//    std::cout << "** " << fLocalDirName << " hits, # per plane: " << nHitsPerView[0] << " / " << nHitsPerView[1] << " / " << nHitsPerView[2] << std::endl;
    
    return;
}
    
// Useful for normalizing histograms
void HitAnalysisAlg::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
    double normFactor(1./numEvents);
    
    for(size_t idx = 0; idx < 3; idx++) fHitsByWire[idx]->Scale(normFactor);
    
    return;
}

}
