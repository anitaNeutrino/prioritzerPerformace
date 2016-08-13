// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to reconstruct HPol pulses from Wais Divide.
********************************************************************************************************* */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "THnSparse.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"
#include "CrossCorrelator.h"
#include "OutputConvention.h"


int main(int argc, char *argv[])
{

  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }
  
  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = atoi(argv[2]);; //argc==3 ? atoi(argv[2]) : firstRun;

  
  TChain* headChain = new TChain("headTree");
  TChain* angResChain = new TChain("angResTree");
  for(Int_t run=firstRun; run<=lastRun; run++){

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/antennaPositionCalib/final/filterOffs/phaseCenter/generateAngularResolutionTreePlots_%d_*.root", run);
    angResChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  UInt_t eventNumber = 0;
  angResChain->SetBranchAddress("eventNumber", &eventNumber);
  Double_t zoomPhiDeg = 0;
  angResChain->SetBranchAddress("zoomPhiDeg", &zoomPhiDeg);
  Double_t zoomThetaDeg = 0;
  angResChain->SetBranchAddress("zoomThetaDeg", &zoomThetaDeg);
  Double_t phiExpected = 0;
  angResChain->SetBranchAddress("phiExpected", &phiExpected);
  Double_t thetaExpected = 0;
  angResChain->SetBranchAddress("thetaExpected", &thetaExpected);
  
  Double_t zoomPeak = 0;
  angResChain->SetBranchAddress("zoomPeak", &zoomPeak);
  
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  headChain->BuildIndex("eventNumber");

  Long64_t nEntries = angResChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  TH2D* hPhiTheta = new TH2D("hPhiTheta", "#phi_{GPU} - #phi_{reco} vs. #theta_{GPU} - #theta_{reco}; #delta#phi (Degrees); #delta#theta (Degrees)", 256, -180, 180, 256, -90, 90);
  TH2D* hPhiThetaClose = new TH2D("hPhiThetaClose", "#phi_{GPU} - #phi_{reco} vs. #theta_{GPU} - #theta_{reco}; #delta#phi (Degrees); #delta#theta (Degrees); #delta#phi (Degrees)", 64, -15, 15, 64, -15, 15);
  TH1D* hPhi = new TH1D("hPhi", "#phi_{GPU} - #phi_{reco}; #delta#phi (Degrees); Events per bin", 256, -180, 180);
  TH1D* hTheta = new TH1D("hTheta", "#theta_{GPU} - #theta_{reco}; #delta#theta (Degrees); Events per bin", 256, -90, 90);
  TH1D* hPeak = new TH1D("hPeak", "P_{GPU} - P_{reco}; #deltaP (no units); Events per bin", 256, -1, 1);


  TH2D* hPhiTheta2 = new TH2D("hPhiTheta2", "#phi_{GPU} - #phi_{WAIS} vs. #theta_{GPU} - #theta_{WAIS}; #delta#phi (Degrees); #delta#theta (Degrees)", 256, -180, 180, 256, -90, 90);
  TH2D* hPhiThetaClose2 = new TH2D("hPhiThetaClose2", "#phi_{GPU} - #phi_{WAIS} vs. #theta_{GPU} - #theta_{WAIS}; #delta#phi (Degrees); #delta#theta (Degrees); #delta#phi (Degrees)", 64, -15, 15, 64, -15, 15);
  TH1D* hPhi2 = new TH1D("hPhi2", "#phi_{GPU} - #phi_{WAIS}; #delta#phi (Degrees); Events per bin", 256, -180, 180);
  TH1D* hTheta2 = new TH1D("hTheta2", "#theta_{GPU} - #theta_{WAIS}; #delta#theta (Degrees); Events per bin", 256, -90, 90);
  
  const int maxPriority = 9;
  const int minPri = 1;

  TH1D* hPriority = new TH1D("hPriority", "priority", maxPriority-minPri+1, minPri, maxPriority+1);
  std::vector<TH2D*> hImagePeakHilbertPeaks;
  for(int pri=minPri; pri <= maxPriority; pri++){

    TString title = TString::Format("Priority %d", pri);
    TString name = TString::Format("hImagePeakHilbertPeak_%d", pri);
    TH2D* hImagePeakHilbertPeak = new TH2D(name, title, 1024, 0, 1, 256, 0, 256);

    hImagePeakHilbertPeak->SetLineColor(pri+1);
    hImagePeakHilbertPeak->SetFillColor(pri+1);

    hImagePeakHilbertPeaks.push_back(hImagePeakHilbertPeak);
  }
  auto geom = AnitaGeomTool::Instance();
  double aftFore = geom->aftForeOffsetAngleVerticalKurtAnita3*TMath::RadToDeg();  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    
    angResChain->GetEntry(entry);
    headChain->GetEntryWithIndex(eventNumber);

    Float_t priPeak = header->getImagePeak();
    Float_t priPhiDeg = header->getPeakPhiDeg();
    Float_t priThetaDeg = header->getPeakThetaDeg();
    Int_t priority = header->priority & 0xf;
    Float_t coherentPeak = header->getCoherentSumPeak();
    if(priority > maxPriority || priority < minPri){
      std::cerr << priority << "?" << std::endl;
    }

    Double_t deltaPhi = RootTools::getDeltaAngleDeg((Double_t)(priPhiDeg)-aftFore, zoomPhiDeg);
    Double_t deltaTheta = (Double_t)priThetaDeg - zoomThetaDeg;

    Double_t deltaPhi2 = RootTools::getDeltaAngleDeg((Double_t)(priPhiDeg)-aftFore, phiExpected);
    Double_t deltaTheta2 = (Double_t)priThetaDeg - thetaExpected;

    Double_t deltaPeak = (Double_t)priPeak - zoomPeak;    
    hImagePeakHilbertPeaks.at(priority-minPri)->Fill(priPeak, (Double_t)coherentPeak);

    hPriority->Fill(priority);
    
    hPhi->Fill(deltaPhi);
    hTheta->Fill(deltaTheta);

    if(TMath::Abs(RootTools::getDeltaAngleDeg(zoomPhiDeg, phiExpected)) < 5){
      hPeak->Fill(deltaPeak);
    }
    hPhiTheta->Fill(deltaPhi, deltaTheta);
    hPhiThetaClose->Fill(deltaPhi, deltaTheta); 

    hPhi2->Fill(deltaPhi2);
    hTheta2->Fill(deltaTheta2);
    hPhiTheta2->Fill(deltaPhi2, deltaTheta2);
    hPhiThetaClose2->Fill(deltaPhi2, deltaTheta2); 
    
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}








