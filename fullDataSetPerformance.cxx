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
  // TChain* gpsChain = new TChain("adu5PatTree");
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


  TH1D* hPhi = new TH1D("hPhi", "#phi", 1024, -180, 180);
  TH1D* hTheta = new TH1D("hTheta", "#theta", 1024, -90, 90);
  TH1D* hPeak = new TH1D("hPeak", "peak", 1024, -1, 1);
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    
    angResChain->GetEntry(entry);

    headChain->GetEntryWithIndex(eventNumber);

    Float_t priPeak = header->getImagePeak();
    Float_t priPhiDeg = header->getPeakPhiDeg();
    Float_t priThetaDeg = header->getPeakThetaDeg();
    Int_t priority = header->priority;
    Float_t hilbertPeak = header->getPeakPhiRad();
    
    // std::cerr << entry << "\t" << eventNumber << "\t" << priPeak << "\t" << priPhiDeg << "\t" << priThetaDeg << "\t" << std::endl;

    // std::cerr << entry << "\t" << eventNumber << "\t" << zoomPeak << "\t" << zoomPhiDeg << "\t" << zoomThetaDeg << "\t" << std::endl;    
    
    hPhi->Fill(RootTools::getDeltaAngleDeg((Double_t)priPhiDeg, zoomPhiDeg));
    hTheta->Fill((Double_t)priThetaDeg - zoomThetaDeg);
    hPeak->Fill((Double_t)priPeak - zoomPeak);    
        
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
