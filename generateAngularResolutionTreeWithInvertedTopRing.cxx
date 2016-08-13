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
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;
  const Double_t maxDeltaTriggerTimeNs = 1200;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;

  CrossCorrelator* cc = new CrossCorrelator();  
  // CrossCorrelator::SimpleNotch notchWideLow("nWideLow", "260MHz Satellite And 200MHz Notch Notch", 0, 286);
  // CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch", 344, 390);
  // cc->addNotch(notchWideLow);
  // cc->addNotch(notch370);
  cc->multiplyTopRingByMinusOne = 1;
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);

  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);
  
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", lindaFileName.Data());
  // TNamed* lindaFileNameReference = new TNamed("lindaFileNameReference", "photogrammetry");  

  // lindaFileNameReference->Write();

  TTree* angResTree = new TTree("angResTree", "angResTree");

  Double_t globalPeak = 0;
  Double_t globalPhiDeg = 0;
  Double_t globalThetaDeg = 0;
  Double_t triggeredPeak = 0;
  Double_t triggeredPhiDeg = 0;
  Double_t triggeredThetaDeg = 0;
  Double_t zoomPeak = 0;
  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t deltaThetaDeg = 0;
  Double_t deltaPhiDeg = 0;

  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  UInt_t eventNumber = 0;
  UInt_t triggerTimeNs = 0;
  UInt_t triggerTimeNsExpected = 0;
  Double_t heading = 0;
  UInt_t l3TrigPattern = 0;
  UInt_t l3TrigPatternH = 0;
  Int_t run = 0;
  Double_t mrms = 0;
  Double_t brms = 0;
  UInt_t realTime = 0;

  
  angResTree->Branch("globalPeak", &globalPeak);
  angResTree->Branch("globalPhiDeg", &globalPhiDeg);
  angResTree->Branch("globalThetaDeg", &globalThetaDeg);

  angResTree->Branch("triggeredPeak", &triggeredPeak);
  angResTree->Branch("triggeredPhiDeg", &triggeredPhiDeg);
  angResTree->Branch("triggeredThetaDeg", &triggeredThetaDeg);

  angResTree->Branch("zoomPeak", &zoomPeak);
  angResTree->Branch("zoomPhiDeg", &zoomPhiDeg);
  angResTree->Branch("zoomThetaDeg", &zoomThetaDeg);
  angResTree->Branch("deltaPhiDeg", &deltaPhiDeg);
  angResTree->Branch("deltaThetaDeg", &deltaThetaDeg);
  
  angResTree->Branch("thetaExpected", &thetaExpected);
  angResTree->Branch("phiExpected", &phiExpected);
  angResTree->Branch("triggerTimeNs", &triggerTimeNs);
  angResTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);
  angResTree->Branch("heading", &heading);
  angResTree->Branch("l3TrigPattern", &l3TrigPattern);
  angResTree->Branch("l3TrigPatternH", &l3TrigPatternH);
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("run", &run);

  angResTree->Branch("mrms", &mrms);
  angResTree->Branch("brms", &brms);
  angResTree->Branch("realTime", &realTime);

  Double_t maxDPhiDeg;
  angResTree->Branch("maxDPhiDeg", &maxDPhiDeg);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  Int_t numSaved = 0;
  Int_t maxToSave = 5;
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    

    headChain->GetEntry(entry);
    gpsChain->GetEntry(entry);        
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
      // if(true){
	eventNumber = header->eventNumber;
	run = header->run;
	
	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	mrms = usefulPat.mrms;
	brms = usefulPat.brms;
	realTime = header->realTime;
	l3TrigPattern = header->l3TrigPattern;
	l3TrigPatternH = header->l3TrigPatternH;
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected*=TMath::RadToDeg();
	thetaExpected*=-1*TMath::RadToDeg();
	
	phiExpected = phiExpected < 0 ? phiExpected + 360 : phiExpected;
	phiExpected = phiExpected >= 360 ? phiExpected - 360 : phiExpected;

	cc->correlateEvent(usefulEvent, pol);
	
	TH2D* hGlobalImageH = cc->makeGlobalImage(pol, globalPeak, globalPhiDeg, globalThetaDeg);
	hGlobalImageH->SetName(TString::Format("hGlobalImageH_%u", eventNumber));

	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						  zoomThetaDeg,
						  globalPhiDeg, globalThetaDeg);
						  // triggeredPhiDeg, triggeredThetaDeg);

	zoomPhiDeg = zoomPhiDeg < 0 ? zoomPhiDeg + 360 : zoomPhiDeg;
	zoomPhiDeg = zoomPhiDeg >= 360 ? zoomPhiDeg - 360 : zoomPhiDeg;

	deltaPhiDeg = RootTools::getDeltaAngleDeg(phiExpected, zoomPhiDeg);
	deltaThetaDeg = RootTools::getDeltaAngleDeg(thetaExpected, zoomThetaDeg);

	TString name1 = TString::Format("hTriggeredZoomImageH%u_1", eventNumber);
	hZoomedImageH->SetName(name1);

	
	if(numSaved < maxToSave){
	  numSaved++;
	}
	else{
	  delete hGlobalImageH;
	  delete hZoomedImageH;
	}

	angResTree->Fill();

	delete usefulEvent;
      }
    }
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
