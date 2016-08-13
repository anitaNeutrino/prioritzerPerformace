// A simple script to hunt for Seavey pulses.
void pulseHuntingTriggerTime(){

  Double_t sourceLat = - (79 + (27.93728/60));
  Double_t sourceLon = -(112 + (6.74974/60));
  Double_t sourceAlt = 1813.42;

  Double_t cutTimeNs = 60e3;
  // Double_t cutTimeLowNs = -1e3;
  // Double_t cutTimeHighNs = -8e2;
  //  Double_t cutTimeNs = 2e3;

  gSystem->Load("libMinuit.so");  
  gSystem->Load("libPhysics.so");  
  gSystem->Load("libAnitaEvent.so");
  gSystem->Load("libAnitaCorrelator.so");

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  for(int run=331; run<355; run++){
  //  for(int run=340; run<341; run++){
    char fileName[1024];
    //sprintf(fileName, "root/run%d/eventHeadFile%d.root", run, run);
    sprintf(fileName, "root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    //sprintf(fileName, "root/run%d/gpsEvent%d.root", run, run);
    sprintf(fileName, "root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("pat->realTime");

  Long64_t numHeaderEntries = headChain->GetEntries();
  cout << numHeaderEntries << " " << gpsChain->GetEntries() << endl;

  TH2D* h0 = new TH2D("h0", "Heading as a function of time, runs 331-354; Entry; Heading (degrees)", 1024, 1420019546, 1420192158, 360, 0, 360);
  TH2D* h1 = new TH2D("h1", "triggerTimeNs vs TriggerTimeNsExpected WAIS pulser runs 331-354; Time since second start (ns); Expected time since second start (ns)", 64, 1e6, 2e6, 64, 1e6, 2e6);
  TH1D* h2 = new TH1D("h2", "TriggerTimeNsExpected - triggerTimeNs WAIS pulser runs 331-354; #deltat (ns); Fraction of events", 10000, -1e9, 1e9);
  TH1D* h3 = new TH1D("h3", "TriggerTimeNsExpected - triggerTimeNs WAIS pulser runs 331-354; #deltat (ns); Fraction of events", 10000, -2e5, 2e5);
  TH1D* h4 = new TH1D("h4", "Priority; Priority; Fraction of events", 10, 0, 10);
  TH1D* h5 = new TH1D("h5", "Priorities near WAIS divide, runs 331-354; Priority; Fraction of events", 10, 0, 10);
  TH2D* h6 = new TH2D("h6", "Phi-sector trigger HPOL vs phi-expected near WAIS divide, runs 331-354; Phi-Sector; Phi expected (degrees)", 16, 0, 16, 16, 0, 360);
  TH2D* h7 = new TH2D("h7", "GPU phi vs Expected phi near WAIS divide, runs 331-354; GPU peak phi (Degrees); Phi expected (degrees)", 32, 0, 360, 32, 0, 360);
  //  TH2D* h8 = new TH2D("h8", "GPU phi - Expected phi near WAIS divide, runs 331-354; GPU peak phi - Phi expected (degrees); Fraction of events", 256, 0, numHeaderEntries, 360, -180, 180);
  TH1D* h8 = new TH1D("h8", "GPU phi - Expected phi (+45) near WAIS divide, runs 331-354; GPU peak phi - Phi expected (degrees); Fraction of events", 360, -180, 180);
  TH2D* h9 = new TH2D("h9", "Heading - Expected phi near WAIS divide, runs 331-354; Entry; Heading - Phi expected (degrees)", 1024, 0, numHeaderEntries, 360, -180, 180);
  TH2D* h10 = new TH2D("h10", "L3 trigger HPOL vs GPU Peak PhiSector, runs 331-354; l3trigger; GPU Peak Phi", 16, 0, 16, 720, 0, 360);
  TH1D* h11 = new TH1D("h11", "GPU theta - Expected theta, near WAIS divide, runs 331-354; GPU Peak Theta - Theta expected (degrees); Fraction of events", 512, -150, 150);
  TH2D* h12 = new TH2D("h12", "GPU theta vs Expected theta, near WAIS divide, runs 331-354; GPU Peak Theta (degrees); Expected theta (degrees)", 150, -75, 75, 150, -100, 100);
  for(Long64_t entry = 0; entry < numHeaderEntries; entry++){
    headChain->GetEntry(entry);

    if(entry==0 || entry == numHeaderEntries-1){
      cout << header->realTime << endl;
    }

    //    cout << entry << "\t" << numHeaderEntries << " " << "hello! this is text" << endl;

    if(header->trigType & 1){

      Long64_t gpsEntry = gpsChain->GetEntryNumberWithIndex(header->realTime);
      if(gpsEntry < 0 ) continue;
      // gpsChain->GetEntry(gpsEntry);
      if(gpsEntry >= 0 ){
	gpsChain->GetEntry(gpsEntry);
      }

      UsefulAdu5Pat usefulPat(pat);

      // cout << sourceLat << " " << sourceLon << endl;
      // cout << pat->latitude << " " << pat->longitude << endl;

      UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      //99756.6;
      // cout << triggerTimeNsExpected << "\t";
      // triggerTimeNsExpected += -99756.6;
      // cout << triggerTimeNsExpected << endl;
      h0->Fill(header->realTime, pat->heading);
      h1->Fill(header->triggerTimeNs, triggerTimeNsExpected);
      h2->Fill(Double_t(triggerTimeNsExpected) - Double_t(header->triggerTimeNs));
      h3->Fill(Double_t(triggerTimeNsExpected) - Double_t(header->triggerTimeNs));

      if(TMath::Abs(Double_t(header->triggerTimeNs) - Double_t(triggerTimeNsExpected)) < cutTimeNs){

	h4->Fill(header->priority&0xf);
	// cout << header->l3TrigPatternH << ": ";
	Double_t thetaWave = -1;
	Double_t phiWave = -1;
	usefulPat.getThetaAndPhiWaveAnita3(sourceLon, sourceLat, sourceAlt, thetaWave, phiWave);
	//usefulPat.getThetaAndPhiWave(pat->longitude, pat->latitude, 0, thetaWave, phiWave);
	Double_t thetaDeg = thetaWave*TMath::RadToDeg(); 
	Double_t phiDeg = phiWave*TMath::RadToDeg();
	phiDeg = phiDeg < 0 ? phiDeg + 360 : phiDeg;
	phiDeg = phiDeg >= 360 ? phiDeg - 360 : phiDeg;
	for(int phi=0; phi<16; phi++){
	  Int_t val = ((header->l3TrigPatternH >> phi) & 1);
	  if(val){
	    h6->Fill(phi, phiDeg);

	    //Int_t phiInd2=(header->prioritizerStuff&0x8fff)>>1;
	    Int_t phiInd2=(header->prioritizerStuff&0x7ff)>>1;
	    Int_t phiInd=(phiInd2&0x3f);
	    Int_t phiSector=(phiInd2>>6);
	    if(phiSector > 15){
	      cout << phiSector << " " << phiInd << " " << phiInd2 << endl;
	    }

	    Double_t gpuPhiDeg = 22.5*phiSector - 11.5 + 22.5*phiInd/64;
	    gpuPhiDeg = gpuPhiDeg < 0 ? gpuPhiDeg + 360 : gpuPhiDeg;

	    Double_t gpuPhiDeg2 = header->getPeakPhiDeg();
	    gpuPhiDeg2 = gpuPhiDeg2 < 0 ? gpuPhiDeg2 + 360 : gpuPhiDeg2;

	    h10->Fill(phi, gpuPhiDeg2);



	  }
	}
	Double_t gpuPhi = header->getPeakPhiDeg();
	if(gpuPhi < 0){
	  // cout << gpuPhi << endl;
	  gpuPhi += 360;
	}

	h7->Fill(header->getPeakPhiDeg(), phiDeg);

	Double_t deltaPhiGPUvsExpected = gpuPhi - phiDeg;
	deltaPhiGPUvsExpected = deltaPhiGPUvsExpected < -180 ? deltaPhiGPUvsExpected + 360 : deltaPhiGPUvsExpected;
	deltaPhiGPUvsExpected = deltaPhiGPUvsExpected > 180 ? deltaPhiGPUvsExpected - 360 : deltaPhiGPUvsExpected;
	//	h8->Fill(entry, deltaPhiGPUvsExpected);
	h8->Fill(deltaPhiGPUvsExpected);
	h11->Fill(header->getPeakThetaDeg() - thetaDeg);
	h12->Fill(header->getPeakThetaDeg(), thetaDeg);

	Double_t deltaHeadingvsPhiExpected = pat->heading - phiDeg;
	deltaHeadingvsPhiExpected = deltaHeadingvsPhiExpected < -180 ? deltaHeadingvsPhiExpected + 360 : deltaHeadingvsPhiExpected;
	deltaHeadingvsPhiExpected = deltaHeadingvsPhiExpected > 180 ? deltaHeadingvsPhiExpected - 360 : deltaHeadingvsPhiExpected;
	h9->Fill(entry, deltaHeadingvsPhiExpected);
      }
      else{
	h5->Fill(header->priority&0xf);
      }
    }
  }

  // TCanvas* c0 = new TCanvas();
  // h0->Draw("colz");

  // TCanvas* c1 = new TCanvas();
  //  gpsChain->Draw("pat->heading:pat->realTime","pat->heading > -500","colz");
  // htemp->SetLineColor(kBlue);
  // headChain->Draw("realTime", "", "same");
  // htemp->SetLineColor(kRed);

  TCanvas* c2 = new TCanvas();
  h1->Draw("colz");
  c1->SetLogz(1);

  TCanvas* c3 = new TCanvas();
  // h3->SetLineColor(kBlue);
  // h3->DrawNormalized();
  h2->SetLineColor(kRed);
  h2->DrawNormalized(); //"same");
  c3->SetLogy(1);

  TCanvas* c4 = new TCanvas();
  h3->SetLineColor(kBlue);
  h3->DrawNormalized();
  c4->SetLogy(1);


  const Int_t nx = 10;
  char *pris[nx] = {"0","1","2","3","4","5","6", "7","8","9"};
  for (int i=1;i<=nx;i++) {
    h5->GetXaxis()->SetBinLabel(i,pris[i-1]);
  }
  h5->GetXaxis()->SetLabelOffset(0.005);
  h5->GetXaxis()->SetLabelSize(0.04);

  TCanvas* c5 = new TCanvas();
  TLegend* l1 = new TLegend(0.7, 0.8, 0.975, 0.975);
  h5->SetLineColor(kRed);
  h5->DrawNormalized("e hist");


  h4->SetLineColor(kBlue);
  h4->DrawNormalized("same e hist");
  l1->AddEntry(h4, "WAIS pulses", "l");
  l1->AddEntry(h5, "Other events", "l");
  l1->Draw();



  TCanvas* c6 = new TCanvas();
  h6->Draw("colz");

  TCanvas* c7 = new TCanvas();
  h7->Draw("colz");

  TCanvas* c8 = new TCanvas();
  h8->DrawNormalized();
  //h8->Draw("colz");

  TCanvas* c9 = new TCanvas();
  h9->Draw("colz");

  TCanvas* c10 = new TCanvas();
  h10->Draw("colz");

  TCanvas* c11 = new TCanvas();
  h11->DrawNormalized();

  TCanvas* c12 = new TCanvas();
  h12->Draw("colz");

}
