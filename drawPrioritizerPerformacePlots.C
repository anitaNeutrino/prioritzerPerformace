{
  /* Draws plots for lovely presenations */

  TFile* f = TFile::Open("prioritizerdPerformancePlots.root");

  gROOT->ProcessLine(".L sc.C");

  TCanvas* c1 = new TCanvas();
  //  gStyle->SetOptStat("emr");
  h5->GetXaxis()->SetLabelOffset(0.007);
  h5->GetXaxis()->SetLabelSize(0.04);
  h5->SetLineColor(kBlue);
  h5->DrawNormalized("ehist");
  h4->SetLineColor(kRed);
  h4->DrawNormalized("sameehist");
  h5->SetTitle("Event Priorities runs 330-356");



  TLegend *l = new TLegend(0.75, 0.75, 0.975, 0.975);
  l->AddEntry(h4, "WAIS pulses", "l");
  l->AddEntry(h5, "Other events", "l");
  l->Draw();

  TCanvas* c2 = new TCanvas();
  THStack* hStack = new THStack("hs","Stacked histogram of events runs 330-356; Priority; Number of events");

  hStack->Add(h5);
  hStack->Add(h4);


  hStack->Draw();

  hStack->GetXaxis()->SetLabelOffset(0.005);
  hStack->GetXaxis()->SetLabelSize(0.04);

  l->Draw();

  c2->SetLogy(1);
  c2->Update();


  TCanvas* c3 = new TCanvas();
  gStyle->SetOptStat("e");
  h14->Draw("colz");
  
  TCanvas* c4 = new TCanvas();
  h13->Divide(h16);
  h13->SetTitle("Fraction of telemetered RF triggered events that were WAIS pulses; Run; Fraction of events per run");
  h13->Draw("ehist");


  TCanvas* c5 = new TCanvas();
  TH1D* hTots = h5->Clone("hTots");
  hTots->Add(h4);

  TH1D* hRatio = h4->Clone("hRatio");
  hRatio->SetTitle("Fraction of events in each Priority bin that were WAIS pulses runs 330-356; Priority; Fraction of events");
  hRatio->Divide(hTots);
  hRatio->Draw();

  const Int_t nx = 10;
  const char *pris[nx] = {"0","1","2","3","4","5","6","7","8","9"};
  for (int i=1;i<=nx;i++) {
    hRatio->GetXaxis()->SetBinLabel(i,pris[i-1]);
  }
  hRatio->GetXaxis()->SetLabelOffset(0.005);
  hRatio->GetXaxis()->SetLabelSize(0.04);
  hRatio->SetMaximum(1);


}
