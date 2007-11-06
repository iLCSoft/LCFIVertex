//
// example ROOT macro to plot LCFIvertex efficiency vs purity curves
//  using PlotProcessor in ROOT mode
//  (i.e. LCFIvertex compiled with USEROOT flag).
//
// you need one or two PlotProcessor output root files. set their names
// (and a brief explanation to be put into the legend)
// in void MakePurityVsEfficiencyRootPlot() close to the end of this macro file, then run
// root -x MakePurityVsEfficiencyRootPlot.C
// and save the plot that it creates.


void drawplot(bool overlay, char* legend_entry) {
  TMultiGraph* FlavourTag
    = dynamic_cast<TMultiGraph*>(gDirectory->Get("FlavourTag"));
  TGraph* bplot=FlavourTag->GetListOfGraphs()->At(0);
  TGraph* cplot=FlavourTag->GetListOfGraphs()->At(1);
  TGraph* bcplot=FlavourTag->GetListOfGraphs()->At(2);

  bplot->SetMarkerStyle(overlay ? 21 : 25);  // 21 for solid square, 25 for open
  bplot->SetMarkerColor(2);
  cplot->SetMarkerStyle(overlay ? 23 : 27);  // 23 for solid triangle, 27 for open
  cplot->SetMarkerColor(3);
  bcplot->SetMarkerStyle(overlay ? 20 : 24); // 20 for solid circle, 24 for open 
  bcplot->SetMarkerColor(4);

  if (!overlay) {
    TCanvas* c1 = new TCanvas();
    TH1F* frame = c1->DrawFrame(0,0,1,1);
  }
  bcplot->Draw("P");
  cplot->Draw("P");
  bplot->Draw("P");
  TText* text = new TText();
  text->SetTextAlign(11);
  text->DrawText(0.8,-0.1,"efficiency");
  text->SetTextAngle(90);
  text->DrawText(-0.08,0.8,"purity");

  TMarker* bmarker = new TMarker(0.1,overlay ? 0.1 : 0.2,overlay ? 21 : 25);
  bmarker->SetMarkerColor(2);
  bmarker->Draw();
  TMarker* cmarker = new TMarker(0.15,overlay ? 0.1 : 0.2,overlay ? 23 : 27);
  cmarker->SetMarkerColor(3);
  cmarker->Draw();
  TMarker* bcmarker = new TMarker(0.2,overlay ? 0.1 : 0.2,overlay ? 20 : 24);
  cmarker->SetMarkerColor(3);
  bcmarker->SetMarkerColor(4);
  bcmarker->Draw();

  text->SetTextAngle(0);
  if (!overlay) {
    text->SetTextAlign(22);
    text->DrawText(0.1,0.28,"b");
    text->DrawText(0.15,0.28,"c");
    text->DrawText(0.2,0.28,"bc");
  }
  text->SetTextAlign(12);
  text->DrawText(0.25,overlay ? 0.1 : 0.2,legend_entry);
  
}


void MakePurityVsEfficiencyRootPlot() {
  TFile* refplotfile = new TFile("PlotProcessorOutput1.root","READ");
  refplotfile->cd();
  drawplot(0,"data contained in file 1");

  // use the following code if you want to overlay another set of curves
  TFile* plotfile = new TFile("PlotProcessorOutput2.root","READ");
  plotfile->cd();
  drawplot(1,"data contained in file 2");
}
