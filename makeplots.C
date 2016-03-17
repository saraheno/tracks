void makeplots(){

  TFile *f1 = new TFile("SCE.root");
  gStyle->SetOptStat(111111);


  TCanvas *c1 = new TCanvas("c1","gerrors2",200,10,700,500);
  c1->SetFillColor(42);
  c1->SetGrid();


  TH1F *hrgen2 = (TH1F*)f1->Get("SCETracks/hrgen2");
  TH1F *hrreco = (TH1F*)f1->Get("SCETracks/hrreco");
  TH1F *hptgen = (TH1F*)f1->Get("SCETracks/hptgen");
  TH1F *hptreco = (TH1F*)f1->Get("SCETracks/hptreco");
  TH1F *hnchgen = (TH1F*)f1->Get("SCETracks/hnchgen");
  TH1F *hnchreco = (TH1F*)f1->Get("SCETracks/hnchreco");
  TH1F *hegen = (TH1F*)f1->Get("SCETracks/hegen");
  TH1F *hereco = (TH1F*)f1->Get("SCETracks/hereco");

  

  hrgen2->SetLineWidth(2);
  hrgen2->SetLineColor(2);
  hrreco->SetLineWidth(2);
  hrreco->SetLineColor(3);


 hrgen2->Draw();
 hrreco->Draw("same");
 c1->SaveAs("plot1.jpg");

 auto effr = new TH1F(*hrreco);
 effr->Divide(hrgen2);
 effr->Draw();
 c1->SaveAs("plot2.jpg");


 auto effpt = new TH1F(*hptreco);
 effpt->Divide(hptgen);
 effpt->Draw();
 c1->SaveAs("plot3.jpg");


 auto effnch = new TH1F(*hnchreco);
 effnch->Divide(hnchgen);
 effnch->Draw();
 c1->SaveAs("plot4.jpg");


 auto effe = new TH1F(*hereco);
 effe->Divide(hegen);
 effe->Draw();
 c1->SaveAs("plot5.jpg");

}
