void makeplots(){

  TFile *f1 = new TFile("SCE.root");
  gStyle->SetOptStat(111111);

  TH1F *hrgen2 = (TH1F*)f1->Get("SCETracks/hrgen2");
  TH1F *hrreco = (TH1F*)f1->Get("SCETracks/hrreco");

  hrgen2->SetLineWidth(2);
  hrgen2->SetLineColor(2);
  hrreco->SetLineWidth(2);
  hrreco->SetLineColor(3);


 hrgen2->Draw();
 hrreco->Draw("same");
}
