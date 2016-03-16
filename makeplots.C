void makeplots(){

  TFile *f1 = new TFile("SCE.root");
  gStyle->SetOptStat(111111);

  TH1F *hrgen2 = (TH1F*)f1->Get("SCETracks/hrgen2");
  TH1F *hrreco = (TH1F*)f1->Get("SCETracks/hrreco");
  TH1F *hptgen = (TH1F*)f1->Get("SCETracks/hptgen");
  TH1F *hptreco = (TH1F*)f1->Get("SCETracks/hptreco");
  TH1F *hnchgen = (TH1F*)f1->Get("SCETracks/hnchgen");
  TH1F *hnchreco = (TH1F*)f1->Get("SCETracks/hnchreco");

  

  hrgen2->SetLineWidth(2);
  hrgen2->SetLineColor(2);
  hrreco->SetLineWidth(2);
  hrreco->SetLineColor(3);


 hrgen2->Draw();
 hrreco->Draw("same");

 auto effr = new TH1F(*hrreco);
 effr->Divide(hrgen2);
 effr->Draw();


 auto effpt = new TH1F(*hptreco);
 effpt->Divide(hptgen);
 effpt->Draw();


 auto effnch = new TH1F(*hnchreco);
 effnch->Divide(hnchgen);
 effnch->Draw();

}
