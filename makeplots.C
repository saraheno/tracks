void makeplots(){

  TFile *f1 = new TFile("SCE.root");
  gStyle->SetOptStat(111111);

  TH1F *hrgen2 = (TH1F*)f1->Get("SCETracks/hrgen2");
  //TH1F *hratvr = (TH1F*)f1->Get("hratvr");



 hrgen2->Draw();
 // hratvr->Draw("same");
}
