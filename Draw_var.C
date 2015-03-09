void Draw_var(char var[100]="recoPt", char cut[300]="recoPt>15.", int drawlog=0){
  TFile *f1 = new TFile("skim.root");
  TTree *t1 = (TTree*)f1->Get("t");

  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->Draw();
  if(drawlog==1) gPad->SetLogy(1);

  char cutn[100];
  char txt[100];
  //TH1F *h1;
  sprintf(cutn,"isMatched==1&&%s",cut);
  sprintf(txt,"%s>>h1",var);
  t1->Draw(txt,cutn);
  TH1F* hs = (TH1F*)h1->Clone();
  hs->Draw();

  sprintf(cutn,"isMatched<1&&%s",cut);
  t1->Draw(txt,cutn);
  TH1F* hf = (TH1F*)h1->Clone();  

  hs->Scale(1./hs->Integral());
  hf->Scale(1./hf->Integral());
   
  float ymax = hs->GetMaximum();
  if(hf->GetMaximum()>hs->GetMaximum()) ymax = hf->GetMaximum();
  hs->SetMaximum(ymax*1.2);
  if(drawlog==1)    hs->SetMaximum(ymax*10.);
  hs->SetTitle("");
  hs->SetNdivisions(505,"XY");
  hs->SetXTitle(var);

  hs->Draw();
  hf->SetLineColor(2);
  hf->Draw("same");

  sprintf(txt,"plots/%s.pdf",var);
  c1->SaveAs(txt);
   
  TLegend *tleg = new TLegend(0.15,0.78,0.4,0.88);
  //tleg->SetHeader(txt);
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  tleg->AddEntry(hs,"true photon","l");
  tleg->AddEntry(hf,"fake photon","l");
  tleg->Draw();

  hs->Draw("same");
  hf->Draw("same");
  gPad->RedrawAxis();



}
  

void Draw_varES(char var[100]="recoPt", char cut[100]="recoPt>15."){
  TFile *f1 = new TFile("ES_pt20.root");
  TFile *f2 = new TFile("ES_pt40.root");
  
  TTree *t1 = (TTree*)f1->Get("t");
  TTree *t2 = (TTree*)f2->Get("t");


  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(2,1);
  c1->cd(1);

  char cutn[100];
  sprintf(cutn,"isMatched==1&&%s",cut);

  t1->Draw(var,cutn);
  c1->cd(2);
  t2->Draw(var,cutn);

}
  
