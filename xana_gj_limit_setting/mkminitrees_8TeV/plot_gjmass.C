Double_t plot_gjmass(char EBEE[10]="EB", int rebin=1, int mass=1500, char dirname[30]="nominal", int closefile=0){  
  
  gROOT->Reset();
  gStyle->SetOptFit(0);

  Double_t* yields = new Double_t[4];

  float Lint = 2.0447.;
  float XS = 5.645; //1500GeV

//   if(mass==700) XS=465.7;
//   if(mass==1000) XS=70.05;
//   if(mass==1200) XS=23.94;
//   if(mass==1500) XS=5.645;
//   if(mass==1700) XS=2.317;
//   if(mass==2000) XS=0.6531;
//   if(mass==2500) XS=0.08782;
//   if(mass==3000) XS=0.01234;
//   if(mass==3500) XS=0.001698;

   if(mass==700) XS=24850.;
  if(mass==1000) XS=4180.;
  if(mass==1200) XS=1552.;
  if(mass==1500) XS=412.4;
  if(mass==1700) XS=184.3;
  if(mass==2000) XS=58.581;
  if(mass==2500) XS=9.768;
  if(mass==3000) XS=1.755;
  if(mass==3500) XS=0.3241;

  char fname[100];
  sprintf(fname,"%s_files/ana_gj_out_mc.root", dirname);
  TFile *fmc = new TFile(fname);

  char hname[100];
  sprintf(fname,"%s_files/ana_gj_out_job_qstar_%d.root", dirname,mass);
    
  TFile *fsignal = new TFile(fname);
  sprintf(fname,"%s_files/ana_gj_out_data.root", dirname);
  TFile *fdata = new TFile(fname);

  TH1F* hsig_evt= (TH1F*)fsignal->Get("h_nVtx_w");

  TH2F* h2_EB_sig    = (TH2F*)fsignal->Get("h_EB_mgj_ptg_phojet");
  TH2F* h2_EB_EGdata = (TH2F*)fdata->Get("h_EB_mgj_ptg_phojet");
  TH2F *h2_EB_gj = (TH2F*)fmc->Get("h_EB_mgj_ptg_phojet");
  TH2F *h2_EB_qcd = (TH2F*)fmc->Get("h_EB_mgj_ptg_dijet");
  TH2F* h2_EE_sig    = (TH2F*)fsignal->Get("h_EE_mgj_ptg_phojet");
  TH2F* h2_EE_EGdata = (TH2F*)fdata->Get("h_EE_mgj_ptg_phojet");
  TH2F *h2_EE_gj = (TH2F*)fmc->Get("h_EE_mgj_ptg_phojet");
  TH2F *h2_EE_qcd = (TH2F*)fmc->Get("h_EE_mgj_ptg_dijet");

  TH2F* h2sig = new TH2F();
  TH2F* h2EGdata = new TH2F();
  TH2F *h2gj = new TH2F();
  TH2F *h2qcd = new TH2F();

  if(strcmp(EBEE,"EB") == 0 ) {
    h2sig = (TH2F*)h2_EB_sig->Clone();
    h2EGdata = (TH2F*)h2_EB_EGdata->Clone();
    h2gj = (TH2F*)h2_EB_gj->Clone();
    h2qcd = (TH2F*)h2_EB_qcd->Clone();
  }else if(strcmp(EBEE,"EE") == 0 ) {
    h2sig = (TH2F*)h2_EE_sig->Clone();
    h2EGdata = (TH2F*)h2_EE_EGdata->Clone();
    h2gj = (TH2F*)h2_EE_gj->Clone();
    h2qcd = (TH2F*)h2_EE_qcd->Clone();
  }else if(strcmp(EBEE,"ALL")==0) {
    h2sig = (TH2F*)h2_EB_sig->Clone();
    h2EGdata = (TH2F*)h2_EB_EGdata->Clone();
    h2gj = (TH2F*)h2_EB_gj->Clone();
    h2qcd = (TH2F*)h2_EB_qcd->Clone();

    h2sig->Add(h2_EE_sig);
    h2EGdata->Add(h2_EE_EGdata);
    h2gj->Add(h2_EE_gj);
    h2qcd->Add(h2_EE_qcd);   
  }


  h2gj->Scale(Lint);
  h2qcd->Scale(Lint);

  TH1F *hsig = new TH1F();    
  TH1F *hgj = new TH1F();    
  TH1F *hqcd = new TH1F();   
  TH1F *hEGdata = new TH1F(); 

  hsig = (TH1F*)h2sig->ProjectionX("_px0");
  hgj = (TH1F*)h2gj->ProjectionX("_px1");
  hqcd = (TH1F*)h2qcd->ProjectionX("_px2");
  hEGdata = (TH1F*)h2EGdata->ProjectionX("_px3");

  Int_t *binrange;
  binrange = find_sig_windows(EBEE,mass);
  printf("integrate bin %d - %d \n", binrange[0], binrange[1]);

  int data_obs = hEGdata->Integral(binrange[0], binrange[1]);
  printf(" signal events scaled by %.2f \n", Lint*XS/hsig_evt->Integral());
  hsig->Scale(Lint*XS/hsig_evt->Integral());
  float nsig = hsig->Integral(binrange[0], binrange[1]);

  hEGdata->Rebin(rebin);
  hgj->Rebin(rebin);
  hqcd->Rebin(rebin);

  TCanvas *c1 = new TCanvas("c1","",800,800);
  c1->Draw();
  gPad->SetLogy();

  float min_mass=500;
  if(strcmp(EBEE,"EE")==0) min_mass=500;
  float max_mass=3000;
  if(strcmp(EBEE,"EE")==0) max_mass=2500;
  int bin1 = (int)(min_mass/(10.*rebin))+1;
  int bin2 = (int)(max_mass/(10.*rebin));

  printf("mass range %f - %f (bin %d - %d) \n", min_mass, max_mass, bin1, bin2);

  float data = hEGdata->Integral(bin1,bin2);
  float purity=0.927;
  if(strcmp(EBEE,"EE")==0) purity = 0.854;

  printf("data %f, gj %f, bkg %f\n", 
	 data, hgj->Integral(bin1,bin2), hqcd->Integral(bin1,bin2));


  //normalized by data and purity
  hgj->Scale(data*purity/hgj->Integral(bin1,bin2));
  hqcd->Scale(data*(1-purity)/hqcd->Integral(bin1,bin2));
  //normalized by MC and XS
//   hgj->Scale(1146.*77100.*1.3*0.0064/6.757937e+06);
//   hqcd->Scale(1146.*1.87e+07*0.00216/4.0187e+07);



  TH1F *hsum = (TH1F*)hgj->Clone();
  hsum->Add(hqcd);
  //hsum->Rebin(rebin);

  printf("data %f, gj %f, qcd %f, sum %f \n", 
	 data, hgj->Integral(bin1,bin2), hqcd->Integral(bin1,bin2),
	 hsum->Integral(bin1,bin2));

  hEGdata->SetMarkerStyle(8);
  hEGdata->SetNdivisions(505,"XY");
  hEGdata->SetTitle();
  hEGdata->SetXTitle("M_{#gammajet} (GeV)"); 
  char text[50]; 
  sprintf(text,"Entries / %d GeV", 10*rebin);
  hEGdata->SetYTitle(text);

//   TF1 *fexp = new TF1("fexp","expo",min_mass, max_mass);
  TF1 *fexp = new TF1("fexp",dsigmadm,min_mass, max_mass,4);
  fexp->SetParameters(1.0e-7, 9., 4., -0.2);
  
  hEGdata->Fit(fexp,"0","", min_mass-20., max_mass);
  Double_t par[10];
  par[0] = fexp->GetParameter(0);
  par[1] = fexp->GetParameter(1);
  par[2] = fexp->GetParameter(2);
  par[3] = fexp->GetParameter(3);

  hEGdata->GetXaxis()->SetRangeUser(min_mass, max_mass);

  hEGdata->SetMaximum(hEGdata->GetBinContent(bin1)*100.);
  hEGdata->SetMinimum(0.1);
  hEGdata->Draw("pe");
  
  hsum->SetLineColor(kCyan+2);
  //hsum->SetFillColor(kCyan+1);
  
  hsum->Draw("hist same");
  hqcd->SetLineColor(kRed-3);
  hqcd->SetFillColor(kRed-10);
  hqcd->Draw("hist same");

  fexp->SetLineStyle(2);
  fexp->Draw("same");

  hsig->SetLineStyle(2);
  hsig->SetLineColor(4);
  hsig->Rebin(rebin);
  printf(" nSignal plotted %f (reweight %f) \n", hsig->Integral(), Lint*XS/hsig_evt->Integral());
  hsig->Draw("hist same");

  hEGdata->Draw("pe same");

  gPad->Update();

  TLegend *tleg = new TLegend(0.6, 0.6, 0.9, 0.88);

//   sprintf(text,"Photon Pt %d GeV",ptbin);
  tleg->SetHeader("");
  tleg->SetFillColor(0);
  tleg->SetShadowColor(0);
  tleg->SetBorderSize(0);
  sprintf(text,"Data %.3f fb^{-1}",Lint);
  tleg->AddEntry(hEGdata,text,"lpe");
  //tleg->AddEntry(hEGdata,"(HLT Photon75+90)","");
  tleg->AddEntry(fexp,"Fit (4 par)","l");
  tleg->AddEntry(hsum,"#gamma+jet","fb");
  tleg->AddEntry(hqcd,"dijet","fb");
  tleg->Draw();

  gPad->RedrawAxis();
  
  TCanvas *c2 = new TCanvas("c2","",600,600);
  c2->Draw();
  gPad->SetGrid(1,1);
//   fexp->Draw();

  TH1F *h_data_fit = new TH1F("h_data_fit","data sub fit", 300, 0., 3000.);
  h_data_fit->Rebin(rebin);
  for(int ibin=bin1; ibin<=bin2; ibin++){
    float fit_val = fexp(hEGdata->GetBinCenter(ibin),par[0]);
//     printf("bin %d, data %f fit val %f \n", ibin, hEGdata->GetBinContent(ibin),fit_val);
    h_data_fit->SetBinContent(ibin,hEGdata->GetBinContent(ibin) - fit_val);
    h_data_fit->SetBinError(ibin,hEGdata->GetBinError(ibin));
    
  }
  h_data_fit->SetMarkerStyle(8);
  h_data_fit->SetNdivisions(505,"XY");
  h_data_fit->SetXTitle("M_{#gammajet} (GeV)"); 
  char text[50]; 
  sprintf(text,"(Data - fit) / %d GeV", 10*rebin);
  h_data_fit->SetYTitle(text);
  h_data_fit->SetTitle("");
  h_data_fit->GetXaxis()->SetRangeUser(min_mass-100, 2500.);
  h_data_fit->Draw("p e");
  hsig->Draw("hist same");

  TLegend *tleg2 = new TLegend(0.6, 0.71, 0.9, 0.88);
  tleg2->SetHeader("");
  tleg2->SetFillColor(0);
  tleg2->SetShadowColor(0);
  tleg2->SetBorderSize(0);
  sprintf(text,"Data %.3f fb^{-1}",Lint);
  tleg2->AddEntry(hEGdata,text,"lpe");
  sprintf(text,"q*(%dGeV)",mass);
  tleg2->AddEntry(hsig,text,"fb");
  tleg2->Draw();

  h_data_fit->Draw("p e same");

  float exp_bg = fexp->Integral((binrange[0]-1.)*10.,binrange[1]*10.)/(10.*rebin);
  printf(" %f, %f \n", (binrange[0]-1.)*10.,binrange[1]*10.);
  printf("observed data %.2f \n", data_obs);
  printf("expected signal %.2f \n", nsig);
  printf("expected background %.2f \n", exp_bg);

  gPad->RedrawAxis();

  yields[0] = nsig;
  yields[1] = exp_bg;
  yields[2] = data_obs;

  if(closefile==1){
    fmc->Close();
    fsignal->Close();
    fdata->Close();
  }

  return yields;

}

Double_t dsigmadm (Double_t *v, Double_t *par)
{
  Double_t s = 7000.;
  Double_t mm = v[0]/s;
  Double_t power = par[2]+par[3]*TMath::Log(mm);

  Double_t num = par[0]*TMath::Power(1-mm, par[1]);
  Double_t den = TMath::Power(mm,power);
  // Double_t den = TMath::Power(v[0],par[2]);

  Double_t func = num/den;

  return func;

}

void plot_fun(){

  TF1 *fexp = new TF1("fexp",dsigmadm,300., 2000.,4);
  fexp->SetParameters(0.00011, 8.6, 4.6, -0.1739);

  fexp->Draw();
}


Double_t bifg(Double_t *v, Double_t *par)
{
  Double_t arg = 0;
  if ( v[0] > par[1] ) arg = (v[0] - par[1]) / par[2];
  else arg = (v[0] - par[1]) / par[3];
  
  //Double_t area = par[1]/(par[3]*sqrt(2*3.1415926));
  Double_t norm_area = 1./((par[2]+par[3])*0.5*sqrt(2*3.1415926));
  Double_t gaus = par[0]*norm_area*TMath::Exp(-0.5*arg*arg);

  if (gaus<=0) gaus=1e-10;
  return gaus;
}


Double_t *find_sig_windows(char EBEE[100]="EB", int mass=1000, char dirname[30]="nominal"){
  Int_t* bin = new Int_t[4];

  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->Draw();
  char fname[100];
  sprintf(fname,"%s_files/ana_gj_out_job_qstar_%d.root",dirname,mass);
  TFile *fsig = new TFile(fname);

  char hname[100];
  sprintf(hname, "h_%s_mgj_ptg_phojet", EBEE);  
  TH2F* h2sig = (TH2F*)fsig->Get(hname);

  TH1F *hsig = new TH1F(); 
  hsig = (TH1F*)h2sig->ProjectionX("_px10");
  TH1F *hsig0 = (TH1F*)hsig->Clone();

  hsig->Rebin(5);

  TF1 *func = new TF1("func",bifg, 500., 3000., 4);
  float peak_pos = hsig->GetBinCenter(hsig->GetMaximumBin());
  printf("peak position %f \n", peak_pos);
  func->SetParameters(hsig->Integral(), peak_pos, 80, 100);
//   func->FixParameter(1, peak_pos);
  func->SetParLimits(1, peak_pos-hsig->GetBinWidth(1), peak_pos+hsig->GetBinWidth(1));
  func->SetParLimits(2, mass*0.01, mass*0.5);
  func->SetParLimits(3, mass*0.01, mass*0.5);
  
  if(mass<17000)
    hsig->Fit(func,"","", peak_pos-150, peak_pos+150);
  else
    hsig->Fit(func,"","", peak_pos-250, peak_pos+250);

  float mlow = func->GetParameter(1)-1.5*func->GetParameter(2);
  float mhigh = func->GetParameter(1)+1.5*func->GetParameter(2);
  printf("fit range %f - %f \n", mlow, mhigh);

//   return bin;

  func->SetParLimits(1, func->GetParameter(1)-func->GetParError(1), func->GetParameter(1)-func->GetParError(1));
  func->SetParLimits(2, func->GetParameter(2)-func->GetParError(2), func->GetParameter(2)-func->GetParError(2));
  func->SetParLimits(3, func->GetParameter(3)-func->GetParError(3), func->GetParameter(3)-func->GetParError(3));
  hsig->Fit(func,"","",mlow,mhigh);

  mlow = func->GetParameter(1)-2.*func->GetParameter(2);
  if(mass==2000)  
    mlow = func->GetParameter(1)-5.*func->GetParameter(2);
  mhigh = func->GetParameter(1)+1.5*func->GetParameter(2);
  func->SetParLimits(1, func->GetParameter(1)-func->GetParError(1), func->GetParameter(1)-func->GetParError(1));
  func->SetParLimits(2, func->GetParameter(2)-func->GetParError(2), func->GetParameter(2)-func->GetParError(2));
  func->SetParLimits(3, func->GetParameter(3)-func->GetParError(3), func->GetParameter(3)-func->GetParError(3));
  hsig->Fit(func,"","",mlow,mhigh);


  mlow = func->GetParameter(1)-1.*func->GetParameter(2);
  mhigh = func->GetParameter(1)+1.*func->GetParameter(2);
  printf("mass windows %f - %f \n", mlow, mhigh);
  printf("mass windows in 10GeV bin %f - %f \n", (int)(mlow/10.)*10., (int)(mhigh/10.+0.5)*10.);
  int bin1 =  (int)(mlow/10.)+1;
  int bin2 = (int)(mhigh/10.+0.5);
  printf("Nsig between bin (%d, %d), %.2f \n", bin1,bin2, hsig0->Integral(bin1,bin2));
  

  hsig->Draw();
  //   func->Draw();
  sprintf(hname, "plots/%s_mgjfit_%d.pdf", EBEE,mass);  
  c3->SaveAs(hname);

  bin[0]=bin1;
  bin[1]=bin2;
  return bin;


}

void get_JES_sys(int mass = 1700){
  Double_t *yields_n;
  Double_t *yields_p;
  Double_t *yields_m;
  
  yields_n = plot_gjmass("EB",5,mass,"nominal");
  yields_p = plot_gjmass("EB",5,mass,"JES_plus1");
  yields_m = plot_gjmass("EB",5,mass,"JES_minus1");

  printf("yields nominal, JES+1, JES-1 \n");
  printf("signal %5.2f, %5.2f, %5.2f\n",yields_n[0], yields_p[0], yields_m[0]);
  printf("exp bg %5.2f, %5.2f, %5.2f\n",yields_n[1], yields_p[1], yields_m[1]);
  printf("data   %5.2f, %5.2f, %5.2f\n",yields_n[2], yields_p[2], yields_m[2]);

  printf("JES on signal %.4f, %.4f \n", yields_p[0]/yields_n[0], yields_m[0]/yields_n[0]);
  printf("JES on exp bg %.4f, %.4f \n", yields_p[1]/yields_n[1], yields_m[1]/yields_n[1]);
  printf("JES on data   %.4f, %.4f \n", yields_p[2]/yields_n[2], yields_m[2]/yields_n[2]);


}
