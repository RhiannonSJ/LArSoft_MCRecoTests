#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"

void plots() {

  // Read in the sample files
  // Single interaction
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/primary_vtx_metrics.root" );

  //TTree *t_single_trk = (TTree*)f_single->Get("fNt_tracks");
  TTree *t_single_pri = (TTree*)f_single->Get("fNt_primary");
  TTree *t_single_oth = (TTree*)f_single->Get("fNt_other");

  // Xr:Yr:Zr:L:thetaxz:thetayz:thetatoz:P

  std::vector< float > Xr_p;
  std::vector< float > Yr_p;
  std::vector< float > Zr_p;
  std::vector< float > L_p;
  std::vector< float > thetaXZ_p;
  std::vector< float > thetaYZ_p;
  std::vector< float > thetatoZ_p;

  std::vector< float > Xr_o;
  std::vector< float > Yr_o;
  std::vector< float > Zr_o;
  std::vector< float > L_o;
  std::vector< float > thetaXZ_o;
  std::vector< float > thetaYZ_o;
  std::vector< float > thetatoZ_o;

  int n_entries_pri = t_single_pri->GetEntries();
  int n_entries_oth = t_single_oth->GetEntries();

  for( unsigned int i = 0; i < n_entries_pri; ++i ){
  
    t_single_pri->GetEntry(i);

    Xr_p.push_back( t_single_pri->GetLeaf("Xr")->GetValue() );
    Yr_p.push_back( t_single_pri->GetLeaf("Yr")->GetValue() );
    Zr_p.push_back( t_single_pri->GetLeaf("Zr")->GetValue() );
    L_p.push_back( t_single_pri->GetLeaf("P")->GetValue() );
    thetaXZ_p.push_back( t_single_pri->GetLeaf("thetaxz")->GetValue() );
    thetaYZ_p.push_back( t_single_pri->GetLeaf("thetayz")->GetValue() );
    thetatoZ_p.push_back( t_single_pri->GetLeaf("thetatoz")->GetValue() );
  }
   
  for( unsigned int i = 0; i < n_entries_oth; ++i ){
  
    t_single_oth->GetEntry(i);

    Xr_o.push_back( t_single_oth->GetLeaf("Xr")->GetValue() );
    Yr_o.push_back( t_single_oth->GetLeaf("Yr")->GetValue() );
    Zr_o.push_back( t_single_oth->GetLeaf("Zr")->GetValue() );
    L_o.push_back( t_single_oth->GetLeaf("P")->GetValue() );
    thetaXZ_o.push_back( t_single_oth->GetLeaf("thetaxz")->GetValue() );
    thetaYZ_o.push_back( t_single_oth->GetLeaf("thetayz")->GetValue() );
    thetatoZ_o.push_back( t_single_oth->GetLeaf("thetatoz")->GetValue() );
  }
  
  TH1D *h_xr_p = new TH1D( "h_xr_p", "Reconstructed x position", 200, -200, 200 );
  TH1D *h_xr_o = new TH1D( "h_xr_o", "Reconstructed x position", 200, -200, 200 );
  
  TH1D *h_yr_p = new TH1D( "h_yr_p", "Reconstructed y position", 200, -200, 200 );
  TH1D *h_yr_o = new TH1D( "h_yr_o", "Reconstructed y position", 200, -200, 200 );
  
  TH1D *h_zr_p = new TH1D( "h_zr_p", "Reconstructed z position", 200, 0, 500 );
  TH1D *h_zr_o = new TH1D( "h_zr_o", "Reconstructed z position", 200, 0, 500 );
 
  TH1D *h_z_dif = new TH1D( "h_z_dif", "Non-primary - primary z position", 200, -250, 250 );

  TH1D *h_L_p  = new TH1D( "h_L_p",  "Reconstructed track length", 200, 0, 500 );
  TH1D *h_L_o  = new TH1D( "h_L_o",  "Reconstructed track length", 200, 0, 500 );
  
  TH1D *h_xz_p = new TH1D( "h_xz_p", "Reconstructed #theta_{XZ}", 200, -50, 50 );
  TH1D *h_xz_o = new TH1D( "h_xz_o", "Reconstructed #theta_{XZ} ", 200, -50, 50 );
  
  TH1D *h_yz_p = new TH1D( "h_yz_p", "Reconstructed #theta_{YZ}", 200, -50, 50 );
  TH1D *h_yz_o = new TH1D( "h_yz_o", "Reconstructed #theta_{YZ} ", 200, -50, 50 );
  
  TH1D *h_z_p  = new TH1D( "h_z_p",  "Reconstructed #theta_{Z}", 200, -50, 50 );
  TH1D *h_z_o  = new TH1D( "h_z_o",  "Reconstructed #theta_{Z} ", 200, -50, 50 );
  
  // Fill histrograms
  for( int i = 0; i < n_entries_pri; ++i ){

    h_xr_p->Fill(Xr_p[i]);
    h_yr_p->Fill(Yr_p[i]);
    h_zr_p->Fill(Zr_p[i]);
    
    h_L_p->Fill(L_p[i]);
  
    h_xz_p->Fill(thetaXZ_p[i]);
    h_yz_p->Fill(thetaYZ_p[i]);
    
    h_z_p->Fill(thetatoZ_p[i]);

  }


  for( int i = 0; i < n_entries_oth; ++i ){

    h_xr_o->Fill(Xr_o[i]);
    h_yr_o->Fill(Yr_o[i]);
    h_zr_o->Fill(Zr_o[i]);
    
    h_L_o->Fill(L_o[i]);
  
    h_xz_o->Fill(thetaXZ_o[i]);
    h_yz_o->Fill(thetaYZ_o[i]);
    
    h_z_o->Fill(thetatoZ_o[i]);
  }

  gStyle->SetPalette(57);
  gStyle->SetNumberContours(250);

  h_xr_p->SetStats(kFALSE);
  h_xr_p->SetLineColor(2);
  h_xr_p->GetXaxis()->SetTitle( "Distance from X plane [cm]" );
  h_xr_p->GetXaxis()->SetTitleOffset(1.2);

  h_xr_o->SetStats(kFALSE);
  h_xr_o->SetLineColor(4);
  h_xr_o->GetXaxis()->SetTitle( "Distance from X plane [cm]" );
  h_xr_o->GetXaxis()->SetTitleOffset(1.2);

  h_yr_p->SetStats(kFALSE);
  h_yr_p->SetLineColor(2);
  h_yr_p->GetXaxis()->SetTitle( "Distance from Y plane [cm]" );
  h_yr_p->GetXaxis()->SetTitleOffset(1.2);

  h_yr_o->SetStats(kFALSE);
  h_yr_o->SetLineColor(4);
  h_yr_o->GetXaxis()->SetTitle( "Distance from Y plane [cm]" );
  h_yr_o->GetXaxis()->SetTitleOffset(1.2);

  h_zr_p->SetStats(kFALSE);
  h_zr_p->SetLineColor(2);
  h_zr_p->GetXaxis()->SetTitle( "Distance from Z plane [cm]" );
  h_zr_p->GetXaxis()->SetTitleOffset(1.2);

  h_zr_o->SetStats(kFALSE);
  h_zr_o->SetLineColor(4);
  h_zr_o->GetXaxis()->SetTitle( "Distance from Z plane [cm]" );
  h_zr_o->GetXaxis()->SetTitleOffset(1.2);

  h_L_p->SetStats(kFALSE);
  h_L_p->SetLineColor(2);
  h_L_p->GetXaxis()->SetTitle( "Track length [cm]" );
  h_L_p->GetXaxis()->SetTitleOffset(1.2);

  h_L_o->SetStats(kFALSE);
  h_L_o->SetLineColor(4);
  h_L_o->GetXaxis()->SetTitle( "Track length [cm]" );
  h_L_o->GetXaxis()->SetTitleOffset(1.2);

  h_xz_p->SetStats(kFALSE);
  h_xz_p->SetLineColor(2);
  h_xz_p->GetXaxis()->SetTitle( "#theta_{XZ}" );
  h_xz_p->GetXaxis()->SetTitleOffset(1.2);

  h_xz_o->SetStats(kFALSE);
  h_xz_o->SetLineColor(4);
  h_xz_o->GetXaxis()->SetTitle( "#theta_{XZ}" );
  h_xz_o->GetXaxis()->SetTitleOffset(1.2);

  h_yz_p->SetStats(kFALSE);
  h_yz_p->SetLineColor(2);
  h_yz_p->GetXaxis()->SetTitle( "#theta_{YZ}" );
  h_yz_p->GetXaxis()->SetTitleOffset(1.2);

  h_yz_o->SetStats(kFALSE);
  h_yz_o->SetLineColor(4);
  h_yz_o->GetXaxis()->SetTitle( "#theta_{YZ}" );
  h_yz_o->GetXaxis()->SetTitleOffset(1.2);

  h_z_p->SetStats(kFALSE);
  h_z_p->SetLineColor(2);
  h_z_p->GetXaxis()->SetTitle( "#theta_{Z}" );
  h_z_p->GetXaxis()->SetTitleOffset(1.2);

  h_z_o->SetStats(kFALSE);
  h_z_o->SetLineColor(4);
  h_z_o->GetXaxis()->SetTitle( "#theta_{Z}" );
  h_z_o->GetXaxis()->SetTitleOffset(1.2);


  TLegend *l_xr = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_yr = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_zr = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_L  = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_xz = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_yz = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_z  = new TLegend(0.7, 0.6, 0.85, 0.85);

  l_xr->AddEntry( h_xr_p, "Primary", "l");
  l_xr->AddEntry( h_xr_o, "Non-primary", "l");

  l_yr->AddEntry( h_yr_p, "Primary", "l");
  l_yr->AddEntry( h_yr_o, "Non-primary", "l");

  l_zr->AddEntry( h_zr_p, "Primary", "l");
  l_zr->AddEntry( h_zr_o, "Non-primary", "l");

  l_L->AddEntry( h_L_p, "Primary", "l");
  l_L->AddEntry( h_L_o, "Non-primary", "l");

  l_xz->AddEntry( h_xz_p, "Primary", "l");
  l_xz->AddEntry( h_xz_o, "Non-primary", "l");

  l_yz->AddEntry( h_yz_p, "Primary", "l");
  l_yz->AddEntry( h_yz_o, "Non-primary", "l");

  l_z->AddEntry( h_z_p, "Primary", "l");
  l_z->AddEntry( h_z_o, "Non-primary", "l");

  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  h_xr_o->Draw();
  h_xr_p->Draw("same");

  l_xr->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/x_distance.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/x_distance.pdf" );
  
  c->Clear();

  h_yr_o->Draw();
  h_yr_p->Draw("same");

  l_yr->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/y_distance.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/y_distance.pdf" );
  
  c->Clear();
  
  h_zr_o->Draw();
  h_zr_p->Draw("same");

  l_zr->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/z_distance.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/z_distance.pdf" );
  
  c->Clear();

  h_L_o->Draw();
  h_L_p->Draw("same");

  l_L->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/track_length.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/track_length.pdf" );
  
  c->Clear();

  h_xz_o->Draw();
  h_xz_p->Draw("same");

  l_xz->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/thetaXZ.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/thetaXZ.pdf" );
  
  c->Clear();

  h_yz_o->Draw();
  h_yz_p->Draw("same");

  l_yz->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/thetaYZ.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/thetaYZ.pdf" );
  
  c->Clear();

  h_z_o->Draw();
  h_z_p->Draw("same");

  l_z->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/thetaZ.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/thetaZ.pdf" );
  
  c->Clear();


}
