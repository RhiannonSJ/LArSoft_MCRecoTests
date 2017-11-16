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
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/sbn_single_interaction_nt.root" );

  TTree *t_single_ev = (TTree*)f_single->Get("fNt_ev");
  TTree *t_single_vtx = (TTree*)f_single->Get("fNt_vtx");

  std::vector< float > primary;
  std::vector< float > dX;
  std::vector< float > dY;
  std::vector< float > dZ;
  std::vector< float > dR;
  std::vector< float > dr;
  std::vector< float > Xr;
  std::vector< float > Yr;
  std::vector< float > Zr;
  std::vector< float > Xt;
  std::vector< float > Yt;
  std::vector< float > Zt;
  std::vector< float > nVtx;
  std::vector< float > nTrk;
  std::vector< float > event;
  std::vector< float > nVtx_above_thresh;

  int n_entries_vtx = t_single_vtx->GetEntries();
  int n_entries_ev  = t_single_ev->GetEntries();

  for( unsigned int i = 0; i < n_entries_vtx; ++i ){
  
    t_single_vtx->GetEntry(i);

    primary.push_back( t_single_vtx->GetLeaf("primary")->GetValue() );
    dX.push_back( t_single_vtx->GetLeaf("dX")->GetValue() );
    dY.push_back( t_single_vtx->GetLeaf("dY")->GetValue() );
    dZ.push_back( t_single_vtx->GetLeaf("dZ")->GetValue() );
    dR.push_back( t_single_vtx->GetLeaf("dR")->GetValue() );
    Xr.push_back( t_single_vtx->GetLeaf("Xr")->GetValue() );
    Yr.push_back( t_single_vtx->GetLeaf("Yr")->GetValue() );
    Zr.push_back( t_single_vtx->GetLeaf("Zr")->GetValue() );

  }
   
  for( unsigned int i = 0; i < n_entries_ev; ++i ){
  
    t_single_ev->GetEntry(i);

    dr.push_back( t_single_ev->GetLeaf("dr")->GetValue() );
    Xt.push_back( t_single_ev->GetLeaf("Xt")->GetValue() );
    Yt.push_back( t_single_ev->GetLeaf("Yt")->GetValue() );
    Zt.push_back( t_single_ev->GetLeaf("Zt")->GetValue() );
    nVtx.push_back( t_single_ev->GetLeaf("nVtx")->GetValue() );
    nTrk.push_back( t_single_ev->GetLeaf("nTrk")->GetValue() );
    event.push_back( t_single_ev->GetLeaf("evt")->GetValue() );
    nVtx_above_thresh.push_back( t_single_ev->GetLeaf("nVtxThresh")->GetValue() );

  }
  TH1D *h_nVtx         = new TH1D( "h_nVtx",         "Number of reconstructed vertices in an event", 71, 0, 70 );
  TH1D *h_nTrk         = new TH1D( "h_nTrk",         "Number of reconstructed tracks in an event", 71, 0, 70 );
  TH2D *h_trk_vtx      = new TH2D( "h_trk_vtx",      "Number of track vs number of vertices", 71, 0, 70, 71, 0, 70 );

  TH1D *h_zr_p         = new TH1D( "h_zr_p",         "z position of the vertex", 200, 0, 500 );
  TH1D *h_zr_n         = new TH1D( "h_zr_n",         "z position of the vertex", 200, 0, 500 );

  TH1D *h_dR           = new TH1D( "h_dR",           "Distance between true and reconstructed vertices", 200, 0, 20);
  TH1D *h_dX           = new TH1D( "h_dX",           "X distance between true and closest reconstructed vertex", 200, -10, 10);
  TH1D *h_dY           = new TH1D( "h_dY",           "Y distance between true and closest reconstructed vertex", 200, -10, 10);
  TH1D *h_dZ           = new TH1D( "h_dZ",           "Z distance between true and closest reconstructed vertex", 200, -10, 10);

  TH1D *h_dr           = new TH1D( "h_dr",           "Distance between true and closest reconstructed vertex", 200, 0, 20);
  
  TH1D *h_xt           = new TH1D( "h_xt",           "True x positions ", 200, -200, 200);
  TH1D *h_xr           = new TH1D( "h_xr",           "Reco x positions ", 200, -200, 200);

  TH1D *h_nVtx_3cm     = new TH1D( "h_nVtx_3cm",     "Number of reconstructed vertices less than 3cm from true primary", 71, 0, 70 );
  //TH1D *h_nTrk_3cm     = new TH1D( "h_nTrk_3cm",     "Number of reconstructed tracks less than 3cm from true primary", 71, 0, 70 );
  //TH2D *h_trk_vtx_3cm  = new TH2D( "h_trk_vtx_3cm",  "Number of track vs number of vertices less than 3cm from true primary", 71, 0, 70, 71, 0, 70 );

  // Counter for number of events with more than 0 vertices below a certain distance
  int n_ev = 0;

  for( int i = 0; i < n_entries_vtx; ++i ){

    h_xr->Fill( Xr[i] );

    h_dR->Fill( dR[i] );
    h_dX->Fill( dX[i] );
    h_dY->Fill( dY[i] );
    h_dZ->Fill( dZ[i] );

    if( primary[i] == 1 ){
      h_zr_p->Fill(Zr[i]);
    }
    else{ 
      h_zr_n->Fill(Zr[i]);
    }
    
  }

  for( int i = 0; i < n_entries_ev; ++i ){
  
    h_nVtx->Fill( nVtx[i] );
    h_nTrk->Fill( nTrk[i] );

    h_trk_vtx->Fill( nTrk[i], nVtx[i] );

    h_xt->Fill( Xt[i]);

    h_dr->Fill( dr[i] );

    h_nVtx_3cm->Fill(nVtx_above_thresh[i]);
    //h_nTrk_3cm->Fill(nTrk[i]);
    
    //h_trk_vtx_3cm->Fill( nTrk[i], nVtx_above_thresh[i] );
    
  }

  gStyle->SetPalette(57);
  gStyle->SetNumberContours(250);

  h_nVtx->SetStats(kFALSE);
  h_nVtx->SetLineColor(12);
  h_nVtx->GetXaxis()->SetTitle( "Number of reconstructed vertices" );
  h_nVtx->GetXaxis()->SetTitleOffset(1.2);

  h_nTrk->SetStats(kFALSE);
  h_nTrk->SetLineColor(12);
  h_nTrk->GetXaxis()->SetTitle( "Number of reconstructed tracks" );
  h_nTrk->GetXaxis()->SetTitleOffset(1.2);

  h_trk_vtx->SetStats(kFALSE);
  h_trk_vtx->GetXaxis()->SetTitle( "Number of reconstructed tracks" );
  h_trk_vtx->GetYaxis()->SetTitle( "Number of reconstructed vertices" );
  
  h_nVtx_3cm->SetStats(kFALSE);
  h_nVtx_3cm->SetLineColor(12);
  h_nVtx_3cm->GetXaxis()->SetTitle( "Number of reconstructed vertices" );
  h_nVtx_3cm->GetXaxis()->SetTitleOffset(1.2);

  /*
  h_nTrk_3cm->SetStats(kFALSE);
  h_nTrk_3cm->SetLineColor(12);
  h_nTrk_3cm->GetXaxis()->SetTitle( "Number of reconstructed tracks" );
  h_nTrk_3cm->GetXaxis()->SetTitleOffset(1.2);

  h_trk_vtx_3cm->SetStats(kFALSE);
  h_trk_vtx_3cm->GetXaxis()->SetTitle( "Number of reconstructed tracks" );
  h_trk_vtx_3cm->GetYaxis()->SetTitle( "Number of reconstructed vertices" );
  */
  h_dR->SetLineColor(12);
  h_dR->GetXaxis()->SetTitle( "Distance between the true neutrino vertex and all reconstructed vertices [cm]" );
  h_dR->GetXaxis()->SetTitleOffset(1.2);

  h_dX->SetLineColor(12);
  h_dX->GetXaxis()->SetTitle( "X distance between the true neutrino vertex and all reconstructed vertices [cm]" );
  h_dX->GetXaxis()->SetTitleOffset(1.2);

  h_dY->SetLineColor(12);
  h_dY->GetXaxis()->SetTitle( "Y distance between the true neutrino vertex and all reconstructed vertices [cm]" );
  h_dY->GetXaxis()->SetTitleOffset(1.2);

  h_dZ->SetLineColor(12);
  h_dZ->GetXaxis()->SetTitle( "Z distance between the true neutrino vertex and all reconstructed vertices [cm]" );
  h_dZ->GetXaxis()->SetTitleOffset(1.2);

  h_dr->SetLineColor(12);
  h_dr->GetXaxis()->SetTitle( "Distance between the true neutrino vertex and closest reconstructed vertex [cm]" );
  h_dr->GetXaxis()->SetTitleOffset(1.2);

  h_zr_p->SetStats(kFALSE);
  h_zr_p->SetLineColor(2);
  h_zr_p->GetXaxis()->SetTitle( "Distance from Z plane [cm]" );
  h_zr_p->GetXaxis()->SetTitleOffset(1.2);

  h_zr_n->SetStats(kFALSE);
  h_zr_n->SetLineColor(4);
  h_zr_n->GetXaxis()->SetTitle( "Distance from Z plane [cm]" );
  h_zr_n->GetXaxis()->SetTitleOffset(1.2);

  h_xt->SetStats(kFALSE);
  h_xt->SetLineColor(4);
  h_xt->GetXaxis()->SetTitle( "X position [cm]" );
  h_xt->GetXaxis()->SetTitleOffset(1.2);

  h_xr->SetStats(kFALSE);
  h_xr->SetLineColor(2);
  h_xr->GetXaxis()->SetTitle( "X position [cm]" );
  h_xr->GetXaxis()->SetTitleOffset(1.2);

  TLegend *l_z  = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_x  = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_dr = new TLegend(0.7, 0.6, 0.85, 0.85);

  l_z->AddEntry( h_zr_p, "Primary", "l");
  l_z->AddEntry( h_zr_n, "Non-primary", "l");

  l_x->AddEntry( h_xt, "Truth", "l");
  l_x->AddEntry( h_xr, "Reconstructed", "l");
  
  l_dr->AddEntry( h_dR, "All vertices", "l");
  l_dr->AddEntry( h_dr, "Closeset vertex", "l");
  
  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  h_nVtx->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nVtx_single.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nVtx_single.pdf" );

  c->Clear();

  h_nTrk->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_single.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_single.pdf" );

  c->Clear();

  h_trk_vtx->Draw( "colz" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_nVtx_single.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_nVtx_single.pdf" );

  c->Clear();

  h_zr_p->Draw();
  h_zr_n->Draw("same");

  l_z->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/z_distance.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/z_distance.pdf" );
  
  c->Clear();

  h_xt->Draw();
  h_xr->Draw("same");

  l_x->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/x_distance.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/x_distance.pdf" );
  
  c->Clear();

  h_dr->Draw();

  //c1->SetLogy();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_closest_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_closest_single.pdf" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_closest_single.png" );
  c->Clear();


  h_dR->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_all_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_all_single.pdf" );

  c->Clear();

  h_dX->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dX_all_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dX_all_single.pdf" );

  c->Clear();

  h_dY->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dY_all_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dY_all_single.pdf" );

  c->Clear();

  h_dZ->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dZ_all_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dZ_all_single.pdf" );

  c->Clear();

  h_dR->SetLineColor(2);
  h_dr->SetLineColor(4);

  h_dR->Draw();
  h_dr->Draw("same");

  l_dr->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_comp_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/dr_comp_single.pdf" );

  c->Clear();


  c->SetLogy();

  h_nVtx_3cm->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nVtx_single_3cm.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nVtx_single_3cm.pdf" );

  c->Clear();

  
  
  
  
  /*  
  h_nTrk_3cm->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_single_3cm.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_single_3cm.pdf" );

  c->Clear();

  h_trk_vtx_3cm->Draw( "colz" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_nVtx_single_3cm.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/single_plots/nTrk_nVtx_single_3cm.pdf" );

  c->Clear();
  */

}
