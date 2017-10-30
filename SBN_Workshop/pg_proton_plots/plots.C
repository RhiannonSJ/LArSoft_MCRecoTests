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
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/sbn_pg_proton_nt.root" );

  TTree *t_single = (TTree*)f_single->Get("fNt");

  std::vector< float > primary;
  std::vector< float > dX;
  std::vector< float > dY;
  std::vector< float > dZ;
  std::vector< float > dR;
  std::vector< float > Xr;
  std::vector< float > Yr;
  std::vector< float > Zr;
  std::vector< float > Xt;
  std::vector< float > Yt;
  std::vector< float > Zt;
  std::vector< float > nVtx;
  std::vector< float > nTrk;
  std::vector< float > pdg;

  int n_entries = t_single->GetEntries();

  for( unsigned int i = 0; i < n_entries; ++i ){
  
    t_single->GetEntry(i);

    primary.push_back( t_single->GetLeaf("primary")->GetValue() );
    dX.push_back( t_single->GetLeaf("dX")->GetValue() );
    dY.push_back( t_single->GetLeaf("dY")->GetValue() );
    dZ.push_back( t_single->GetLeaf("dZ")->GetValue() );
    dR.push_back( t_single->GetLeaf("dR")->GetValue() );
    Xr.push_back( t_single->GetLeaf("Xr")->GetValue() );
    Yr.push_back( t_single->GetLeaf("Yr")->GetValue() );
    Zr.push_back( t_single->GetLeaf("Zr")->GetValue() );
    Xt.push_back( t_single->GetLeaf("Xt")->GetValue() );
    Yt.push_back( t_single->GetLeaf("Yt")->GetValue() );
    Zt.push_back( t_single->GetLeaf("Zt")->GetValue() );
    nVtx.push_back( t_single->GetLeaf("nVtx")->GetValue() );
    nTrk.push_back( t_single->GetLeaf("nTrk")->GetValue() );
    pdg.push_back( t_single->GetLeaf("pdg")->GetValue() );

  }
   

  TH1D *h_nVtx         = new TH1D( "h_nVtx",         "Number of reconstructed vertices in an event", 71, 0, 70 );
  TH1D *h_nTrk         = new TH1D( "h_nTrk",         "Number of reconstructed tracks in an event", 71, 0, 70 );
  TH2D *h_trk_vtx      = new TH2D( "h_trk_vtx",      "Number of track vs number of vertices", 71, 0, 70, 71, 0, 70 );

  TH1D *h_zr_p         = new TH1D( "h_zr_p",         "z position of the vertex", 200, 0, 500 );
  TH1D *h_zr_n         = new TH1D( "h_zr_n",         "z position of the vertex", 200, 0, 500 );

  TH1D *h_dr           = new TH1D( "h_dr",           "Distance between true and reconstructed vertices", 200, 0, 20);

  TH1D *h_xt           = new TH1D( "h_xt",           "True x positions ", 200, -200, 200);
  TH1D *h_xr           = new TH1D( "h_xr",           "Reco x positions ", 200, -200, 200);

  TH1D *h_nVtx_3mm     = new TH1D( "h_nVtx_3mm",     "Number of reconstructed vertices less than 3mm from true primary", 71, 0, 70 );
  TH1D *h_nTrk_3mm     = new TH1D( "h_nTrk_3mm",     "Number of reconstructed tracks less than 3mm from true primary", 71, 0, 70 );
  TH2D *h_trk_vtx_3mm  = new TH2D( "h_trk_vtx_3mm",  "Number of track vs number of vertices less than 3mm from true primary", 71, 0, 70, 71, 0, 70 );

  // Counter for number of events with more than 0 vertices below a certain distance
  int n_ev = 0;

  for( int i = 0; i < n_entries; ++i ){
  
    if( pdg[i] == 2211 ){
  
      h_nVtx->Fill( nVtx[i] );
      h_nTrk->Fill( nTrk[i] );

      h_trk_vtx->Fill( nTrk[i], nVtx[i] );

      h_xt->Fill( Xt[i]);

      h_xr->Fill( Xr[i] );

      h_dr->Fill( dR[i] );

      if( primary[i] == 1 ){
        h_zr_p->Fill(Zr[i]);
      }
      else{ 
      h_zr_n->Fill(Zr[i]);
      }

      if( dR[i] < 0.3 ){

        n_ev++;

        h_nVtx_3mm->Fill(nVtx[i]);
        h_nTrk_3mm->Fill(nTrk[i]);
    
        h_trk_vtx_3mm->Fill( nTrk[i], nVtx[i] );
    
      }
    }
  }
  std::cout << " N vertices : " << n_ev << std::endl;

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
  
  h_nVtx_3mm->SetStats(kFALSE);
  h_nVtx_3mm->SetLineColor(12);
  h_nVtx_3mm->GetXaxis()->SetTitle( "Number of reconstructed vertices" );
  h_nVtx_3mm->GetXaxis()->SetTitleOffset(1.2);

  h_nTrk_3mm->SetStats(kFALSE);
  h_nTrk_3mm->SetLineColor(12);
  h_nTrk_3mm->GetXaxis()->SetTitle( "Number of reconstructed tracks" );
  h_nTrk_3mm->GetXaxis()->SetTitleOffset(1.2);

  h_trk_vtx_3mm->SetStats(kFALSE);
  h_trk_vtx_3mm->GetXaxis()->SetTitle( "Number of reconstructed tracks" );
  h_trk_vtx_3mm->GetYaxis()->SetTitle( "Number of reconstructed vertices" );
  
  h_dr->SetLineColor(12);
  h_dr->GetXaxis()->SetTitle( "Distance between the true neutrino vertex and all reconstructed vertices [cm]" );
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

  TLegend *l_z = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_x = new TLegend(0.7, 0.6, 0.85, 0.85);

  l_z->AddEntry( h_zr_p, "Primary", "l");
  l_z->AddEntry( h_zr_n, "Non-primary", "l");

  l_x->AddEntry( h_xt, "Truth", "l");
  l_x->AddEntry( h_xr, "Reconstructed", "l");
  
  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  h_nVtx->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nVtx_single.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nVtx_single.pdf" );

  c->Clear();

  h_nTrk->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_single.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_single.pdf" );

  c->Clear();

  h_trk_vtx->Draw( "colz" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_nVtx_single.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_nVtx_single.pdf" );

  c->Clear();

  h_zr_p->Draw();
  h_zr_n->Draw("same");

  l_z->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/z_distance.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/z_distance.pdf" );
  
  c->Clear();

  h_xt->Draw();
  h_xr->Draw("same");

  l_x->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/x_distance.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/x_distance.pdf" );
  
  c->Clear();

  h_dr->Draw();

  //c1->SetLogy();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/dr_all_single.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/dr_all_single.pdf" );

  c->Clear();

  c->SetLogy();

  h_nVtx_3mm->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nVtx_single_3mm.root" );
  
  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nVtx_single_3mm.pdf" );

  c->Clear();
  
  h_nTrk_3mm->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_single_3mm.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_single_3mm.pdf" );

  c->Clear();

  h_trk_vtx_3mm->Draw( "colz" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_nVtx_single_3mm.root" );

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/nTrk_nVtx_single_3mm.pdf" );

  c->Clear();
}
