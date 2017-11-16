#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"

void calo_plots() {

  // Read in the sample files
  // Single interaction
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_55_00/LArSoft-v06_55_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/primary_vtx_metrics.root" );

  TTree *t_single_trk = (TTree*)f_single->Get("fNt_tracks");
  TTree *t_single_pri = (TTree*)f_single->Get("fNt_primary");

  std::vector< float > L;
  std::vector< float > KE;
  std::vector< float > nVtx;
  std::vector< float > r_vtx;
  std::vector< float > r_end;
  std::vector< float > t_vtx;
  std::vector< float > t_end;

  // Number of tracks out of the primary or other vertices
  std::vector< float > n_primary;
  std::vector< float > n_other;
  
  int n_entries_trk = t_single_trk->GetEntries();
  int n_entries_pri = t_single_pri->GetEntries();

  for( unsigned int i = 0; i < n_entries_trk; ++i ){
  
    t_single_trk->GetEntry(i);

    L.push_back( t_single_trk->GetLeaf("L")->GetValue() );
    KE.push_back( t_single_trk->GetLeaf("KE")->GetValue() );
    nVtx.push_back( t_single_trk->GetLeaf("nVtx")->GetValue() );
    r_vtx.push_back( t_single_trk->GetLeaf("r_dist_vtx")->GetValue() );
    r_end.push_back( t_single_trk->GetLeaf("r_dist_end")->GetValue() );
    t_vtx.push_back( t_single_trk->GetLeaf("t_dist_vtx")->GetValue() );
    t_end.push_back( t_single_trk->GetLeaf("t_dist_end")->GetValue() );
  
  }
   
  for( unsigned int i = 0; i < n_entries_pri; ++i ){
  
    t_single_pri->GetEntry(i);

    n_primary.push_back( t_single_pri->GetLeaf("nPrimary")->GetValue() );
    n_other.push_back( t_single_pri->GetLeaf("nOther")->GetValue() );
  
  }


  TH1D *h_L     = new TH1D( "h_L", "Longest track's length", 200, 0, 500 );
  
  TH1D *h_KE    = new TH1D( "h_KE", "Longest track kinetic energy", 200, 0, 2000 );
  
  TH1D *h_nVtx  = new TH1D( "h_nVtx", "Number of vertices", 10, 0, 10 );
  
  TH1D *h_r_vtx = new TH1D( "h_r_vtx", "Distance of reconstructed vertex from track vertex", 200, 0, 4 );
  TH1D *h_r_end = new TH1D( "h_r_end", "Distance of reconstructed vertex from track end", 200, 0, 500 );
  
  TH1D *h_t_vtx = new TH1D( "h_t_vtx", "Distance of true vertex from track vertex", 200, 0, 4 );
  TH1D *h_t_end = new TH1D( "h_t_end", "Distance of true vertex from track end", 200, 0, 500 );
  
  TH1D *h_pri = new TH1D( "h_pri", "Number of tracks out of a vertex", 8, 0, 8 );
  TH1D *h_oth = new TH1D( "h_oth", "Number of tracks out of a vertex", 8, 0, 8 );
  
  // Fill histrograms
  for( int i = 0; i < n_entries_trk; ++i ){

    h_L->Fill( L[i] );
    h_KE->Fill( KE[i] );
    h_nVtx->Fill( nVtx[i] );
  
    h_r_vtx->Fill( r_vtx[i] );
    h_r_end->Fill( r_end[i] );

    h_t_vtx->Fill( t_vtx[i] );
    h_t_end->Fill( t_end[i] );

  }

  for( int i = 0; i < n_entries_pri; ++i ){
  
    h_pri->Fill( n_primary[i] );
    h_oth->Fill( n_other[i] );

  }


  gStyle->SetPalette(57);
  gStyle->SetNumberContours(250);

  h_L->SetStats(kFALSE);
  h_L->GetXaxis()->SetTitle( "Longest track length [cm]" );
  h_L->GetXaxis()->SetTitleOffset(1.2);

  h_KE->SetStats(kFALSE);
  h_KE->GetXaxis()->SetTitle( "Longest track kinetic energy [MeV]" );
  h_KE->GetXaxis()->SetTitleOffset(1.2);

  h_r_vtx->SetStats(kFALSE);
  h_r_vtx->SetLineColor(2);
  h_r_vtx->GetXaxis()->SetTitle( "Distance of track vertex from vertex [cm]" );
  h_r_vtx->GetXaxis()->SetTitleOffset(1.2);

  h_t_vtx->SetStats(kFALSE);
  h_t_vtx->SetLineColor(4);
  h_t_vtx->GetXaxis()->SetTitle( "Distance of track vertex from vertex [cm]" );
  h_t_vtx->GetXaxis()->SetTitleOffset(1.2);

  h_r_end->SetStats(kFALSE);
  h_r_end->SetLineColor(2);
  h_r_end->GetXaxis()->SetTitle( "Distance of track end from vertex [cm]" );
  h_r_end->GetXaxis()->SetTitleOffset(1.2);

  h_t_end->SetStats(kFALSE);
  h_t_end->SetLineColor(4);
  h_t_end->GetXaxis()->SetTitle( "Distance of track end from vertex [cm]" );
  h_t_end->GetXaxis()->SetTitleOffset(1.2);

  h_pri->SetStats(kFALSE);
  h_pri->SetLineColor(2);
  h_pri->GetXaxis()->SetTitle( "Number of tracks from a vertex" );
  h_pri->GetXaxis()->SetTitleOffset(1.2);

  h_oth->SetStats(kFALSE);
  h_oth->SetLineColor(4);
  h_oth->GetXaxis()->SetTitle( "Number of tracks from a vertex" );
  h_oth->GetXaxis()->SetTitleOffset(1.2);

  TLegend *l_vtx = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_end = new TLegend(0.7, 0.6, 0.85, 0.85);
  TLegend *l_trk = new TLegend(0.7, 0.6, 0.85, 0.85);
  
  l_vtx->AddEntry( h_r_vtx, "Reconstructed primary vertex", "l");
  l_vtx->AddEntry( h_t_vtx, "True primary vertex", "l");
  
  l_end->AddEntry( h_r_end, "Reconstructed primary vertex", "l");
  l_end->AddEntry( h_t_end, "True primary vertex", "l");

  l_trk->AddEntry( h_pri, "Primary vertex", "l");
  l_trk->AddEntry( h_oth, "Non-primary vertex", "l");
  
  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  h_L->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_55_00/LArSoft-v06_55_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/longest_track_length.root" );

  c->Clear();

  h_KE->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_55_00/LArSoft-v06_55_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/longest_track_KE.root" );

  c->Clear();

  h_r_vtx->Draw();
  h_t_vtx->Draw("same");

  l_vtx->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_55_00/LArSoft-v06_55_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/track_vertex_dist.root" );

  c->Clear();

  h_t_end->Draw();
  h_r_end->Draw("same");

  l_end->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_55_00/LArSoft-v06_55_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/track_end_dist.root" );

  c->Clear();

  h_pri->Draw();
  h_oth->Draw("same");

  l_trk->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_55_00/LArSoft-v06_55_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/n_tracks_from_vertex.root" );

  c->Clear();

}
