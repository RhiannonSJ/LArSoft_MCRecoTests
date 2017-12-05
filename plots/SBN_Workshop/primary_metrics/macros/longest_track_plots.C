#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"

void longest_track_plots() {

  // Read in the sample files
  // Single interaction
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/longest_track_info.root" );

  TTree *t_single = (TTree*)f_single->Get("fNt_length");

  // From track nTuple
  std::vector< float > primary;
  std::vector< float > dR;
  std::vector< float > L;
  std::vector< float > KE;

  int n_entries = t_single->GetEntries();

  for( unsigned int i = 0; i < n_entries; ++i ){
  
    t_single->GetEntry(i);

    primary.push_back( t_single->GetLeaf("primary")->GetValue() );
    L.push_back(       t_single->GetLeaf(      "L")->GetValue() );
    KE.push_back(      t_single->GetLeaf(     "KE")->GetValue() );
    dR.push_back(      t_single->GetLeaf(     "dR")->GetValue() );
  
  }

  TH1D *h_L_p        = new TH1D( "h_L_p",  "Track length", 100, 0, 500 );
  TH1D *h_L_n        = new TH1D( "h_L_n",  "Track length", 100, 0, 500 );
  
  TH1D *h_KE_p       = new TH1D( "h_KE_p", "Kinetic energy", 100, 0, 2000 );
  TH1D *h_KE_n       = new TH1D( "h_KE_n", "Kinetic energy", 100, 0, 2000 );
  
  TH1D *h_dR_p       = new TH1D( "h_dR_p", "Distance between reconstructed and MC vertex", 150, 0, 150 );
  TH1D *h_dR_n       = new TH1D( "h_dR_n", "Distance between reconstructed and MC vertex", 150, 0, 150 );
  
  // Fill histrograms
  for( int i = 0; i < n_entries; ++i ){
    if( primary[i] == 1 ){
    
      h_L_p->Fill( L[i] );
      h_KE_p->Fill( KE[i] );
      h_dR_p->Fill( dR[i] );

    }
    else{
    
      h_L_n->Fill( L[i] );
      h_KE_n->Fill( KE[i] );
      h_dR_n->Fill( dR[i] );
    }
  }

  gStyle->SetPalette(57);
  gStyle->SetNumberContours(250);

  // Set the style for histograms, normalise for shape dependencies 
  h_L_p->SetStats(kFALSE);
  h_L_p->SetLineColor(2);
  h_L_p->GetXaxis()->SetTitle( "Track length [cm]" );
  h_L_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_L_p = h_L_p->Integral();
  h_L_p->Scale(1/scale_L_p);

  h_L_n->SetStats(kFALSE);
  h_L_n->SetLineColor(4);
  h_L_n->GetXaxis()->SetTitle( "Track length [cm]" );
  h_L_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_L_n = h_L_n->Integral();
  h_L_n->Scale(1/scale_L_n);

  h_KE_p->SetStats(kFALSE);
  h_KE_p->SetLineColor(2);
  h_KE_p->GetXaxis()->SetTitle( "Track kinetic energy [MeV]" );
  h_KE_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_KE_p = h_KE_p->Integral();
  h_KE_p->Scale(1/scale_KE_p);

  h_KE_n->SetStats(kFALSE);
  h_KE_n->SetLineColor(4);
  h_KE_n->GetXaxis()->SetTitle( "Track kinetic energy [MeV]" );
  h_KE_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_KE_n = h_KE_n->Integral();
  h_KE_n->Scale(1/scale_KE_n);

  h_dR_p->SetStats(kFALSE);
  h_dR_p->SetLineColor(2);
  h_dR_p->GetXaxis()->SetTitle( "True - reco vertex distance [cm]" );
  h_dR_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_dR_p = h_dR_p->Integral();
  h_dR_p->Scale(1/scale_dR_p);

  h_dR_n->SetStats(kFALSE);
  h_dR_n->SetLineColor(4);
  h_dR_n->GetXaxis()->SetTitle( "True - reco vertex distance [cm]" );
  h_dR_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_dR_n = h_dR_n->Integral();
  h_dR_n->Scale(1/scale_dR_n);

  TLegend *l_L  = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_KE = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_dR = new TLegend(0.5, 0.65, 0.85, 0.85);
 
  l_L->AddEntry(  h_L_p,  "Primary",     "l" );
  l_L->AddEntry(  h_L_n,  "Non-primary", "l" );
  
  l_KE->AddEntry( h_KE_p, "Primary",     "l" );
  l_KE->AddEntry( h_KE_n, "Non-primary", "l" );
  
  l_dR->AddEntry( h_dR_p, "Primary",     "l" );
  l_dR->AddEntry( h_dR_n, "Non-primary", "l" );
  
  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  h_L_n->Draw("hist");
  h_L_p->Draw("hist same");
  l_L->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/longest_track_plots/track_length.root" );

  c->Clear();

  h_KE_n->Draw("hist");
  h_KE_p->Draw("hist same");
  l_KE->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/longest_track_plots/track_KE.root" );

  c->Clear();

  h_dR_p->Draw("hist");
  h_dR_n->Draw("hist same");
  l_dR->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/longest_track_plots/true_reco_vertex_distance.root" );

  c->Clear();
}
