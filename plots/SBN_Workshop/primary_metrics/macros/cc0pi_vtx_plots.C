#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"

void cc0pi_vtx_plots() {

  // Read in the sample files
  // Single interaction
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/cc0pi_primary_vtx_metrics.root" );

  TTree *t_single_cal = (TTree*)f_single->Get("fNt_calo");
  TTree *t_single_pri = (TTree*)f_single->Get("fNt_primary");

  // From track nTuple
  std::vector< float > t_primary;
  std::vector< float > L;
  std::vector< float > range;
  std::vector< float > hits;
  std::vector< float > ke;

  // From vertex nTuple
  std::vector< float > v_primary;
  std::vector< float > nTrk;
  std::vector< float > keSum;
  std::vector< float > maxL;
  std::vector< float > maxKE;
  std::vector< float > maxHits;
  std::vector< float > maxRange;

  int n_entries_cal = t_single_cal->GetEntries();
  int n_entries_pri = t_single_pri->GetEntries();

  for( unsigned int i = 0; i < n_entries_cal; ++i ){
  
    t_single_cal->GetEntry(i);

    t_primary.push_back( t_single_cal->GetLeaf("primary")->GetValue() );
    L.push_back( t_single_cal->GetLeaf("L")->GetValue() );
    ke.push_back( t_single_cal->GetLeaf("ke")->GetValue() );
    hits.push_back( t_single_cal->GetLeaf("hits")->GetValue() );
    range.push_back( t_single_cal->GetLeaf("range")->GetValue() );
  
  }
   
  for( unsigned int i = 0; i < n_entries_pri; ++i ){
  
    t_single_pri->GetEntry(i);

    v_primary.push_back( t_single_pri->GetLeaf("primary")->GetValue() );
    nTrk.push_back( t_single_pri->GetLeaf("nTrk")->GetValue() );
    keSum.push_back( t_single_pri->GetLeaf("keSum")->GetValue() );
    maxL.push_back( t_single_pri->GetLeaf("maxL")->GetValue() );
    maxKE.push_back( t_single_pri->GetLeaf("maxKE")->GetValue() );
    maxHits.push_back( t_single_pri->GetLeaf("maxHits")->GetValue() );
    maxRange.push_back( t_single_pri->GetLeaf("maxRange")->GetValue() );
  
  }


  TH1D *h_L_p     = new TH1D( "h_L_p",     "Track length", 100, 0, 500 );
  TH1D *h_L_n     = new TH1D( "h_L_n",     "Track length", 100, 0, 500 );
  
  TH1D *h_ke_p    = new TH1D( "h_ke_p",    "Kinetic energy", 100, 0, 2000 );
  TH1D *h_ke_n    = new TH1D( "h_ke_n",    "Kinetic energy", 100, 0, 2000 );
  
  TH1D *h_hits_p  = new TH1D( "h_hits_p",  "Number of hits per track", 2000, 0, 2000 );
  TH1D *h_hits_n  = new TH1D( "h_hits_n",  "Number of hits per track", 2000, 0, 2000 );
  
  TH1D *h_range_p = new TH1D( "h_range_p", "Track range", 100, 0, 500 );
  TH1D *h_range_n = new TH1D( "h_range_n", "Track range", 100, 0, 500 );
  
  TH1D *h_nTrk_p  = new TH1D( "h_nTrk_p",  "Number of tracks from a vertex", 8, 0, 8 );
  TH1D *h_nTrk_n  = new TH1D( "h_nTrk_n",  "Number of tracks from a vertex", 8, 0, 8 );
  
  TH1D *h_keSum_p = new TH1D( "h_keSum_p", "Sum of the track kinetic energy out of a vertex", 100, 0, 2000 );
  TH1D *h_keSum_n = new TH1D( "h_keSum_n", "Sum of the track kinetic energy out of a vertex", 100, 0, 2000 );
  
  TH1D *h_maxL_p = new TH1D( "h_maxL_p", "Length of the longest track out of a vertex", 100, 0, 500 );
  TH1D *h_maxL_n = new TH1D( "h_maxL_n", "Length of the longest track out of a vertex", 100, 0, 500 );
  
  TH1D *h_maxKE_p = new TH1D( "h_maxKE_p", "Kinetic energy of the longest track out of a vertex", 100, 0, 2000 );
  TH1D *h_maxKE_n = new TH1D( "h_maxKE_n", "Kinetic energy of the longest track out of a vertex", 100, 0, 2000 );
  
  TH1D *h_maxHits_p = new TH1D( "h_maxHits_p", "Number of hits in the longest track out of a vertex", 200, 0, 200 );
  TH1D *h_maxHits_n = new TH1D( "h_maxHits_n", "Number of hits in the longest track out of a vertex", 200, 0, 200 );
  
  TH1D *h_maxRange_p = new TH1D( "h_maxRange_p", "Range of the longest track out of a vertex", 100, 0, 500 );
  TH1D *h_maxRange_n = new TH1D( "h_maxRange_n", "Range of the longest track out of a vertex", 100, 0, 500 );
  
  // Fill histrograms
  for( int i = 0; i < n_entries_cal; ++i ){
    if( t_primary[i] == 1 ){
    
      h_L_p->Fill( L[i] );
      h_ke_p->Fill( ke[i] );
      h_range_p->Fill( range[i] );
      h_hits_p->Fill( hits[i] );

    }
    else{
    
      h_L_n->Fill( L[i] );
      h_ke_n->Fill( ke[i] );
      h_range_n->Fill( range[i] );
      h_hits_n->Fill( hits[i] );
    }
  }

  for( int i = 0; i < n_entries_pri; ++i ){
    if( v_primary[i] == 1 ){
    
      h_nTrk_p->Fill( nTrk[i] );
      h_keSum_p->Fill( keSum[i] );
      h_maxL_p->Fill( maxL[i] );
      h_maxKE_p->Fill( maxKE[i] );
      h_maxRange_p->Fill( maxRange[i] );
      h_maxHits_p->Fill( maxHits[i] );

    }
    else{
    
      h_nTrk_n->Fill( nTrk[i] );
      h_keSum_n->Fill( keSum[i] );
      h_maxL_n->Fill( maxL[i] );
      h_maxKE_n->Fill( maxKE[i] );
      h_maxRange_n->Fill( maxRange[i] );
      h_maxHits_n->Fill( maxHits[i] );

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

  h_ke_p->SetStats(kFALSE);
  h_ke_p->SetLineColor(2);
  h_ke_p->GetXaxis()->SetTitle( "Track kinetic energy [MeV]" );
  h_ke_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_ke_p = h_ke_p->Integral();
  h_ke_p->Scale(1/scale_ke_p);

  h_ke_n->SetStats(kFALSE);
  h_ke_n->SetLineColor(4);
  h_ke_n->GetXaxis()->SetTitle( "Track kinetic energy [MeV]" );
  h_ke_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_ke_n = h_ke_n->Integral();
  h_ke_n->Scale(1/scale_ke_n);

  h_range_p->SetStats(kFALSE);
  h_range_p->SetLineColor(2);
  h_range_p->GetXaxis()->SetTitle( "Track range [cm]" );
  h_range_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_range_p = h_range_p->Integral();
  h_range_p->Scale(1/scale_range_p);

  h_range_n->SetStats(kFALSE);
  h_range_n->SetLineColor(4);
  h_range_n->GetXaxis()->SetTitle( "Track range [cm]" );
  h_range_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_range_n = h_range_n->Integral();
  h_range_n->Scale(1/scale_range_n);

  h_hits_p->SetStats(kFALSE);
  h_hits_p->SetLineColor(2);
  h_hits_p->GetXaxis()->SetTitle( "Number of hits per track" );
  h_hits_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_hits_p = h_hits_p->Integral();
  h_hits_p->Scale(1/scale_hits_p);

  h_hits_n->SetStats(kFALSE);
  h_hits_n->SetLineColor(4);
  h_hits_n->GetXaxis()->SetTitle( "Number of hits per track" );
  h_hits_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_hits_n = h_hits_n->Integral();
  h_hits_n->Scale(1/scale_hits_n);

  h_nTrk_p->SetStats(kFALSE);
  h_nTrk_p->SetLineColor(2);
  h_nTrk_p->GetXaxis()->SetTitle( "Number of tracks per vertex" );
  h_nTrk_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_nTrk_p = h_nTrk_p->Integral();
  h_nTrk_p->Scale(1/scale_nTrk_p);

  h_nTrk_n->SetStats(kFALSE);
  h_nTrk_n->SetLineColor(4);
  h_nTrk_n->GetXaxis()->SetTitle( "Number of tracks per vertex" );
  h_nTrk_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_nTrk_n = h_nTrk_n->Integral();
  h_nTrk_n->Scale(1/scale_nTrk_n);

  h_keSum_p->SetStats(kFALSE);
  h_keSum_p->SetLineColor(2);
  h_keSum_p->GetXaxis()->SetTitle( "Sum of track kinetic energy from a vertex [MeV]" );
  h_keSum_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_keSum_p = h_keSum_p->Integral();
  h_keSum_p->Scale(1/scale_keSum_p);

  h_keSum_n->SetStats(kFALSE);
  h_keSum_n->SetLineColor(4);
  h_keSum_n->GetXaxis()->SetTitle( "Sum of track kinetic energy from a vertex [MeV]" );
  h_keSum_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_keSum_n = h_keSum_n->Integral();
  h_keSum_n->Scale(1/scale_keSum_n);

  h_maxL_p->SetStats(kFALSE);
  h_maxL_p->SetLineColor(2);
  h_maxL_p->GetXaxis()->SetTitle( "Track length [cm]" );
  h_maxL_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxL_p = h_maxL_p->Integral();
  h_maxL_p->Scale(1/scale_maxL_p);

  h_maxL_n->SetStats(kFALSE);
  h_maxL_n->SetLineColor(4);
  h_maxL_n->GetXaxis()->SetTitle( "Track length [cm]" );
  h_maxL_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxL_n = h_maxL_n->Integral();
  h_maxL_n->Scale(1/scale_maxL_n);

  h_maxKE_p->SetStats(kFALSE);
  h_maxKE_p->SetLineColor(2);
  h_maxKE_p->GetXaxis()->SetTitle( "Track kinetic energy [MeV]" );
  h_maxKE_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxKE_p = h_maxKE_p->Integral();
  h_maxKE_p->Scale(1/scale_maxKE_p);

  h_maxKE_n->SetStats(kFALSE);
  h_maxKE_n->SetLineColor(4);
  h_maxKE_n->GetXaxis()->SetTitle( "Track kinetic energy [MeV]" );
  h_maxKE_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxKE_n = h_maxKE_n->Integral();
  h_maxKE_n->Scale(1/scale_maxKE_n);

  h_maxRange_p->SetStats(kFALSE);
  h_maxRange_p->SetLineColor(2);
  h_maxRange_p->GetXaxis()->SetTitle( "Track range [cm]" );
  h_maxRange_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxRange_p = h_maxRange_p->Integral();
  h_maxRange_p->Scale(1/scale_maxRange_p);

  h_maxRange_n->SetStats(kFALSE);
  h_maxRange_n->SetLineColor(4);
  h_maxRange_n->GetXaxis()->SetTitle( "Track range [cm]" );
  h_maxRange_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxRange_n = h_maxRange_n->Integral();
  h_maxRange_n->Scale(1/scale_maxRange_n);

  h_maxHits_p->SetStats(kFALSE);
  h_maxHits_p->SetLineColor(2);
  h_maxHits_p->GetXaxis()->SetTitle( "Number of hits per track" );
  h_maxHits_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxHits_p = h_maxHits_p->Integral();
  h_maxHits_p->Scale(1/scale_maxHits_p);

  h_maxHits_n->SetStats(kFALSE);
  h_maxHits_n->SetLineColor(4);
  h_maxHits_n->GetXaxis()->SetTitle( "Number of hits per track" );
  h_maxHits_n->GetXaxis()->SetTitleOffset(1.2);
  double scale_maxHits_n = h_maxHits_n->Integral();
  h_maxHits_n->Scale(1/scale_maxHits_n);

  TLegend *l_L        = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_ke       = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_range    = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_hits     = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_keSum    = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_nTrk     = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_maxL     = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_maxKE    = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_maxRange = new TLegend(0.5, 0.65, 0.85, 0.85);
  TLegend *l_maxHits  = new TLegend(0.5, 0.65, 0.85, 0.85);
 
  l_L->AddEntry(        h_L_p,        "Primary",     "l" );
  l_L->AddEntry(        h_L_n,        "Non-primary", "l" );
  
  l_ke->AddEntry(       h_ke_p,       "Primary",     "l" );
  l_ke->AddEntry(       h_ke_n,       "Non-primary", "l" );
  
  l_range->AddEntry(    h_range_p,    "Primary",     "l" );
  l_range->AddEntry(    h_range_n,    "Non-primary", "l" );
  
  l_hits->AddEntry(     h_hits_p,     "Primary",     "l" );
  l_hits->AddEntry(     h_hits_n,     "Non-primary", "l" );
  
  l_keSum->AddEntry(    h_keSum_p,    "Primary",     "l" );
  l_keSum->AddEntry(    h_keSum_n,    "Non-primary", "l" );
  
  l_nTrk->AddEntry(     h_nTrk_p,     "Primary",     "l" );
  l_nTrk->AddEntry(     h_nTrk_n,     "Non-primary", "l" );
  
  l_maxL->AddEntry(     h_maxL_p,     "Primary",     "l" );
  l_maxL->AddEntry(     h_maxL_n,     "Non-primary", "l" );
  
  l_maxKE->AddEntry(    h_maxKE_p,    "Primary",     "l" );
  l_maxKE->AddEntry(    h_maxKE_n,    "Non-primary", "l" );
  
  l_maxRange->AddEntry( h_maxRange_p, "Primary",     "l" );
  l_maxRange->AddEntry( h_maxRange_n, "Non-primary", "l" );
  
  l_maxHits->AddEntry(  h_maxHits_p,  "Primary",     "l" );
  l_maxHits->AddEntry(  h_maxHits_n,  "Non-primary", "l" );
  
  TCanvas *c  = new TCanvas( "c", "Canvas", 800, 600 );

  h_L_p->Draw("hist");
  h_L_n->Draw("hist same");
  l_L->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/track_length.root" );

  c->Clear();

  h_ke_p->Draw("hist");
  h_ke_n->Draw("hist same");
  l_ke->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/track_KE.root" );

  c->Clear();

  h_range_p->Draw("hist");
  h_range_n->Draw("hist same");
  l_range->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/track_range.root" );

  c->Clear();

  h_hits_p->Draw("hist");
  h_hits_n->Draw("hist same");
  l_hits->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/track_nHits.root" );

  c->Clear();

  h_keSum_n->Draw("hist");
  h_keSum_p->Draw("hist same");
  l_keSum->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/vertex_keSum.root" );

  c->Clear();

  h_nTrk_n->Draw("hist");
  h_nTrk_p->Draw("hist same");
  l_nTrk->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/vertex_nTracks.root" );

  c->Clear();

  h_maxL_p->Draw("hist");
  h_maxL_n->Draw("hist same");
  l_maxL->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/max_track_length.root" );

  c->Clear();

  h_maxKE_p->Draw("hist");
  h_maxKE_n->Draw("hist same");
  l_maxKE->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/max_track_KE.root" );

  c->Clear();

  h_maxRange_p->Draw("hist");
  h_maxRange_n->Draw("hist same");
  l_maxRange->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/max_track_range.root" );

  c->Clear();

  h_maxHits_p->Draw("hist");
  h_maxHits_n->Draw("hist same");
  l_maxHits->Draw();

  c->SaveAs( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/plots/max_track_nHits.root" );

  c->Clear();



}
