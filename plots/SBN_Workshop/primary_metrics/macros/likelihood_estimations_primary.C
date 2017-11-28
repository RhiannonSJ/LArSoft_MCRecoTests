#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"

void likelihood_estimations_primary() {

  // Read in the sample files
  // Single interaction
  TFile *f_single = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/cc0pi_primary_vtx_metrics.root" );

  gStyle->SetPalette(57);
  gStyle->SetNumberContours(250);

  /*
  // Set the style for histograms, normalise for shape dependencies 
  h_L_p->SetStats(kFALSE);
  h_L_p->SetLineColor(2);
  h_L_p->GetXaxis()->SetTitle( "Track length [cm]" );
  h_L_p->GetXaxis()->SetTitleOffset(1.2);
  double scale_L_p = h_L_p->Integral();
  h_L_p->Scale(1/scale_L_p);

  */

}
