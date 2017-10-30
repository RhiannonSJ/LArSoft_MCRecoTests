////////////////////////////////////////////////////////////////////////
// Class:       VtxParameters
// Plugin Type: analyzer (art v2_08_03)
// File:        VtxParameters_module.cc
//
// Generated at Mon Oct  9 11:01:52 2017 by Rhiannon Jones using cetskelgen
// from cetlib version v3_01_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"

namespace recotests {
  class VtxParameters;
}


class recotests::VtxParameters : public art::EDAnalyzer {
public:
  explicit VtxParameters(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  VtxParameters(VtxParameters const &) = delete;
  VtxParameters(VtxParameters &&) = delete;
  VtxParameters & operator = (VtxParameters const &) = delete;
  VtxParameters & operator = (VtxParameters &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Initiate nTuple to hold the interesting parameters

  // Counter to find out how many times we reconstruct the primary to within an amount
  int primary_found, not_found, event_number, n_tracks, n_flip, no_vtx, one_vtx, more_vtx;

  double cut, eff;

  TNtuple *fNt_tracks;
  //TNtuple *fNt_primary;
  //TNtuple *fNt_other;

};


recotests::VtxParameters::VtxParameters(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void recotests::VtxParameters::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  event_number++;

  // Get the truth handle
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel( "generator", mct_handle );

  // Get the vertex positions from each event
  // Get the Vertex handle
  art::Handle< std::vector< recob::Vertex > > vtx_handle;
  e.getByLabel( "pmalgtrackmaker", vtx_handle );

  // Get the Track handle
  art::Handle< std::vector< recob::Track > > trk_handle;
  e.getByLabel( "pmalgtrackmaker", trk_handle );

  // Get something out of the vertices and print
  // Check the validity of the handle
  int trk_size = trk_handle->size();
  int vtx_size = vtx_handle->size();
  int mct_size = mct_handle->size();

  // If at least one vertex has been reconstructed
  if( trk_size && vtx_size && mct_size && mct_handle.isValid() ){
  
    // Loop over truth
    for( auto & mct : (*mct_handle) ){

      // True neutrino vertex ( primary vertex position)
      double nu_x, nu_y, nu_z;
      int nu_pdg, nu_CCNC;

      //int nu_int = mct.GetNeutrino().InteractionType();

      nu_CCNC = mct.GetNeutrino().CCNC(); 
      
      nu_pdg  = mct.GetNeutrino().Nu().PdgCode();
      
      nu_x = mct.GetNeutrino().Lepton().Vx();
      nu_y = mct.GetNeutrino().Lepton().Vy();
      nu_z = mct.GetNeutrino().Lepton().Vz();

      // Only print the events which occur within the active volume and only CC muon neutrino interactions
      if( nu_x < 200 && nu_x > -200 && nu_y < 200 && nu_y > -200 && nu_z < 500 && nu_z > 0 && nu_pdg == 14 && nu_CCNC == 0 ){

        if( vtx_size && vtx_handle.isValid() ){
  
          // Vector to hold the square-distances between the true primary vertex and the reconstructed vertices
          // Find which is the smallest and assume this was reconstructed as the primary
          // Vectors to hold vertex location information for determining the true primary
          std::vector< double > dR, dX, dY, dZ, X, Y, Z;

          dR.clear();
          dX.clear();
          dY.clear();
          dZ.clear();
  
          int vtx_cut_count = 0;

          for( auto & vtx : (*vtx_handle) ){
        
            // Access 3D reconstructed vertex position using:  
            // Array to hold points
            double xyz[3];
  
            // Set array to be current vertex position
            vtx.XYZ(xyz);
  
            // x,y,z of current vertex
            double x, y, z, R;
  
            x = xyz[0];
            y = xyz[1];
            z = xyz[2];
  
            // Find square distances between true vertex and reconstructed vertex in y-z plane
            R =  sqrt( pow( ( x - nu_x ), 2 ) + pow( ( y - nu_y ), 2 ) + pow( ( z - nu_z ), 2 ) ); 
            dR.push_back( R );

            X.push_back( x );
            Y.push_back( y );
            Z.push_back( z );

            // If the distance to the true vertex is within the cut value of the true 
            if( R < cut ){
              
              vtx_cut_count++;

            } 
          }
         
          if( vtx_cut_count == 0 ){
          
            no_vtx++;

          }
          else if( vtx_cut_count == 1 ){
            
            one_vtx++;
          
          }
          else{
          
            more_vtx++;

          }

          // Find the minimum value in the vector of distances and append txt file
          std::vector< double >::iterator result  = std::min_element( std::begin( dR ), std::end( dR ) );

          // Primary vertex position to compare with track vertices 
          double x_r = X[ std::distance( std::begin( dR ), result ) ];
          double y_r = Y[ std::distance( std::begin( dR ), result ) ];
          double z_r = Z[ std::distance( std::begin( dR ), result ) ];

          // Find the track vertex distance and corresponding track lengths
          std::vector< double > reco_trk_vtx_dists;
          std::vector< double > reco_trk_end_dists;
          std::vector< double > true_trk_vtx_dists;
          std::vector< double > true_trk_end_dists;
          std::vector< double > trk_lengths;

          for( auto & trk : (*trk_handle) ){

            n_tracks++;

            TVector3 trk_vtx = trk.Vertex();
            TVector3 trk_end = trk.End();
            
            double reco_trk_vtx_dist = sqrt( pow( ( x_r - trk_vtx[0] ), 2 ) + pow( ( y_r - trk_vtx[1] ), 2 ) + pow( ( z_r - trk_vtx[2] ), 2 ) );
            double reco_trk_end_dist = sqrt( pow( ( x_r - trk_end[0] ), 2 ) + pow( ( y_r - trk_end[1] ), 2 ) + pow( ( z_r - trk_end[2] ), 2 ) );
            double true_trk_vtx_dist = sqrt( pow( ( trk_vtx[0] - nu_x ), 2 ) + pow( ( trk_vtx[1] - nu_y ), 2 ) + pow( ( trk_vtx[2] - nu_z ), 2 ) );
            double true_trk_end_dist = sqrt( pow( ( trk_end[0] - nu_x ), 2 ) + pow( ( trk_end[1] - nu_y ), 2 ) + pow( ( trk_end[2] - nu_z ), 2 ) );

            if( true_trk_vtx_dist >  true_trk_end_dist ){
            
              n_flip++;
            
            }

            reco_trk_vtx_dists.push_back( reco_trk_vtx_dist );
            reco_trk_end_dists.push_back( reco_trk_end_dist );
            true_trk_vtx_dists.push_back( true_trk_vtx_dist );
            true_trk_end_dists.push_back( true_trk_end_dist );
            
            trk_lengths.push_back( trk.Length() );

          }

          // Find the closest track to the tre primary and call it the primary
          std::vector< double >::iterator max = std::max_element( std::begin( trk_lengths ), std::end( trk_lengths ) );

          double max_length                   = trk_lengths[ std::distance( std::begin( trk_lengths ), max ) ]; ;
          double max_length_dist_reco_vtx     = reco_trk_vtx_dists[ std::distance( std::begin( trk_lengths ), max ) ]; ;
          double max_length_dist_reco_end     = reco_trk_end_dists[ std::distance( std::begin( trk_lengths ), max ) ]; ;
          double max_length_dist_true_vtx     = true_trk_vtx_dists[ std::distance( std::begin( trk_lengths ), max ) ]; ;
          double max_length_dist_true_end     = true_trk_end_dists[ std::distance( std::begin( trk_lengths ), max ) ]; ;

          fNt_tracks->Fill( event_number, max_length, max_length_dist_reco_vtx, max_length_dist_reco_end, max_length_dist_true_vtx, max_length_dist_true_end, vtx_cut_count );

          // Count how often we find the true vertex to within 5 cm
          if( ( max_length_dist_true_vtx > max_length_dist_true_end && max_length_dist_true_end < cut ) || ( max_length_dist_true_vtx < max_length_dist_true_end && max_length_dist_true_vtx < cut ) ){
            primary_found++;
          }
          else{
            not_found++;
          }
        }
      }
    }
  }
}
          /*
            for( auto & trk : (*trk_handle) ){

            TVector3 trk_vtx = trk.Vertex();
           
            double reco_trk_vtx_dist = sqrt( pow( ( x_r - trk_vtx[0] ), 2 ) + pow( ( y_r - trk_vtx[1] ), 2 ) + pow( (z_r - trk_vtx[2] ), 2 ) );

            // The current track is at the position of a true primary track - it is primary
            if( reco_trk_vtx_dist == trk_vtx_dist ){
            
              // Vertex is the primary
              // p_"" refers to primary
              double p_trk_L, p_thetaxz, p_thetayz, p_theta_to_z, p_dir; 

              ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > p_start_dir;  
              p_start_dir  = trk.StartDirection();
              p_dir        = sqrt( pow( p_start_dir.X(), 2 ) + pow( p_start_dir.Y(), 2 ) + pow( p_start_dir.Z(), 2 ) );
             
              p_trk_L      = trk.Length();
              p_theta_to_z = acos( ( p_start_dir.Z() ) / ( p_dir * sqrt( 2 ) ) );
              p_thetaxz    = atan( p_start_dir.X() / p_start_dir.Z() );
              p_thetayz    = atan( p_start_dir.Z() / p_start_dir.Y() );

              fNt_primary->Fill( trk_vtx[0], trk_vtx[1], trk_vtx[2], p_trk_L, p_thetaxz, p_thetayz, p_theta_to_z );
            
            }
            
            else{
              
              // Vertex is not the primary
              // n_"" refers to non-primary
              double n_trk_L, n_thetaxz, n_thetayz, n_theta_to_z, n_dir;
    
              ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > n_start_dir;  
              n_start_dir  = trk.StartMomentumVector();
              n_dir        = sqrt( pow( n_start_dir.X(), 2 ) + pow( n_start_dir.Y(), 2 ) + pow( n_start_dir.Z(), 2 ) );

              n_trk_L      = trk.Length();
              n_theta_to_z = acos( ( n_start_dir.Z() ) / ( n_dir * sqrt( 2 ) ) );
              n_thetaxz    = atan( n_start_dir.X() / n_start_dir.Z() );
              n_thetayz    = atan( n_start_dir.Z() / n_start_dir.Y() );
  
              double z_diff = trk_vtx[2] - z_r;

              fNt_other->Fill( trk_vtx[0], trk_vtx[1], trk_vtx[2], n_trk_L, n_thetaxz, n_thetayz, n_theta_to_z, z_diff );

            }
          }
        }
      }
    }
  }
}*/

void recotests::VtxParameters::beginJob()
{

  primary_found = 0;
  not_found     = 0;
  event_number  = 0;
  n_tracks      = 0;
  n_flip        = 0;
  no_vtx        = 0;
  one_vtx       = 0;
  more_vtx      = 0;

  cut = 7; // 7 cm

  // Implementation of optional member function here.
  // Define nTuple
  //  reconstructed vertex position, angle in xz plane, 
  //  angle in yz plane, vertex momentum, length
  fNt_tracks       = new TNtuple( "fNt_tracks",  "Track length and distance from a vertex", "evt:L:r_dist_vtx:r_dist_end:t_dist_vtx:t_dist_end:nVtx" );
  //fNt_primary      = new TNtuple( "fNt_primary", "Primary vertex metrics",         "Xr:Yr:Zr:L:thetaxz:thetayz:thetatoz" );
  //fNt_other        = new TNtuple( "fNt_other",   "Non-primary vertex metrics",     "Xr:Yr:Zr:L:thetaxz:thetayz:thetatoz:zdiff" );

  fNt_tracks->SetDirectory(0);
  //fNt_primary->SetDirectory(0);
  //fNt_other->SetDirectory(0);

}

void recotests::VtxParameters::endJob()
{

  // Effiency of finding the true primary to within some distance
  eff = primary_found / double(primary_found + not_found );

  std::cout << std::setw(5) << " Primary found "     << primary_found << " times " << std::endl;
  std::cout << std::setw(5) << " Primary not found " << not_found     << " times " << std::endl;
  std::cout << std::setw(5) << " Efficiency of finding the primary vertex : " << 100 * eff << std::endl;

  std::cout << " ------------------------------------------------- " << std::endl;

  std::cout << std::setw(5) << " Number of tracks : "   << n_tracks << std::endl;
  std::cout << std::setw(5) << " Number to flip : "     << n_flip   << std::endl;
  std::cout << std::setw(5) << " Percentage to flip : " << 100 * ( n_flip / double( n_tracks ) ) << std::endl;
  
  std::cout << " ------------------------------------------------- " << std::endl;
 
  std::cout << std::setw(5) << " No vertices found "       << no_vtx << " times " << std::endl; 
  std::cout << std::setw(5) << " 1 vertex found "          << one_vtx << " times " << std::endl; 
  std::cout << std::setw(5) << " > 1 vertices found "      << more_vtx << " times " << std::endl; 
  std::cout << std::setw(5) << " > 1 vertices found "      << more_vtx << " times " << std::endl; 
  std::cout << std::setw(5) << " Percentage of 1 found : " << 100 * ( one_vtx / double( no_vtx + one_vtx + more_vtx ) ) << std::endl; 

  // Implementation of optional member function here.
  // Initiate file and write the nTuples
  TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/primary_metrics/primary_vtx_metrics.root", "RECREATE" );

  fNt_tracks->Write();
  //fNt_primary->Write();
  //fNt_other->Write();

  f->Close();

  delete f;
  delete fNt_tracks;
  //delete fNt_primary;
  //delete fNt_other;
}

void recotests::VtxParameters::reconfigure(fhicl::ParameterSet const & /*p*/)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::VtxParameters)
