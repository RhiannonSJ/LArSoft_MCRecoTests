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
#include "canvas/Persistency/Common/FindMany.h"
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
#include "lardataobj/AnalysisBase/Calorimetry.h"

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
  int primary_found, not_found, event_number, n_tracks, n_flip, no_vtx, one_vtx, more_vtx, n_calo, n_vtx;

  double cut, eff;

  TNtuple *fNt_tracks;
  TNtuple *fNt_primary;
  TNtuple *fNt_calo;

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
  if( trk_size && vtx_size && mct_size && mct_handle.isValid() && vtx_handle.isValid() && trk_handle.isValid() ){
    // Loop over truth
    for( auto & mct : (*mct_handle) ){

      // True neutrino vertex ( primary vertex position)
      double nu_x, nu_y, nu_z;
      int nu_pdg, nu_CCNC;

      nu_CCNC = mct.GetNeutrino().CCNC(); 
      nu_pdg  = mct.GetNeutrino().Nu().PdgCode();
      nu_x = mct.GetNeutrino().Lepton().Vx();
      nu_y = mct.GetNeutrino().Lepton().Vy();
      nu_z = mct.GetNeutrino().Lepton().Vz();

      // Only print the events which occur within the active volume and only CC muon neutrino interactions
      if( nu_x < 200 && nu_x > -200 && nu_y < 200 && nu_y > -200 && nu_z < 500 && nu_z > 0 && nu_pdg == 14 && nu_CCNC == 0 ){
        if( vtx_size ){
  
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
          std::vector< double >   reco_trk_vtx_dists;
          std::vector< double >   reco_trk_end_dists;
          std::vector< double >   true_trk_vtx_dists;
          std::vector< double >   true_trk_end_dists;
          std::vector< double >   trk_lengths;
          std::vector< double >   KE;

          std::vector< TVector3 > vertices;

          // Use track <=> vertex associations to find how many tracks leave a vertex
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
          
          // Initialise vertex and calo handles and the primary and non-primary counters
          art::FindMany< recob::Track > ftrk( vtx_handle, e, "pmalgtrackmaker" );
          art::FindMany< recob::Vertex > fvtx( trk_handle, e, "pmalgtrackmaker" );
          art::FindMany<anab::Calorimetry> fmcal(trk_handle, e, "pmatrackcalo");
            
          // Loop over tracks
          for( int i = 0; i < trk_size; ++i ){

            KE.push_back( -std::numeric_limits<double>::max() );

            // Define calorimetry handle
            std::vector<const anab::Calorimetry*> cal_assn = fmcal.at(i);
     
            // Loop over calorimetry association
            for ( size_t j = 0; j < cal_assn.size(); ++j ){

              if (!cal_assn[j]) continue;
              if (!cal_assn[j]->PlaneID().isValid) continue;
                
              // Get the plane number
              int planenum = cal_assn[j]->PlaneID().Plane;
    
              // Should only be 1,2 or 3
              if (planenum<0||planenum>2) continue;

              if( planenum == 2 ){
                  
                n_calo++;
                  
                // For the collection plane only (ICARUS' best defined plane)
                const size_t NHits = cal_assn[j]->dEdx().size();
                KE[i] = cal_assn[j]->KineticEnergy();
                fNt_calo->Fill(cal_assn[j]->KineticEnergy(), cal_assn[j]->Range(), cal_assn[j]->TrkPitchC(), (int) NHits );

              }
            } 
          } 
            
          // Find the closest track to the tre primary and call it the primary
          std::vector< double >::iterator max = std::max_element( std::begin( trk_lengths ), std::end( trk_lengths ) );

          double max_length                   = trk_lengths[ std::distance( std::begin( trk_lengths ), max ) ];
          double max_length_dist_reco_vtx     = reco_trk_vtx_dists[ std::distance( std::begin( trk_lengths ), max ) ];
          double max_length_dist_reco_end     = reco_trk_end_dists[ std::distance( std::begin( trk_lengths ), max ) ];
          double max_length_dist_true_vtx     = true_trk_vtx_dists[ std::distance( std::begin( trk_lengths ), max ) ];
          double max_length_dist_true_end     = true_trk_end_dists[ std::distance( std::begin( trk_lengths ), max ) ];
          double max_KE                       = KE[ std::distance( std::begin( trk_lengths ), max ) ];

          fNt_tracks->Fill( event_number, max_length, max_length_dist_reco_vtx, max_length_dist_reco_end, max_length_dist_true_vtx, max_length_dist_true_end, vtx_cut_count, max_KE );

          // Count how often we find the true vertex to within 5 cm
          // Also, if primary has been found - count number of tracks corresponding to that reconstructed vertex
          if( ( max_length_dist_true_vtx > max_length_dist_true_end && max_length_dist_true_end < cut ) || ( max_length_dist_true_vtx < max_length_dist_true_end && max_length_dist_true_vtx < cut ) ){
            primary_found++;
          }
          else{
            not_found++;
          }
         
          // Loop over vertices, if the current vertex is within 7 cm of the true primary
          // count the number of tracks as from the primary
          // if not, count them as the number of tracks from other vertices
          for( int i = 0; i < vtx_size; ++i ){

            art::Ptr< recob::Vertex > vtx( vtx_handle, i );

            bool primary = false; 
            int track_counter = 0;
            double KE_sum = 0.;

            // Vertex location
            double xyz[3];
  
            // Set array to be current vertex position
            vtx->XYZ(xyz);
  
            double x, y, z, r;

            x = xyz[0];
            y = xyz[1];
            z = xyz[2];
            r = sqrt( pow( x - x_r, 2 ) + pow( y - y_r, 2 ) + pow( z - z_r, 2 ) );
                
            // Vertexing infomation
            // Define track association 
            std::vector<const recob::Track*> trk_assn = ftrk.at(i);
     
            for( size_t j = 0; j < trk_assn.size(); ++j ){
              
              // x,y,z of current vertex
              double x_t_v, y_t_v, z_t_v, r_t_v;
              TVector3 trk_vtx = trk_assn[j]->Vertex();
  
              x_t_v = trk_vtx[0];
              y_t_v = trk_vtx[1];
              z_t_v = trk_vtx[2];
              r_t_v = sqrt( pow( x_t_v - x, 2 ) + pow( y_t_v - y, 2 ) + pow( z_t_v - z, 2 ) );  
                
              // Find the number of tracks coming out of every vertex
              if( r_t_v < 0.5 ){
                track_counter++;
              }
              
              // Get the calorimetry associated with the tracks associated with each vertex
              std::vector<const anab::Calorimetry*> cal_assn_vtx = fmcal.at(j);
     
              // Loop over calorimetry association
              for ( size_t k = 0; k < cal_assn_vtx.size(); ++k ){

                if (!cal_assn_vtx[k]) continue;
                if (!cal_assn_vtx[k]->PlaneID().isValid) continue;
                
                // Get the plane number
                int planenum = cal_assn_vtx[k]->PlaneID().Plane;
    
                // Should only be 1,2 or 3
                if (planenum<0||planenum>2) continue;

                if( planenum == 2 ){
                  
                  // For the collection plane only (ICARUS' best defined plane)
                  KE_sum += cal_assn_vtx[k]->KineticEnergy();

                }
              } 
            }

            if( r < 0.5 ){
              primary = true;
            }

            // nTuple of counters to plot
            fNt_primary->Fill( primary, track_counter, KE_sum );
          
          }
        }
      }
    } 
  }
}

void recotests::VtxParameters::beginJob()
{

  primary_found = 0;
  not_found     = 0;
  event_number  = 0;
  n_tracks      = 0;
  n_vtx         = 0;
  n_calo        = 0;
  n_flip        = 0;
  no_vtx        = 0;
  one_vtx       = 0;
  more_vtx      = 0;

  cut = 7; // 7 cm

  // Implementation of optional member function here.
  // Define nTuple
  //  reconstructed vertex position, angle in xz plane, 
  //  angle in yz plane, vertex momentum, length
  fNt_tracks       = new TNtuple( "fNt_tracks",  "Track length and distance from a vertex", "evt:L:r_dist_vtx:r_dist_end:t_dist_vtx:t_dist_end:nVtx:KE" );
  fNt_primary      = new TNtuple( "fNt_primary", "Primary vertex metrics",                  "primary:nTrk:keSum" );
  fNt_calo         = new TNtuple( "fNt_calo",    "Track calorimetry information",           "ke:range:pitch:hits" );

  fNt_tracks->SetDirectory(0);
  fNt_primary->SetDirectory(0);
  fNt_calo->SetDirectory(0);


}

void recotests::VtxParameters::endJob()
{

  // Effiency of finding the true primary to within some distance
  eff = primary_found / double(primary_found + not_found );

  std::cout << " ------------------------------------------------- " << std::endl;

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
  std::cout << std::setw(5) << " Percentage of 1 found : " << 100 * ( one_vtx / double( no_vtx + one_vtx + more_vtx ) ) << std::endl; 

  std::cout << " ------------------------------------------------- " << std::endl;

  std::cout << " N, tracks : " << n_tracks << ", vertices : " << n_vtx << ", calorimetry : " << n_calo << std::endl;
  
  std::cout << " ------------------------------------------------- " << std::endl;

  // Implementation of optional member function here.
  // Initiate file and write the nTuples
  TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/primary_vtx_metrics.root", "RECREATE" );

  fNt_tracks->Write();
  fNt_primary->Write();
  fNt_calo->Write();

  f->Close();

  delete f;
  delete fNt_tracks;
  delete fNt_primary;
  delete fNt_calo;
}

void recotests::VtxParameters::reconfigure(fhicl::ParameterSet const & /*p*/)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::VtxParameters)
