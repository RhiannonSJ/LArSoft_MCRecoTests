////////////////////////////////////////////////////////////////////////
// Class:       TopologyPrimaryVtx
// Plugin Type: analyzer (art v2_08_04)
// File:        TopologyPrimaryVtx_module.cc
//
// Generated at Fri Nov 17 10:44:05 2017 by Rhiannon Jones using cetskelgen
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
#include "lardataobj/RecoBase/PFParticle.h"
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
  class TopologyPrimaryVtx;
}


class recotests::TopologyPrimaryVtx : public art::EDAnalyzer {
public:
  explicit TopologyPrimaryVtx(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TopologyPrimaryVtx(TopologyPrimaryVtx const &) = delete;
  TopologyPrimaryVtx(TopologyPrimaryVtx &&) = delete;
  TopologyPrimaryVtx & operator = (TopologyPrimaryVtx const &) = delete;
  TopologyPrimaryVtx & operator = (TopologyPrimaryVtx &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  std::map< std::vector< int >, int > m_selection;
  float m_detectorHalfLengthX;
  float m_detectorHalfLengthY;
  float m_detectorHalfLengthZ;
  float m_coordinateOffsetX;
  float m_coordinateOffsetY;
  float m_coordinateOffsetZ;
  float m_selectedBorderX;
  float m_selectedBorderY;
  float m_selectedBorderZ;

  // Initialisers
  int cc0pi, other, event_number, n_tracks, no_vtx, one_vtx, more_vtx, n_calo, n_vtx, flip_yes, flip_no, both, primary_reco;

  double cut, eff, primary_reco_cut;

  TNtuple *fNt_primary;
  TNtuple *fNt_calo;

};


recotests::TopologyPrimaryVtx::TopologyPrimaryVtx(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  this->reconfigure(p);

}

void recotests::TopologyPrimaryVtx::analyze(art::Event const & e)
{
  // Implementation of required member function here
  //
  // Boolean for whether or not we have the correct topology within
  // the fiducial volume of the detector.

  event_number++;

  // Initialise to be true, and if any of the counters don't match, 
  // change to false
  bool correct = true;

  typedef::std::map< std::vector< int >, int > topology_map;
 
  /*
  // Loop over everything and check it's working  
  // Show the chosen topology to filter on
  std::cout << "Filtering on topology:" << std::endl;
 
  for ( topology_map::iterator it = m_selection.begin(); 
        it != m_selection.end(); 
        ++it ) {
    
    std::vector< int > pdgVect = it->first;
    int                count   = it->second;
 
    std::cout << "  " << count << " particles with PDG in [";
    
    for ( int pdg : pdgVect ) {
        std::cout << " " << pdg << " ";
    }
    
    std::cout << "]" << std::endl;
 
  }
 */
  // Get the MCTruth information 
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel("generator", mct_handle );
  int size = mct_handle->size();
 
  if(mct_handle.isValid() && size) {
  
    // Loop over the truth info
    for(auto const& mct : (*mct_handle)) {
 
      // Check the neutrino came from the beam
      if(mct.Origin() != simb::kBeamNeutrino) continue;
 
      //-----------------------------------------------------------
      // Check the neutrino vertex is within the fiducial volume
      //-----------------------------------------------------------
      float nuVtxX = mct.GetNeutrino().Nu().Vx();
      float nuVtxY = mct.GetNeutrino().Nu().Vy();
      float nuVtxZ = mct.GetNeutrino().Nu().Vz();
 
 
      if (    (nuVtxX > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
           || (nuVtxX < (-m_coordinateOffsetX + m_selectedBorderX)) 
           || (nuVtxY > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
           || (nuVtxY < (-m_coordinateOffsetY + m_selectedBorderY)) 
           || (nuVtxZ > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
           || (nuVtxZ < (-m_coordinateOffsetZ + m_selectedBorderZ))){
 
        correct = false;
        //continue;
    
      }

 
      //-----------------------------------------------------------
      //            Check the interaction
      //-----------------------------------------------------------
 
      // Get the number of particles
      const int n_particles = mct.NParticles();
      //std::cout << " Number of particles in the event : " << n_particles << std::endl;
 
      // Loop over the map and get the elements
      for( topology_map::iterator it = m_selection.begin(); it != m_selection.end(); ++it ){

        std::vector< int > pdgVect = it->first;
        int                count   = it->second;

        // Initialise a counter for the number within the mct vector
        int particle_counter = 0;

        // Loop over the particles in the event
        for( int i = 0; i < n_particles; ++i ){

          // Find if the pdg code of the current particle is one of the ones in the map
          if( std::find( pdgVect.begin(), pdgVect.end(), mct.GetParticle(i).PdgCode() ) != pdgVect.end() ) ++particle_counter;

        }

        // If the counters don't match, return false
        if( particle_counter != count ){

          // Not cc0pi event
          correct =  false;

        }
        /*
        else{
          std::cout << " PDG codes of particles before filtering : " << std::endl;
          // Loop over the particles in the event
          for( int i = 0; i < n_particles; ++i ){
            std::cout << "           " << mct.GetParticle(i).PdgCode() << std::endl;
        }*/
      }
    }
  }

  // Make sure we are definining the filtered events correctly
  if( correct ){
  
    // Looking at a cc 0pi event so do analysis
    cc0pi++;
  
    // Get the vertex positions from each event
    // Get the Vertex handle
    art::Handle< std::vector< recob::Vertex > > vtx_handle;
    e.getByLabel( "pmalgtrackmaker", vtx_handle );

    // Get the Track handle
    art::Handle< std::vector< recob::Track > > trk_handle;
    e.getByLabel( "pmalgtrackmaker", trk_handle );

    // Get the Vertex handle
    art::Handle< std::vector< recob::PFParticle > > pfp_handle;
    e.getByLabel( "pmalgtrackmaker", pfp_handle );
    
    // Get something out of the vertices and print
    // Check the validity of the handle
    int trk_size = trk_handle->size();
    int vtx_size = vtx_handle->size();
    int pfp_size = pfp_handle->size();

    std::cout << " Event number of CC 0pi : " << event_number << std::endl;
    std::cout << " It has : " << vtx_size << " reconstructed vertices and : " << trk_size << " reconstructed tracks " << std::endl;
    std::cin.get();

    // If at least one vertex has been reconstructed
    if( pfp_size && trk_size && vtx_size && size && mct_handle.isValid() && vtx_handle.isValid() && trk_handle.isValid() && pfp_handle.isValid() ){
      // Loop over truth
      for( auto & mct : (*mct_handle) ){

        // True neutrino vertex ( primary vertex position)
        double nu_x, nu_y, nu_z;

        nu_x = mct.GetNeutrino().Lepton().Vx();
        nu_y = mct.GetNeutrino().Lepton().Vy();
        nu_z = mct.GetNeutrino().Lepton().Vz();

        // Only print the events which occur within the active volume and only CC muon neutrino interactions
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

          // 3D distance between true neutrino vertex and closest reconstructed
          // vertex
          double r_r = dR[ std::distance( std::begin( dR ), result ) ];

          // Continue if r_r < 2.523, a cut value determined by fitting 
          // the distribution of these distances and finding the mean + 3sigma

          if( r_r < primary_reco_cut ){
            
            ++primary_reco;

            // Initialise vertex and calo handles and the primary and non-primary counters
            art::FindMany< recob::Track >      ftrk(  vtx_handle, e, "pmalgtrackmaker" );
            art::FindMany< recob::Vertex >     fvtx(  trk_handle, e, "pmalgtrackmaker" );
            art::FindMany< recob::PFParticle > fpfp(  vtx_handle, e, "pmalgtrackmaker" );
            art::FindMany< recob::PFParticle > fpfpt( trk_handle, e, "pmalgtrackmaker" );
            art::FindMany< anab::Calorimetry > fmcal( trk_handle, e, "pmatrackcalo" );
            

            // Is the track a primary
            bool trk_primary = false;
            std::vector< double > KE;

            // Loop over tracks
            for( int i = 0; i < trk_size; ++i ){

              // Push back large number onto KE vector so that we can cut on 
              // these later and all elements are filled
              KE.push_back( -std::numeric_limits<double>::max() );

              // Get a pointer to the current track
              art::Ptr< recob::Track > trk( trk_handle, i );
              
              // x,y,z of current track's vertex
              double x_t_v, y_t_v, z_t_v, r_t_v, x_t_e, y_t_e, z_t_e, r_t_e;
              TVector3 trk_vtx = trk->Vertex();
              TVector3 trk_end = trk->End();

              x_t_v = trk_vtx[0];
              y_t_v = trk_vtx[1];
              z_t_v = trk_vtx[2];
              x_t_e = trk_end[0];
              y_t_e = trk_end[1];
              z_t_e = trk_end[2];
              r_t_v = sqrt( pow( x_t_v - x_r, 2 ) + pow( y_t_v - y_r, 2 ) + pow( z_t_v - z_r, 2 ) );  
              r_t_e = sqrt( pow( x_t_e - x_r, 2 ) + pow( y_t_e - y_r, 2 ) + pow( z_t_e - z_r, 2 ) );  
                
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
                    
                  // Find if track is associated to reconstructed primary
                  if( r_t_v < 0.5 ){
                   
                    flip_no++;
                    trk_primary  = true;

                  }
                  else if( r_t_e < 0.5 ){
                  
                    flip_yes++;
                    trk_primary = true;
                  
                  }
   
                  else if( r_t_e < 0.5 && r_t_v < 0.5 ){
                    both++;
                  }

                  // For the collection plane only (ICARUS' best defined plane)
                  const size_t NHits = cal_assn[j]->dEdx().size();
                  KE[i] = cal_assn[j]->KineticEnergy();
                  fNt_calo->Fill(trk_primary, cal_assn[j]->KineticEnergy(), cal_assn[j]->Range(), cal_assn[j]->TrkPitchC(), trk->Length(), (int) NHits );

                }
              } 
            } 
              
            // Loop over vertices, if the current vertex is within 7 cm of the true primary
            // count the number of tracks as from the primary
            // if not, count them as the number of tracks from other vertices
            for( int i = 0; i < vtx_size; ++i ){

              art::Ptr< recob::Vertex > vtx( vtx_handle, i );

              bool primary = false; 
              int track_counter = 0;
              double KE_sum = 0.;

              ++n_vtx;

              // Vertex location
              double xyz[3];
    
              // Set array to be current vertex position
              vtx->XYZ(xyz);
    
              double x, y, z, r; // s;

              x = xyz[0];
              y = xyz[1];
              z = xyz[2];
              r = sqrt( pow( x - x_r, 2 ) + pow( y - y_r, 2 ) + pow( z - z_r, 2 ) );
              //s = sqrt( pow( x, 2 ) + pow( y, 2 ) + pow( z, 2 ) );  

              // Vertexing infomation
              // Define track association 
              std::vector<const recob::Track*> trk_assn = ftrk.at(i);
      
              std::vector< double > lengths;
              std::vector< int > js;
              lengths.clear();

              // Find longest track associated with each vertex
              for( size_t j = 0; j < trk_assn.size(); ++j ){
              
                lengths.push_back( trk_assn[j]->Length() );
                js.push_back( j );

              }
             
              // For the current vertex, get the maximum track length
              std::vector< double >::iterator max_l = std::max_element( std::begin( lengths ), std::end( lengths ) );
              double m_length                       = lengths[          std::distance( std::begin( lengths ), max_l ) ];
              unsigned int m_j                      = js[               std::distance( std::begin( lengths ), max_l ) ];

              // Variables for kinetic energy, range and hits
              double longest_track_KE    = std::numeric_limits< double >::lowest();
              double longest_track_range = std::numeric_limits< double >::lowest();
              int longest_track_hits     = std::numeric_limits< int >::lowest();
              
              for( size_t j = 0; j < trk_assn.size(); ++j ){
              
                ++n_tracks;
              
                // x,y,z of current vertex
                double x_t_v, y_t_v, z_t_v, r_t_v, x_t_e, y_t_e, z_t_e, r_t_e;
                TVector3 trk_vtx = trk_assn[j]->Vertex();
                TVector3 trk_end = trk_assn[j]->End();
                
                x_t_v = trk_vtx[0];
                y_t_v = trk_vtx[1];
                z_t_v = trk_vtx[2];
                x_t_e = trk_end[0];
                y_t_e = trk_end[1];
                z_t_e = trk_end[2];
                r_t_v = sqrt( pow( x_t_v - x, 2 ) + pow( y_t_v - y, 2 ) + pow( z_t_v - z, 2 ) );  
                r_t_e = sqrt( pow( x_t_e - x, 2 ) + pow( y_t_e - y, 2 ) + pow( z_t_e - z, 2 ) );  
                
                // Find the number of tracks coming out of every vertex
                if( r_t_v < 0.5 || r_t_e < 0.5 ){
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
                    
                    // If we are looking at the longest track
                    
                    if( j == m_j ){
                 
                      longest_track_hits  = cal_assn_vtx[k]->dEdx().size();
                      longest_track_range = cal_assn_vtx[k]->Range();
                      longest_track_KE    = cal_assn_vtx[k]->KineticEnergy();

                    }
                  }
                } 
              
                std::vector<const recob::PFParticle*> pfp_assn = fpfpt.at(j);
                
                // Loop over pfparticle association
                for ( size_t k = 0; k < pfp_assn.size(); ++k ){

                  std::cout << " Parents: " << pfp_assn[k]->Parent() << std::endl;
                }
              }

              if( r < 0.5 ){
                primary = true;
              }

              // If the doubles exist 
              // nTuple of counters to plot
              fNt_primary->Fill( primary, track_counter, KE_sum, m_length, longest_track_KE, longest_track_hits, longest_track_range, z );
            }
          }
        }
      }
    }
  }
  else{
  
    // Not looking at a cc 0pi event
    other++;
 
  }
}

void recotests::TopologyPrimaryVtx::beginJob()
{
  // Implementation of optional member function here.

  // Counters for topology check
  cc0pi             = 0;
  other             = 0;
  event_number      = 0;
  n_tracks          = 0;
  n_vtx             = 0;
  n_calo            = 0;
  no_vtx            = 0;
  one_vtx           = 0;
  more_vtx          = 0;
  flip_yes          = 0;
  flip_no           = 0;
  both              = 0;
  primary_reco      = 0;

  primary_reco_cut = 2.523;
  cut = 10; // 10 cm

  // Implementation of optional member function here.
  // Define nTuple
  //  reconstructed vertex position, angle in xz plane, 
  //  angle in yz plane, vertex momentum, length
  fNt_primary      = new TNtuple( "fNt_primary", "Primary vertex metrics",                  "primary:nTrk:keSum:maxL:maxKE:maxHits:maxRange:zPos" );
  fNt_calo         = new TNtuple( "fNt_calo",    "Track calorimetry information",           "primary:ke:range:pitch:L:hits" );

  fNt_primary->SetDirectory(0);
  fNt_calo->SetDirectory(0);


}

void recotests::TopologyPrimaryVtx::endJob()
{
  // Implementation of optional member function here.
  
  std::cout << "-------------------------------------" << std::endl;
  std::cout << " Number of CC 0pi events   : " << cc0pi << std::endl;
  std::cout << " Number of other events    : " << other << std::endl;
  std::cout << " Fraction of CC 0pi events : " << cc0pi / double( other + cc0pi ) << std::endl;
  std::cout << "-------------------------------------" << std::endl;

  // Effiency of finding the true primary to within some distance
  eff = primary_reco / double( event_number );

  std::cout << " ------------------------------------------------- " << std::endl;

  std::cout << std::setw(5) << " Primary found "     << primary_reco << " times "   << std::endl;
  std::cout << std::setw(5) << " Primary not found " << event_number - primary_reco << " times " << std::endl;
  std::cout << std::setw(5) << " Efficiency of finding the primary vertex : " << 100 * eff << std::endl;

  std::cout << " ------------------------------------------------- " << std::endl;

  std::cout << std::setw(5) << " Number of tracks : "   << n_tracks   << std::endl;
  std::cout << std::setw(5) << " Number to flip : "     << flip_yes   << std::endl;
  std::cout << std::setw(5) << " Percentage to flip : " << 100 * ( flip_yes / double( n_tracks ) ) << std::endl;
  
  std::cout << " ------------------------------------------------- " << std::endl;
 
  std::cout << std::setw(5) << " No vertices found "       << no_vtx   << " times " << std::endl; 
  std::cout << std::setw(5) << " 1 vertex found "          << one_vtx  << " times " << std::endl; 
  std::cout << std::setw(5) << " > 1 vertices found "      << more_vtx << " times " << std::endl; 
  std::cout << std::setw(5) << " Percentage of 1 found : " << 100 * ( one_vtx / double( no_vtx + one_vtx + more_vtx ) ) << std::endl; 

  std::cout << " ------------------------------------------------- " << std::endl;

  std::cout << " N, tracks : " << n_tracks << ", vertices : " << n_vtx << ", calorimetry : " << n_calo << std::endl;
  
  std::cout << " ------------------------------------------------- " << std::endl;

  std::cout << " Track end closest " << flip_yes << " times, track start closest " << flip_no << " times, both : " << both << std::endl;
  
  std::cout << " ------------------------------------------------- " << std::endl;

  // Implementation of optional member function here.
  // Initiate file and write the nTuples
  TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/cc0pi_primary_vtx_metrics.root", "RECREATE" );

  fNt_primary->Write();
  fNt_calo->Write();

  f->Close();

  delete f;
  delete fNt_primary;
  delete fNt_calo;
}

void recotests::TopologyPrimaryVtx::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.

   std::vector< int > blankVect;
   std::vector< std::vector< int > > input;
 
   std::vector< int > selection1 = p.get< std::vector< int > >("Selection1",        blankVect);
   if ( selection1.size() != 0 ) input.push_back(selection1);
 
   std::vector< int > selection2 = p.get< std::vector< int > >("Selection2",        blankVect);
   if ( selection2.size() != 0 ) input.push_back(selection2);
 
   std::vector< int > selection3 = p.get< std::vector< int > >("Selection3",        blankVect);
   if ( selection3.size() != 0 ) input.push_back(selection3);
 
 
   for ( auto & inputVect : input ) {
     if ( inputVect.size() < 2 ) {
       std::cerr << " Error: Selection vector must have at least 2 elements " <<    std::endl;
       std::cerr << "        First element:     Number of particles of PDG code(s)  specified " << std::endl;
       std::cerr << "        Remaining element: PDG codes to filter on " << std::   endl;
       exit(1);
     }
 
     int count = inputVect[0];
     inputVect.erase( inputVect.begin() );
 
     m_selection.insert( std::make_pair( inputVect, count ) );
   }
 
   m_detectorHalfLengthX = p.get<float>("DetectorHalfLengthX");
   m_detectorHalfLengthY = p.get<float>("DetectorHalfLengthY");
   m_detectorHalfLengthZ = p.get<float>("DetectorHalfLengthZ");
   m_coordinateOffsetX   = p.get<float>("CoordinateOffsetX");
   m_coordinateOffsetY   = p.get<float>("CoordinateOffsetY");
   m_coordinateOffsetZ   = p.get<float>("CoordinateOffsetZ");
   m_selectedBorderX     = p.get<float>("SelectedBorderX");
   m_selectedBorderY     = p.get<float>("SelectedBorderY");
   m_selectedBorderZ     = p.get<float>("SelectedBorderZ");
 

}

DEFINE_ART_MODULE(recotests::TopologyPrimaryVtx)
