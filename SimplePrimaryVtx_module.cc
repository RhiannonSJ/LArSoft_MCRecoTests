////////////////////////////////////////////////////////////////////////
// Class:       SimplePrimaryVtx
// Plugin Type: analyzer (art v2_08_04)
// File:        SimplePrimaryVtx_module.cc
//
// Generated at Sun Dec  3 14:09:26 2017 by Rhiannon Jones using cetskelgen
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
  class SimplePrimaryVtx;
}


class recotests::SimplePrimaryVtx : public art::EDAnalyzer {
public:
  explicit SimplePrimaryVtx(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimplePrimaryVtx(SimplePrimaryVtx const &) = delete;
  SimplePrimaryVtx(SimplePrimaryVtx &&) = delete;
  SimplePrimaryVtx & operator = (SimplePrimaryVtx const &) = delete;
  SimplePrimaryVtx & operator = (SimplePrimaryVtx &&) = delete;

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

  // Counter initialisation
  int not_fiducial, fiducial_event, numu_topology, reco_correct, reco_wrong, event;

  // nTuple initialisation
  TNtuple *fNt_length;

};


recotests::SimplePrimaryVtx::SimplePrimaryVtx(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

  this->reconfigure(p);

}

void recotests::SimplePrimaryVtx::analyze(art::Event const & e)
{
  event++;
  // Implementation of required member function here.
  // Check the interaction takes place within the fiducial volume
  bool contained_topology = true;

  typedef::std::map< std::vector< int >, int > topology_map;
  
  // Get the MCTruth information 
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel("generator", mct_handle );
  int size = mct_handle->size();
 
  if(mct_handle.isValid() && size) {
  
    // Loop over the truth info
    for(auto const& mct : (*mct_handle)) {
 
      // Check the neutrino came from the beam
      if(mct.Origin() != simb::kBeamNeutrino) continue;
 
      //-------------------------------------------------------------------------
      //   Check the neutrino interaction vertex is within the fiducial volume
      //-------------------------------------------------------------------------
      float nuVtxX = mct.GetNeutrino().Lepton().Vx();
      float nuVtxY = mct.GetNeutrino().Lepton().Vy();
      float nuVtxZ = mct.GetNeutrino().Lepton().Vz();
 
 
      if (    (nuVtxX > (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
           || (nuVtxX < (-m_coordinateOffsetX + m_selectedBorderX)) 
           || (nuVtxY > (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
           || (nuVtxY < (-m_coordinateOffsetY + m_selectedBorderY)) 
           || (nuVtxZ > (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
           || (nuVtxZ < (-m_coordinateOffsetZ + m_selectedBorderZ))){
 
        contained_topology = false;
   
        not_fiducial++;

      }
      else{
      
        fiducial_event++;

      }
      
      //-------------------------------------------------------------------------
      //                           Check the interaction
      //-------------------------------------------------------------------------
 
      // Get the number of particles
      const int n_particles = mct.NParticles();
 
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

          // Not topology event
          contained_topology =  false;

        }
      }
    }
  }

  if( contained_topology ){
  
    numu_topology++;

    // Get handles and make relevant validity checks
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

    // If at least one vertex has been reconstructed
    if( trk_size && vtx_size && size && 
        mct_handle.isValid() && 
        vtx_handle.isValid() && 
        trk_handle.isValid() ){
      
      // Loop over truth
      for( auto & mct : (*mct_handle) ){

        // True neutrino vertex ( primary vertex position)
        double nu_x, nu_y, nu_z;

        // Neutrino vertex
        nu_x = mct.GetNeutrino().Lepton().Vx();
        nu_y = mct.GetNeutrino().Lepton().Vy();
        nu_z = mct.GetNeutrino().Lepton().Vz();
        
        // Vector to hold the square-distances between the true primary vertex 
        // and the reconstructed vertices
        // Find which is the smallest and assume this was reconstructed as the 
        // primary
        // Vectors to hold vertex location information for determining the 
        // true primary
        // Vectors to hold the length of the longest track in the event
        std::vector< double > dR, dX, dY, dZ, X, Y, Z, lengths, KE;
        std::vector< TVector3 > vertices, ends;
        std::vector< int > true_it_vect;

        int reco_it = -1;
        int temp_it = -1;
        
        double current_max = -1;
        double current_KE  = -1;
        double current_dR  = -1;
        double current_dX  = -1;
        double current_dY  = -1;
        double current_dZ  = -1;

        TVector3 current_vtx, current_end;

        dR.clear();
        dX.clear();
        dY.clear();
        dZ.clear();
        KE.clear();
        lengths.clear();
        true_it_vect.clear();
        vertices.clear();
        ends.clear();

        art::FindMany< recob::Track >      ftrk(  vtx_handle, e, "pmalgtrackmaker" );
        art::FindMany< anab::Calorimetry > fcal( trk_handle, e, "pmatrackcalo" );
        
        // Boolean to check whether the vertex is the true reconstructed 
        // primary
 
        bool is_primary = true;

        // Loop over vertices, get track associations and find longest track
        for( int i = 0; i < vtx_size; ++i ){

          art::Ptr< recob::Vertex > vtx( vtx_handle, i );

          // Access 3D reconstructed vertex position using:  
          // Array to hold points
          double xyz[3];

          // Set array to be current vertex position
          vtx->XYZ(xyz);

          // x,y,z of current vertex
          double x, y, z, R;

          x = xyz[0];
          y = xyz[1];
          z = xyz[2];
          R =  sqrt( pow( ( x - nu_x ), 2 ) + pow( ( y - nu_y ), 2 ) + pow( ( z - nu_z ), 2 ) ); 
          
          dR.push_back( R );
          dX.push_back( x - nu_x );
          dY.push_back( y - nu_y );
          dZ.push_back( z - nu_z );
          X.push_back( x );
          Y.push_back( y );
          Z.push_back( z );
          true_it_vect.push_back( i );
       
          // Track - vertex association
          std::vector<const recob::Track*> trk_assn = ftrk.at(i);
          
          // Find longest track associated with each vertex
          for( size_t j = 0; j < trk_assn.size(); ++j ){
          
            TVector3 trk_vtx = trk_assn[j]->Vertex();
            TVector3 trk_end = trk_assn[j]->End();
        
            // Push back large number onto KE vector so that we can cut on 
            // these later and all elements are filled
            KE.push_back( -std::numeric_limits<double>::max() );
            lengths.push_back( trk_assn[j]->Length() );
            vertices.push_back( trk_assn[j]->Vertex() );
            ends.push_back( trk_assn[j]->End() );

            // Calorimetry - track association
            std::vector<const anab::Calorimetry*> cal_assn = fcal.at(j);
          
            // Loop over calorimetry association
            for ( size_t k = 0; k < cal_assn.size(); ++k ){

              if (!cal_assn[k]) continue;
              if (!cal_assn[k]->PlaneID().isValid) continue;
              
              // Get the plane number
              int planenum = cal_assn[k]->PlaneID().Plane;
  
              // Should only be 1,2 or 3
              if (planenum<0||planenum>2) continue;

              if( planenum == 2 ){
                
                KE[j] = cal_assn[k]->KineticEnergy();
               
              }
            }
          }
          
          std::vector< double >::iterator temp_max = std::max_element(       std::begin( lengths ), std::end( lengths ) );
          double temp_max_length                   = lengths[  std::distance( std::begin( lengths ), temp_max ) ];
          double temp_max_length_KE                = KE[       std::distance( std::begin( lengths ), temp_max ) ];
          TVector3 temp_max_vertex                 = vertices[ std::distance( std::begin( lengths ), temp_max ) ];
          TVector3 temp_max_end                    = ends[     std::distance( std::begin( lengths ), temp_max ) ];
        
          if( temp_max_length > current_max ){
        
            current_dR  = R; 
            current_dX  = x - nu_x; 
            current_dY  = y - nu_y; 
            current_dZ  = z - nu_z; 
            current_max = temp_max_length;
            current_KE  = temp_max_length_KE;
            current_vtx = temp_max_vertex;
            current_end = temp_max_end;
            temp_it     = i;

          }
        }
        
        // -------------------------------------------------------------------
        //                The reconstructed primary
        // -------------------------------------------------------------------
        // For the current vertex, get the maximum track length
        double   max_length     = current_max;
        double   max_length_KE  = current_KE;
        double   max_length_dR  = current_dR;
        double   max_length_dX  = current_dX;
        double   max_length_dY  = current_dY;
        double   max_length_dZ  = current_dZ;
        TVector3 max_length_vtx = current_vtx;
        TVector3 max_length_end = current_end;
        reco_it = temp_it;

        double x_t_v, y_t_v, z_t_v, x_t_e, y_t_e, z_t_e;
        
        x_t_v = max_length_vtx[0];
        x_t_e = max_length_end[0];
        y_t_v = max_length_vtx[1];
        y_t_e = max_length_end[1];
        z_t_v = max_length_vtx[2];
        z_t_e = max_length_end[2];
        
        // -------------------------------------------------------------------
        //                The true reconstructed primary
        // -------------------------------------------------------------------
        // 3D distance between true neutrino vertex and closest reconstructed
        // vertex
        std::vector< double >::iterator closest = std::min_element(            std::begin( dR ), std::end( dR ) );
        //double primary_dR                       = dR[           std::distance( std::begin( dR ), closest ) ];
        int true_it                             = true_it_vect[ std::distance( std::begin( dR ), closest ) ];           

        // If the 
        if (    (x_t_v < (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
             && (x_t_v > (-m_coordinateOffsetX + m_selectedBorderX)) 
             && (x_t_e < (m_detectorHalfLengthX - m_coordinateOffsetX - m_selectedBorderX)) 
             && (x_t_e > (-m_coordinateOffsetX + m_selectedBorderX)) 
             && (y_t_v < (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
             && (y_t_v > (-m_coordinateOffsetY + m_selectedBorderY)) 
             && (y_t_e < (m_detectorHalfLengthY - m_coordinateOffsetY - m_selectedBorderY)) 
             && (y_t_e > (-m_coordinateOffsetY + m_selectedBorderY)) 
             && (z_t_v < (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
             && (z_t_v > (-m_coordinateOffsetZ + m_selectedBorderZ))
             && (z_t_e < (m_detectorHalfLengthZ - m_coordinateOffsetZ - m_selectedBorderZ)) 
             && (z_t_e > (-m_coordinateOffsetZ + m_selectedBorderZ))){
              
          // If the iterators match, we have found the primary!
          if( reco_it == true_it ){
          
            reco_correct++;

          } 
          else{
         
            is_primary = false;
            reco_wrong++;

          }

          fNt_length->Fill( is_primary, max_length, max_length_dR, max_length_KE, max_length_dX, max_length_dY, max_length_dZ ); 
        
        } 
      }
    }
  }
}

void recotests::SimplePrimaryVtx::beginJob()
{
  // Implementation of optional member function here.
  // Initialise counters
  event          = 0;
  fiducial_event = 0;
  not_fiducial   = 0;
  numu_topology  = 0;
  reco_correct   = 0;
  reco_wrong     = 0;

  fNt_length = new TNtuple( "fNt_length", "Longest track vertex information", "primary:L:dR:KE:dX:dY:dZ" );
  fNt_length->SetDirectory(0);

}

void recotests::SimplePrimaryVtx::endJob()
{
  // Implementation of optional member function here.
  // Read out some interesting quantities

  std::cout << " ====================================================================== ";
  std::cout << std::endl;

  std::cout << " Number of events                                         : ";
  std::cout << event;
  std::cout << std::endl;

  std::cout << " Number of events within the fiducial volume              : ";
  std::cout << fiducial_event;
  std::cout << std::endl;

  std::cout << " Number of events not within the fiducial volume          : ";
  std::cout << not_fiducial;
  std::cout << std::endl;

  std::cout << " Number of fiducial events of chosen topology             : ";
  std::cout << numu_topology;
  std::cout << std::endl;

  std::cout << " Number of fiducial events of other topologies            : ";
  std::cout << fiducial_event - numu_topology;
  std::cout << std::endl;

  std::cout << " Fraction of fiducial events due to chosen interactions   : ";
  std::cout << numu_topology / double( fiducial_event );
  std::cout << std::endl;

  std::cout << " ---------------------------------------------------------------------- ";
  std::cout << std::endl;

  std::cout << " Longest track corresponds to the vertex        : ";
  std::cout << reco_correct;
  std::cout << " times ";
  std::cout << std::endl;

  std::cout << " Longest track doesn't correspond to the vertex : ";
  std::cout << reco_wrong;
  std::cout << " times ";
  std::cout << std::endl;

  std::cout << " Fraction of correct reconstruction attempts    : ";
  std::cout << reco_correct / double( reco_correct + reco_wrong );
  std::cout << std::endl;

  std::cout << " ---------------------------------------------------------------------- ";
  std::cout << std::endl;

  std::cout << " ====================================================================== ";
  std::cout << std::endl;

  // Initiate file and write the nTuples
  TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_56_00/LArSoft-v06_56_00/srcs/recoperformance/recoperformance/plots/SBN_Workshop/primary_metrics/root/longest_track_info.root", "RECREATE" );

  fNt_length->Write();

  f->Close();

  delete f;
  delete fNt_length;

}

void recotests::SimplePrimaryVtx::reconfigure(fhicl::ParameterSet const & p)
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

DEFINE_ART_MODULE(recotests::SimplePrimaryVtx)
