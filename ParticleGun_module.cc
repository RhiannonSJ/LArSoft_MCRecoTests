////////////////////////////////////////////////////////////////////////
// Class:	     ParticleGun
// Plugin Type: analyzer (art v2_07_03)
// File:        ParticleGun_module.cc
//
// Generated at Mon Aug 14 10:29:45 2017 by Rhiannon Jones using cetskelgen
// from cetlib version v3_00_01.
//
// Module to look at truth and reconstructed information to try and guage
// the reconstruction performance in SBND
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
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"

namespace recotests {
	class ParticleGun;
}


class recotests::ParticleGun : public art::EDAnalyzer {
public:
	explicit ParticleGun(fhicl::ParameterSet const & p);
	// The compiler-generated destructor is fine for non-base
	// classes without bare pointers or other resource use.

	// Plugins should not be copied or assigned.
	ParticleGun(ParticleGun const &) = delete;
	ParticleGun(ParticleGun &&) = delete;
	ParticleGun & operator = (ParticleGun const &) = delete;
	ParticleGun & operator = (ParticleGun &&) = delete;

	// Required functions.
	void analyze(art::Event const & e) override;

	// Selected optional functions.
	void reconfigure(fhicl::ParameterSet const & p) override;
    void beginJob() override;
    void endJob() override;

private:

	// Declare member data here.
    TNtuple *fNt;
    int event_number;
    std::vector< int > events, tracks, vertices;

};


recotests::ParticleGun::ParticleGun(fhicl::ParameterSet const & p):
    EDAnalyzer(p)// ,
    // More initializers here.
{}


void recotests::ParticleGun::analyze(art::Event const & e)
{

    event_number++;

    
	  // Get the truth handle
	  art::Handle< std::vector< simb::MCTruth > > mct_handle;
	  e.getByLabel( "generator", mct_handle );

    // Try and get the vertex positions from each event
    // Get the Vertex handle and print something to check it's working
    art::Handle< std::vector< recob::Vertex > > vtx_handle;
	  e.getByLabel( "pmalgtrackmaker", vtx_handle );

    // Get the Track handle and print something to check it's working
    art::Handle< std::vector< recob::Track > > trk_handle;
	  e.getByLabel( "pmalgtrackmaker", trk_handle );

    // Get something out of the vertices and print
    // Check the validity of the handle
    int vtx_size = vtx_handle->size();
    int trk_size = trk_handle->size();
    int mct_size = mct_handle->size();

    if( vtx_size > 50 ){
     
      events.push_back( event_number );

      tracks.push_back( trk_size );
      
      vertices.push_back( vtx_size );

    }

    fNt->Fill( vtx_size, trk_size );

    /*
    // Open file to append
    std::ofstream file;
    file.open( "/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/distances.txt", std::ios_base::app );
    */

    
    int total_events  = 0;
    int reco_vertices = 0;

    
    if( mct_size && mct_handle.isValid() ){
   
        for( auto & mct : (*mct_handle) ){

            // True neutrino vertex ( primary vertex position)
            double nu_x, nu_y, nu_z;

            int primary_pdg = mct.GetParticle(0).PdgCode();

            nu_x = mct.GetParticle(0).Vx();
            nu_y = mct.GetParticle(0).Vy();
            nu_z = mct.GetParticle(0).Vz();
            
            // Only print the events which occur within the active volume and if they are nu_mu
            if( nu_x < 200 && nu_x > -200 && nu_y < 200 && nu_y > -200 && nu_z < 500 && nu_z > 0 ){

              total_events++;


                if( vtx_size && vtx_handle.isValid() ){

                  reco_vertices++;
    
                    // Vector to hold the square-distances between the true primary vertex and the reconstructed vertices
                    // Find which is the smallest and assume this was reconstructed as the primary
                    std::vector< double > dR, dX, dY, dZ, X, Y, Z;
    
                    dR.clear();
                    dX.clear();
                    dY.clear();
                    dZ.clear();
    
                    for( auto & vtx : (*vtx_handle) ){
                
                        // Access 3D reconstructed vertex position using:  
                        // Array to hold points
                        double xyz[3];
    
                        // Set array to be current vertex position
                        vtx.XYZ(xyz);
    
                        // x,y,z of current vertex
                        double x, y, z;
    
                        x = xyz[0];
                        y = xyz[1];
                        z = xyz[2];
    

                        // Find square distances between true vertex and reconstructed vertex in y-z plane
                        //      x-axis is unreliable
                        double r = sqrt( pow( ( nu_x - x ), 2 ) + pow( ( nu_y - y ), 2 ) + pow( ( nu_z - z ), 2 ) );
                        dR.push_back( r );
                        dX.push_back( x - nu_x ); // Reco - true vertex
                        dY.push_back( y - nu_y );
                        dZ.push_back( z - nu_z );
    
                        X.push_back( x );
                        Y.push_back( y );
                        Z.push_back( z );

                    }
         
                    // Find the minimum value in the vector of distances and append txt file
                    std::vector< double >::iterator result = std::min_element( std::begin( dR ), std::end( dR ) );

                    /*
                    // Define 2D distance and 1D distances for the 'primary' reconstructed vertex
                    double dr = dR[ std::distance( std::begin( dR ), result ) ]; 
                    double dx = dX[ std::distance( std::begin( dR ), result ) ];
                    double dy = dY[ std::distance( std::begin( dR ), result ) ];
                    double dz = dZ[ std::distance( std::begin( dR ), result ) ];
    
                    double x_r = X[ std::distance( std::begin( dR ), result ) ];
                    double y_r = Y[ std::distance( std::begin( dR ), result ) ];
                    double z_r = Z[ std::distance( std::begin( dR ), result ) ];
                    */

                    // loop over vectors, if primary set bool = 1
                    bool primary = false;

                    for( unsigned int i = 0; i < dR.size(); ++i ){
                      
                      if( i == std::distance( std::begin( dR ), result ) ){
                        primary = true;
                      }
            
                      fNt->Fill( primary, dX[i], dY[i], dZ[i], dR[i], X[i], Y[i], Z[i], nu_x, nu_y, nu_z, vtx_size, trk_size, primary_pdg );
                    
                    }
                }
            }
        }
        
  /*      double eff = ( reco_vertices / total_events ) * 100;
        std::cout << " Efficiency : " << eff << std::endl; */
    }

}

void recotests::ParticleGun::beginJob()
{
    
  // Counter to find out event number with > 50 vtx
  event_number = 0;

  // Vector to hold events with this number of vertices
  events.clear();
  tracks.clear();
  vertices.clear();

  // Implementation of optional member function here.
  fNt = new TNtuple( "fNt", "True and reconstructed primary vertex position comparison", "primary:dX:dY:dZ:dR:Xr:Yr:Zr:Xt:Yt:Zt:nVtx:nTrk:pdg");
  //fNt = new TNtuple( "fNt", "True and reconstructed primary vertex position comparison", "nVtx:nTrk");

  fNt->SetDirectory(0);
}

void recotests::ParticleGun::endJob()
{
 
  // Implementation of optional member function here.
    TFile *f = new TFile( "/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/SBN_Workshop/pg_proton_plots/sbn_pg_proton_nt.root", "UPDATE" );
   
    fNt->Write();
    f->Close();

    delete f;
    delete fNt;
}

void recotests::ParticleGun::reconfigure(fhicl::ParameterSet const & /*p*/)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::ParticleGun)
