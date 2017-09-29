////////////////////////////////////////////////////////////////////////
// Class:	     MCRecoComp
// Plugin Type: analyzer (art v2_07_03)
// File:        MCRecoComp_module.cc
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
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace recotests {
	class MCRecoComp;
}


class recotests::MCRecoComp : public art::EDAnalyzer {
public:
	explicit MCRecoComp(fhicl::ParameterSet const & p);
	// The compiler-generated destructor is fine for non-base
	// classes without bare pointers or other resource use.

	// Plugins should not be copied or assigned.
	MCRecoComp(MCRecoComp const &) = delete;
	MCRecoComp(MCRecoComp &&) = delete;
	MCRecoComp & operator = (MCRecoComp const &) = delete;
	MCRecoComp & operator = (MCRecoComp &&) = delete;

	// Required functions.
	void analyze(art::Event const & e) override;

	// Selected optional functions.
	void reconfigure(fhicl::ParameterSet const & p) override;

private:

	// Declare member data here.

};


recotests::MCRecoComp::MCRecoComp(fhicl::ParameterSet const & p):
    EDAnalyzer(p)// ,
    // More initializers here.
{}


void recotests::MCRecoComp::analyze(art::Event const & e)
{

	// Get the truth handle
	art::Handle< std::vector< simb::MCTruth > > mct_handle;
	e.getByLabel( "generator", mct_handle );

    // Try and get the vertex positions from each event
    // Get the Vertex handle and print something to check it's working
    art::Handle< std::vector< recob::Vertex > > vtx_handle;
	e.getByLabel( "pmalgtrackmaker", vtx_handle );

    // Get something out of the vertices and print
    // Check the validity of the handle
    int vtx_size = vtx_handle->size();
    int mct_size = mct_handle->size();

    // Open file to append
    std::ofstream file;
    file.open( "/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/distances.txt", std::ios_base::app );

    if( vtx_size && mct_size && mct_handle.isValid() ){
   
        for( auto & mct : (*mct_handle) ){

            // True neutrino vertex ( primary vertex position)
            double nu_x, nu_y, nu_z;

            nu_x = mct.GetNeutrino().Lepton().Vx();
            nu_y = mct.GetNeutrino().Lepton().Vy();
            nu_z = mct.GetNeutrino().Lepton().Vz();

            std::cout << "---------------------------------------------" << std::endl;
            std::cout << " Location of neutrino vertex      : (" << nu_x << ", " << nu_y << ", " << nu_z << ")" << std::endl;

            std::cout << "---------------------------------------------" << std::endl;
            std::cout << " Number of reconstructed vertices : " << vtx_size << std::endl;

            if( vtx_size && vtx_handle.isValid() ){

                // Vector to hold the square-distances between the true primary vertex and the reconstructed vertices
                // Find which is the smallest and assume this was reconstructed as the primary
                std::vector< double > sqdistances; 
               
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

                    std::cout << " Location of reconstructed vertex  : (" << x << ", " << y << ", " << z << ")" << std::endl;

                    // Find square distances between true vertex and reconstructed vertex in y-z plane
                    //      x-axis is unreliable
                    sqdistances.push_back( pow( ( nu_y - y ), 2 ) + pow( ( nu_z - z ), 2 ) );
                }
     
                std::cout << "---------------------------------------------" << std::endl;
                for( unsigned int i = 0; i < sqdistances.size(); ++i ){
                    std::cout << sqdistances[i] << std::endl;
                }
                // Find the minimum value in the vector of distances and append txt file
                std::vector< double >::iterator result = std::min_element( std::begin( sqdistances ), std::end( sqdistances ) );

                double min_dist_2D = sqdistances[ std::distance( std::begin( sqdistances ), result ) ]; 

                // Write the minimum distance in the y-z plane to a .txt file for reading in a ROOT macro
                file << sqrt( min_dist_2D ) << std::endl; 

            }
            std::cout << "---------------------------------------------" << std::endl;
        }

    }



	/*
    //
	// Comparing truth and reco information
	// 
	// Start by printing event info to a txt file
	
	// Get the truth information
	// Get the hit handle
	art::Handle< std::vector< recob::Hit > > hit_handle;
	e.getByLabel("gaushit", hit_handle );
 
	// Get the truth handle
	art::Handle< std::vector< simb::MCTruth > > mct_handle;
	e.getByLabel( "generator", mct_handle );

	// Get the truth track handle
	art::Handle< std::vector< sim::MCTrack > > mctrack_handle;
	e.getByLabel( "pmalgtrackmaker", mctrack_handle );

	// Get the 3D track information
	art::Handle< std::vector< pma::Track3D > > track3d_handle;
	e.getByLabel("pmalgtrackmaker", track3d_handle );

    */
}

void recotests::MCRecoComp::reconfigure(fhicl::ParameterSet const & /*p*/)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::MCRecoComp)
