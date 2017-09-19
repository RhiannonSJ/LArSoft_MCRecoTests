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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"

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

	// Check validity of the handle
	int mctrack_size = mctrack_handle->size();
	int mct_size     = mct_handle->size();

    // MCTrack
	if( mctrack_handle.isValid() && mctrack_size ){

        // Let's start by simply printing out some information
        for ( auto const& mctrack : (*mctrack_handle)){
    
            std::cout << "'Reco?' PDG Codes : " << mctrack.PdgCode() << std::endl; 
            std::cout << "'Reco?' Track ID  : " << mctrack.TrackID() << std::endl; 

        }

    }

    // Truth stuff
    if( mct_handle.isValid() && mct_size ){

        // Let's start by simply printing out some information
        for ( auto const& mct : (*mct_handle)){
   
            int n_particles = mct.NParticles();
            for ( int i = 0; i < n_particles; ++i ){
                std::cout << "Truth PDG Codes  : " << mct.GetParticle(i).PdgCode() << std::endl; 
                std::cout << "Truth Track ID   : " << mct.GetParticle(i).TrackId() << std::endl; 

            }
        }
    }
}

void recotests::MCRecoComp::reconfigure(fhicl::ParameterSet const & /*p*/)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::MCRecoComp)
