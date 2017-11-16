////////////////////////////////////////////////////////////////////////
// Class:       CalorimetryVtx
// Plugin Type: analyzer (art v2_08_04)
// File:        CalorimetryVtx_module.cc
//
// Generated at Mon Oct 30 12:13:16 2017 by Rhiannon Jones using cetskelgen
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
#include "lardataobj/AnalysisBase/Calorimetry.h"
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
  class CalorimetryVtx;
}


class recotests::CalorimetryVtx : public art::EDAnalyzer {
public:
  explicit CalorimetryVtx(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CalorimetryVtx(CalorimetryVtx const &) = delete;
  CalorimetryVtx(CalorimetryVtx &&) = delete;
  CalorimetryVtx & operator = (CalorimetryVtx const &) = delete;
  CalorimetryVtx & operator = (CalorimetryVtx &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


recotests::CalorimetryVtx::CalorimetryVtx(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void recotests::CalorimetryVtx::analyze(art::Event const & e)
{
  // Implementation of required member function here.
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

  // Vector sizes
  int trk_size = trk_handle->size();
  int vtx_size = vtx_handle->size();
  int mct_size = mct_handle->size();

  // If at least one vertex has been reconstructed
  if( trk_size && vtx_size && mct_size && mct_handle.isValid() ){

    std::cout << " Good to go!" << std::endl;
    /*
    // Loop over truth
    for( auto & mct : (*mct_handle) ){
    
    }*/
  }

  /*
  // Attempt to access Calorimetry information
  art::FindMany<anab::Calorimetry> fmcal(trk_handle, e, "pmatrackcalo");

  // Access things from the calorimetry
  // Loop over the tracks and access calo information for each track
  for( int i = 0; i < trk_size; ++i ){ 
    if (fmcal.isValid()){
         
      // Define calorimetry handle
      std::vector<const anab::Calorimetry*> cal_handle = fmcal.at(i);
      
      // Loop over calorimetry handle
      for ( size_t j = 0; j < cal_handle.size(); ++j ){
        
        if (!cal_handle[j]) continue;
        
        if (!cal_handle[j]->PlaneID().isValid) continue;
        
        int planenum = cal_handle[j]->PlaneID().Plane;
        
        if (planenum<0||planenum>2) continue;
        
        TrackerData.trkke[i][planenum]    = cal_handle[j]->KineticEnergy();
          TrackerData.trkrange[i][planenum] = cal_handle[j]->Range();
           //For now make the second argument as 13 for muons. 
           TrackerData.trkpitchc[i][planenum]= cal_handle[j] -> TrkPitchC();
           const size_t NHits = cal_handle[j] -> dEdx().size();
           TrackerData.ntrkhits[i][planenum] = (int) NHits;
           if (NHits > TrackerData.GetMaxHitsPerTrack(i, planenum)) {
             // if you get this error, you'll have to increase kMaxTrackHits
             mf::LogError("AnalysisTree:limits")
               << "the " << fTrackModuleLabel[iTracker] << " track #" << i
               << " has " << NHits << " hits on calorimetry plane #" << planenum
               <<", only "
               << TrackerData.GetMaxHitsPerTrack(i, planenum) << " stored in     tree";
           }

*/
}

void recotests::CalorimetryVtx::beginJob()
{
  // Implementation of optional member function here.
}

void recotests::CalorimetryVtx::endJob()
{
  // Implementation of optional member function here.
}

void recotests::CalorimetryVtx::reconfigure(fhicl::ParameterSet const & /*p*/)
{
   // Implementation of optional member function here.
}
DEFINE_ART_MODULE(recotests::CalorimetryVtx)
