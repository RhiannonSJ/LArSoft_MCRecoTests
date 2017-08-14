////////////////////////////////////////////////////////////////////////
// Class:       MCRecoComp
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


recotests::MCRecoComp::MCRecoComp(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void recotests::MCRecoComp::analyze(art::Event const & e)
{
  //
  // Comparing
}

void recotests::MCRecoComp::reconfigure(fhicl::ParameterSet const & /*p*/)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::MCRecoComp)
