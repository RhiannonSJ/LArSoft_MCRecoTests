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
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Initiate nTuple to hold the interesting parameters


};


recotests::VtxParameters::VtxParameters(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void recotests::VtxParameters::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

void recotests::VtxParameters::beginJob()
{
  // Implementation of optional member function here.
}

void recotests::VtxParameters::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recotests::VtxParameters)
