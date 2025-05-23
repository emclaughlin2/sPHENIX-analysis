/// ---------------------------------------------------------------------------
/*! \file   SLambdaJetHunter.cc
 *  \author Derek Anderson
 *  \date   01.25.2024
 *
 *  A minimal analysis module to find lambda-tagged
 *  jets in pythia events.
 */
/// ---------------------------------------------------------------------------

#define SLAMBDAJETHUNTER_CC

// analysis utilities
#include "SLambdaJetHunter.h"
#include "SLambdaJetHunter.sys.h"
#include "SLambdaJetHunter.ana.h"
#include "SLambdaJetHunterConfig.h"

// make common namespaces implicit
using namespace std;
using namespace fastjet;



namespace SColdQcdCorrelatorAnalysis {

  // ctor/dtor ================================================================

  // --------------------------------------------------------------------------
  //! Module ctor accepting name and debug on/off flag
  // --------------------------------------------------------------------------
  SLambdaJetHunter::SLambdaJetHunter(const string &name, const bool debug) : SubsysReco(name) {

    if (debug) {
      cout << "SLambdaJetHunter::SLambdaJetHunter(string&, bool) Calling ctor" << endl;
    }

  }  // end ctor(string&, bool)



  // --------------------------------------------------------------------------
  //! Module ctor accepting config struct
  // --------------------------------------------------------------------------
  SLambdaJetHunter::SLambdaJetHunter(SLambdaJetHunterConfig& config) : SubsysReco(config.moduleName) {

    m_config = config;
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::SLambdaJetHunter(SLambdaJetHunterConfig&) Calling ctor" << endl;
    }

  }  // end ctor(SLambdaJetHunter&)



  // --------------------------------------------------------------------------
  //! Module dtor
  // --------------------------------------------------------------------------
  SLambdaJetHunter::~SLambdaJetHunter() {

    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::~SLambdaJetHunter() Calling dtor" << endl;
    }

  }  // end dtor



  // f4a methods ==============================================================

  // --------------------------------------------------------------------------
  //! F4A module initialization
  // --------------------------------------------------------------------------
  int SLambdaJetHunter::Init(PHCompositeNode *topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::Init(PHCompositeNode *topNode) Initializing" << endl;
    }

    /* TODO check user fastjet input */

    InitOutput();
    InitTree();
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'Init(PHCompositeNode*)'



  // --------------------------------------------------------------------------
  //! Process an event inside F4A
  // --------------------------------------------------------------------------
  int SLambdaJetHunter::process_event(PHCompositeNode *topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::process_event(PHCompositeNode *topNode) Processing Event" << endl;
    }

    // make sure output containers are emtpy
    ResetOutput();

    // run analysis routines
    GrabEventInfo(topNode);
    FindLambdas(topNode);
    MakeJets(topNode);
    CollectJetOutput(topNode);
    AssociateLambdasToJets(topNode);

    // fill tree and return
    FillOutputTree();
    return Fun4AllReturnCodes::EVENT_OK;

  }  // end 'process_event(PHCompositeNode*)'



  // --------------------------------------------------------------------------
  //! F4A module wind-down
  // --------------------------------------------------------------------------
  int SLambdaJetHunter::End(PHCompositeNode *topNode) {

    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::End(PHCompositeNode *topNode) This is the End..." << endl;
    }

    SaveAndCloseOutput();
    return Fun4AllReturnCodes::EVENT_OK;

  }

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
