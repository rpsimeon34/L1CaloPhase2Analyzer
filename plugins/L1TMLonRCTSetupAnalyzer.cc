/*
 *  \file L1TMLonRCTSetupAnalyzer.cc
 *  Author R. Simeon
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
// #include <cstdint>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <TLorentzVector.h>
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#endif

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/RCT_output.h"


// ECAL TPs
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// HCAL TPs
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

// Output tower collection
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloPFCluster.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1CaloPhase2Analyzer/interface/L1TMLonRCTSetupAnalyzer.h"
#include "DataFormats/Math/interface/deltaR.h"


// ECAL propagation
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"

//////////////////////////////////////////////////////////////////////////////////////

using namespace edm;

L1TMLonRCTSetupAnalyzer::L1TMLonRCTSetupAnalyzer( const ParameterSet & cfg ) :
  decoderToken_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
  caloGeometryToken_(esConsumes<CaloGeometry, CaloGeometryRecord>(edm::ESInputTag("", ""))),
  hbTopologyToken_(esConsumes<HcalTopology, HcalRecNumberingRecord>(edm::ESInputTag("", ""))),
  EGammaSLR3Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("EGammaSLR3"))),
  EGammaSLR2Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("EGammaSLR2"))),
  EGammaSLR1Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("EGammaSLR1"))),
  EGammaSLR0Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("EGammaSLR0"))),
  ECALUnclusteredSLR3Src_ (consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("ECALUnclusteredSLR3"))),
  ECALUnclusteredSLR2Src_ (consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("ECALUnclusteredSLR2"))),
  ECALUnclusteredSLR1Src_ (consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("ECALUnclusteredSLR1"))),
  ECALUnclusteredSLR0Src_ (consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("ECALUnclusteredSLR0"))),
  HCAL8Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("HCAL8"))),
  HCAL7Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("HCAL7"))),
  HCAL6Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("HCAL6"))),
  HCAL5Src_(consumes<l1tp2::rctOutputLinkCollection>(cfg.getParameter<edm::InputTag>("HCAL5")))
{
    folderName_          = cfg.getUntrackedParameter<std::string>("folderName");

    linkTree = tfs_->make<TTree>("linkTree", "ML on RCT Link Tree");

    linkTree->Branch("run",    &run,     "run/I");
    linkTree->Branch("lumi",   &lumi,    "lumi/I");
    linkTree->Branch("event",  &event,   "event/I");
    
    ////putting bufsize at 32000 and changing split level to 0 so that the branch isn't split into multiple branches

    // Cluster branches SLR3
    linkTree->Branch("SLR3_cluster_seed_energy", "vector<vector<int>>", &SLR3_cluster_seed_energy, 32000, 0);
    linkTree->Branch("SLR3_cluster_energy", "vector<vector<int>>", &SLR3_cluster_energy, 32000, 0);
    linkTree->Branch("SLR3_cluster_eta", "vector<vector<int>>", &SLR3_cluster_eta, 32000, 0);
    linkTree->Branch("SLR3_cluster_phi", "vector<vector<int>>", &SLR3_cluster_phi, 32000, 0);
    linkTree->Branch("SLR3_cluster_et5x5", "vector<vector<int>>", &SLR3_cluster_et5x5, 32000, 0);
    linkTree->Branch("SLR3_cluster_et2x5", "vector<vector<int>>", &SLR3_cluster_et2x5, 32000, 0);
    linkTree->Branch("SLR3_cluster_timing", "vector<vector<int>>", &SLR3_cluster_timing, 32000, 0);
    linkTree->Branch("SLR3_cluster_spike", "vector<vector<int>>", &SLR3_cluster_spike, 32000, 0);
    linkTree->Branch("SLR3_cluster_satur", "vector<vector<int>>", &SLR3_cluster_satur, 32000, 0);
    linkTree->Branch("SLR3_cluster_brems", "vector<vector<int>>", &SLR3_cluster_brems, 32000, 0);
    linkTree->Branch("SLR3_cluster_spare", "vector<vector<int>>", &SLR3_cluster_spare, 32000, 0);

    // Cluster branches SLR2
    linkTree->Branch("SLR2_cluster_seed_energy", "vector<vector<int>>", &SLR2_cluster_seed_energy, 32000, 0);
    linkTree->Branch("SLR2_cluster_energy", "vector<vector<int>>", &SLR2_cluster_energy, 32000, 0);
    linkTree->Branch("SLR2_cluster_eta", "vector<vector<int>>", &SLR2_cluster_eta, 32000, 0);
    linkTree->Branch("SLR2_cluster_phi", "vector<vector<int>>", &SLR2_cluster_phi, 32000, 0);
    linkTree->Branch("SLR2_cluster_et5x5", "vector<vector<int>>", &SLR2_cluster_et5x5, 32000, 0);
    linkTree->Branch("SLR2_cluster_et2x5", "vector<vector<int>>", &SLR2_cluster_et2x5, 32000, 0);
    linkTree->Branch("SLR2_cluster_timing", "vector<vector<int>>", &SLR2_cluster_timing, 32000, 0);
    linkTree->Branch("SLR2_cluster_spike", "vector<vector<int>>", &SLR2_cluster_spike, 32000, 0);
    linkTree->Branch("SLR2_cluster_satur", "vector<vector<int>>", &SLR2_cluster_satur, 32000, 0);
    linkTree->Branch("SLR2_cluster_brems", "vector<vector<int>>", &SLR2_cluster_brems, 32000, 0);
    linkTree->Branch("SLR2_cluster_spare", "vector<vector<int>>", &SLR2_cluster_spare, 32000, 0);

    // Cluster branches SLR1
    linkTree->Branch("SLR1_cluster_seed_energy", "vector<vector<int>>", &SLR1_cluster_seed_energy, 32000, 0);
    linkTree->Branch("SLR1_cluster_energy", "vector<vector<int>>", &SLR1_cluster_energy, 32000, 0);
    linkTree->Branch("SLR1_cluster_eta", "vector<vector<int>>", &SLR1_cluster_eta, 32000, 0);
    linkTree->Branch("SLR1_cluster_phi", "vector<vector<int>>", &SLR1_cluster_phi, 32000, 0);
    linkTree->Branch("SLR1_cluster_et5x5", "vector<vector<int>>", &SLR1_cluster_et5x5, 32000, 0);
    linkTree->Branch("SLR1_cluster_et2x5", "vector<vector<int>>", &SLR1_cluster_et2x5, 32000, 0);
    linkTree->Branch("SLR1_cluster_timing", "vector<vector<int>>", &SLR1_cluster_timing, 32000, 0);
    linkTree->Branch("SLR1_cluster_spike", "vector<vector<int>>", &SLR1_cluster_spike, 32000, 0);
    linkTree->Branch("SLR1_cluster_satur", "vector<vector<int>>", &SLR1_cluster_satur, 32000, 0);
    linkTree->Branch("SLR1_cluster_brems", "vector<vector<int>>", &SLR1_cluster_brems, 32000, 0);
    linkTree->Branch("SLR1_cluster_spare", "vector<vector<int>>", &SLR1_cluster_spare, 32000, 0);

    // Cluster branches SLR0
    linkTree->Branch("SLR0_cluster_seed_energy", "vector<vector<int>>", &SLR0_cluster_seed_energy, 32000, 0);
    linkTree->Branch("SLR0_cluster_energy", "vector<vector<int>>", &SLR0_cluster_energy, 32000, 0);
    linkTree->Branch("SLR0_cluster_eta", "vector<vector<int>>", &SLR0_cluster_eta, 32000, 0);
    linkTree->Branch("SLR0_cluster_phi", "vector<vector<int>>", &SLR0_cluster_phi, 32000, 0);
    linkTree->Branch("SLR0_cluster_et5x5", "vector<vector<int>>", &SLR0_cluster_et5x5, 32000, 0);
    linkTree->Branch("SLR0_cluster_et2x5", "vector<vector<int>>", &SLR0_cluster_et2x5, 32000, 0);
    linkTree->Branch("SLR0_cluster_timing", "vector<vector<int>>", &SLR0_cluster_timing, 32000, 0);
    linkTree->Branch("SLR0_cluster_spike", "vector<vector<int>>", &SLR0_cluster_spike, 32000, 0);
    linkTree->Branch("SLR0_cluster_satur", "vector<vector<int>>", &SLR0_cluster_satur, 32000, 0);
    linkTree->Branch("SLR0_cluster_brems", "vector<vector<int>>", &SLR0_cluster_brems, 32000, 0);
    linkTree->Branch("SLR0_cluster_spare", "vector<vector<int>>", &SLR0_cluster_spare, 32000, 0);

    // ECAL unclustered energy tower branches SLR3
    linkTree->Branch("ECALUnclusteredSLR3_tower_et", "vector<vector<int>>", &ECALUnclusteredSLR3_tower_et, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR3_tower_eta", "vector<vector<int>>", &ECALUnclusteredSLR3_tower_eta, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR3_tower_phi", "vector<vector<int>>", &ECALUnclusteredSLR3_tower_phi, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR3_tower_timing", "vector<vector<int>>", &ECALUnclusteredSLR3_tower_timing, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR3_tower_spike", "vector<vector<int>>", &ECALUnclusteredSLR3_tower_spike, 32000, 0);

    // ECAL unclustered energy tower branches SLR2
    linkTree->Branch("ECALUnclusteredSLR2_tower_et", "vector<vector<int>>", &ECALUnclusteredSLR2_tower_et, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR2_tower_eta", "vector<vector<int>>", &ECALUnclusteredSLR2_tower_eta, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR2_tower_phi", "vector<vector<int>>", &ECALUnclusteredSLR2_tower_phi, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR2_tower_timing", "vector<vector<int>>", &ECALUnclusteredSLR2_tower_timing, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR2_tower_spike", "vector<vector<int>>", &ECALUnclusteredSLR2_tower_spike, 32000, 0);

    // ECAL unclustered energy tower branches SLR1
    linkTree->Branch("ECALUnclusteredSLR1_tower_et", "vector<vector<int>>", &ECALUnclusteredSLR1_tower_et, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR1_tower_eta", "vector<vector<int>>", &ECALUnclusteredSLR1_tower_eta, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR1_tower_phi", "vector<vector<int>>", &ECALUnclusteredSLR1_tower_phi, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR1_tower_timing", "vector<vector<int>>", &ECALUnclusteredSLR1_tower_timing, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR1_tower_spike", "vector<vector<int>>", &ECALUnclusteredSLR1_tower_spike, 32000, 0);

    // ECAL unclustered energy tower branches SLR0
    linkTree->Branch("ECALUnclusteredSLR0_tower_et", "vector<vector<int>>", &ECALUnclusteredSLR0_tower_et, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR0_tower_eta", "vector<vector<int>>", &ECALUnclusteredSLR0_tower_eta, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR0_tower_phi", "vector<vector<int>>", &ECALUnclusteredSLR0_tower_phi, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR0_tower_timing", "vector<vector<int>>", &ECALUnclusteredSLR0_tower_timing, 32000, 0);
    linkTree->Branch("ECALUnclusteredSLR0_tower_spike", "vector<vector<int>>", &ECALUnclusteredSLR0_tower_spike, 32000, 0);

    // HCAL link 8 tower branches
    linkTree->Branch("HCAL8_tower_et", "vector<vector<int>>", &HCAL8_tower_et, 32000, 0);
    linkTree->Branch("HCAL8_tower_eta", "vector<vector<int>>", &HCAL8_tower_eta, 32000, 0);
    linkTree->Branch("HCAL8_tower_phi", "vector<vector<int>>", &HCAL8_tower_phi, 32000, 0);
    linkTree->Branch("HCAL8_tower_fb", "vector<vector<int>>", &HCAL8_tower_fb, 32000, 0);

    // HCAL link 7 tower branches
    linkTree->Branch("HCAL7_tower_et", "vector<vector<int>>", &HCAL7_tower_et, 32000, 0);
    linkTree->Branch("HCAL7_tower_eta", "vector<vector<int>>", &HCAL7_tower_eta, 32000, 0);
    linkTree->Branch("HCAL7_tower_phi", "vector<vector<int>>", &HCAL7_tower_phi, 32000, 0);
    linkTree->Branch("HCAL7_tower_fb", "vector<vector<int>>", &HCAL7_tower_fb, 32000, 0);

    // HCAL link 6 tower branches
    linkTree->Branch("HCAL6_tower_et", "vector<vector<int>>", &HCAL6_tower_et, 32000, 0);
    linkTree->Branch("HCAL6_tower_eta", "vector<vector<int>>", &HCAL6_tower_eta, 32000, 0);
    linkTree->Branch("HCAL6_tower_phi", "vector<vector<int>>", &HCAL6_tower_phi, 32000, 0);
    linkTree->Branch("HCAL6_tower_fb", "vector<vector<int>>", &HCAL6_tower_fb, 32000, 0);

    // HCAL link 5 tower branches
    linkTree->Branch("HCAL5_tower_et", "vector<vector<int>>", &HCAL5_tower_et, 32000, 0);
    linkTree->Branch("HCAL5_tower_eta", "vector<vector<int>>", &HCAL5_tower_eta, 32000, 0);
    linkTree->Branch("HCAL5_tower_phi", "vector<vector<int>>", &HCAL5_tower_phi, 32000, 0);
    linkTree->Branch("HCAL5_tower_fb", "vector<vector<int>>", &HCAL5_tower_fb, 32000, 0);

  }

void L1TMLonRCTSetupAnalyzer::beginJob( const EventSetup & es) {
}

void L1TMLonRCTSetupAnalyzer::analyze( const Event& evt, const EventSetup& es )
 {

  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();

  // Set up handles to read in links
  edm::Handle<l1tp2::rctOutputLinkCollection> EGammaSLR3_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> EGammaSLR2_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> EGammaSLR1_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> EGammaSLR0_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> ECALUnclusteredSLR3_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> ECALUnclusteredSLR2_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> ECALUnclusteredSLR1_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> ECALUnclusteredSLR0_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> HCAL8_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> HCAL7_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> HCAL6_link;
  edm::Handle<l1tp2::rctOutputLinkCollection> HCAL5_link;

  // Clear cluster vectors
  // SLR3
  SLR3_cluster_seed_energy->clear();
  SLR3_cluster_energy->clear();
  SLR3_cluster_eta->clear();
  SLR3_cluster_phi->clear();
  SLR3_cluster_et5x5->clear();
  SLR3_cluster_et2x5->clear();
  SLR3_cluster_timing->clear();
  SLR3_cluster_spike->clear();
  SLR3_cluster_satur->clear();
  SLR3_cluster_brems->clear();
  SLR3_cluster_spare->clear();
  // SLR2
  SLR2_cluster_seed_energy->clear();
  SLR2_cluster_energy->clear();
  SLR2_cluster_eta->clear();
  SLR2_cluster_phi->clear();
  SLR2_cluster_et5x5->clear();
  SLR2_cluster_et2x5->clear();
  SLR2_cluster_timing->clear();
  SLR2_cluster_spike->clear();
  SLR2_cluster_satur->clear();
  SLR2_cluster_brems->clear();
  SLR2_cluster_spare->clear();
  // SLR1
  SLR1_cluster_seed_energy->clear();
  SLR1_cluster_energy->clear();
  SLR1_cluster_eta->clear();
  SLR1_cluster_phi->clear();
  SLR1_cluster_et5x5->clear();
  SLR1_cluster_et2x5->clear();
  SLR1_cluster_timing->clear();
  SLR1_cluster_spike->clear();
  SLR1_cluster_satur->clear();
  SLR1_cluster_brems->clear();
  SLR1_cluster_spare->clear();
  // SLR0
  SLR0_cluster_seed_energy->clear();
  SLR0_cluster_energy->clear();
  SLR0_cluster_eta->clear();
  SLR0_cluster_phi->clear();
  SLR0_cluster_et5x5->clear();
  SLR0_cluster_et2x5->clear();
  SLR0_cluster_timing->clear();
  SLR0_cluster_spike->clear();
  SLR0_cluster_satur->clear();
  SLR0_cluster_brems->clear();
  SLR0_cluster_spare->clear();

  // Clear ECAL unclustered energy vectors
  // SLR3
  ECALUnclusteredSLR3_tower_et->clear();
  ECALUnclusteredSLR3_tower_eta->clear();
  ECALUnclusteredSLR3_tower_phi->clear();
  ECALUnclusteredSLR3_tower_timing->clear();
  ECALUnclusteredSLR3_tower_spike->clear();
  // SLR2
  ECALUnclusteredSLR2_tower_et->clear();
  ECALUnclusteredSLR2_tower_eta->clear();
  ECALUnclusteredSLR2_tower_phi->clear();
  ECALUnclusteredSLR2_tower_timing->clear();
  ECALUnclusteredSLR2_tower_spike->clear();
  // SLR1
  ECALUnclusteredSLR1_tower_et->clear();
  ECALUnclusteredSLR1_tower_eta->clear();
  ECALUnclusteredSLR1_tower_phi->clear();
  ECALUnclusteredSLR1_tower_timing->clear();
  ECALUnclusteredSLR1_tower_spike->clear();
  // SLR0
  ECALUnclusteredSLR0_tower_et->clear();
  ECALUnclusteredSLR0_tower_eta->clear();
  ECALUnclusteredSLR0_tower_phi->clear();
  ECALUnclusteredSLR0_tower_timing->clear();
  ECALUnclusteredSLR0_tower_spike->clear();

  // Clear HCAL vectors
  // HCAL link 8
  HCAL8_tower_et->clear();
  HCAL8_tower_eta->clear();
  HCAL8_tower_phi->clear();
  HCAL8_tower_fb->clear();
  // HCAL link 7
  HCAL7_tower_et->clear();
  HCAL7_tower_eta->clear();
  HCAL7_tower_phi->clear();
  HCAL7_tower_fb->clear();
  // HCAL link 6
  HCAL6_tower_et->clear();
  HCAL6_tower_eta->clear();
  HCAL6_tower_phi->clear();
  HCAL6_tower_fb->clear();
  // HCAL link 5
  HCAL5_tower_et->clear();
  HCAL5_tower_eta->clear();
  HCAL5_tower_phi->clear();
  HCAL5_tower_fb->clear();

  // Read out SLR3 EGamma clusters
  if(evt.getByToken(EGammaSLR3Src_, EGammaSLR3_link)){
    for(const auto & link : *EGammaSLR3_link){

      getEGammaClusters(
        link.data(),
        3,
        RCT_cluster_seed_energy,
        RCT_cluster_energy,
        RCT_cluster_eta,
        RCT_cluster_phi,
        RCT_cluster_et5x5,
        RCT_cluster_et2x5,
        RCT_cluster_timing,
        RCT_cluster_spike,
        RCT_cluster_satur,
        RCT_cluster_brems,
        RCT_cluster_spare
      );
      SLR3_cluster_seed_energy->push_back(*RCT_cluster_seed_energy);
      SLR3_cluster_energy->push_back(*RCT_cluster_energy);
      SLR3_cluster_eta->push_back(*RCT_cluster_eta);
      SLR3_cluster_phi->push_back(*RCT_cluster_phi);
      SLR3_cluster_et5x5->push_back(*RCT_cluster_et5x5);
      SLR3_cluster_et2x5->push_back(*RCT_cluster_et2x5);
      SLR3_cluster_timing->push_back(*RCT_cluster_timing);
      SLR3_cluster_spike->push_back(*RCT_cluster_spike);
      SLR3_cluster_satur->push_back(*RCT_cluster_satur);
      SLR3_cluster_brems->push_back(*RCT_cluster_brems);
      SLR3_cluster_spare->push_back(*RCT_cluster_spare);
    }
  }
  // Read out SLR2 EGamma clusters
  if(evt.getByToken(EGammaSLR2Src_, EGammaSLR2_link)){
    for(const auto & link : *EGammaSLR2_link){

      getEGammaClusters(
        link.data(),
        2,
        RCT_cluster_seed_energy,
        RCT_cluster_energy,
        RCT_cluster_eta,
        RCT_cluster_phi,
        RCT_cluster_et5x5,
        RCT_cluster_et2x5,
        RCT_cluster_timing,
        RCT_cluster_spike,
        RCT_cluster_satur,
        RCT_cluster_brems,
        RCT_cluster_spare
      );
      SLR2_cluster_seed_energy->push_back(*RCT_cluster_seed_energy);
      SLR2_cluster_energy->push_back(*RCT_cluster_energy);
      SLR2_cluster_eta->push_back(*RCT_cluster_eta);
      SLR2_cluster_phi->push_back(*RCT_cluster_phi);
      SLR2_cluster_et5x5->push_back(*RCT_cluster_et5x5);
      SLR2_cluster_et2x5->push_back(*RCT_cluster_et2x5);
      SLR2_cluster_timing->push_back(*RCT_cluster_timing);
      SLR2_cluster_spike->push_back(*RCT_cluster_spike);
      SLR2_cluster_satur->push_back(*RCT_cluster_satur);
      SLR2_cluster_brems->push_back(*RCT_cluster_brems);
      SLR2_cluster_spare->push_back(*RCT_cluster_spare);
    }
  }
  // Read out SLR1 EGamma clusters
  if(evt.getByToken(EGammaSLR1Src_, EGammaSLR1_link)){
    for(const auto & link : *EGammaSLR1_link){

      getEGammaClusters(
        link.data(),
        1,
        RCT_cluster_seed_energy,
        RCT_cluster_energy,
        RCT_cluster_eta,
        RCT_cluster_phi,
        RCT_cluster_et5x5,
        RCT_cluster_et2x5,
        RCT_cluster_timing,
        RCT_cluster_spike,
        RCT_cluster_satur,
        RCT_cluster_brems,
        RCT_cluster_spare
      );
      SLR1_cluster_seed_energy->push_back(*RCT_cluster_seed_energy);
      SLR1_cluster_energy->push_back(*RCT_cluster_energy);
      SLR1_cluster_eta->push_back(*RCT_cluster_eta);
      SLR1_cluster_phi->push_back(*RCT_cluster_phi);
      SLR1_cluster_et5x5->push_back(*RCT_cluster_et5x5);
      SLR1_cluster_et2x5->push_back(*RCT_cluster_et2x5);
      SLR1_cluster_timing->push_back(*RCT_cluster_timing);
      SLR1_cluster_spike->push_back(*RCT_cluster_spike);
      SLR1_cluster_satur->push_back(*RCT_cluster_satur);
      SLR1_cluster_brems->push_back(*RCT_cluster_brems);
      SLR1_cluster_spare->push_back(*RCT_cluster_spare);
    }
  }
  // Read out SLR0 EGamma clusters
  if(evt.getByToken(EGammaSLR0Src_, EGammaSLR0_link)){
    for(const auto & link : *EGammaSLR0_link){

      getEGammaClusters(
        link.data(),
        0,
        RCT_cluster_seed_energy,
        RCT_cluster_energy,
        RCT_cluster_eta,
        RCT_cluster_phi,
        RCT_cluster_et5x5,
        RCT_cluster_et2x5,
        RCT_cluster_timing,
        RCT_cluster_spike,
        RCT_cluster_satur,
        RCT_cluster_brems,
        RCT_cluster_spare
      );
      SLR0_cluster_seed_energy->push_back(*RCT_cluster_seed_energy);
      SLR0_cluster_energy->push_back(*RCT_cluster_energy);
      SLR0_cluster_eta->push_back(*RCT_cluster_eta);
      SLR0_cluster_phi->push_back(*RCT_cluster_phi);
      SLR0_cluster_et5x5->push_back(*RCT_cluster_et5x5);
      SLR0_cluster_et2x5->push_back(*RCT_cluster_et2x5);
      SLR0_cluster_timing->push_back(*RCT_cluster_timing);
      SLR0_cluster_spike->push_back(*RCT_cluster_spike);
      SLR0_cluster_satur->push_back(*RCT_cluster_satur);
      SLR0_cluster_brems->push_back(*RCT_cluster_brems);
      SLR0_cluster_spare->push_back(*RCT_cluster_spare);
    }
  }

  // Read out SLR3 ECAL unclustered energy
  if(evt.getByToken(ECALUnclusteredSLR3Src_, ECALUnclusteredSLR3_link)){
    for(const auto & link : *ECALUnclusteredSLR3_link){

      getECALUnclusteredEnergy(
        link.data(),
        3,
        RCT_ECAL_tower_et,
        RCT_ECAL_tower_eta,
        RCT_ECAL_tower_phi,
        RCT_ECAL_tower_timing,
        RCT_ECAL_tower_spike
      );
      ECALUnclusteredSLR3_tower_et->push_back(*RCT_ECAL_tower_et);
      ECALUnclusteredSLR3_tower_eta->push_back(*RCT_ECAL_tower_eta);
      ECALUnclusteredSLR3_tower_phi->push_back(*RCT_ECAL_tower_phi);
      ECALUnclusteredSLR3_tower_timing->push_back(*RCT_ECAL_tower_timing);
      ECALUnclusteredSLR3_tower_spike->push_back(*RCT_ECAL_tower_spike);
    }
  }
  // Read out SLR2 ECAL unclustered energy
  if(evt.getByToken(ECALUnclusteredSLR2Src_, ECALUnclusteredSLR2_link)){
    for(const auto & link : *ECALUnclusteredSLR2_link){

      getECALUnclusteredEnergy(
        link.data(),
        2,
        RCT_ECAL_tower_et,
        RCT_ECAL_tower_eta,
        RCT_ECAL_tower_phi,
        RCT_ECAL_tower_timing,
        RCT_ECAL_tower_spike
      );
      ECALUnclusteredSLR2_tower_et->push_back(*RCT_ECAL_tower_et);
      ECALUnclusteredSLR2_tower_eta->push_back(*RCT_ECAL_tower_eta);
      ECALUnclusteredSLR2_tower_phi->push_back(*RCT_ECAL_tower_phi);
      ECALUnclusteredSLR2_tower_timing->push_back(*RCT_ECAL_tower_timing);
      ECALUnclusteredSLR2_tower_spike->push_back(*RCT_ECAL_tower_spike);
    }
  }
  // Read out SLR1 ECAL unclustered energy
  if(evt.getByToken(ECALUnclusteredSLR1Src_, ECALUnclusteredSLR1_link)){
    for(const auto & link : *ECALUnclusteredSLR1_link){

      getECALUnclusteredEnergy(
        link.data(),
        1,
        RCT_ECAL_tower_et,
        RCT_ECAL_tower_eta,
        RCT_ECAL_tower_phi,
        RCT_ECAL_tower_timing,
        RCT_ECAL_tower_spike
      );
      ECALUnclusteredSLR1_tower_et->push_back(*RCT_ECAL_tower_et);
      ECALUnclusteredSLR1_tower_eta->push_back(*RCT_ECAL_tower_eta);
      ECALUnclusteredSLR1_tower_phi->push_back(*RCT_ECAL_tower_phi);
      ECALUnclusteredSLR1_tower_timing->push_back(*RCT_ECAL_tower_timing);
      ECALUnclusteredSLR1_tower_spike->push_back(*RCT_ECAL_tower_spike);
    }
  }
  // Read out SLR0 ECAL unclustered energy
  if(evt.getByToken(ECALUnclusteredSLR0Src_, ECALUnclusteredSLR0_link)){
    for(const auto & link : *ECALUnclusteredSLR0_link){

      getECALUnclusteredEnergy(
        link.data(),
        0,
        RCT_ECAL_tower_et,
        RCT_ECAL_tower_eta,
        RCT_ECAL_tower_phi,
        RCT_ECAL_tower_timing,
        RCT_ECAL_tower_spike
      );
      ECALUnclusteredSLR0_tower_et->push_back(*RCT_ECAL_tower_et);
      ECALUnclusteredSLR0_tower_eta->push_back(*RCT_ECAL_tower_eta);
      ECALUnclusteredSLR0_tower_phi->push_back(*RCT_ECAL_tower_phi);
      ECALUnclusteredSLR0_tower_timing->push_back(*RCT_ECAL_tower_timing);
      ECALUnclusteredSLR0_tower_spike->push_back(*RCT_ECAL_tower_spike);
    }
  }

  // Write out HCAL TPs
  int cc; // Keep track of which RCT card we are on for secondhalfstarts purposes
  bool secondhalfstarts;
  if(evt.getByToken(HCAL8Src_, HCAL8_link)){
    cc = 0;
    for(const auto & link : *HCAL8_link){

      secondhalfstarts = (((cc + 3) % 4) > 1); //True for cards 0,3,4,7,etc.

      getHCALTowers(
        link.data(),
        secondhalfstarts,
        8,
        RCT_HCAL_tower_et,
        RCT_HCAL_tower_eta,
        RCT_HCAL_tower_phi,
        RCT_HCAL_tower_fb
      );
      HCAL8_tower_et->push_back(*RCT_HCAL_tower_et);
      HCAL8_tower_eta->push_back(*RCT_HCAL_tower_eta);
      HCAL8_tower_phi->push_back(*RCT_HCAL_tower_phi);
      HCAL8_tower_fb->push_back(*RCT_HCAL_tower_fb);

      cc = cc + 1;
    }
  }
  if(evt.getByToken(HCAL7Src_, HCAL7_link)){
    cc = 0;
    for(const auto & link : *HCAL7_link){

      secondhalfstarts = (((cc + 3) % 4) > 1); //True for cards 0,3,4,7,etc.

      getHCALTowers(
        link.data(),
        secondhalfstarts,
        7,
        RCT_HCAL_tower_et,
        RCT_HCAL_tower_eta,
        RCT_HCAL_tower_phi,
        RCT_HCAL_tower_fb
      );
      HCAL7_tower_et->push_back(*RCT_HCAL_tower_et);
      HCAL7_tower_eta->push_back(*RCT_HCAL_tower_eta);
      HCAL7_tower_phi->push_back(*RCT_HCAL_tower_phi);
      HCAL7_tower_fb->push_back(*RCT_HCAL_tower_fb);

      cc = cc + 1;
    }
  }
  if(evt.getByToken(HCAL6Src_, HCAL6_link)){
    cc = 0;
    for(const auto & link : *HCAL6_link){

      secondhalfstarts = (((cc + 3) % 4) > 1); //True for cards 0,3,4,7,etc.

      getHCALTowers(
        link.data(),
        secondhalfstarts,
        6,
        RCT_HCAL_tower_et,
        RCT_HCAL_tower_eta,
        RCT_HCAL_tower_phi,
        RCT_HCAL_tower_fb
      );
      HCAL6_tower_et->push_back(*RCT_HCAL_tower_et);
      HCAL6_tower_eta->push_back(*RCT_HCAL_tower_eta);
      HCAL6_tower_phi->push_back(*RCT_HCAL_tower_phi);
      HCAL6_tower_fb->push_back(*RCT_HCAL_tower_fb);

      cc = cc + 1;
    }
  }
  if(evt.getByToken(HCAL5Src_, HCAL5_link)){
    cc = 0;
    for(const auto & link : *HCAL5_link){

      secondhalfstarts = (((cc + 3) % 4) > 1); //True for cards 0,3,4,7,etc.

      getHCALTowers(
        link.data(),
        secondhalfstarts,
        5,
        RCT_HCAL_tower_et,
        RCT_HCAL_tower_eta,
        RCT_HCAL_tower_phi,
        RCT_HCAL_tower_fb
      );
      HCAL5_tower_et->push_back(*RCT_HCAL_tower_et);
      HCAL5_tower_eta->push_back(*RCT_HCAL_tower_eta);
      HCAL5_tower_phi->push_back(*RCT_HCAL_tower_phi);
      HCAL5_tower_fb->push_back(*RCT_HCAL_tower_fb);

      cc = cc + 1;
    }
  }

  linkTree->Fill();
 
 }


void L1TMLonRCTSetupAnalyzer::endJob() {
}

L1TMLonRCTSetupAnalyzer::~L1TMLonRCTSetupAnalyzer(){
}

//////////////////////////// Utility functions ///////////////////////////////
void getEGammaClusters(
  ap_uint<576> Data,
  int SLR,
  std::vector<int>* RCT_seed_energy,
  std::vector<int>* RCT_energy,
  std::vector<int>* RCT_eta,
  std::vector<int>* RCT_phi,
  std::vector<int>* RCT_et5x5,
  std::vector<int>* RCT_et2x5,
  std::vector<int>* RCT_timing,
  std::vector<int>* RCT_spike,
  std::vector<int>* RCT_satur,
  std::vector<int>* RCT_brems,
  std::vector<int>* RCT_spare
) {

  // If cluster is from SLR other than 0, shift iEta up appropriately
  int maybe_iEta_offset = 75 - SLR*25;
  int iEta_offset = std::max(maybe_iEta_offset,0);
  
  RCT_seed_energy->clear();
  RCT_energy->clear();
  RCT_eta->clear();
  RCT_phi->clear();
  RCT_et5x5->clear();
  RCT_et2x5->clear();
  RCT_timing->clear();
  RCT_spike->clear();
  RCT_satur->clear();
  RCT_brems->clear();
  RCT_spare->clear();

  for(int i=0; i<9; i++){
    int start = i*64;
    
    int this_seed_energy = (int)Data.range(start+9,start);
    RCT_seed_energy->push_back(this_seed_energy);
    int this_energy = (int)Data.range(start+21,start+10);
    RCT_energy->push_back(this_energy);
    int this_eta = (int)Data.range(start+26,start+22) + iEta_offset;
    RCT_eta->push_back(this_eta);
    int this_phi = (int)Data.range(start+31,start+27);
    RCT_phi->push_back(this_phi);
    int this_et5x5 = (int)Data.range(start+41,start+32);
    RCT_et5x5->push_back(this_et5x5);
    int this_wps = (int)Data.range(start+51,start+42);
    RCT_et2x5->push_back(this_wps);
    int this_timing = (int)Data.range(start+56,start+52);
    RCT_timing->push_back(this_timing);
    int this_spike = (int)Data.range(start+57,start+57);
    RCT_spare->push_back(this_spike);
    int this_satur = (int)Data.range(start+58,start+58);
    RCT_spare->push_back(this_satur);
    int this_brems = (int)Data.range(start+60,start+59);
    RCT_brems->push_back(this_brems);
    int this_spare = (int)Data.range(start+63,start+61);
    RCT_spare->push_back(this_spare);
  }

}

void getECALUnclusteredEnergy(
  ap_uint<576> Data,
  int SLR,
  std::vector<int>* RCT_et,
  std::vector<int>* RCT_eta,
  std::vector<int>* RCT_phi,
  std::vector<int>* RCT_timing,
  std::vector<int>* RCT_spike
) {

  // If tower is from SLR other than 0, shift iEta up appropriately
  int maybe_iEta_offset = SLR*5 - 3;
  int iEta_offset = std::max(maybe_iEta_offset,0);

  // If SLR is not 0, then nTowers is 30. Otherwise, it is 12
  int nTowers;
  if(SLR > 0){
    nTowers = 30;
  }
  else {
    nTowers = 12;
  }

  int this_et;
  int this_eta;
  int this_phi;
  int this_timing;
  int this_spike;

  RCT_et->clear();
  RCT_eta->clear();
  RCT_phi->clear();
  RCT_timing->clear();
  RCT_spike->clear();

  for(int i=0; i<nTowers; i++) {
    int start = i*18;

    this_et = (int)Data.range(start+11,start);
    RCT_et->push_back(this_et);
    this_eta = i/6 + iEta_offset;
    RCT_eta->push_back(this_eta);
    this_phi = i%6;
    RCT_phi->push_back(this_phi);
    this_timing = (int)Data.range(start+16,start+12);
    RCT_timing->push_back(this_timing);
    this_spike = (int)Data.range(start+17,start+17);
    RCT_spike->push_back(this_spike);
  }

}

void getHCALTowers(
  ap_uint<576> Data,
  bool secondhalfstarts,
  int nLink,
  std::vector<int>* RCT_et,
  std::vector<int>* RCT_eta,
  std::vector<int>* RCT_phi,
  std::vector<int>* RCT_fb
) {

  int this_et;
  int this_eta;
  int this_phi;
  int this_fb;

  RCT_et->clear();
  RCT_eta->clear();
  RCT_phi->clear();
  RCT_fb->clear();

  for(int i=0; i<32; i++) {
    // Lower iPhi in this link
    int start = i*16;

    this_et = (int)Data.range(start+9,start);
    RCT_et->push_back(this_et);
    this_fb = (int)Data.range(start+15,start+10);
    RCT_fb->push_back(this_fb);

    this_eta = ((nLink-5)%2)*8 + (i/4);
    if(!secondhalfstarts){
      this_phi = ((nLink-5)/2)*4 + (i%4);
    }
    else {
      this_phi = 2 - ((nLink-5)/2)*4 + (i%4);
    }
    RCT_eta->push_back(this_eta);
    RCT_phi->push_back(this_phi);
  }

}

DEFINE_FWK_MODULE(L1TMLonRCTSetupAnalyzer);
