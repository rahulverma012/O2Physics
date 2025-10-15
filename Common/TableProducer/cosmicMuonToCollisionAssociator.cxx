// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief  Cosmic Muon tracking for Alice Data
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Marian I Ivanov (mivanov@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Core/HfHelper.h"

// #include "DerivedPIDTables.h"
#include "Common/DataModel/OccupancyTables.h"

#include <typeinfo>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <chrono>

#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"
#include "TPDGCode.h"

#include "TableHelper.h"
#include "Common/Tools/TrackTuner.h"

namespace o2
{
namespace aod{
using Tracks = aod::Tracks;

namespace csmu
{
DECLARE_SOA_INDEX_COLUMN_FULL(UTrack, uTrack, int64_t, Tracks, "_UTrk"); //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(LTrack, lTrack, int64_t, Tracks, "_LTrk"); //! Negative track

DECLARE_SOA_COLUMN(UTrkGlobalBcStart, uTrkGlobalBcStart, int64_t);
DECLARE_SOA_COLUMN(UTrkGlobalBcEnd, uTrkGlobalBcEnd, int64_t);

DECLARE_SOA_COLUMN(LTrkGlobalBcStart, lTrkGlobalBcStart, int64_t);
DECLARE_SOA_COLUMN(LTrkGlobalBcEnd, lTrkGlobalBcEnd, int64_t);

DECLARE_SOA_COLUMN(UTimestamp, uTimestamp, int64_t);
DECLARE_SOA_COLUMN(LTimestamp, lTimestamp, int64_t);

DECLARE_SOA_COLUMN(UTfId, uTfId, int64_t);
DECLARE_SOA_COLUMN(LTfId, lTfId, int64_t);

DECLARE_SOA_COLUMN(UBCinTF, uBCinTF, int64_t);
DECLARE_SOA_COLUMN(LBCinTF, lBCinTF, int64_t);

DECLARE_SOA_COLUMN(UTrkX, uTrkX, float);
DECLARE_SOA_COLUMN(UTrkAlpha, uTrkAlpha, float);
DECLARE_SOA_COLUMN(UTrkY, uTrkY, float);
DECLARE_SOA_COLUMN(UTrkZ, uTrkZ, float);
DECLARE_SOA_COLUMN(UTrkSnp, uTrkSnp, float);
DECLARE_SOA_COLUMN(UTrkTgl, uTrkTgl, float);
DECLARE_SOA_COLUMN(UTrkSigned1Pt, uTrkSigned1Pt, float);
DECLARE_SOA_COLUMN(UTrkPt, uTrkPt, float);
DECLARE_SOA_COLUMN(UTrkP, uTrkP, float);
DECLARE_SOA_COLUMN(UTrkEta, uTrkEta, float);
DECLARE_SOA_COLUMN(UTrkPhi, uTrkPhi, float);

DECLARE_SOA_COLUMN(UTrkTpcNClsFindable, uTrkTpcNClsFindable, uint8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsFindableMinusFound, uTrkTpcNClsFindableMinusFound, int8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsFindableMinusCrossedRows, uTrkTpcNClsFindableMinusCrossedRows, int8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsShared, uTrkTpcNClsShared, uint8_t);

DECLARE_SOA_COLUMN(UTrkTime, uTrkTime, float);
DECLARE_SOA_COLUMN(UTrkTimeRes, uTrkTimeRes, float);

DECLARE_SOA_COLUMN(UTrkTpcNClsFound, uTrkTPCNClsFound, uint8_t);
DECLARE_SOA_COLUMN(UTrkTpcNClsCrossedRows, uTrkTPCNClsCrossedRows, uint8_t);

DECLARE_SOA_COLUMN(UTrkTpcCrossedRowsOverFindableCls, uTrkTPCCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(UTrkTpcFoundOverFindableCls, uTrkTPCFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(UTrkTpcFractionSharedCls, uTrkTPCFractionSharedCls, float);

DECLARE_SOA_COLUMN(UTrkDcaCalc, uTrkDcaCalc, float);
DECLARE_SOA_COLUMN(UTrkDcaXY, uTrkDcaXY, float);
DECLARE_SOA_COLUMN(UTrkDcaZ, uTrkDcaZ, float);
DECLARE_SOA_COLUMN(UTrkSigmaDcaXY2, uTrkSigmaDcaXY2, float);
DECLARE_SOA_COLUMN(UTrkSigmaDcaZ2, uTrkSigmaDcaZ2, float);

DECLARE_SOA_COLUMN(UTrkTpcSignal, uTrkTpcSignal, float);
DECLARE_SOA_COLUMN(UTrkTrdSignal, uTrkTrdSignal, float);

DECLARE_SOA_COLUMN(LTrkX, lTrkX, float);
DECLARE_SOA_COLUMN(LTrkAlpha, lTrkAlpha, float);
DECLARE_SOA_COLUMN(LTrkY, lTrkY, float);
DECLARE_SOA_COLUMN(LTrkZ, lTrkZ, float);
DECLARE_SOA_COLUMN(LTrkSnp, lTrkSnp, float);
DECLARE_SOA_COLUMN(LTrkTgl, lTrkTgl, float);
DECLARE_SOA_COLUMN(LTrkSigned1Pt, lTrkSigned1Pt, float);
DECLARE_SOA_COLUMN(LTrkPt, lTrkPt, float);
DECLARE_SOA_COLUMN(LTrkP, lTrkP, float);
DECLARE_SOA_COLUMN(LTrkEta, lTrkEta, float);
DECLARE_SOA_COLUMN(LTrkPhi, lTrkPhi, float);

DECLARE_SOA_COLUMN(LTrkTpcNClsFindable, lTrkTpcNClsFindable, uint8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsFindableMinusFound, lTrkTpcNClsFindableMinusFound, int8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsFindableMinusCrossedRows, lTrkTpcNClsFindableMinusCrossedRows, int8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsShared, lTrkTpcNClsShared, uint8_t);

DECLARE_SOA_COLUMN(LTrkTime, lTrkTime, float);
DECLARE_SOA_COLUMN(LTrkTimeRes, lTrkTimeRes, float);

DECLARE_SOA_COLUMN(LTrkTpcNClsFound, lTrkTPCNClsFound, uint8_t);
DECLARE_SOA_COLUMN(LTrkTpcNClsCrossedRows, lTrkTPCNClsCrossedRows, uint8_t);

DECLARE_SOA_COLUMN(LTrkTpcCrossedRowsOverFindableCls, lTrkTPCCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(LTrkTpcFoundOverFindableCls, lTrkTPCFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(LTrkTpcFractionSharedCls, lTrkTPCFractionSharedCls, float);

DECLARE_SOA_COLUMN(LTrkDcaCalc, lTrkDcaCalc, float);
DECLARE_SOA_COLUMN(LTrkDcaXY, lTrkDcaXY, float);
DECLARE_SOA_COLUMN(LTrkDcaZ, lTrkDcaZ, float);
DECLARE_SOA_COLUMN(LTrkSigmaDcaXY2, lTrkSigmaDcaXY2, float);
DECLARE_SOA_COLUMN(LTrkSigmaDcaZ2, lTrkSigmaDcaZ2, float);

DECLARE_SOA_COLUMN(LTrkTpcSignal, lTrkTpcSignal, float);
DECLARE_SOA_COLUMN(LTrkTrdSignal, lTrkTrdSignal, float);

DECLARE_SOA_COLUMN(B_kG, b_kG, int8_t);

DECLARE_SOA_COLUMN(PtCut15, ptCut15, bool);
DECLARE_SOA_COLUMN(IsOKSum6, isOKSum6, bool);
DECLARE_SOA_COLUMN(IsOK6, isOK6, bool);
DECLARE_SOA_COLUMN(TOK0, tOK0, bool);
DECLARE_SOA_COLUMN(IsSelected, isSelected, bool);


DECLARE_SOA_COLUMN(HasCosmicPair, hasCosmicPair, bool);
DECLARE_SOA_COLUMN(HasIsSelectedCosmic, hasIsSelectedCosmic, bool);
DECLARE_SOA_COLUMN(CosmicPairIndexList, cosmicPairIndexList, std::vector<int>);

}; // namespace csmu

DECLARE_SOA_TABLE(CosmicPairs, "AOD", "COSMICPAIRS", o2::soa::Index<>,
o2::aod::csmu::UTrackId,
o2::aod::csmu::LTrackId,

o2::aod::csmu::UTrkGlobalBcStart,
o2::aod::csmu::UTrkGlobalBcEnd,
// o2::aod::csmu::UTimestamp,
o2::aod::csmu::UTfId,
o2::aod::csmu::UBCinTF,

o2::aod::csmu::LTrkGlobalBcStart,
o2::aod::csmu::LTrkGlobalBcEnd,
// o2::aod::csmu::LTimestamp,
o2::aod::csmu::LTfId,
o2::aod::csmu::LBCinTF,

o2::aod::csmu::B_kG,

o2::aod::csmu::UTrkX,
o2::aod::csmu::UTrkAlpha,
o2::aod::csmu::UTrkY,
o2::aod::csmu::UTrkZ,
o2::aod::csmu::UTrkSnp,
o2::aod::csmu::UTrkTgl,
o2::aod::csmu::UTrkSigned1Pt,
o2::aod::csmu::UTrkPt,
o2::aod::csmu::UTrkP,
o2::aod::csmu::UTrkEta,
o2::aod::csmu::UTrkPhi,

o2::aod::csmu::UTrkTpcNClsFindable,
o2::aod::csmu::UTrkTpcNClsFindableMinusFound,
o2::aod::csmu::UTrkTpcNClsFindableMinusCrossedRows,
o2::aod::csmu::UTrkTpcNClsShared,
o2::aod::csmu::UTrkTpcNClsFound,
o2::aod::csmu::UTrkTpcNClsCrossedRows,
o2::aod::csmu::UTrkTpcCrossedRowsOverFindableCls,
o2::aod::csmu::UTrkTpcFoundOverFindableCls,
o2::aod::csmu::UTrkTpcFractionSharedCls,
o2::aod::csmu::UTrkDcaCalc,
o2::aod::csmu::UTrkDcaXY,
o2::aod::csmu::UTrkDcaZ,
o2::aod::csmu::UTrkSigmaDcaXY2,
o2::aod::csmu::UTrkSigmaDcaZ2,

o2::aod::csmu::UTrkTpcSignal,                       // dE/dx signal in the TPC
o2::aod::csmu::UTrkTrdSignal,                       // PID signal in the TRD

o2::aod::csmu::LTrkX,
o2::aod::csmu::LTrkAlpha,
o2::aod::csmu::LTrkY,
o2::aod::csmu::LTrkZ,
o2::aod::csmu::LTrkSnp,
o2::aod::csmu::LTrkTgl,
o2::aod::csmu::LTrkSigned1Pt,
o2::aod::csmu::LTrkPt,
o2::aod::csmu::LTrkP,
o2::aod::csmu::LTrkEta,
o2::aod::csmu::LTrkPhi,

o2::aod::csmu::LTrkTpcNClsFindable,
o2::aod::csmu::LTrkTpcNClsFindableMinusFound,
o2::aod::csmu::LTrkTpcNClsFindableMinusCrossedRows,
o2::aod::csmu::LTrkTpcNClsShared,
o2::aod::csmu::LTrkTpcNClsFound,
o2::aod::csmu::LTrkTpcNClsCrossedRows,
o2::aod::csmu::LTrkTpcCrossedRowsOverFindableCls,
o2::aod::csmu::LTrkTpcFoundOverFindableCls,
o2::aod::csmu::LTrkTpcFractionSharedCls,
o2::aod::csmu::LTrkDcaCalc,
o2::aod::csmu::LTrkDcaXY,
o2::aod::csmu::LTrkDcaZ,
o2::aod::csmu::LTrkSigmaDcaXY2,
o2::aod::csmu::LTrkSigmaDcaZ2,

o2::aod::csmu::LTrkTpcSignal,                       // dE/dx signal in the TPC
o2::aod::csmu::LTrkTrdSignal,                       // PID signal in the TRD

o2::aod::csmu::PtCut15,
o2::aod::csmu::IsOKSum6,
o2::aod::csmu::IsOK6,
o2::aod::csmu::TOK0,
o2::aod::csmu::IsSelected
);

DECLARE_SOA_TABLE(CollisionCosmicFlags, "AOD", "COLLCOSMICFLAGS", o2::soa::Index<>,
o2::aod::csmu::HasCosmicPair, o2::aod::csmu::HasIsSelectedCosmic, o2::aod::csmu::CosmicPairIndexList);

} //namespace aod
} // namespace o2

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

template<typename T>
void PrintTime(T Start, std::string String){
  auto Stop = std::chrono::high_resolution_clock::now();
  auto Duration = duration_cast<std::chrono::microseconds>(Stop-Start);
  LOG(info)<<String<<float(Duration.count())/float(1000000)<<" seconds";//<<endl;
}

int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

struct cosmicMuonToCollisionAssociator{

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  //declare production of tables
  Produces<aod::CosmicPairs> genCosmicPairs;
  Produces<aod::CollisionCosmicFlags> genCollisionCosmicFlags;

  //Histogram registry;
  HistogramRegistry histReg   {"histReg"   , {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  //Configurables
  
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<int> nBCinTF{"nBCinTF", 114048, "nBCinTF"};

  struct : ConfigurableGroup {
  Configurable<float> csmMuThrPt        {"csmMuThrPt"        , 0.5, "csmMuThrPt"};
  Configurable<float> csmMuThrDCAMin    {"csmMuThrDCAMin"    , 0.5, "csmMuThrDCAMin"};
  Configurable<float> csmMuThrDCAMax    {"csmMuThrDCAMax"    , 120, "csmMuThrDCAMax"};

  Configurable<float> csmMuSumPtPair    {"csmMuSumPtPair"    , 2.0 , "csmMuSumPtPair"};
  Configurable<float> csmMuSumQPtPair   {"csmMuSumQPtPair"   , 0.2 , "csmMuSumQPtPair"};
  Configurable<float> csmMuSumTglPair   {"csmMuSumTglPair"   , 0.01, "csmMuSumTglPair"};
  Configurable<float> csmMuSumDcaXY     {"csmMuSumDcaXY"     , 6.0 , "csmMuSumDcaXY"};
  Configurable<float> csmMuDiffDcaXY    {"csmMuDiffDcaXY"    , 0.6 , "csmMuDiffDcaXY"};
  Configurable<float> csmMuDiffAlphaPair{"csmMuDiffAlphaPair", 0.1, "csmMuDiffAlphaPair"};

  Configurable<bool> csmMuCheckSumDcaXY     {"csmMuCheckSumDcaXY"     , true, "csmMuCheckSumDcaXY"};
  Configurable<bool> csmMuCheckDiffDcaXY    {"csmMuCheckDiffDcaXY"    , true, "csmMuCheckDiffDcaXY"};
  Configurable<bool> csmMuCheckDiffAlphaPair{"csmMuCheckDiffAlphaPair", true, "csmMuCheckDiffAlphaPair"};
  } cfgCM;

  Configurable<std::vector<double>> countBins{
    "countBins",
    {
      // 0 to 10 (step 1)
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
      // 10 to 100 (step 10)
      20, 30, 40, 50, 60, 70, 80, 90, 100,
      // 100 to 1,000 (step 100)
      200, 300, 400, 500, 600, 700, 800, 900, 1000,
      // 1,000 to 10,000 (step 1,000)
      2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
      // 10,000 to 100,000 (step 10,000)
      20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000,
      // 100,000 to 1,000,000 (step 100,000)
      200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000,
      // 1,000,000 to 10,000,000 (step 1,000,000)
      2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000,
      // 10,000,000 to 100,000,000 (step 10,000,000)
      20000000, 30000000, 40000000, 50000000, 60000000, 70000000, 80000000, 90000000, 100000000,
      // 100,000,000 to 1,000,000,000 (step 100,000,000)
      200000000, 300000000, 400000000, 500000000, 600000000, 700000000, 800000000, 900000000, 1000000000
      },
      "Particle Count Bins"
  };


  enum CosmicPairRejectionType{
    kPairPassed = 0,
    kFailSumPtPair,    
    kFailSumQPtPair,   
    kFailSumTglPair,   
    kFailSumDcaXY,     
    kFailDiffDcaXY,    
    kFailDiffAlphaPair
  };

  static constexpr std::string_view CMRejectionTag[]{
    "kPairPassed",
    "kFailSumPtPair",    
    "kFailSumQPtPair",   
    "kFailSumTglPair",   
    "kFailSumDcaXY",     
    "kFailDiffDcaXY",    
    "kFailDiffAlphaPair"
  };
  
  void init(InitContext const&){
    const AxisSpec axisCount{countBins, "Counts"};

    histReg.add("hTiming_nDFsProcessed", "hTiming_nDFsProcessed", kTH1F, {{1, 0, 1}});
    histReg.add("hTiming_nTFsPerDataFrame", "hTiming_nTFsPerDataFrame", kTH1F, {{5000, 0, 5000}});
    histReg.add("hTiming_nTFsProcessed", "hTiming_nTFsProcessed", kTH1F, {{1, 0, 1}});

    histReg.add("hCosmicPairs_nCosmicsPerDF", "hCosmicPairs_nCosmicsPerDF", kTH1F, {{5000, 0, 5000}});
    histReg.add("hCosmicPairs_nIsSelCosmicsPerDF", "hCosmicPairs_nIsSelCosmicsPerDF", kTH1F, {{5000, 0, 5000}});

    histReg.add("hRatio",           "Ratio;Ratio;Entries"       , kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioVtx8",       "Ratio Vtx<8;Ratio;Entries" , kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioVtx10",      "Ratio Vtx<10;Ratio;Entries", kTH1F, {{1000, 0, 1}});

    histReg.add("hRatioIsSel",      "Ratio IsSelected;Ratio;Entries",        kTH1F,  {{1000, 0, 1}});
    histReg.add("hRatioIsSelVtx8",  "Ratio IsSelected Vtx<8;Ratio;Entries",  kTH1F, {{1000, 0, 1}});
    histReg.add("hRatioIsSelVtx10", "Ratio IsSelected Vtx<10;Ratio;Entries", kTH1F,  {{1000, 0, 1}});

    histReg.add("hCosmicFlag",      "Cosmic Flag;Count;Entries",        kTH1F, {axisCount});
    histReg.add("hCosmicFlagVtx8",  "Cosmic Flag Vtx<8;Count;Entries",  kTH1F, {axisCount});
    histReg.add("hCosmicFlagVtx10", "Cosmic Flag Vtx<10;Count;Entries", kTH1F, {axisCount});

    histReg.add("hCosmicFlagIsSel",      "IsSelected Cosmic Flag;Count;Entries",        kTH1F, {axisCount});
    histReg.add("hCosmicFlagIsSelVtx8",  "IsSelected Cosmic Flag Vtx<8;Count;Entries",  kTH1F, {axisCount});
    histReg.add("hCosmicFlagIsSelVtx10", "IsSelected Cosmic Flag Vtx<10;Count;Entries", kTH1F, {axisCount});

    histReg.add("hColls",      "Total Colls;Count;Entries",  kTH1F, {axisCount});
    histReg.add("hCollsVtx8",  "Colls Vtx<8;Count;Entries",  kTH1F, {axisCount});
    histReg.add("hCollsVtx10", "Colls Vtx<10;Count;Entries", kTH1F,  {axisCount});
  }

  void GetRunInfo(const int& run, int& nBCsPerTF, int64_t& bcSOR)
  {
    auto runDuration = ccdb->getRunDuration(run, true);
    int64_t tsSOR = runDuration.first;
    auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
    int64_t tsOrbitReset = (*ctpx)[0];
    uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
    int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
    orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
    bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit; // customOrbitOffset is a Configurable
    nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
  }

  template <typename T>
  void GetTimingInfo(const T& bc, int& lastRun, int32_t& nBCsPerTF, int64_t& bcSOR, uint64_t& time, int64_t& tfIdThis, int64_t& bcInTF)
  {
    int run = bc.runNumber();
    if (run != lastRun) { // update run info
      lastRun = run;
      GetRunInfo(run, nBCsPerTF, bcSOR); // update nBCsPerTF && bcSOR
    }
    // update the information
    time = bc.timestamp();
    tfIdThis = (bc.globalBC() - bcSOR) / nBCsPerTF;
    bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
  }

  int runNumber = -1;
  o2::parameters::GRPMagField* grpmag = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  int8_t B_kG; 

  struct : ConfigurableGroup {
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } cfgB;

  template<typename T>
  void initCCDB(const T& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(cfgB.lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgB.grpmagPath, bc.timestamp());
    B_kG = grpmag->getNominalL3Field();
    runNumber = bc.runNumber();
  }

  static constexpr double kB2C = 0.299792458;
  static constexpr double kAlmost0Field = 1e-4;
  double getCurvature(double signed1Pt, double b){ return signed1Pt * b * kB2C; }
  double getDCA(double trackX, double trackY, double alpha, double snp, double signed1Pt, 
                  double xRef, double yRef, double b)
  {
      // Handle zero or near-zero magnetic field by simple linear distance
      if (std::abs(b) < kAlmost0Field) { return std::hypot(trackX - xRef, trackY - yRef);}

      double curvature = getCurvature(signed1Pt, b);

      // Rotate reference point into track frame (rotation by -alpha)
      double sinAlpha = std::sin(alpha);
      double cosAlpha = std::cos(alpha);

      double xRot = xRef * cosAlpha + yRef * sinAlpha;
      double yRot = -xRef * sinAlpha + yRef * cosAlpha;

      // Compute difference vector in local track frame
      double dx = trackX - xRot;
      double dy = trackY - yRot;

      double sqrtTerm = std::sqrt((1.0 - snp) * (1.0 + snp));

      double sn = curvature * dx - snp;
      double cs = curvature * dy + sqrtTerm;

      double numerator = 2.0 * (dx * snp - dy * sqrtTerm) - curvature * (dx * dx + dy * dy);
      double denominator = 1.0 + std::sqrt(sn * sn + cs * cs);

      return -numerator / denominator;
  }

  template<typename T>
  std::vector<int> getCosmicTriggerList(int64_t globalBC, const T& list){
    std::vector<int> triggeredCosmicIndexList;
    for(int i = 0; i < list.size(); i++){
      if(std::abs(globalBC - list[i]) < 6000) {
        triggeredCosmicIndexList.push_back(i); //it had a cosmic entry;
      }
    }
    return triggeredCosmicIndexList;
  }

  enum VariableType{
   kTrkGI, //track.globalIndex()
   kTrkFirstGBC, //first global bunch crossing of Track
   kTrkLastGBC, //second global bunch crossing of Track
   kTrkTfIdThis, //Timeframe if of the Track  
   kTrkBCinTF  //Bunch crossing of track in time frame
  };

  int dfCount = 0;
  int32_t nBCsPerTF = -999;
  int64_t bcSOR = -999;
  uint64_t time = -1;
  int64_t tfIdThis = -1;
  int64_t bcInTF = -1;
  int lastRun = -999;  
  std::chrono::high_resolution_clock::time_point Start0 = std::chrono::high_resolution_clock::now();

  template<typename B, typename C, typename T, typename A>
  void executeProcess(const B& BCs, const C& colls, const T& tracks, const A& ambgTracks,const auto& Origins ){
    dfCount++; auto Start1 = std::chrono::high_resolution_clock::now();
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: DF_"<<Origins.iteratorAt(0).dataframeID()
             <<" :: BCs.size() = "<<BCs.size()
             <<" :: colls.size() = "<<colls.size()
             <<" :: tracks.size() = "<<tracks.size()
             <<" :: ambgTracks.size() = "<<ambgTracks.size()
             ;

    if(colls.size() == 0) return;

    //Step 1: loop over tracks and get good tracks.
    float trackDcaCalc;
    int64_t lastCollId = -999;
    std::vector<int64_t> upperTrkGIlist;
    std::vector<int64_t> lowerTrkGIlist;

    std::vector<int64_t> trackRunNumber;
    std::vector<int64_t> trackGlobalBC;
    std::vector<int64_t> trackTriggerMask;
    std::vector<int64_t> trackTimestamp;
    std::vector<int64_t> trackTFidThis;
    std::vector<int64_t> trackBcInTF;
    std::vector<int64_t> fullTimeFrameIdList;

    std::vector<std::array<int64_t, 5>> upperTrkGIListWithGlobalBCInfo;
    std::vector<std::array<int64_t, 5>> lowerTrkGIListWithGlobalBCInfo;

    std::vector<float> upperTrkDCAInfo;
    std::vector<float> lowerTrkDCAInfo;

    initCCDB(BCs.begin());
    auto coll = colls.begin();
    auto bc = BCs.begin();
    int64_t firstGlobalBC = 0;
    int64_t lastGlobalBC =0; 
    auto trackLoopStart = std::chrono::high_resolution_clock::now();

    for(const auto& track: tracks){
      if(track.pt() < cfgCM.csmMuThrPt ) {continue;} 
      if(track.tpcNClsFound() < 100 ) {continue;} 
      trackDcaCalc = getDCA(track.x(), track.y(), track.alpha(), track.snp(), track.signed1Pt(), 0, 0, B_kG);
      if(std::abs(trackDcaCalc) < cfgCM.csmMuThrDCAMin ) {continue;}
      if(std::abs(trackDcaCalc) > cfgCM.csmMuThrDCAMax ) {continue;}

      //Get Timing info of the track
      if(track.collisionId()< 0 ) {
        lastCollId = -999;
        trackRunNumber.clear();
        trackGlobalBC.clear();
        trackTriggerMask.clear();        
        trackTimestamp.clear();
        trackTFidThis.clear();
        trackBcInTF.clear();
        auto bcs = track.template ambgTrack_as<A>().template bc_as<B>();
        if( bcs.size() == 0 ){
          trackRunNumber.push_back(0);
          trackGlobalBC.push_back(0);
          trackTriggerMask.push_back(0);
          trackTimestamp.push_back(0);
          trackTFidThis.push_back(-1);
          trackBcInTF.push_back(-1);
          firstGlobalBC = -10000; 
          lastGlobalBC = -10000;
        } else{
          for (const auto& bc: bcs){
            GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
            trackRunNumber.push_back(bc.runNumber());
            trackGlobalBC.push_back(bc.globalBC());
            trackTriggerMask.push_back(bc.triggerMask());
            trackTimestamp.push_back(time);
            trackTFidThis.push_back(tfIdThis);
            trackBcInTF.push_back(bcInTF);
          }
          if(bcs.size() == 1){
            firstGlobalBC = trackGlobalBC[0];
            lastGlobalBC  = -10000;
          } else {
            firstGlobalBC = trackGlobalBC[0];
            lastGlobalBC  = trackGlobalBC[bcs.size()-1];
          }
        }
      }
      else {
        coll = track.template collision_as<C>();
        if (lastCollId != coll.globalIndex()) {
          trackRunNumber.clear();
          trackGlobalBC.clear();
          trackTriggerMask.clear();
          trackTimestamp.clear();
          trackTFidThis.clear();
          trackBcInTF.clear();
          lastCollId = coll.globalIndex();
          bc = coll.template bc_as<B>();
          GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);
          trackRunNumber.push_back(bc.runNumber());
          trackGlobalBC.push_back(bc.globalBC());
          trackTriggerMask.push_back(bc.triggerMask());
          trackTimestamp.push_back(time);
          trackTFidThis.push_back(tfIdThis);
          trackBcInTF.push_back(bcInTF);
          firstGlobalBC = trackGlobalBC[0];
          lastGlobalBC  = -10000;
        }
      }

      fullTimeFrameIdList.push_back(trackTFidThis[0]);
      if(track.alpha() < 0 ) {
        lowerTrkGIListWithGlobalBCInfo.push_back({track.globalIndex(), firstGlobalBC, lastGlobalBC, trackTFidThis[0], trackBcInTF[0]}); lowerTrkDCAInfo.push_back(trackDcaCalc);}
      if(track.alpha() > 0 ) {
        upperTrkGIListWithGlobalBCInfo.push_back({track.globalIndex(), firstGlobalBC, lastGlobalBC, trackTFidThis[0], trackBcInTF[0]}); upperTrkDCAInfo.push_back(trackDcaCalc);//} //continue; Fill upper track
      }
    }

    std::sort(fullTimeFrameIdList.begin(), fullTimeFrameIdList.end());
    auto last = std::unique(fullTimeFrameIdList.begin(), fullTimeFrameIdList.end());
    fullTimeFrameIdList.erase(last, fullTimeFrameIdList.end());

    histReg.fill(HIST("hTiming_nDFsProcessed"), 0.5);
    histReg.fill(HIST("hTiming_nTFsProcessed"), 0.5, fullTimeFrameIdList.size());
    histReg.fill(HIST("hTiming_nTFsPerDataFrame"), fullTimeFrameIdList.size());

    LOG(info)<<"DEBUG :: upperTrkGIListWithGlobalBCInfo = "<<upperTrkGIListWithGlobalBCInfo.size();
    LOG(info)<<"DEBUG :: lowerTrkGIListWithGlobalBCInfo = "<<lowerTrkGIListWithGlobalBCInfo.size();
    PrintTime(trackLoopStart, Form("DEBUG :: df_%d :: DF Reading :: Track loop Time :: ",dfCount));

    auto cosmicMuonLoopStart = std::chrono::high_resolution_clock::now();
    std::vector<int64_t> cosmicMuUTrkBC0;
    std::vector<int64_t> cosmicMuUTrkBC1;
    std::vector<int64_t> cosmicMuLTrkBC0;
    std::vector<int64_t> cosmicMuLTrkBC1;
    std::vector<bool> cosmicMuIsSelectedFlag;

    int cosmicPairCounter = 0;
    float fSumPtPair = 0, fSumQPtPair = 0,  uTrkDcaCalc = 0, lTrkDcaCalc = 0;

    constexpr float sP0 = 1.f;
    constexpr float sP3 = 0.006f;
    constexpr float sP4 = 0.03f;
    constexpr float sAlpha = 0.01f;
    constexpr float sixSquared = 36.f;  // 6*6 //to avoid sqrt compare with square

    float dP0, dP3, dP4, dAlpha, ptMin;
    float resP0, resP3, resP4, resAlpha, sumSq;
    bool ptCut15, isOKSum6, isOK6, tOK0, isSelected;

    for(uint idxUTrk = 0 ; idxUTrk < upperTrkGIListWithGlobalBCInfo.size(); idxUTrk++){
      const auto& uTrk = upperTrkGIListWithGlobalBCInfo[idxUTrk];
      const auto& upperTrk = tracks.iteratorAt(uTrk[kTrkGI]);
 
      for(uint idxLTrk = 0 ; idxLTrk < lowerTrkGIListWithGlobalBCInfo.size(); idxLTrk++){
        const auto& lTrk = lowerTrkGIListWithGlobalBCInfo[idxLTrk];        
        const auto& lowerTrk = tracks.iteratorAt(lTrk[kTrkGI]);

        if(uTrk[kTrkTfIdThis] != lTrk[kTrkTfIdThis]) {continue;}

        fSumPtPair      = upperTrk.pt()        + lowerTrk.pt();
        fSumQPtPair     = upperTrk.signed1Pt() + lowerTrk.signed1Pt();

        uTrkDcaCalc = upperTrkDCAInfo[idxUTrk];
        lTrkDcaCalc = lowerTrkDCAInfo[idxLTrk];

        if ( std::abs(fSumPtPair) < cfgCM.csmMuSumPtPair                                     ) {continue;}
        if ( std::abs(fSumQPtPair) > cfgCM.csmMuSumQPtPair                                   ) {continue;}

        if( !(std::abs(std::abs(upperTrk.signed1Pt())-std::abs(lowerTrk.signed1Pt())) < 0.2 )) {continue;}
        if( !(std::abs(std::abs(upperTrk.tgl())      -std::abs(lowerTrk.tgl())) < 0.2 )      ) {continue;}
        if( !(std::abs(std::abs(uTrkDcaCalc)         -std::abs(lTrkDcaCalc)) < 5 )           ) {continue;}

        cosmicMuUTrkBC0.push_back(uTrk[kTrkFirstGBC]);
        cosmicMuUTrkBC1.push_back(uTrk[kTrkLastGBC ]);
        cosmicMuLTrkBC0.push_back(lTrk[kTrkFirstGBC]);
        cosmicMuLTrkBC1.push_back(lTrk[kTrkLastGBC ]);

        dP0 = lowerTrk.dcaXY() + upperTrk.dcaXY();
        dP3 = lowerTrk.tgl() + upperTrk.tgl();
        dP4 = lowerTrk.signed1Pt() + upperTrk.signed1Pt();
        dAlpha = upperTrk.alpha()-lowerTrk.alpha()-TMath::Pi();

        ptMin = (lowerTrk.pt() < upperTrk.pt()) ? lowerTrk.pt() : upperTrk.pt();
        ptCut15 = ptMin > 1.5f; //this is residual distortion dependent

        //Compute square of normalised residuals for fast calculation
        resP0 = dP0 / sP0;
        resP3 = dP3 / sP3;
        resP4 = dP4 / sP4;
        resAlpha = dAlpha / sAlpha;
        sumSq = resP0 * resP0 + resP3 * resP3 + resP4 * resP4 + resAlpha * resAlpha;
        isOKSum6 = sumSq < sixSquared;
        isOK6 = abs(resP0)<6 && abs(resP3)<6&&abs(resP4)<6&&abs(dAlpha/sAlpha)<6;
        tOK0 = abs(uTrk[kTrkBCinTF]-lTrk[kTrkBCinTF])<4000;
        isSelected = isOK6 && tOK0 && ptCut15;

        cosmicMuIsSelectedFlag.push_back(isSelected);
        cosmicPairCounter++;
        genCosmicPairs(
           upperTrk.globalIndex()
          ,lowerTrk.globalIndex()

          ,uTrk[kTrkFirstGBC]
          ,uTrk[kTrkLastGBC ]
          ,uTrk[kTrkTfIdThis]
          ,uTrk[kTrkBCinTF]

          ,lTrk[kTrkFirstGBC]
          ,lTrk[kTrkLastGBC ]
          ,lTrk[kTrkTfIdThis]
          ,lTrk[kTrkBCinTF]

          ,B_kG

          ,upperTrk.x()     	  
          ,upperTrk.alpha() 	  
          ,upperTrk.y() 	      
          ,upperTrk.z() 	      
          ,upperTrk.snp() 	    
          ,upperTrk.tgl() 	    
          ,upperTrk.signed1Pt() 
          ,upperTrk.pt() 	      
          ,upperTrk.p() 	      
          ,upperTrk.eta() 	    
          ,upperTrk.phi() 	    

          // // TracksExtra
          ,upperTrk.tpcNClsFindable() 	               
          ,upperTrk.tpcNClsFindableMinusFound() 	     
          ,upperTrk.tpcNClsFindableMinusCrossedRows() 
          ,upperTrk.tpcNClsShared() 	                 
                    
          ,upperTrk.tpcNClsFound() 	                 
          ,upperTrk.tpcNClsCrossedRows() 	           
          
          ,upperTrk.tpcCrossedRowsOverFindableCls() 	
          ,upperTrk.tpcFoundOverFindableCls()         
          ,upperTrk.tpcFractionSharedCls() 	          

          ,uTrkDcaCalc
          ,upperTrk.dcaXY()  
          ,upperTrk.dcaZ()	  
          ,upperTrk.sigmaDcaXY2()
          ,upperTrk.sigmaDcaZ2()  

          ,upperTrk.tpcSignal()
          ,upperTrk.trdSignal()  

          ////////////////Lower Track Variables ////////////////////////////////////////////////////////////////////////////////////////////////
          ,lowerTrk.x() 	        
          ,lowerTrk.alpha() 	    
          ,lowerTrk.y() 	        
          ,lowerTrk.z() 	        
          ,lowerTrk.snp() 	      
          ,lowerTrk.tgl() 	      
          ,lowerTrk.signed1Pt() 	
          ,lowerTrk.pt() 	        
          ,lowerTrk.p() 	        
          ,lowerTrk.eta() 	      
          ,lowerTrk.phi() 	      

          ,lowerTrk.tpcNClsFindable() 	               
          ,lowerTrk.tpcNClsFindableMinusFound() 	     
          ,lowerTrk.tpcNClsFindableMinusCrossedRows() 
          ,lowerTrk.tpcNClsShared() 	                 
                    
          ,lowerTrk.tpcNClsFound() 	                
          ,lowerTrk.tpcNClsCrossedRows() 	          
          
          ,lowerTrk.tpcCrossedRowsOverFindableCls() 	 
          ,lowerTrk.tpcFoundOverFindableCls()          
          ,lowerTrk.tpcFractionSharedCls() 	           

          ,lTrkDcaCalc
          ,lowerTrk.dcaXY()   
          ,lowerTrk.dcaZ()	   
          ,lowerTrk.sigmaDcaXY2() 
          ,lowerTrk.sigmaDcaZ2()
          ,lowerTrk.tpcSignal()
          ,lowerTrk.trdSignal()
          ,ptCut15
          ,isOKSum6
          ,isOK6
          ,tOK0
          ,isSelected
        );
      }
    }
    
    int isSelCosmicPairCounter = std::count(cosmicMuIsSelectedFlag.begin(), cosmicMuIsSelectedFlag.end(), true);

    histReg.fill(HIST("hCosmicPairs_nCosmicsPerDF"), cosmicPairCounter);
    histReg.fill(HIST("hCosmicPairs_nIsSelCosmicsPerDF"), isSelCosmicPairCounter);

    PrintTime(cosmicMuonLoopStart, Form("DEBUG :: df_%d :: DF Reading :: cosmicMuonLoop Time :: ", dfCount));    
    LOG(info)<<"DEBUG :: df_"<<dfCount<<" :: cosmic counted = "<<cosmicPairCounter;

    // Collision looping and flag table creator
    auto collLoopStart = std::chrono::high_resolution_clock::now();
    std::unordered_set<int> finalCosmicIndexSet; //to keep only unique entries.
    std::vector<int> tempCosmicIndexList;
    std::vector<int> finalCosmicIndexList;
    bool flagHasCosmicMuon = false;
    bool flagHasIsSelectedCosmic = false;

    int counterCollVtx10 = 0;
    int counterCollVtx8 = 0;
    int counterFlagHasCosmic = 0;
    int counterFlagHasIsSelectedCosmic = 0 ;
    int counterVtx10FlagHasCosmic = 0;
    int counterVtx10FlagHasIsSelectedCosmic = 0 ;
    int counterVtx8FlagHasCosmic = 0;
    int counterVtx8FlagHasIsSelectedCosmic = 0 ;

    auto addIndices = [](const std::vector<int>& indices, std::unordered_set<int>& indexSet, bool& flagHasCosmic) {
      if (!indices.empty()) {
        flagHasCosmic = true;
        for (int idx : indices) indexSet.insert(idx);
      }
    };

    for(const auto& coll : colls){
      const auto& bc = coll.template bc_as<B>();
      GetTimingInfo(bc, lastRun, nBCsPerTF, bcSOR, time, tfIdThis, bcInTF);

      if (nBCsPerTF > nBCinTF) {
        LOG(fatal) << "DEBUG :: FATAL ERROR :: nBCsPerTF > nBCinTF i.e " << nBCsPerTF << " > " << nBCinTF << " will cause crash in further process";
        return;
      }
    
      flagHasCosmicMuon = false;
      finalCosmicIndexSet.clear();

      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuUTrkBC0), finalCosmicIndexSet, flagHasCosmicMuon);
      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuUTrkBC1), finalCosmicIndexSet, flagHasCosmicMuon);
      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuLTrkBC0), finalCosmicIndexSet, flagHasCosmicMuon);
      addIndices(getCosmicTriggerList(bc.globalBC(), cosmicMuLTrkBC1), finalCosmicIndexSet, flagHasCosmicMuon);

      finalCosmicIndexList.clear();
      finalCosmicIndexList.assign(finalCosmicIndexSet.begin(), finalCosmicIndexSet.end());

      flagHasIsSelectedCosmic = false;
      for (int idx : finalCosmicIndexSet) {
        if (cosmicMuIsSelectedFlag[idx]) {
          flagHasIsSelectedCosmic = true;
          break;
        }
      }

      const float absZ = std::abs(coll.posZ());
      const bool z8  = absZ < 8.0f;
      const bool z10 = absZ < 10.0f;

      counterCollVtx8  += z8;
      counterCollVtx10 += z10;

      if (flagHasCosmicMuon) {
        counterFlagHasCosmic++;
        counterVtx8FlagHasCosmic  += z8;
        counterVtx10FlagHasCosmic += z10;
      }

      if (flagHasIsSelectedCosmic) {
        counterFlagHasIsSelectedCosmic++;
        counterVtx8FlagHasIsSelectedCosmic  += z8;
        counterVtx10FlagHasIsSelectedCosmic += z10;
      }

      genCollisionCosmicFlags(flagHasCosmicMuon, flagHasIsSelectedCosmic, finalCosmicIndexList);
    }

    histReg.fill(HIST("hRatio"), static_cast<double>(counterFlagHasCosmic) / static_cast<double>(colls.size()));
    histReg.fill(HIST("hRatioVtx8"), static_cast<double>(counterVtx8FlagHasCosmic) / static_cast<double>(counterCollVtx8));
    histReg.fill(HIST("hRatioVtx10"), static_cast<double>(counterVtx10FlagHasCosmic) / static_cast<double>(counterCollVtx10));

    histReg.fill(HIST("hRatioIsSel"), static_cast<double>(counterFlagHasIsSelectedCosmic) /  static_cast<double>(colls.size()));
    histReg.fill(HIST("hRatioIsSelVtx8"), static_cast<double>(counterVtx8FlagHasIsSelectedCosmic) / static_cast<double>(counterCollVtx8));
    histReg.fill(HIST("hRatioIsSelVtx10"), static_cast<double>(counterVtx10FlagHasIsSelectedCosmic) / static_cast<double>(counterCollVtx10));

    histReg.fill(HIST("hCosmicFlag"), static_cast<double>(counterFlagHasCosmic));
    histReg.fill(HIST("hCosmicFlagVtx8"), static_cast<double>(counterVtx8FlagHasCosmic));
    histReg.fill(HIST("hCosmicFlagVtx10"), static_cast<double>(counterVtx10FlagHasCosmic));

    histReg.fill(HIST("hCosmicFlagIsSel"), static_cast<double>(counterFlagHasIsSelectedCosmic));
    histReg.fill(HIST("hCosmicFlagIsSelVtx8"), static_cast<double>(counterVtx8FlagHasIsSelectedCosmic));
    histReg.fill(HIST("hCosmicFlagIsSelVtx10"), static_cast<double>(counterVtx10FlagHasIsSelectedCosmic));

    histReg.fill(HIST("hColls"), static_cast<double>(colls.size()));
    histReg.fill(HIST("hCollsVtx8"), static_cast<double>(counterCollVtx8));
    histReg.fill(HIST("hCollsVtx10"), static_cast<double>(counterCollVtx10));

    PrintTime(collLoopStart, Form("DEBUG :: df_%d :: DF Reading :: coll loop Time :: ",dfCount));    
    PrintTime(Start1, Form("DEBUG :: df_%d :: DF End    :: DF Read Time :: ",dfCount));
    PrintTime(Start0, Form("DEBUG :: df_%d :: DF End    :: Elapsed Time :: ",dfCount));
    LOG(info)<<"DEBUG ::";

  }//Process function ends

  void process(aod::BCsWithTimestamps const& BCs,
               aod::Collisions const& collisions
               ,soa::Join<aod::Tracks, aod::TrackToAmbgTrk, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov> const& tracks
               ,o2::aod::AmbiguousTracks const& ambgTracks
               ,o2::aod::Origins const& Origins
              ){
    executeProcess(BCs, collisions, tracks, ambgTracks, Origins);
  }
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
  return WorkflowSpec{ 
                      adaptAnalysisTask<cosmicMuonToCollisionAssociator>(cfgc)
  };
}
