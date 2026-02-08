#ifndef NANOAOD_WRITER_HPP
#define NANOAOD_WRITER_HPP

#include <filesystem>
#include <vector>

#include "skelana_analysis.hpp"
#include "skelana.hpp"

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/SMatrixFfwd.h>
#include <Math/SVector.h>

#include <TTree.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

//DataProcessing dependencies
#include "DataProcessing/include/particleData.h"
#include "DataProcessing/include/eventSelection.h"
#include "DataProcessing/include/eventData.h"
#include "DataProcessing/include/thrustTools.h"
#include "DataProcessing/include/sphericityTools.h"

using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleModel = ROOT::Experimental::RNTupleModel;

namespace sk = skelana;

typedef ROOT::Math::XYZTVectorF XYZTVectorF;
typedef ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag > XYZVectorF;
typedef ROOT::Math::PositionVector3D< ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag > XYZPointF;

typedef ROOT::Math::SMatrixSym3F SMatrixSym3F;
typedef ROOT::Math::SMatrixSym5F SMatrixSym5F;
typedef ROOT::Math::SVector<float, 5> SVector5F;

enum class DataKind {
    data,
    gen,
    sim
};

class NanoAODWriter : public sk::Analysis
{
private:
    typedef sk::Analysis super;

public:
    NanoAODWriter(const NanoAODWriter &) = delete;
    NanoAODWriter &operator=(const NanoAODWriter &) = delete;
    virtual ~NanoAODWriter();
    static NanoAODWriter *getInstance();
    void setOutput(const std::filesystem::path &output);
    void setMC();

protected:
    NanoAODWriter();
    virtual void user00();
    virtual int user01();
    virtual void user02();
    virtual void user99();

private:
    void defineEvent(std::unique_ptr<RNTupleModel> &model);
    void fillEvent();
    void definePart(std::unique_ptr<RNTupleModel> &model);
    void fillPart();
    void defineJet(std::unique_ptr<RNTupleModel> &model);
    void fillJet();
    void defineVtx(std::unique_ptr<RNTupleModel> &model);
    void fillVtx();
    void defineSimPart(std::unique_ptr<RNTupleModel> &model);
    void fillSimPart();
    void defineGenPart(std::unique_ptr<RNTupleModel> &model);
    void fillGenPart();
    void defineSimVtx(std::unique_ptr<RNTupleModel> &model);
    void fillSimVtx();
    void defineBeamSpot(std::unique_ptr<RNTupleModel> &model);
    void fillBeamSpot();
    void defineTrac(std::unique_ptr<RNTupleModel> &model);
    void fillTrac();
    void defineMuid(std::unique_ptr<RNTupleModel> &model);
    void fillMuid();
    void defineElid(std::unique_ptr<RNTupleModel> &model);
    void fillElid();
    void defineHadid(std::unique_ptr<RNTupleModel> &model);
    void fillHadid();
    void defineBtag(std::unique_ptr<RNTupleModel> &model);
    void fillBtag();
    void fillPartLoop(particleData& pData, eventData& eData, DataKind cat = DataKind::data);
    void fillSelection(particleData& pData, eventData& eData);

    std::filesystem::path output_;
    std::unique_ptr<RNTupleWriter> writer_;
    bool mc_;

    // Create a TTree to store event-based data
    TFile* file_ = nullptr;
    TTree* out_t = nullptr;
    particleData    out_pData;
    eventData       out_eData;

    TTree* out_tgen = nullptr;
    particleData    out_pData_gen;
    eventData       out_eData_gen;

    TTree* out_tsim = nullptr;
    particleData    out_pData_sim;
    eventData       out_eData_sim;

    // Variables for event information and particle data
    float emf[particleData::nMaxPart];
    float hpc[particleData::nMaxPart];
    float hac[particleData::nMaxPart];
    float stic[particleData::nMaxPart];
    int lock[particleData::nMaxPart];
    TDatabasePDG* pdgDatabase = nullptr;

    std::shared_ptr<int> Event_runNumber_;
    std::shared_ptr<int> Event_eventNumber_;
    std::shared_ptr<int> Event_fillNumber_;
    std::shared_ptr<int> Event_date_;
    std::shared_ptr<int> Event_time_;
    std::shared_ptr<float> Event_magField_;
    std::shared_ptr<float> Event_cmEnergy_;
    std::shared_ptr<int8_t> Event_shortDstVersion_;
    std::shared_ptr<bool> Event_hadronTagT4_;
    std::shared_ptr<int16_t> Event_chargedMultT4_;
    std::shared_ptr<int16_t> Event_chargedMult_;
    std::shared_ptr<int16_t> Event_neutralMult_;
    std::shared_ptr<float> Event_totalChargedEnergy_;
    std::shared_ptr<float> Event_totalEMEnergy_;
    std::shared_ptr<float> Event_totalHadronicEnergy_;
    std::shared_ptr<std::string> Event_DSTType_;

    std::shared_ptr<int16_t> nPart_;
    std::shared_ptr<std::vector<XYZTVectorF>> Part_fourMomentum_;
    std::shared_ptr<std::vector<int8_t>> Part_charge_;
    std::shared_ptr<std::vector<int16_t>> Part_pdgId_;
    std::shared_ptr<std::vector<int>> Part_massId_;
    std::shared_ptr<std::vector<int8_t>> Part_jetIdx_;
    std::shared_ptr<std::vector<int8_t>> Part_hemisphereIdx_;
    std::shared_ptr<std::vector<int8_t>> Part_vtxCode_;
    std::shared_ptr<std::vector<int8_t>> Part_vtxIdx_;
    std::shared_ptr<std::vector<int>> Part_lock_;
    std::shared_ptr<std::vector<float>> Part_hpcShowerEnergy_;
    std::shared_ptr<std::vector<float>> Part_hpcShowerTheta_;
    std::shared_ptr<std::vector<float>> Part_hpcShowerPhi_;
    std::shared_ptr<std::vector<int>> Part_hpcParticleCode_;
    std::shared_ptr<std::vector<int>> Part_hpcNumLayers_;
    std::shared_ptr<std::vector<int>> Part_hpcLayerHitPattern_;
    std::shared_ptr<std::vector<int>> Part_hpcNumAssociatedShowers_;
    std::shared_ptr<std::vector<float>> Part_hpcTotalShowerEnergy_;
    std::shared_ptr<std::vector<float>> Part_hacShowerEnergy_;
    std::shared_ptr<std::vector<float>> Part_hacShowerTheta_;
    std::shared_ptr<std::vector<float>> Part_hacShowerPhi_;
    std::shared_ptr<std::vector<int>> Part_hacParticleCode_;
    std::shared_ptr<std::vector<int>> Part_hacNumTowers_;
    std::shared_ptr<std::vector<int>> Part_hacTowerHitPattern_;
    std::shared_ptr<std::vector<int>> Part_hacNumAssociatedShowers_;
    std::shared_ptr<std::vector<float>> Part_hacTotalShowerEnergy_;
    std::shared_ptr<std::vector<float>> Part_sticShowerEnergy_;
    std::shared_ptr<std::vector<float>> Part_sticShowerTheta_;
    std::shared_ptr<std::vector<float>> Part_sticShowerPhi_;
    std::shared_ptr<std::vector<int>> Part_sticNumTowers_;
    std::shared_ptr<std::vector<int>> Part_sticChargedTag_;
    std::shared_ptr<std::vector<int>> Part_sticSiliconVertexPos_;
    std::shared_ptr<std::vector<int16_t>> Part_simIdx_;
    std::shared_ptr<std::vector<int16_t>> Part_originVtxIdx_;
    std::shared_ptr<std::vector<int16_t>> Part_decayVtxIdx_;

    std::shared_ptr<int16_t> nJet_;
    std::shared_ptr<std::vector<XYZTVectorF>> Jet_fourMomentum_;
    std::shared_ptr<std::vector<int8_t>> Jet_charge_;

    std::shared_ptr<float> Jet_thrust_;
    std::shared_ptr<float> Jet_oblatness_;
    std::shared_ptr<std::vector<XYZVectorF>> Jet_thrustVector_;
    std::shared_ptr<float> Jet_sphericity_;
    std::shared_ptr<float> Jet_aplanarity_;
    std::shared_ptr<std::vector<XYZVectorF>> Jet_sphericityVector_;

    std::shared_ptr<int16_t> nSimPart_;
    std::shared_ptr<std::vector<XYZTVectorF>> SimPart_fourMomentum_;
    std::shared_ptr<std::vector<int16_t>> SimPart_charge_;
    std::shared_ptr<std::vector<int16_t>> SimPart_pdgId_;
    std::shared_ptr<std::vector<int16_t>> SimPart_partIdx_;
    std::shared_ptr<std::vector<int16_t>> SimPart_genIdx_;
    std::shared_ptr<std::vector<int16_t>> SimPart_originVtxIdx_;
    std::shared_ptr<std::vector<int16_t>> SimPart_decayVtxIdx_;

    std::shared_ptr<int16_t> nGenPart_;
    std ::shared_ptr<std::vector<int16_t>> GenPart_status_;
    std::shared_ptr<std::vector<int16_t>> GenPart_pdgId_;
    std::shared_ptr<std::vector<int16_t>> GenPart_parentIdx_;
    std::shared_ptr<std::vector<int16_t>> GenPart_firstChildIdx_;
    std::shared_ptr<std::vector<int16_t>> GenPart_lastChildIdx_;
    std::shared_ptr<std::vector<XYZTVectorF>> GenPart_fourMomentum_;
    std::shared_ptr<std::vector<XYZTVectorF>> GenPart_fourPosition_;
    std::shared_ptr<std::vector<float>> GenPart_tau_;
    std::shared_ptr<std::vector<int16_t>> GenPart_simIdx_;
    std::shared_ptr<std::vector<float>> GenPart_mass_;

    std::shared_ptr<int16_t> nVtx_;
    std::shared_ptr<std::vector<int16_t>> Vtx_outgoingIdx_;
    std::shared_ptr<std::vector<int16_t>> Vtx_incomingIdx_;
    std::shared_ptr<std::vector<int16_t>> Vtx_nOutgoing_;
    std::shared_ptr<std::vector<int16_t>> Vtx_ndf_;
    std::shared_ptr<std::vector<XYZPointF>> Vtx_position_;
    std::shared_ptr<std::vector<float>> Vtx_chi2_;
    std::shared_ptr<std::vector<SMatrixSym3F>> Vtx_errorMatrix_;
    std::shared_ptr<std::vector<int16_t>> Vtx_errorFlag_;
    std::shared_ptr<std::vector<int16_t>> Vtx_status_;

    std::shared_ptr<int16_t> nSimVtx_;
    std::shared_ptr<std::vector<int16_t>> SimVtx_outgoingIdx_;
    std::shared_ptr<std::vector<int16_t>> SimVtx_incomingIdx_;
    std::shared_ptr<std::vector<int16_t>> SimVtx_nOutgoing_;
    std::shared_ptr<std::vector<int16_t>> SimVtx_masscode_;
    std::shared_ptr<std::vector<XYZPointF>> SimVtx_position_;

    std::shared_ptr<XYZVectorF> Beam_position_;
    std::shared_ptr<XYZVectorF> Beam_size_;

    std::shared_ptr<std::vector<int16_t>> Trac_partIdx_;
    std::shared_ptr<std::vector<int16_t>> Trac_originVtxIdx_;
    std::shared_ptr<std::vector<int16_t>> Trac_decayVtxIdx_;
    std::shared_ptr<std::vector<SVector5F>> Trac_perigee_;
    std::shared_ptr<std::vector<ROOT::Math::SMatrixSym5F>> Trac_weightMatrix_;
    std::shared_ptr<std::vector<float>> Trac_length_;
    std::shared_ptr<std::vector<int16_t>> Trac_detectors_;
    std::shared_ptr<std::vector<float>> Trac_firstPointR_;
    std::shared_ptr<std::vector<float>> Trac_firstPointZ_;
    std::shared_ptr<std::vector<float>> Trac_chi2NoVD_;
    std::shared_ptr<std::vector<float>> Trac_chi2VD_;
    std::shared_ptr<std::vector<int16_t>> Trac_ndfNoVD_;
    std::shared_ptr<std::vector<int16_t>> Trac_ndfVD_;
    std::shared_ptr<std::vector<int16_t>> Trac_vdHitsRPhi_;
    std::shared_ptr<std::vector<int16_t>> Trac_vdHitsZ_;
    std::shared_ptr<std::vector<float>> Trac_resRPhiFirstPoint_;
    std::shared_ptr<std::vector<float>> Trac_errorResRPhiFirstPoint_;
    std::shared_ptr<std::vector<float>> Trac_resZFirstPoint_;
    std::shared_ptr<std::vector<float>> Trac_errorResZFirstPoint_;
    std::shared_ptr<std::vector<float>> Trac_impParRPhi_;
    std::shared_ptr<std::vector<float>> Trac_impParZ_;
    std::shared_ptr<std::vector<float>> Trac_impParToVertexRPhi_;
    std::shared_ptr<std::vector<float>> Trac_impParToVertexZ_;
    std::shared_ptr<std::vector<float>> Trac_impParToBeamSpotRPhi_;
    std::shared_ptr<std::vector<float>> Trac_chi2VDHits_;

    std::shared_ptr<std::vector<int>> Muid_partIdx_;
    std::shared_ptr<std::vector<int>> Muid_tag_;
    std::shared_ptr<std::vector<float>> Muid_looseChi2_;
    std::shared_ptr<std::vector<int16_t>> Muid_hitPattern_;

    std::shared_ptr<std::vector<int>> Elid_partIdx_;
    std::shared_ptr<std::vector<int>> Elid_tag_;
    std::shared_ptr<std::vector<int16_t>> Elid_gammaConversion_;
    std::shared_ptr<std::vector<float>> Elid_px_;
    std::shared_ptr<std::vector<float>> Elid_py_;
    std::shared_ptr<std::vector<float>> Elid_pz_;
    std::shared_ptr<std::vector<float>> Elou_eOverP_;
    std::shared_ptr<std::vector<float>> Elou_deltaZHpc_;
    std::shared_ptr<std::vector<float>> Elou_deltaPhiDirectionHpc_;
    std::shared_ptr<std::vector<float>> Elou_deltaPhiPositionHpc_;
    std::shared_ptr<std::vector<float>> Elou_dedxMeasurement_;
    std::shared_ptr<std::vector<float>> Elou_dedxMeasurementError_;
    std::shared_ptr<std::vector<float>> Elou_dedxNumTpcWires_;
    std::shared_ptr<std::vector<float>> Elou_probFromEOverP_;
    std::shared_ptr<std::vector<float>> Elou_probFromShowerFit_;
    std::shared_ptr<std::vector<float>> Elou_probFromDeltaZHpc_;
    std::shared_ptr<std::vector<float>> Elou_probFromDeltaPhiDirectionHpc_;
    std::shared_ptr<std::vector<float>> Elou_probFromDedxElectron_;
    std::shared_ptr<std::vector<float>> Elou_probFromDedxPion_;
    std::shared_ptr<std::vector<float>> Elou_hpcElectronProbability_;

    std::shared_ptr<std::vector<int>> Haid_sign_;
    std::shared_ptr<std::vector<int>> Haid_kaonDedx_;
    std::shared_ptr<std::vector<int>> Haid_protonDedx_;
    std::shared_ptr<std::vector<int>> Haid_kaonRich_;
    std::shared_ptr<std::vector<int>> Haid_protonRich_;
    std::shared_ptr<std::vector<int>> Haid_pionRich_;
    std::shared_ptr<std::vector<float>> Haid_kaonCombined_;
    std::shared_ptr<std::vector<float>> Haid_protonCombined_;
    std::shared_ptr<std::vector<int>>Haid_richQuality_;

    std::shared_ptr<std::vector<int8_t>> Haidn_pionTag_;
    std::shared_ptr<std::vector<int8_t>> Haidn_kaonTag_;
    std::shared_ptr<std::vector<int8_t>> Haidn_protonTag_;
    std::shared_ptr<std::vector<int8_t>> Haidn_heavyTag_;
    std::shared_ptr<std::vector<int8_t>> Haidn_pionTrackSelection_;
    std::shared_ptr<std::vector<int8_t>> Haidn_kaonTrackSelection_;
    std::shared_ptr<std::vector<int8_t>> Haidn_protonTrackSelection_;
    std::shared_ptr<std::vector<int8_t>> Haidn_heavyTrackSelection_;

    std::shared_ptr<std::vector<int8_t>> Haidr_pionTag_;
    std::shared_ptr<std::vector<int8_t>> Haidr_kaonTag_;
    std::shared_ptr<std::vector<int8_t>> Haidr_protonTag_;
    std::shared_ptr<std::vector<int8_t>> Haidr_heavyTag_;
    std::shared_ptr<std::vector<int8_t>> Haidr_electronTag_;
    std::shared_ptr<std::vector<int8_t>> Haidr_selectionFlag_;

    std::shared_ptr<std::vector<int8_t>> Haide_pionTag_;
    std::shared_ptr<std::vector<int8_t>> Haide_kaonTag_;
    std::shared_ptr<std::vector<int8_t>> Haide_protonTag_;
    std::shared_ptr<std::vector<int8_t>> Haide_heavyTag_;
    std::shared_ptr<std::vector<int8_t>> Haide_electronTag_;
    std::shared_ptr<std::vector<int8_t>> Haide_selectionFlag_;

    std::shared_ptr<std::vector<int8_t>> Haidc_pionTag_;
    std::shared_ptr<std::vector<int8_t>> Haidc_kaonTag_;
    std::shared_ptr<std::vector<int8_t>> Haidc_protonTag_;
    std::shared_ptr<std::vector<int8_t>> Haidc_heavyTag_;
    std::shared_ptr<std::vector<int8_t>> Haidc_electronTag_;
    std::shared_ptr<std::vector<int8_t>> Haidc_selectionFlag_;

    std::shared_ptr<std::vector<float>>  Dedx_value_;
    std::shared_ptr<std::vector<float>>  Dedx_width_;
    std::shared_ptr<std::vector<int16_t>> Dedx_nrWires_; 
    std::shared_ptr<std::vector<float>> Dedx_gapWires_; 
    std::shared_ptr<std::vector<float>>  Dedx_error_;
    std::shared_ptr<std::vector<float>>  Dedx_valueVD_; 
    std::shared_ptr<std::vector<int16_t>> Dedx_nrVDHits_;

    std::shared_ptr<std::vector<float>> Rich_theg_;
    std::shared_ptr<std::vector<float>> Rich_sigg_;
    std::shared_ptr<std::vector<int16_t>> Rich_nphg_;
    std::shared_ptr<std::vector<float>> Rich_nepg_;
    std::shared_ptr<std::vector<int16_t>> Rich_flagg_;
    std::shared_ptr<std::vector<float>> Rich_thel_;
    std::shared_ptr<std::vector<float>> Rich_sigl_;
    std::shared_ptr<std::vector<int16_t>> Rich_nphl_;
    std::shared_ptr<std::vector<float>> Rich_nepl_;
    std::shared_ptr<std::vector<int16_t>> Rich_flagl_;

    std::shared_ptr<std::vector<float>> Btag_probNegIP_;
    std::shared_ptr<std::vector<float>> Btag_probPosIP_;
    std::shared_ptr<std::vector<float>> Btag_probAllIP_;
    std::shared_ptr<XYZVectorF> Btag_thrustVector_; 

};

#endif // NANOAOD_WRITER_HPP
