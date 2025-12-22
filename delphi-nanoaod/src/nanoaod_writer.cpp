#include "nanoaod_writer.hpp"
#include "phdst.hpp"
#include <Math/PositionVector3D.h>
#include <Math/DisplacementVector3D.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/SMatrix.h>

#include <TTree.h>

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

//Additional common block to extract jet releateted variables
extern "C" struct
{
    float thr;
    float obl;
    float sph;
    float apl;
} psjval_;

float &SPH = psjval_.sph;
float &APL = psjval_.apl;
float &THR = psjval_.thr;
float &OBL = psjval_.obl;

// MakeField with type deduction
template <typename T>
void MakeField(std::unique_ptr<RNTupleModel> &model, const std::string &name, const std::string &description, std::shared_ptr<T> &ptr)
{
    ptr = model->MakeField<T>({name, description});
}

// Fill a vector with n elements using a function f
template <typename T, typename Func>
void fillVector(std::shared_ptr<std::vector<T>> &ptr, int i, int n, Func f)
{
    ptr->clear();
    std::generate_n(std::back_inserter(*ptr), n, [i, f]() mutable
                    { return f(i++); });
}

NanoAODWriter::NanoAODWriter() : output_("nanoaod.root"), mc_(false) {};

NanoAODWriter::~NanoAODWriter() {};

NanoAODWriter *NanoAODWriter::getInstance()
{
    if (instance_ == nullptr)
    {
        instance_ = new NanoAODWriter();
    }
    return static_cast<NanoAODWriter *>(instance_);
};

void NanoAODWriter::setOutput(const std::filesystem::path &output)
{
    output_ = output;
};

void NanoAODWriter::setMC()
{
    mc_ = true;
};

void NanoAODWriter::user00()
{
    // std::cout << "NanoAODWriter::user00: Initialising" << std::endl;

    super::user00();

    std::unique_ptr<RNTupleModel> model = RNTupleModel::Create();

    defineEvent(model);
    definePart(model);
    defineVtx(model);

    if (sk::IFLJET > 0)
    {
        defineJet(model);
    }
    if (mc_ && sk::IFLSIM > 0)
    {
        defineSimPart(model);
        defineGenPart(model);
        defineSimVtx(model);
    }
    if (sk::IFLTRA > 0)
    {
        defineTrac(model);
    }
    if (sk::IFLMUO > 0)
    {
        defineMuid(model);
    }
    if (sk::IFLELE > 0)
    {
        defineElid(model);
    }
    if (sk::IFLHAD > 0)
    {
        defineHadid(model);
    }
    if (sk::IFLBTG > 0)
    {
        defineBtag(model);
    }
    // writer_ = RNTupleWriter::Recreate(std::move(model), "Events", output_.string());

    file_ = TFile::Open(output_.string().c_str(), "recreate");
    writer_ = RNTupleWriter::Append(std::move(model), "Events", *file_);
    // file_ = TFile::CurrentFile();
    // file_->cd();
    out_t = new TTree("t", "t");
    out_tgen = new TTree("tgen", "tgen");
    out_tsim = new TTree("tgenBefore", "tgenBefore");
    out_t->SetDirectory(file_);
    out_tgen->SetDirectory(file_);
    out_tsim->SetDirectory(file_);

    // register branches
    do_chThrust           = false;
    do_neuThrust          = false;
    do_thrustCorr         = false;
    do_thrustCorrInverse  = false;
    do_thrustMissP        = true;
    out_pData.SetBranchWrite(out_t, 1);
    out_eData.SetBranchWrite(out_t, 0);

    out_pData_gen.SetBranchWrite(out_tgen, 1);
    out_eData_gen.SetBranchWrite(out_tgen, 0);

    out_pData_sim.SetBranchWrite(out_tsim, 1);
    out_eData_sim.SetBranchWrite(out_tsim, 0);

    out_t->Branch("EMF", emf, "EMF[nParticle]/F");
    out_t->Branch("HPC", hpc, "HPC[nParticle]/F");
    out_t->Branch("HAC", hac, "HAC[nParticle]/F");
    out_t->Branch("STIC", stic, "STIC[nParticle]/F");
    out_t->Branch("LOCK", lock, "LOCK[nParticle]/I");

    pdgDatabase = TDatabasePDG::Instance();
};

int NanoAODWriter::user01()
{
    return super::user01();
};

void NanoAODWriter::user02()
{
    super::user02();

    fillEvent();
    fillPart();
    fillVtx();
    if (sk::IFLJET > 0)
    {
        fillJet();
    }
    if (mc_ && sk::IFLSIM > 0)
    {
        fillSimPart();
        fillGenPart();
        fillSimVtx();
    }
    if (sk::IFLTRA > 0)
    {
        fillTrac();
    }
    if (sk::IFLMUO > 0)
    {
        fillMuid();
    }
    if (sk::IFLELE > 0)
    {
        fillElid();
    }
    if (sk::IFLBTG > 0)
    {
        fillBtag();
    }
    if (sk::IFLHAD > 0)
    {
        fillHadid();
    }

    fillPartLoop(out_pData, out_eData, DataKind::data);
    fillSelection(out_pData, out_eData);

    writer_->Fill();
    out_t->Fill();
    if (mc_) {
        fillPartLoop(out_pData_gen, out_eData_gen, DataKind::gen);
        fillSelection(out_pData_gen, out_eData_gen);
        out_tgen->Fill();

        fillPartLoop(out_pData_sim, out_eData_sim, DataKind::sim);
        fillSelection(out_pData_sim, out_eData_sim);
        out_tsim->Fill();
    }
};

void NanoAODWriter::defineEvent(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Event_runNumber", "Run number", Event_runNumber_);
    MakeField(model, "Event_evtNumber", "Event number", Event_eventNumber_);
    MakeField(model, "Event_fillNumber", "Fill number", Event_fillNumber_);
    MakeField(model, "Event_date", "Date", Event_date_);
    MakeField(model, "Event_time", "Time", Event_time_);
    MakeField(model, "Event_magField", "Magnetic field", Event_magField_);
    MakeField(model, "Event_shortDstVersion", "Short DST version", Event_shortDstVersion_);
    MakeField(model, "Event_hadronTagT4", "Team4 Hadron tag", Event_hadronTagT4_);
    MakeField(model, "Event_chargedMultT4", "Charged Multiplicity (Team 4)", Event_chargedMultT4_);
    MakeField(model, "Event_chargedMult", "Charged Multiplicity", Event_chargedMult_);
    MakeField(model, "Event_neutralMult", "Neutral Multiplicity", Event_neutralMult_);
    MakeField(model, "Event_cmEnergy", "Center of mass energy", Event_cmEnergy_);
    MakeField(model, "Event_totalChargedEnergy", "Total charged energy", Event_totalChargedEnergy_);
    MakeField(model, "Event_totalEMEnergy", "Total electromagnetic energy", Event_totalEMEnergy_);
    MakeField(model, "Event_totalHadronicEnergy", "Total hadronic energy", Event_totalHadronicEnergy_);
    MakeField(model, "Event_DSTType", "DST type", Event_DSTType_);
}

void NanoAODWriter::fillEvent()
{
    *Event_runNumber_ = phdst::IIIRUN;
    *Event_eventNumber_ = phdst::IIIEVT;
    *Event_fillNumber_ = phdst::IIFILL;
    *Event_date_ = phdst::IIIDAT;
    *Event_time_ = phdst::IIITIM;
    *Event_magField_ = sk::BMAG;
    *Event_shortDstVersion_ = sk::ISVER;
    *Event_hadronTagT4_ = sk::IHAD4 > 0;
    *Event_chargedMultT4_ = sk::NCTR4;
    *Event_chargedMult_ = sk::NCTRK;
    *Event_neutralMult_ = sk::NNTRK;
    *Event_cmEnergy_ = sk::ECMAS;
    *Event_totalChargedEnergy_ = sk::ECHAR;
    *Event_totalEMEnergy_ = sk::EMNEU;
    *Event_totalHadronicEnergy_ = sk::EHNEU;
    *Event_DSTType_ = sk::CDTYPE();


    out_pData.RunNo = *Event_runNumber_;
    out_pData.EventNo = *Event_eventNumber_;
    out_pData_gen.RunNo = out_pData.RunNo;
    out_pData_gen.EventNo = out_pData.EventNo;
    out_pData_sim.RunNo = out_pData.RunNo;
    out_pData_sim.EventNo = out_pData.EventNo;

    out_pData.year = 0; // temp //
	  // out_pData.subDir = -999;
	  // out_pData.process = -999;
    out_pData.source = 0; // temp //
    out_pData.isMC = mc_;
    out_pData.isOnres = false;
    out_pData_gen.isMC = true;
    out_pData_gen.isOnres = false;
    out_pData_sim.isMC = true;
    out_pData_sim.isOnres = false;
    out_pData.particleWeight = 1;
    out_pData_gen.particleWeight = 1;
    out_pData_sim.particleWeight = 1;
	  // out_pData.uniqueID = 0;
	  // out_pData.Energy = 0;
	  // out_pData.bFlag = -999;
	  // out_pData.bx = -999;
	  // out_pData.by = -999;
	  // out_pData.ebx = -999;
	  // out_pData.eby = -999;
    out_pData.Energy = *Event_cmEnergy_;
    out_pData_gen.Energy = *Event_cmEnergy_;
    out_pData_sim.Energy = *Event_cmEnergy_;
}

void NanoAODWriter::definePart(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nPart", "Number of particles", nPart_);
    MakeField(model, "Part_fourMomentum", "Particle 4-momentum", Part_fourMomentum_);
    MakeField(model, "Part_charge", "Particle charge", Part_charge_);
    MakeField(model, "Part_pdgId", "Particle PDG ID", Part_pdgId_);
    MakeField(model, "Part_massId", "Particle mass ID", Part_massId_);
    MakeField(model, "Part_jetIdx", "Jet index", Part_jetIdx_);
    MakeField(model, "Part_hemisphereIdx", "Hemisphere index", Part_hemisphereIdx_);
    MakeField(model, "Part_vtxCode", "Vertex code", Part_vtxCode_);
    MakeField(model, "Part_vtxIdx", "Vertex index", Part_vtxIdx_);
    MakeField(model, "Part_lock", "Particle lock", Part_lock_);

    if (mc_ && sk::IFLSIM > 0)
    {
        MakeField(model, "Part_simIdx", "Particle simulation index", Part_simIdx_);
        MakeField(model, "Part_originVtxIdx", "Particle origin vertex index", Part_originVtxIdx_);
        MakeField(model, "Part_decayVtxIdx", "Particle decay vertex index", Part_decayVtxIdx_);
    }
}

void NanoAODWriter::fillPart()
{
    *nPart_ = sk::NVECP;
    fillVector(Part_fourMomentum_, sk::LVPART, sk::NVECP, [](int i)
               { return XYZTVectorF(sk::VECP(1, i), sk::VECP(2, i), sk::VECP(3, i), sk::VECP(4, i)); });
    fillVector(Part_charge_, sk::LVPART, sk::NVECP, [](int i)
               { return int(sk::VECP(7, i)); });
    fillVector(Part_pdgId_, sk::LVPART, sk::NVECP, [](int i)
               { return int(sk::VECP(8, i)); });
    fillVector(Part_massId_, sk::LVPART, sk::NVECP, [](int i)
               { return sk::IVECP(9, i); });
    fillVector(Part_jetIdx_, sk::LVPART, sk::NVECP, [](int i)
               { return (sk::IVECP(10, i) / 1000) % 10 - 1; });
    fillVector(Part_hemisphereIdx_, sk::LVPART, sk::NVECP, [](int i)
               { return (sk::IVECP(10, i) / 100) % 10 - 1; });
    fillVector(Part_vtxCode_, sk::LVPART, sk::NVECP, [](int i)
               { return (sk::IVECP(10, i) / 10) % 10; });
    fillVector(Part_vtxIdx_, sk::LVPART, sk::NVECP, [](int i)
               { return sk::IVECP(10, i) % 10 - 1; });
    fillVector(Part_lock_, sk::LVPART, sk::NVECP, [](int i)
               { return sk::LVLOCK(i); });

    if (mc_ && sk::IFLSIM > 0)
    {
        fillVector(Part_simIdx_, sk::LVPART, sk::NVECP, [](int i)
                   { return sk::IPAST(i) - 1; });
        fillVector(Part_originVtxIdx_, sk::LVPART, sk::NVECP, [](int i)
                   { return sk::IPAPV(1, i) - 1; });
        fillVector(Part_decayVtxIdx_, sk::LVPART, sk::NVECP, [](int i)
                   { return sk::IPAPV(2, i) - 1; });
    }

}

void NanoAODWriter::defineJet(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nJet", "Number of jets", nJet_);
    MakeField(model, "Jet_fourMomentum", "Jet 4-momentum", Jet_fourMomentum_);
    MakeField(model, "Jet_charge", "Jet charge", Jet_charge_);

    if ((sk::IFLJET / 10) % 10 == 0)
    {
        MakeField(model, "Jet_oblatness", "Jet oblatness", Jet_oblatness_);
        MakeField(model, "Jet_thrust", "Jet thrust", Jet_thrust_);
        MakeField(model, "Jet_thrustVector", "Jet thrust vector", Jet_thrustVector_);
    }
    if ((sk::IFLJET / 100) % 10 == 0)
    {
        MakeField(model, "Jet_sphericity", "Jet sphericity", Jet_sphericity_);
        MakeField(model, "Jet_aplanarity", "Jet aplanarity", Jet_aplanarity_);
        MakeField(model, "Jet_sphericityVector", "Jet sphericity vector", Jet_sphericityVector_);
    }
}

void NanoAODWriter::fillJet()
{
    *nJet_ = sk::NJET;
    fillVector(Jet_fourMomentum_, sk::LVJET + 1, sk::NJET, [](int i)
               { return XYZTVectorF(sk::VECP(1, i), sk::VECP(2, i), sk::VECP(3, i), sk::VECP(4, i)); });
    fillVector(Jet_charge_, sk::LVJET + 1, sk::NJET, [](int i)
               { return int(sk::VECP(6, i)); });

    if ((sk::IFLJET / 10) % 10 == 0)
    {
        *Jet_thrust_ = THR;
        *Jet_oblatness_ = OBL;
        fillVector(Jet_thrustVector_, sk::LVTHRU, 3, [](int i)
                   { return XYZVectorF(sk::VECP(1, i), sk::VECP(2, i), sk::VECP(3, i)); });
    }

    if ((sk::IFLJET / 100) % 10 == 0)
    {
        *Jet_sphericity_ = SPH;
        *Jet_aplanarity_ = APL;
        fillVector(Jet_sphericityVector_, sk::LVSPHE, 3, [](int i)
                   { return XYZVectorF(sk::VECP(1, i), sk::VECP(2, i), sk::VECP(3, i)); });
    }
}

void NanoAODWriter::defineVtx(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nVtx", "Number of vertices", nVtx_);
    MakeField(model, "Vtx_outgoingIdx", "First outgoing particle index", Vtx_outgoingIdx_);
    MakeField(model, "Vtx_incomingIdx", "Incoming particle index", Vtx_incomingIdx_);
    MakeField(model, "Vtx_nOutgoing", "Number of outgoing particles", Vtx_nOutgoing_);
    MakeField(model, "Vtx_ndf", "Number of degrees of freedom of the Vertex Fit", Vtx_ndf_);
    MakeField(model, "Vtx_position", "Vertex position", Vtx_position_);
    MakeField(model, "Vtx_chi2", "Vertex fit chi2", Vtx_chi2_);
    MakeField(model, "Vtx_errorMatrix", "Vertex error matrix", Vtx_errorMatrix_);
    MakeField(model, "Vtx_errorFlag", "Vertex error flag", Vtx_errorFlag_);
    MakeField(model, "Vtx_status", "Vertex status", Vtx_status_);
}

void NanoAODWriter::fillVtx()
{
    *nVtx_ = sk::NVTX;
    fillVector(Vtx_outgoingIdx_, 1, sk::NVTX, [](int i)
               { return sk::KVTX(1, i); });
    fillVector(Vtx_incomingIdx_, 1, sk::NVTX, [](int i)
               { return sk::KVTX(2, i); });
    fillVector(Vtx_nOutgoing_, 1, sk::NVTX, [](int i)
               { return sk::KVTX(3, i); });
    fillVector(Vtx_ndf_, 1, sk::NVTX, [](int i)
               { return sk::KVTX(4, i); });
    fillVector(Vtx_position_, 1, sk::NVTX, [](int i)
               { return XYZPointF(sk::QVTX(6, i), sk::QVTX(7, i), sk::QVTX(8, i)); });
    fillVector(Vtx_chi2_, 1, sk::NVTX, [](int i)
               { return sk::QVTX(9, i); });
    fillVector(Vtx_errorMatrix_, 1, sk::NVTX, [](int i)
               { return ROOT::Math::SMatrixSym3F(&sk::QVTX(10, i), 6); });
    fillVector(Vtx_errorFlag_, 1, sk::NVTX, [](int i)
               { return sk::KVTX(16, i); });
    fillVector(Vtx_status_, 1, sk::NVTX, [](int i)
               { return sk::KVTX(17, i); });
}

void NanoAODWriter::defineSimPart(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nSimPart", "Number of simulated particles", nSimPart_);
    MakeField(model, "SimPart_fourMomentum", "Simulated particle 4-momentum", SimPart_fourMomentum_);
    MakeField(model, "SimPart_charge", "Simulated particle charge", SimPart_charge_);
    MakeField(model, "SimPart_pdgId", "Simulated particle PDG ID", SimPart_pdgId_);
    MakeField(model, "SimPart_partIdx", "Simulated particle index", SimPart_partIdx_);
    MakeField(model, "SimPart_genIdx", "Simulated particle generation index", SimPart_genIdx_);
    MakeField(model, "SimPart_originVtxIdx", "Simulated particle origin vertex index", SimPart_originVtxIdx_);
    MakeField(model, "SimPart_decayVtxIdx", "Simulated particle decay vertex index", SimPart_decayVtxIdx_);
}

void NanoAODWriter::fillSimPart()
{
    *nSimPart_ = sk::NVECMC;
    fillVector(SimPart_fourMomentum_, sk::MTRACK + sk::LVPART, sk::NVECMC, [](int i)
    { return XYZTVectorF(sk::VECP(1, i), sk::VECP(2, i), sk::VECP(3, i), sk::VECP(4, i)); });
    fillVector(SimPart_charge_, sk::MTRACK + sk::LVPART, sk::NVECMC, [](int i)
    { return int(sk::VECP(7, i)); });
    fillVector(SimPart_pdgId_, sk::MTRACK + sk::LVPART, sk::NVECMC, [](int i)
    { return int(sk::VECP(8, i)); });
    fillVector(SimPart_partIdx_, sk::LVPART, sk::NVECMC, [](int i)
    { return sk::ISTPA(i) - 1; });
    fillVector(SimPart_genIdx_, sk::LVPART, sk::NVECMC, [](int i)
    { return sk::ISTSH(i) - 1; });
    fillVector(SimPart_originVtxIdx_, sk::LVPART, sk::NVECMC, [](int i)
    { return sk::ISTVX(1, i) - 1; });
    fillVector(SimPart_decayVtxIdx_, sk::LVPART, sk::NVECMC, [](int i)
    { return sk::ISTVX(2, i) - 1; });
}

void NanoAODWriter::defineGenPart(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nGenPart", "Number of generated particles", nGenPart_);
    MakeField(model, "GenPart_status", "Generated particle status", GenPart_status_);
    MakeField(model, "GenPart_pdgId", "Generated particle PDG ID", GenPart_pdgId_);
    MakeField(model, "GenPart_parentIdx", "Generated particle parent index", GenPart_parentIdx_);
    MakeField(model, "GenPart_firstChildIdx", "Generated particle first child index", GenPart_firstChildIdx_);
    MakeField(model, "GenPart_lastChildIdx", "Generated particle last child index", GenPart_lastChildIdx_);
    MakeField(model, "GenPart_vector", "Generated particle 4-momentum", GenPart_fourMomentum_);
    MakeField(model, "GenPart_vertex", "Generated particle 4-position", GenPart_fourPosition_);
    MakeField(model, "GenPart_tau", "Generated particle lifetime", GenPart_tau_);
    MakeField(model, "GenPart_simIdx", "Generated particle simulation index", GenPart_simIdx_);
    MakeField(model, "GenPart_mass", "Generated particle mass", GenPart_mass_);
}

void NanoAODWriter::fillGenPart()
{
    *nGenPart_ = sk::NP;
    fillVector(GenPart_status_, 1, sk::NP, [](int i)
               { return sk::KP(i, 1); });
    fillVector(GenPart_pdgId_, 1, sk::NP, [](int i)
               { return sk::KP(i, 2); });
    fillVector(GenPart_parentIdx_, 1, sk::NP, [](int i)
               { return sk::KP(i, 3) - 1; });
    fillVector(GenPart_firstChildIdx_, 1, sk::NP, [](int i)
               { return sk::KP(i, 4) - 1; });
    fillVector(GenPart_lastChildIdx_, 1, sk::NP, [](int i)
               { return sk::KP(i, 5) - 1; });
    fillVector(GenPart_fourMomentum_, 1, sk::NP, [](int i)
               { return XYZTVectorF(sk::PP(i, 1), sk::PP(i, 2), sk::PP(i, 3), sk::PP(i, 4)); });
    fillVector(GenPart_fourPosition_, 1, sk::NP, [](int i)
               { return XYZTVectorF(sk::VP(i, 1), sk::VP(i, 2), sk::VP(i, 3), sk::VP(i, 4)); });
    fillVector(GenPart_tau_, 1, sk::NP, [](int i)
               { return sk::VP(i, 5); });
    fillVector(GenPart_simIdx_, 1, sk::NP, [](int i)
               { return sk::ISHST(i) - 1; });

    fillVector(GenPart_mass_, 1, sk::NP, [](int i)
    { return sk::PP(i, 5) - 1; });

}

void NanoAODWriter::defineSimVtx(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nSimVtx", "Number of simulated vertices", nSimVtx_);
    MakeField(model, "SimVtx_outgoingIdx", "First outgoing particle index", SimVtx_outgoingIdx_);
    MakeField(model, "SimVtx_incomingIdx", "Incoming particle index", SimVtx_incomingIdx_);
    MakeField(model, "SimVtx_nOut", "Number of outgoing particles", SimVtx_nOutgoing_);
    MakeField(model, "SimVtx_masscode", "Mass code", SimVtx_masscode_);
    MakeField(model, "SimVtx_position", "Position", SimVtx_position_);
}

void NanoAODWriter::fillSimVtx()
{
    *nSimVtx_ = sk::NVTXMC;
    fillVector(SimVtx_outgoingIdx_, sk::NVTXMX + 1, sk::NVTXMC, [](int i)
    { return sk::KVTX(1, i) - 1; });
    fillVector(SimVtx_incomingIdx_, sk::NVTXMX + 1, sk::NVTXMC, [](int i)
    { return sk::KVTX(2, i) - 1; });
    fillVector(SimVtx_nOutgoing_, sk::NVTXMX + 1, sk::NVTXMC, [](int i)
    { return sk::KVTX(3, i); });
    fillVector(SimVtx_masscode_, sk::NVTXMX + 1, sk::NVTXMC, [](int i)
    { return sk::KVTX(4, i); });
    fillVector(SimVtx_position_, sk::NVTXMX + 1, sk::NVTXMC, [](int i)
    { return XYZPointF(sk::QVTX(1, i), sk::QVTX(2, i), sk::QVTX(3, i)); });
}

void NanoAODWriter::defineTrac(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Trac_partIdx", "Particle track index", Trac_partIdx_);
    MakeField(model, "Trac_originVtxIdx", "Origin vertex index", Trac_originVtxIdx_);
    MakeField(model, "Trac_decayVtxIdx", "Decay vertex index", Trac_decayVtxIdx_);
    MakeField(model, "Trac_perigee", "Perigee parameter", Trac_perigee_);
    MakeField(model, "Trac_weightMatrix", "Track weight matrix", Trac_weightMatrix_);
    MakeField(model, "Trac_length", "Track length", Trac_length_);
    MakeField(model, "Trac_detectors", "Track detectors", Trac_detectors_);
    MakeField(model, "Trac_firstPointR", "First point in R", Trac_firstPointR_);
    MakeField(model, "Trac_firstPointZ", "First point in z", Trac_firstPointZ_);
    MakeField(model, "Trac_chi2NoVD", "Track fit chi2 no vertex detector", Trac_chi2NoVD_);
    MakeField(model, "Trac_chi2VD", "Track fit chi2 vertex detector", Trac_chi2VD_);
    MakeField(model, "Trac_ndfNoVD", "Track fit ndf no vertex detector", Trac_ndfNoVD_);
    MakeField(model, "Trac_ndfVD", "Track fit ndf vertex detector", Trac_ndfVD_);
    MakeField(model, "Trac_vdHitsRPhi", "Hits in vertex detector RPhi", Trac_vdHitsRPhi_);
    MakeField(model, "Trac_vdHitsZ", "Hits in vertex detector Z", Trac_vdHitsZ_);
    MakeField(model, "Trac_resRPhiFirstPoint", "Residual in RPhi to first point", Trac_resRPhiFirstPoint_);
    MakeField(model, "Trac_errorResRPhiFirstPoint", "Error on residual RPhi to first point", Trac_errorResRPhiFirstPoint_);
    MakeField(model, "Trac_resZFirstPoint", "Residual in Z to first point", Trac_resZFirstPoint_);
    MakeField(model, "Trac_errorResZFirstPoint", "Error on residual Z to first point", Trac_errorResZFirstPoint_);
    MakeField(model, "Trac_impParToVertexRPhi", "Impact parameter to vertex in RPhi (geometric sign)", Trac_impParToVertexRPhi_);
    MakeField(model, "Trac_impParToVertexZ", "Impact parameter to vertex in Z (geometric sign)", Trac_impParToVertexZ_);
    MakeField(model, "Trac_impParToBeamSpotRPhi", "Impact parameter to beam spot in RPhi (geometric sign)", Trac_impParToBeamSpotRPhi_);
    MakeField(model, "Trac_chi2VDHits", "Chi2 vertex detector", Trac_chi2VDHits_);
}

void NanoAODWriter::fillTrac()
{
    Trac_partIdx_->clear();
    Trac_originVtxIdx_->clear();
    Trac_decayVtxIdx_->clear();
    Trac_perigee_->clear();
    Trac_weightMatrix_->clear();
    Trac_length_->clear();
    Trac_detectors_->clear();
    Trac_firstPointR_->clear();
    Trac_firstPointZ_->clear();
    Trac_chi2NoVD_->clear();
    Trac_chi2VD_->clear();
    Trac_ndfNoVD_->clear();
    Trac_ndfVD_->clear();
    Trac_vdHitsRPhi_->clear();
    Trac_vdHitsZ_->clear();
    Trac_resRPhiFirstPoint_->clear();
    Trac_errorResRPhiFirstPoint_->clear();
    Trac_resZFirstPoint_->clear();
    Trac_errorResZFirstPoint_->clear();
    Trac_impParToVertexZ_->clear();
    Trac_impParToVertexRPhi_->clear();
    Trac_impParToBeamSpotRPhi_->clear();
    Trac_chi2VDHits_->clear();

    for (int i = sk::LVPART; i <= sk::NVECP; i++)
    {
        // if (int(sk::VECP(7, i)) != 0)
        // {
            Trac_partIdx_->push_back(i - 1);
            Trac_originVtxIdx_->push_back(sk::KTRAC(1, i) - 1);
            Trac_decayVtxIdx_->push_back(sk::KTRAC(2, i) - 1);
            Trac_perigee_->push_back(SVector5F(&sk::QTRAC(3, i), 5));
            Trac_weightMatrix_->push_back(SMatrixSym5F(&sk::QTRAC(8, i), 15));
            Trac_length_->push_back(sk::QTRAC(24, i));
            Trac_detectors_->push_back(sk::KTRAC(25, i));
            Trac_firstPointR_->push_back(sk::QTRAC(26, i));
            Trac_firstPointZ_->push_back(sk::QTRAC(27, i));
            Trac_chi2NoVD_->push_back(sk::QTRAC(28, i));
            Trac_chi2VD_->push_back(sk::QTRAC(29, i));
            Trac_ndfNoVD_->push_back(sk::KTRAC(30, i));
            Trac_ndfVD_->push_back(sk::KTRAC(31, i));
            Trac_vdHitsRPhi_->push_back(sk::KTRAC(32, i));
            Trac_vdHitsZ_->push_back(sk::KTRAC(33, i));
            Trac_resRPhiFirstPoint_->push_back(sk::QTRAC(34, i));
            Trac_errorResRPhiFirstPoint_->push_back(sk::QTRAC(35, i));
            Trac_resZFirstPoint_->push_back(sk::QTRAC(36, i));
            Trac_errorResZFirstPoint_->push_back(sk::QTRAC(37, i));
            Trac_impParToVertexRPhi_->push_back(sk::QTRAC(38, i));
            Trac_impParToVertexZ_->push_back(sk::QTRAC(39, i));
            Trac_impParToBeamSpotRPhi_->push_back(sk::QTRAC(40, i));
            Trac_chi2VDHits_->push_back(sk::QTRAC(42, i));
        // }
    }
}

void NanoAODWriter::defineMuid(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Muid_partIdx", "Muon ID link to track ID", Muid_partIdx_);
    MakeField(model, "Muid_tag", "Muon ID tag", Muid_tag_);
    MakeField(model, "Muid_looseChi2", "Muon ID loose chi2", Muid_looseChi2_);
    MakeField(model, "Muid_hitPattern", "Muon ID hit pattern", Muid_hitPattern_);
}

void NanoAODWriter::fillMuid()
{
    Muid_partIdx_->clear();
    Muid_tag_->clear();
    Muid_looseChi2_->clear();
    Muid_hitPattern_->clear();

    for (int i = sk::LVPART; i <= sk::NVECP; i++)
      {
        // if (sk::KMUID(1, i) != 0)
        // {
            Muid_partIdx_->push_back(i - 1);
            Muid_tag_->push_back(sk::KMUID(1, i));
            Muid_looseChi2_->push_back(sk::QMUID(2, i));
            Muid_hitPattern_->push_back(sk::KMUID(3, i));
        // }
    }
}

void NanoAODWriter::defineElid(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Elid_partIdx", "Electron ID link to track ID", Elid_partIdx_);
    MakeField(model, "Elid_tag", "Electron ID tag", Elid_tag_);
    MakeField(model, "Elid_gammaConversion", "Electron ID gamma conversion", Elid_gammaConversion_);
    MakeField(model, "Elid_px", "Best electron px estimation", Elid_px_);
    MakeField(model, "Elid_py", "Best electron py estimation", Elid_py_);
    MakeField(model, "Elid_pz", "Best electron pz estimation", Elid_pz_);
}

void NanoAODWriter::fillElid()
{
    Elid_partIdx_->clear();
    Elid_tag_->clear();
    Elid_gammaConversion_->clear();
    Elid_px_->clear();
    Elid_py_->clear();
    Elid_pz_->clear();
    for (int i = sk::LVPART; i <= sk::NVECP; i++)
    {
        // if (sk::KELID(1, i) != 0)
        // {
            Elid_partIdx_->push_back(i - 1);
            Elid_tag_->push_back(sk::KELID(1, i));
            Elid_gammaConversion_->push_back(sk::KELID(2, i));
            Elid_px_->push_back(sk::QELID(3, i));
            Elid_py_->push_back(sk::QELID(4, i));
            Elid_pz_->push_back(sk::QELID(5, i));
        // }
    }
}

void NanoAODWriter::defineBtag(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Btag_probNegIP", "B-tagging probability for tracks with negative impact parameter, hemisphere 1", Btag_probNegIP_);
    MakeField(model, "Btag_probPosIP", "B-tagging probability for tracks with negative impact parameter, hemisphere 2", Btag_probPosIP_);
    MakeField(model, "Btag_probAllIP", "B-tagging probability for tracks with negative impact parameter, whole event", Btag_probAllIP_);
    MakeField(model, "Btag_thrustVector", "B-tagging probability for tracks with positive impact parameter, hemisphere 1", Btag_thrustVector_);
}

void NanoAODWriter::fillBtag()
{
    fillVector(Btag_probNegIP_, 1, 3, [](int i)
               { return sk::QBTPRN(i); });
    fillVector(Btag_probPosIP_, 1, 3, [](int i)
               { return sk::QBTPRP(i); });
    fillVector(Btag_probAllIP_, 1, 3, [](int i)
               { return sk::QBTPRS(i); });

    *Btag_thrustVector_ = XYZVectorF(sk::QBTTHR(1), sk::QBTTHR(2), sk::QBTTHR(3));
}

void NanoAODWriter::defineHadid(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Haid_sign", "Used for combined tag", Haid_sign_);
    MakeField(model, "Haid_kaonDedx", "Kaon signature with DEDX", Haid_kaonDedx_);
    MakeField(model, "Haid_protonDedx", "Proton signature with DEDX", Haid_protonDedx_);
    MakeField(model, "Haid_kaonRich", "Kaon signature with RICH", Haid_kaonRich_);
    MakeField(model, "Haid_protonRich", "Proton signature with RICH", Haid_protonRich_);
    MakeField(model, "Haid_pionRich", "Pion signature with RICH", Haid_pionRich_);
    MakeField(model, "Haid_kaonCombined", "Kaon signature with combined tag", Haid_kaonCombined_);
    MakeField(model, "Haid_protonCombined", "Proton signature with combined tag", Haid_protonCombined_);
    MakeField(model, "Haid_richQuality", "RICH quality status", Haid_richQuality_);

    MakeField(model, "Haidn_pionTag", "Pion tag", Haidn_pionTag_);
    MakeField(model, "Haidn_kaonTag", "Kaon tag", Haidn_kaonTag_);
    MakeField(model, "Haidn_protonTag", "Proton tag", Haidn_protonTag_);
    MakeField(model, "Haidn_heavyTag", "Heavy particle tag", Haidn_heavyTag_);
    MakeField(model, "Haidn_pionTrackSelection", "Pion track selection", Haidn_pionTrackSelection_);
    MakeField(model, "Haidn_kaonTrackSelection", "Kaon track selection", Haidn_kaonTrackSelection_);
    MakeField(model, "Haidn_protonTrackSelection", "Proton track selection", Haidn_protonTrackSelection_);
    MakeField(model, "Haidn_heavyTrackSelection", "Heavy particle track selection", Haidn_heavyTrackSelection_);

    MakeField(model, "Haidr_pionTag", "Pion tag", Haidr_pionTag_);
    MakeField(model, "Haidr_kaonTag", "Kaon tag", Haidr_kaonTag_);
    MakeField(model, "Haidr_protonTag", "Proton tag", Haidr_protonTag_);
    MakeField(model, "Haidr_heavyTag", "Heavy particle tag", Haidr_heavyTag_);
    MakeField(model, "Haidr_electronTag", "Electron tag", Haidr_electronTag_);
    MakeField(model, "Haidr_selectionFlag", "Selection flag", Haidr_selectionFlag_);

    MakeField(model, "Haide_pionTag", "Pion tag", Haide_pionTag_);
    MakeField(model, "Haide_kaonTag", "Kaon tag", Haide_kaonTag_);
    MakeField(model, "Haide_protonTag", "Proton tag", Haide_protonTag_);
    MakeField(model, "Haide_heavyTag", "Heavy particle tag", Haide_heavyTag_);
    MakeField(model, "Haide_electronTag", "Electron tag", Haide_electronTag_);
    MakeField(model, "Haide_selectionFlag", "Selection flag", Haide_selectionFlag_);

    MakeField(model, "Haidc_pionTag", "Pion tag", Haidc_pionTag_);
    MakeField(model, "Haidc_kaonTag", "Kaon tag", Haidc_kaonTag_);
    MakeField(model, "Haidc_protonTag", "Proton tag", Haidc_protonTag_);
    MakeField(model, "Haidc_heavyTag", "Heavy particle tag", Haidc_heavyTag_);
    MakeField(model, "Haidc_electronTag", "Electron tag", Haidc_electronTag_);
    MakeField(model, "Haidc_selectionFlag", "Selection flag", Haidc_selectionFlag_);

    MakeField(model, "Dedx_value", "Dedx value (1 for mips)", Dedx_value_);
    MakeField(model, "Dedx_width", "Dedx Landau width", Dedx_width_);
    MakeField(model, "Dedx_nrWires", "Dedx number of wires", Dedx_nrWires_);
    MakeField(model, "Dedx_gapWires", "Mean distance between wires", Dedx_gapWires_);
    MakeField(model, "Dedx_error", "Dedx value error", Dedx_error_);
    MakeField(model, "Dedx_valueVD", "Dedx VD value", Dedx_valueVD_);
    MakeField(model, "Dedx_nrVDHits", "Number of VD hits", Dedx_nrVDHits_);

    MakeField(model, "Rich_theg", "Cherenkov angle gas radiator", Rich_theg_);
    MakeField(model, "Rich_sigg", "Sigma of Cherenkov angle gas radiator", Rich_sigg_);
    MakeField(model, "Rich_nphg", "Observed number of photons gas radiator", Rich_nphg_);
    MakeField(model, "Rich_nepg", "Expected number of photons gas radiator", Rich_nepg_);
    MakeField(model, "Rich_flagg", "Flag gas radiator", Rich_flagg_);
    MakeField(model, "Rich_thel", "Cherenkov angle liquid radiator", Rich_thel_);
    MakeField(model, "Rich_sigl", "Sigma of Cherenkov angle liquid radiator", Rich_sigl_);
    MakeField(model, "Rich_nphl", "Observed number of photons liquid radiator", Rich_nphl_);
    MakeField(model, "Rich_nepl", "Expected number of photons liquid radiator", Rich_nepl_);
    MakeField(model, "Rich_flagl", "Flag liquid radiator", Rich_flagl_);
}

void NanoAODWriter::fillHadid()
{
    Haid_sign_->clear();
    Haid_kaonDedx_->clear();
    Haid_protonDedx_->clear();
    Haid_kaonRich_->clear();
    Haid_protonRich_->clear();
    Haid_pionRich_->clear();
    Haid_kaonCombined_->clear();
    Haid_protonCombined_->clear();
    Haid_richQuality_->clear();

    Haidn_pionTag_->clear();
    Haidn_kaonTag_->clear();
    Haidn_protonTag_->clear();
    Haidn_heavyTag_->clear();
    Haidn_pionTrackSelection_->clear();
    Haidn_kaonTrackSelection_->clear();
    Haidn_protonTrackSelection_->clear();
    Haidn_heavyTrackSelection_->clear();

    Haidr_pionTag_->clear();
    Haidr_kaonTag_->clear();
    Haidr_protonTag_->clear();
    Haidr_heavyTag_->clear();
    Haidr_electronTag_->clear();
    Haidr_selectionFlag_->clear();

    Haide_pionTag_->clear();
    Haide_kaonTag_->clear();
    Haide_protonTag_->clear();
    Haide_heavyTag_->clear();
    Haide_electronTag_->clear();
    Haide_selectionFlag_->clear();

    Haidc_pionTag_->clear();
    Haidc_kaonTag_->clear();
    Haidc_protonTag_->clear();
    Haidc_heavyTag_->clear();
    Haidc_electronTag_->clear();
    Haidc_selectionFlag_->clear();

    Dedx_value_->clear();
    Dedx_width_->clear();
    Dedx_nrWires_->clear();
    Dedx_gapWires_->clear();
    Dedx_error_->clear();
    Dedx_valueVD_->clear();
    Dedx_nrVDHits_->clear();

    Rich_theg_->clear();
    Rich_sigg_->clear();
    Rich_nphg_->clear();
    Rich_nepg_->clear();
    Rich_flagg_->clear();
    Rich_thel_->clear();
    Rich_sigl_->clear();
    Rich_nphl_->clear();
    Rich_nepl_->clear();
    Rich_flagl_->clear();

    for (int i = sk::LVPART; i <= sk::NVECP; ++i)
    {
        // if (int(sk::VECP(7, i)) != 0)
        // {
            Haid_sign_->push_back(sk::KHAID(1, i));
            Haid_kaonDedx_->push_back(sk::KHAID(2, i));
            Haid_protonDedx_->push_back(sk::KHAID(3, i));
            Haid_kaonRich_->push_back(sk::KHAID(4, i));
            Haid_protonRich_->push_back(sk::KHAID(5, i));
            Haid_pionRich_->push_back(sk::KHAID(6, i));
            Haid_kaonCombined_->push_back(sk::QHAID(7, i));
            Haid_protonCombined_->push_back(sk::QHAID(8, i));
            Haid_richQuality_->push_back(sk::KHAID(9, i));

            Haidn_pionTag_->push_back(sk::KHAIDN(1, i));
            Haidn_kaonTag_->push_back(sk::KHAIDN(2, i));
            Haidn_protonTag_->push_back(sk::KHAIDN(3, i));
            Haidn_heavyTag_->push_back(sk::KHAIDN(4, i));
            Haidn_pionTrackSelection_->push_back(sk::KHAIDT(1, i));
            Haidn_kaonTrackSelection_->push_back(sk::KHAIDT(2, i));
            Haidn_protonTrackSelection_->push_back(sk::KHAIDT(3, i));
            Haidn_heavyTrackSelection_->push_back(sk::KHAIDT(4, i));

            Haidr_pionTag_->push_back(sk::KHAIDR(1, i));
            Haidr_kaonTag_->push_back(sk::KHAIDR(2, i));
            Haidr_protonTag_->push_back(sk::KHAIDR(3, i));
            Haidr_heavyTag_->push_back(sk::KHAIDR(4, i));
            Haidr_electronTag_->push_back(sk::KHAIDR(5, i));
            Haidr_selectionFlag_->push_back(sk::KHAIDR(6, i));

            Haide_pionTag_->push_back(sk::KHAIDE(1, i));
            Haide_kaonTag_->push_back(sk::KHAIDE(2, i));
            Haide_protonTag_->push_back(sk::KHAIDE(3, i));
            Haide_heavyTag_->push_back(sk::KHAIDE(4, i));
            Haide_electronTag_->push_back(sk::KHAIDE(5, i));
            Haide_selectionFlag_->push_back(sk::KHAIDE(6, i));

            Haidc_pionTag_->push_back(sk::KHAIDC(1, i));
            Haidc_kaonTag_->push_back(sk::KHAIDC(2, i));
            Haidc_protonTag_->push_back(sk::KHAIDC(3, i));
            Haidc_heavyTag_->push_back(sk::KHAIDC(4, i));
            Haidc_electronTag_->push_back(sk::KHAIDC(5, i));
            Haidc_selectionFlag_->push_back(sk::KHAIDC(6, i));

            Dedx_value_->push_back(sk::QDEDX(1, i));
            Dedx_width_->push_back(sk::QDEDX(2, i));
            Dedx_nrWires_->push_back(sk::KDEDX(3, i));
            Dedx_gapWires_->push_back(sk::QDEDX(4, i));
            Dedx_error_->push_back(sk::QDEDX(5, i));
            Dedx_valueVD_->push_back(sk::QDEDX(6, i));
            Dedx_nrVDHits_->push_back(sk::KDEDX(7, i));

            Rich_theg_->push_back(sk::THEG(i));
            Rich_sigg_->push_back(sk::SIGG(i));
            Rich_nphg_->push_back(sk::NPHG(i));
            Rich_nepg_->push_back(sk::NEPG(i));
            Rich_flagg_->push_back(sk::FLAGG(i));
            Rich_thel_->push_back(sk::THEL(i));
            Rich_sigl_->push_back(sk::SIGL(i));
            Rich_nphl_->push_back(sk::NPHL(i));
            Rich_nepl_->push_back(sk::NEPL(i));
            Rich_flagl_->push_back(sk::FLAGL(i));
        // }
    }
}

void NanoAODWriter::fillPartLoop(particleData& pData,
                                 eventData& eData,
                                 DataKind cat) {
    int nParticle = 0;
    int nParticleHP = 0;
    int nChargedParticle = 0;
    int nChargedParticleHP = 0;
    TVector3 netP(0, 0, 0);
    TVector3 netChargedP(0, 0, 0);
    TVector3 netPGen(0, 0, 0);
    TVector3 netChargedPGen(0, 0, 0);
    std::vector<int> llp;
    int nSize = (cat == DataKind::data) ? *nPart_
        : (cat == DataKind::gen) ? *nGenPart_
        : (cat == DataKind::sim) ? *nSimPart_
        : 0;
    for (auto iSize = 0; iSize < nSize; ++iSize) {
        bool pass = false;
        // calculate charge
        float q = 0;
        if (cat == DataKind::gen) {
            int status = sk::KP(iSize + 1, 1);
            pass = (status == 1 || status == 4); // final state particle, decayed or not decayed
        } else if (cat == DataKind::sim) {
            int status = sk::KP(sk::ISTSH(iSize+1), 1);
            int particleID = sk::KP(sk::ISTSH(iSize+1), 2);
            int mom = sk::ISTVX(1, iSize+1);
            int dau = sk::ISTVX(2, iSize+1);
            if (status == 1 || (status == 4 && dau == 0)) {
                pass = 1;
            } else if (status == 4) {
                llp.push_back(dau);
                pass = 0;
            } else if (particleID == 0 && std::find(llp.begin(), llp.end(), mom) != llp.end()) {
                pass = 1;
            } else {
                pass = 0;
            }
        } else {
            pass = !Part_lock_->at(iSize); // bad particle
        }
        if (pass) {
            // FORTRAN index starts with 1
            int i = iSize + 1;
            XYZTVectorF temp;
            if (cat == DataKind::gen) {
                temp = GenPart_fourMomentum_->at(iSize);
                int pdgid = sk::KP(i, 2);
                TParticlePDG* particle = pdgDatabase->GetParticle(pdgid);
                if (particle) {
                    q = particle->Charge() / 3.0;
                } else {
                    particle = pdgDatabase->GetParticle(-pdgid);
                    q = -particle->Charge() / 3.0;
                }
                pData.pid[nParticle] = pdgid;
                // for gen, use this var to store Lund status
                pData.pwflag[nParticle] = (sk::KP(i, 1) == 1) ? 1 : 4; 
                pData.highPurity[nParticle]= 1;
                pData.index[nParticle] = i;
                pData.correspondenceIndex[nParticle] = sk::ILUST(i);

                nParticleHP++;
                if (abs(q) > 0.5) {
                    nChargedParticle++;
                    nChargedParticleHP++;
                }
            } else if (cat == DataKind::sim) {
                temp = SimPart_fourMomentum_->at(iSize);
                q = sk::VECP(7, sk::MTRACK+i);
                pData.pid[nParticle] = sk::KP(sk::ISTSH(iSize+1), 2);
                // store mass code. -- 2 for electron (only thing useful as of June 11 2025)
                pData.pwflag[nParticle] = sk::VECP(8, sk::MTRACK+i);
                pData.highPurity[nParticle]= 1;
                pData.index[nParticle] = i;
                pData.correspondenceIndex[nParticle] = sk::ISTPA(i);
                
                nParticleHP++;
                nChargedParticle++;
                nChargedParticleHP++;
            } else {
                q = Part_charge_->at(iSize);
                temp = Part_fourMomentum_->at(iSize);
                lock[nParticle] = Part_lock_->at(iSize);
                // LEPTON ID with standard MUID and ELID
                if (Part_charge_->at(iSize) == 0) {
                    pData.pwflag[nParticle] = 4;
                    // standard muon selection (3rd bit from the right)
                } else if (sk::KMUID(1, i) & (1 << 2)) {
                    pData.pwflag[nParticle] = 1;
                    // standard electron selection
                } else if (sk::KELID(1, i) >= 4) {
                    pData.pwflag[nParticle] = 2;
                    // loose conversion ele tag
                } else if (sk::KELID(2, i) >= 1) {
                    pData.pwflag[nParticle] = 3;
                } else {
                    pData.pwflag[nParticle] = 0;
                }
                // TODO: use vdHits?
                pData.ntpc[nParticle] = (pData.charge[nParticle]!=0)? 7: 0;
                pData.d0[nParticle] = sk::QTRAC(4, i);
                pData.z0[nParticle] = sk::QTRAC(5, i);
                // use this variable to store track length for DELPHI
                pData.weight[nParticle] = sk::QTRAC(24, i);
                
                // below we follow the same definition in eventSelection.h
                // TODO: Use DELPHI selections
                if (pData.pwflag[nParticle]<=2) {
                    pData.highPurity[nParticle]= pData.pwflag[nParticle]<=2 && temp.Pt() >= 0.2;
                } else if (pData.pwflag[nParticle]==4) {
                    pData.highPurity[nParticle]= pData.pwflag[nParticle]==4 && temp.Pt() >= 0.4;
                }
                
                pData.index[nParticle] = i;
                pData.correspondenceIndex[nParticle] = sk::IPAST(i);
                
                if(pData.pwflag[nParticle]<=2) {
                    nChargedParticle++;
                    if (pData.highPurity[nParticle]) nChargedParticleHP++;
                }
                if (pData.highPurity[nParticle]) {
                    nParticleHP++;
                }
            }
            
            netP -= TVector3(temp.x(), temp.y(), temp.z());
            if (abs(q) > 0.5) {
                netChargedP -= TVector3(temp.x(), temp.y(), temp.z());
            }

            pData.px[nParticle] = temp.x();
            pData.py[nParticle] = temp.y();
            pData.pz[nParticle] = temp.z();
            emf[nParticle] = sk::QEMF(8,i);
            hpc[nParticle] = sk::QHPC(8,i);
            hac[nParticle] = sk::QHAC(8,i);
            stic[nParticle] = sk::QSTIC(1,i);
            pData.pt[nParticle] = temp.Pt();
            pData.pmag[nParticle]   = temp.P();
            pData.rap[nParticle]    = temp.Rapidity();
            pData.eta[nParticle]    = temp.Eta();
            pData.theta[nParticle]  = temp.Theta();
            pData.phi[nParticle]    = temp.Phi();
            pData.mass[nParticle]   = temp.M();
            pData.charge[nParticle] = q;
            
            ++nParticle;
        }
    }
    pData.nParticle         = nParticle;
    eData.nChargedParticle  = nChargedParticle;
    eData.nParticleHP       = nParticleHP;
    eData.nChargedParticleHP= nChargedParticleHP;

    eData.missP = netP.Mag();
    eData.missPt = netP.Perp();
    eData.missTheta = netP.Theta();
    eData.missPhi = netP.Phi();

    eData.missChargedP = netChargedP.Mag();
    eData.missChargedPt = netChargedP.Perp();
    eData.missChargedTheta = netChargedP.Theta();
    eData.missChargedPhi = netChargedP.Phi();
}

void NanoAODWriter::fillSelection(particleData& pData,
                                  eventData& eData) {
    TVector3 thrust             = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL);
    TVector3 thrustWithMissP    = getThrust(pData.nParticle, pData.px, pData.py, pData.pz, THRUST::OPTIMAL,false,false,NULL,true,pData.pwflag);

    eData.Thrust        = thrust.Mag();
    eData.TTheta        = thrust.Theta();
    eData.TPhi          = thrust.Phi();
    eData.ThrustWithMissP   = thrustWithMissP.Mag();
    eData.TThetaWithMissP   = thrustWithMissP.Theta();
    eData.TPhiWithMissP     = thrustWithMissP.Phi();
    if ( do_thrustMissP ) {
        setThrustVariables(&pData, &eData, TVector3(),
                           TVector3(), TVector3(), TVector3(), TVector3(),
                           thrustWithMissP);
    }

    Sphericity spher = Sphericity(pData.nParticle,
                        pData.px,
                        pData.py,
                        pData.pz,
                        pData.pwflag,
                        false);
    spher.setTree(&eData);

    eventSelection eSelection;
    eSelection.setEventSelection(&pData, &eData);

    eData.passesTotalChgEnergyMin = eSelection.getPassesTotalChgEnergyMin();
    eData.passesNeuNch = eSelection.getPassesNeuNch();
    eData.passesNTrkMin = eSelection.getPassesNTrkMin();
    eData.passesSTheta = eSelection.getPassesSTheta();
    eData.passesBELLE = eSelection.getPassesNeuNch() && 
      eSelection.getPassesSTheta() &&	
      eSelection.getPassesTotalChgEnergyMin() &&
      eSelection.getPassesNTrkMin();
    eData.passesISR = eSelection.getPassesISR();
    eData.passesWW = eSelection.getPassesWW();
    eData.Mvis = eSelection.getMvis();
    eData.sPrime = eSelection.getsPrime();
    eData.d2 = eSelection.getd2();
    eData.cW = eSelection.getcW();
}


void NanoAODWriter::user99()
{
    // std::cout << "NanoAODWriter::user99: Finalising" << std::endl;
    file_->cd();
    out_t->Write();
    out_tgen->Write();
    out_tsim->Write();
    writer_.reset();
    file_->Close();
};
