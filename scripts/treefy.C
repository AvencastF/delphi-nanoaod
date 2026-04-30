// see.C
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector3D.h>
#include <TSystem.h>
#include <TString.h>
#include <TRegexp.h>
#include <iostream>
#include <fstream>
#include <set>
#include <regex>
#include <map>
using ROOT::VecOps::RVec;
using Vec3 = ROOT::Math::XYZVector;

std::map<std::string, double> compute_thrust(const RVec<float>& px,
                                          const RVec<float>& py,
                                          const RVec<float>& pz,
                                          const RVec<int>& flag_valid,
                                          double eps,
                                          bool include_met);
ROOT::RDF::RNode define_branches(ROOT::RDF::RNode df, std::vector<std::string>& branches, bool is_data);
int categorize_tau_decay(const RVec<short>& decay_products_pdgId, int idx_start, int idx_end);
int categorize_event(const RVec<short>& truth_pdgId);

// enum for tau truth decay categories
// Event truth category = 10*tau_plust_cat + tau_minus_cat
enum TauDecayCategory {
  NotTau = 0,
  SinglePi = 1,
  Rho = 2,
  // Lep = 3,
  El = 3,
  Mu = 4,
  Others = 5,
};
 

void treefy(const char* infile = "test.root") {
    std::string infile_str(infile);
    std::string outfile = infile_str.substr(0, infile_str.find_last_of(".")) + "_ttree.root";
    int nevt = 0;

    ROOT::RDataFrame raw_df("Events", infile);
    ROOT::RDF::RNode df = raw_df;

    // Collect column names and types
    auto colNames = df.GetColumnNames();
    std::set<std::string> classes;
    std::regex class_pattern(R"(ROOT::[^\n]+)");

    // Inspect column types and collect needed classes
    bool is_data = true;
    for (auto&& col : colNames) {
      if (col.find("GenPart")!=std::string::npos) is_data = false;
      std::string colType = df.GetColumnType(col);
      std::cout << colType << " " << col << std::endl;
      std::smatch match;
      if (std::regex_search(colType, match, class_pattern)) {
          std::string cls = match.str(0);
          if (std::count(cls.begin(), cls.end(), '<') != std::count(cls.begin(), cls.end(), '>')) {
              while (!cls.empty() && cls.back() == '>') cls.pop_back();
          }
          classes.insert(cls);
      }
      if (colType.find("ROOT::VecOps::RVec<std::int8_t>") != std::string::npos){
        // redefine the column to be RVec<short>
        df = df.Redefine(col, "ROOT::VecOps::RVec<short> new_col(" + col + ".begin(), " + col + ".end()); return new_col;");
      }
    }

    // Define new branches
    std::vector<std::string> new_branches;
    df = define_branches(df, new_branches, is_data);
    for (const auto& nb : new_branches) {
      colNames.push_back(nb);
    }

    if (gSystem->AccessPathName("dict.so")) {
      std::cout << "Building dict" << std::endl;
      // Write dict.h
      std::ofstream fout("dict.h");
      fout << "#include \"ROOT/RVec.hxx\"\n";
      fout << "#include \"Math/Vector4D.h\"\n";
      fout << "#include \"Math/GenVector/DisplacementVector3D.h\"\n";
      fout << "#include \"Math/GenVector/Cartesian3D.h\"\n";
      fout << "#include \"Math/GenVector/PositionVector3D.h\"\n";
      fout << "#include \"Math/GenVector/VectorUtil.h\"\n";
      fout << "#include \"Math/SVector.h\"\n";
      fout << "#include \"Math/SMatrix.h\"\n";
      fout << "#ifdef __CLING__\n";
      for (const auto& cls : classes) {
        fout << "#pragma link C++ class " << cls << "+;\n";
      }
      fout << "#endif\n";
      fout.close();

      // Build dictionary
      gSystem->Exec("rootcling -f dict.cxx -c dict.h");
      gSystem->Exec("g++ -shared -fPIC $(root-config --cflags --libs) dict.cxx -o dict.so");
    }
    gSystem->Load("dict.so");

    // Select columns
    std::vector<std::string> filtered;
    // std::regex part_pattern("^Part_fourMomentum");
    std::vector<std::string> bad_types = {"Event_shortDstVersion"};

    for (auto&& col : colNames) {
        std::string colstr = std::string(col);
        bool bad = false;
        for (auto&& badpat : bad_types) {
          if (colstr.find(badpat) != std::string::npos) {
            bad = true;
            break;
          }
        }
        if (!bad) filtered.push_back(colstr);
    }


    // Snapshot
    if (nevt > 0) {
      df.Range(nevt).Snapshot("t", outfile, filtered);
    } else {
      df.Snapshot("t", outfile, filtered);
    }
}


// Define new branches
ROOT::RDF::RNode define_branches(ROOT::RDF::RNode df, std::vector<std::string>& branches, bool is_data) {
  // Define event category based on tau decay modes
  if (is_data){
    df = df.Define("event_category", "-1.");
  } else {
    df = df.Define("event_category", categorize_event, {"GenPart_pdgId"});
  }
  branches.push_back("event_category");
  // Select particle p > 92/2 * 0.07 GeV = 3.22 GeV, so pt^2 > 10.3684 && |cos(theta)| > 0.035 (boundary region of TPC)
  // TODO: use p_beam instead of 92/2?
  // df = df.Define("Part_isGood", "((Part_fourMomentum.fCoordinates.fX*Part_fourMomentum.fCoordinates.fX + Part_fourMomentum.fCoordinates.fY*Part_fourMomentum.fCoordinates.fY + Part_fourMomentum.fCoordinates.fZ*Part_fourMomentum.fCoordinates.fZ) > 10.3684) && (abs(Part_fourMomentum.fCoordinates.fZ / sqrt(Part_fourMomentum.fCoordinates.fX*Part_fourMomentum.fCoordinates.fX + Part_fourMomentum.fCoordinates.fY*Part_fourMomentum.fCoordinates.fY + Part_fourMomentum.fCoordinates.fZ*Part_fourMomentum.fCoordinates.fZ)) > 0.035)");
  df = df.Define("Part_isGood", "Part_lock == 0");
  df = df.Define("nGoodPart", "(int) Part_isGood[Part_isGood].size()");
  df = df.Define("Part_thrust_calc_flag", "(Part_charge != 0) & (Part_isGood)");
  branches.push_back("Part_isGood");
  branches.push_back("nGoodPart");

  // Define thrust
  df = df.Define("thrust_map",
                  [](const RVec<float>& px,
                     const RVec<float>& py,
                     const RVec<float>& pz,
                     const RVec<int>& flag_valid) {
                    return compute_thrust(px, py, pz, flag_valid, 1e-12, false);
                  },
                  {"Part_fourMomentum.fCoordinates.fX",
                   "Part_fourMomentum.fCoordinates.fY",
                   "Part_fourMomentum.fCoordinates.fZ",
                   "Part_thrust_calc_flag"}
                );

  df = df.Define("thrust_Mag", "thrust_map[\"thrust_Mag\"]");
  df = df.Define("thrust_x", "thrust_map[\"thrust_x\"]");
  df = df.Define("thrust_y", "thrust_map[\"thrust_y\"]");
  df = df.Define("thrust_z", "thrust_map[\"thrust_z\"]");
  branches.push_back("thrust_Mag");
  branches.push_back("thrust_x");
  branches.push_back("thrust_y");
  branches.push_back("thrust_z");

  return df;
}



std::map<std::string, double> compute_thrust(const RVec<float>& px,
                                          const RVec<float>& py,
                                          const RVec<float>& pz,
                                          const RVec<int>& flag_valid,
                                          double eps = 1e-12,
                                          bool include_met = false){
  
  std::map<std::string, double> result;
  result["thrust_Mag"] = 0.0;
  result["thrust_x"] = 0.0;
  result["thrust_y"] = 0.0;
  result["thrust_z"] = 0.0;

  // Convert to Vec3
  std::vector<Vec3> p3;
  p3.reserve(px.size());
  for (size_t k = 0; k < px.size(); ++k) {
    if (!flag_valid[k]) continue;
    p3.emplace_back(px[k], py[k], pz[k]);
  }

  if (include_met) {
    Vec3 sum(0,0,0);
    for (auto& p : p3) sum += p;
    p3.push_back(-sum);
  }

  const int n = (int)p3.size();
  if (n==0) return result;

  // Precompute |p_i| and sum |p_i|
  std::vector<double> p_norm(n);
  double p_sum = 0.0;
  for (int i = 0; i < n; ++i) {
    p_norm[i] = std::sqrt(p3[i].Mag2());
    p_sum += p_norm[i];
  }
  if (p_sum <= 0.0) return result;

  // List all "good" pairs (i<j) where |p_i x p_j| >= eps
  std::vector<int> i_idx, j_idx;
  std::vector<Vec3> cross_vec;
  i_idx.reserve(n*(n-1)/2);
  j_idx.reserve(n*(n-1)/2);
  cross_vec.reserve(n*(n-1)/2);

  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      Vec3 c = p3[i].Cross(p3[j]);
      double c_norm = std::sqrt(c.Mag2());
      if (c_norm >= eps) {
        i_idx.push_back(i);
        j_idx.push_back(j);
        cross_vec.push_back(c);
      }
    }
  }

  // Number of pairs
  const int Mp = (int)i_idx.size();
  if (Mp == 0) return result;

  // Sign variants for (+/- p_i +/- p_j)
  constexpr std::array<std::array<int,2>,4> SIGN_VARIANTS {{
    {{+1,+1}},
    {{+1,-1}},
    {{-1,+1}},
    {{-1,-1}}
  }};

  double best_pair2_global = -1.0;
  Vec3   best_vec_global(0,0,0);

  // Loop over surviving pairs
  for (int ell = 0; ell < Mp; ++ell) {
    const int i = i_idx[ell];
    const int j = j_idx[ell];
    const Vec3& c = cross_vec[ell];

    // Determine signs s_k = sign(p_k · c)
    std::vector<int> s(n, -1);
    for (int k = 0; k < n; ++k) {
      double d = p3[k].Dot(c);
      s[k] = (d > 0.0) ? +1 : -1;
    }

    // base_sum = Σ_{k != i,j} s_k p_k
    Vec3 base_sum(0,0,0);
    for (int k = 0; k < n; ++k) {
      if (k == i || k == j) continue;
      base_sum += (double)s[k] * p3[k];
    }

    // Try four (±p_i ±p_j) variants
    double best_pair2_local = -1.0;
    Vec3 best_vec_local(0,0,0);

    for (int m = 0; m < (int)SIGN_VARIANTS.size(); ++m) {
      int si = SIGN_VARIANTS[m][0];
      int sj = SIGN_VARIANTS[m][1];
      Vec3 v = base_sum + (double)si * p3[i] + (double)sj * p3[j];
      double mag2 = v.Mag2();
      if (mag2 > best_pair2_local) {
        best_pair2_local = mag2;
        best_vec_local = v;
      }
    }

    // Update global best
    if (best_pair2_local > best_pair2_global) {
      best_pair2_global = best_pair2_local;
      best_vec_global = best_vec_local;
    }
  }

  if (best_pair2_global <= 0.0) return result;

  double best_mag = std::sqrt(best_pair2_global);
  double T = best_mag / p_sum;

  result["thrust_Mag"] = T;
  result["thrust_x"] = best_vec_global.X() / best_mag * T;
  result["thrust_y"] = best_vec_global.Y() / best_mag * T;
  result["thrust_z"] = best_vec_global.Z() / best_mag * T;
  return result;
}


int categorize_tau_decay(const RVec<short>& decay_products_pdgId, int idx_start, int idx_end) {
  int n_charged_pions = 0;
  int n_neutral_pions = 0;
  // int n_leptons = 0;
  int n_el = 0;
  int n_mu = 0;
  int n_kaons = 0;

  for (int i = idx_start; i < idx_end; ++i) {
    int pdgId = decay_products_pdgId[i];
    if (pdgId == 211 || pdgId == -211) {
      n_charged_pions++;
    } else if (pdgId == 111) {
      n_neutral_pions++;
    } else if (abs(pdgId) == 11) {
      n_el++;
    } else if (abs(pdgId) == 13) {
      n_mu++;
    } else if (abs(pdgId) == 321 || abs(pdgId) == 311 || abs(pdgId) == 130 || abs(pdgId) == 310) {
      n_kaons++;
    }
  }
  if (n_charged_pions == 1 && n_neutral_pions == 0 && n_el == 0 && n_mu == 0 && n_kaons == 0) {
    return TauDecayCategory::SinglePi;
  } else if (n_charged_pions == 1 && n_neutral_pions == 1 && n_el == 0 && n_mu == 0 && n_kaons == 0) {
    return TauDecayCategory::Rho;
  } else if (n_el==1 && n_mu==0 && n_charged_pions==0 && n_neutral_pions==0 && n_kaons==0) {
    return TauDecayCategory::El;
  } else if (n_el==0 && n_mu==1 && n_charged_pions==0 && n_neutral_pions==0 && n_kaons==0) {
    return TauDecayCategory::Mu;
  } else {
    return TauDecayCategory::Others;
  }
}

int categorize_event(const RVec<short>& truth_pdgId) {
  int event_category = 0;
  // truth_pdgId is distributed as: x,x,...,15,-15,(22),16,Products of tau-,-16,Products of tau+
  // Find tau+ and tau-
  int num_tau_plus = 0, num_tau_minus = 0;
  int idx_start_tau_minus = -1, idx_end_tau_minus = -1, idx_start_tau_plus = -1, idx_end_tau_plus = -1;
  for (size_t i = 0; i < truth_pdgId.size(); ++i) {
    if (truth_pdgId[i] == 15) num_tau_minus++;
    else if (truth_pdgId[i] == -15) num_tau_plus++;

    // Find indices of tau decay products if both taus are found
    if (abs(truth_pdgId[i]) == 16 && (num_tau_minus>0) && (num_tau_plus>0)) {
      for (size_t j = i+1; j < truth_pdgId.size(); ++j) {
        if (abs(truth_pdgId[j]) == 16) {
          if (truth_pdgId[j] == 16) {
            idx_start_tau_plus = i;
            idx_end_tau_plus = j;
            idx_start_tau_minus = j;
            idx_end_tau_minus = truth_pdgId.size();
          } else if (truth_pdgId[j] == -16) {
            idx_start_tau_minus = i;
            idx_end_tau_minus = j;
            idx_start_tau_plus = j;
            idx_end_tau_plus = truth_pdgId.size();
          }
          break;
        }
      }
      break;
    }
  }

  if (num_tau_plus > 0 && num_tau_minus > 0 &&
      idx_start_tau_minus != -1 && idx_end_tau_minus != -1 &&
      idx_start_tau_plus != -1 && idx_end_tau_plus != -1) {
    int tau_minus_cat = categorize_tau_decay(truth_pdgId, idx_start_tau_minus, idx_end_tau_minus);
    int tau_plus_cat = categorize_tau_decay(truth_pdgId, idx_start_tau_plus, idx_end_tau_plus);
    event_category = 10 * tau_plus_cat + tau_minus_cat;
  } else {
    event_category = 0; // Not tau
  }
  return event_category;
}

int main(int argc, char** argv) {
  const char* infile = "test.root";
  // bool is_signal_MC = false;

  if (argc > 1) {
    infile = argv[1];
  }

  treefy(infile);

  return 0;
}