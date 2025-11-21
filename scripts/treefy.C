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
                                          const RVec<bool>& flag_valid,
                                          double eps,
                                          bool include_met);
ROOT::RDF::RNode define_branches(ROOT::RDF::RNode df, std::vector<std::string>& branches);

void treefy(const char* infile = "test.root",
            std::string outfile = "",
            int nevt = 0) {
    std::string infile_str(infile);
    if (outfile == "") {
      outfile = infile_str.substr(0, infile_str.find_last_of(".")) + "_ttree.root";
    }

    ROOT::RDataFrame raw_df("Events", infile);
    ROOT::RDF::RNode df = raw_df;

    // Collect column names and types
    auto colNames = df.GetColumnNames();
    std::set<std::string> classes;
    std::regex class_pattern(R"(ROOT::[^\n]+)");

    // Inspect column types and collect needed classes
    for (auto&& col : colNames) {
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
    df = define_branches(df, new_branches);
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

    // Filter tautau->pi pi nu nu events: 
    // no kaon, lambda and Xi0: no status=4 particles
    // exactly one pi+ and one pi- in final status particles
    // no neutral pions, no short-lived particles, no kaons, eta, omega, neutrinos other than nu_tau
    auto df_pipi = df.Filter("\
      GenPart_status[GenPart_status==4].size()==0 && \
      GenPart_pdgId[(GenPart_pdgId==211)&&(GenPart_status==1)].size()==1 && \
      GenPart_pdgId[(GenPart_pdgId==-211)&&(GenPart_status==1)].size()==1 && \
      GenPart_pdgId[\
         (abs(GenPart_pdgId)==111) || (abs(GenPart_pdgId)==321) || (abs(GenPart_pdgId)==221) || \
         (abs(GenPart_pdgId)==223) || (abs(GenPart_pdgId)==12) || (abs(GenPart_pdgId)==14) \
        ].size()==0 \
    ");
    df_pipi.Snapshot("t", outfile.substr(0, outfile.find_last_of(".")) + "_pipi.root", filtered);
}


// Define new branches
ROOT::RDF::RNode define_branches(ROOT::RDF::RNode df, std::vector<std::string>& branches) {
  // Select particle pt > 92/2 * 0.07 GeV = 3.22 GeV, so pt^2 > 10.3684
  // TODO: use p_beam instead of 92/2?
  // df = df.Define("Part_isGood", "(Part_fourMomentum_fCoordinates_fX*Part_fourMomentum_fCoordinates_fX + Part_fourMomentum_fCoordinates_fY*Part_fourMomentum_fCoordinates_fY) > 10.3684");
  df = df.Define("Part_isGood", "(Part_fourMomentum.fCoordinates.fX*Part_fourMomentum.fCoordinates.fX + Part_fourMomentum.fCoordinates.fY*Part_fourMomentum.fCoordinates.fY) > 10.3684");
  df = df.Define("nGoodPart", "Part_isGood[Part_isGood].size()");
  branches.push_back("Part_isGood");
  branches.push_back("nGoodPart");

  // Define thrust
  df = df.Define("thrust_map", "compute_thrust(Part_fourMomentum.fCoordinates.fX, Part_fourMomentum.fCoordinates.fY, Part_fourMomentum.fCoordinates.fZ, Part_isGood, 1e-12, false)");

  df = df.Define("thrust", "thrust_map[\"thrust\"]");
  df = df.Define("thrust_x", "thrust_map[\"thrust_x\"]");
  df = df.Define("thrust_y", "thrust_map[\"thrust_y\"]");
  df = df.Define("thrust_z", "thrust_map[\"thrust_z\"]");
  branches.push_back("thrust");
  branches.push_back("thrust_x");
  branches.push_back("thrust_y");
  branches.push_back("thrust_z");

  return df;
}



std::map<std::string, double> compute_thrust(const RVec<float>& px,
                                          const RVec<float>& py,
                                          const RVec<float>& pz,
                                          const RVec<bool>& flag_valid,
                                          double eps = 1e-12,
                                          bool include_met = false){
  
  std::map<std::string, double> result;
  result["thrust"] = 0.0;
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