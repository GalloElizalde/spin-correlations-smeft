

// ROOT core I/O and histogramming
#include <TFile.h>
#include <TTree.h>
#include <TGraph2DErrors.h>
#include <TF2.h>

// ROOT plotting (optional)
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>

// Standard C++ libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>


// ================================== FUNCTION DECLARATION =====================================
void fill_graph2d(TGraph2DErrors* g2, const char* observable, const char* input_operator1, 
                  const char* input_operator2, const char* input_file, int& n_points, 
                  int& match_order, double norm_factor);



// =================================================================================================================
//                                             M A I N    C O D E
// =================================================================================================================

void get_fit2d_cms_normalized() {

    // Path to ROOT file containing the observable values per SMEFT weight
    const char* input_file = "../results/observables_800GeV_costheta_04.root";

    // Path to output text file where fit coefficients will be saved
    std::string output_path = "../jupyter/2d_fits_800GeV_costheta_04_cms_normalized.txt";
    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        std::cerr << "[ERROR] Cannot open output file." << std::endl;
        return;
    }

    // =========================
    // Spin correlations (central values only)
    // =========================

    // Region: m(ttbar) > 800 GeV
    double Crr_highmass   = -0.13309;
    double Cnn_highmass   =  0.16968;
    double Ckk_highmass   = -0.041208;

    // Region: m(ttbar) > 800 GeV, |cos(theta)| < 0.4
    double Crr_highmass_cuts = -0.69407;
    double Cnn_highmass_cuts =  0.56297;
    double Ckk_highmass_cuts = -0.50611;

    // =========================
    // Choose region
    // =========================
    bool region_with_cuts = true;  // <--- cambia a true para usar la región con cortes

    // Write header of CSV output
    outfile << "observable,operator_x,min_x,max_x,operator_y,min_y,max_y,"
            << "[1],[2],[3],[4],[5],[0]\n";

    // Wilson coefficients to be scanned (as appear in tree branches)
    std::vector<std::string> wilson_names = {
        "tGRe", "tGIm", 
        "td1", "tu1", "tj1", "Qu1", "Qd1", "Qj11", "Qj31",
        "td8", "tu8", "tj8", "Qu8", "Qd8", "Qj18", "Qj38"
    };

    // Observables to be fitted as functions of 2 Wilson coefficients
    std::vector<std::string> observable_names = {
        "Ckk", "Crr", "Cnn",
    };

    // FITTING LOOP OVER ALL OPERATOR PAIRS AND OBSERVABLES
    for (size_t i = 0; i < wilson_names.size(); ++i) {
        for (size_t j = i + 1; j < wilson_names.size(); ++j) {

            std::string operator1 = wilson_names[i];
            std::string operator2 = wilson_names[j];
            for (const auto& observable : observable_names) {

                // Select normalization factor depending on region
                double norm_factor = 1.0;
                if (!region_with_cuts) {
                    if (observable == "Crr") norm_factor = Crr_highmass;
                    else if (observable == "Cnn") norm_factor = Cnn_highmass;
                    else if (observable == "Ckk") norm_factor = Ckk_highmass;
                } else {
                    if (observable == "Crr") norm_factor = Crr_highmass_cuts;
                    else if (observable == "Cnn") norm_factor = Cnn_highmass_cuts;
                    else if (observable == "Ckk") norm_factor = Ckk_highmass_cuts;
                }

                // Fill graph with observable values over 2D grid of Wilson coefficients
                TGraph2DErrors* g2 = new TGraph2DErrors();
                int count, match_order;
                fill_graph2d(g2, observable.c_str(), operator1.c_str(), operator2.c_str(),
                             input_file, count, match_order, norm_factor);

                // Require full grid: SM + 12 EFT points
                if (count != 13) {
                    delete g2;
                    continue;
                }

                // Enforce a consistent operator order
                std::string x_operator = match_order ? operator1 : operator2;
                std::string y_operator = match_order ? operator2 : operator1;

                std::cout << x_operator << " " << y_operator << std::endl;

                // Determine min/max scan ranges for fit
                double min_x = std::numeric_limits<double>::max();
                double max_x = -std::numeric_limits<double>::max();
                double min_y = std::numeric_limits<double>::max();
                double max_y = -std::numeric_limits<double>::max();

                for (int ip = 0; ip < g2->GetN(); ++ip) {
                    double x, y, z;
                    g2->GetPoint(ip, x, y, z);
                    if (x < min_x) min_x = x;
                    if (x > max_x) max_x = x;
                    if (y < min_y) min_y = y;
                    if (y > max_y) max_y = y;
                }

                // Perform quadratic surface fit
                TF2* fit = new TF2("fit", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y", 
                                   min_x, max_x, min_y, max_y);
                g2->Fit(fit, "WQ");

                // Store coefficients
                std::vector<double> p(6);
                for (int k = 0; k < 6; ++k)
                    p[k] = fit->GetParameter(k);

                // Write output line to text file
                outfile << observable << "," << x_operator << "," << min_x << "," << max_x << ","
                        << y_operator << "," << min_y << "," << max_y << ","
                        << p[1] << "," << p[2] << "," << p[3] << "," << p[4] << "," << p[5] << "," << p[0] << "\n";

                delete fit;
                delete g2;
            }
        }
    }

    // FINALIZE
    outfile.close();
    std::cout << "\n[INFO] Fit results written to '" << output_path << "'\n";
}



// =============================== AUXILIARY FUNCTION ===========================================

void fill_graph2d(TGraph2DErrors* g2, const char* observable, const char* input_operator1, 
                  const char* input_operator2, const char* input_file, int& n_points, 
                  int& match_order, double norm_factor) {

    TFile* file = TFile::Open(input_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open file: " << input_file << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("tree_observables");
    if (!tree) {
        std::cerr << "[ERROR] TTree 'tree_observables' not found.\n";
        file->Close();
        return;
    }

    std::string* weight_name = nullptr;
    std::string* operator1 = nullptr;
    std::string* operator2 = nullptr;
    double wilson1, wilson2;

    tree->SetBranchAddress("weight_name", &weight_name);
    tree->SetBranchAddress("wilson1_operator", &operator1);
    tree->SetBranchAddress("wilson1_value", &wilson1);
    tree->SetBranchAddress("wilson2_operator", &operator2);
    tree->SetBranchAddress("wilson2_value", &wilson2);

    double obs_val;
    double obs_err;
    std::string err_name = std::string(observable) + "_err";

    if (tree->SetBranchAddress(observable, &obs_val) != 0 ||
        tree->SetBranchAddress(err_name.c_str(), &obs_err) != 0) {
        std::cerr << "[ERROR] Observable or error branch not found: " << observable << std::endl;
        file->Close();
        return;
    }

    // Get SM point
    double sm_val = 0.0, sm_err = 0.0;
    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") {
            sm_val = obs_val;
            sm_err = obs_err;
            break;
        }
    }

    // Insert SM point (0,0) with normalization
    g2->SetPoint(0, 0.0, 0.0, norm_factor);
    g2->SetPointError(0, 0.0, 0.0, sm_err / sm_val * norm_factor);
    int index = 1;

    // Determine order of operator match
    match_order = false;
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        if (*operator1 == input_operator1 && *operator2 == input_operator2) {
            match_order = true;
            break;
        }
    }

    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") continue;

        bool match1 = (*operator1 == input_operator1 && operator2->empty());
        bool match2 = (*operator1 == input_operator2 && operator2->empty());
        bool match3 = (*operator1 == input_operator1 && *operator2 == input_operator2);
        bool match4 = (*operator1 == input_operator2 && *operator2 == input_operator1);

        if (match_order) {
            if (match1 || match3) {
                g2->SetPoint(index, wilson1, wilson2, (obs_val/sm_val) * norm_factor);
                g2->SetPointError(index, 0.0, 0.0, (obs_err/sm_val) * norm_factor);
                ++index;
            } else if (match2) {
                g2->SetPoint(index, wilson2, wilson1, (obs_val/sm_val) * norm_factor);
                g2->SetPointError(index, 0.0, 0.0, (obs_err/sm_val) * norm_factor);
                ++index;
            }
        } else {
            if (match1) {
                g2->SetPoint(index, wilson2, wilson1, (obs_val/sm_val) * norm_factor);
                g2->SetPointError(index, 0.0, 0.0, (obs_err/sm_val) * norm_factor);
                ++index;
            } else if (match2 || match4) {
                g2->SetPoint(index, wilson1, wilson2, (obs_val/sm_val) * norm_factor);
                g2->SetPointError(index, 0.0, 0.0, (obs_err/sm_val) * norm_factor);
                ++index;
            }
        }
    }

    n_points = index;
    file->Close();
}
