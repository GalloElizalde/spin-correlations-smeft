// =============================================================================================
// SMEFT 2D Quadratic Fit Tool (with errors)
// =============================================================================================
/*
 * Author: Josue Elizalde
 * Date  : 03.07.2025
 * 
 * Description:
 *   Performs 2D quadratic fits of quantum observables as functions of two Wilson coefficients (WC),
 *   using TGraph2DErrors. The observable is modeled as:
 * 
 *       O(C_i, C_j) = [0] + [1] * C_i + [2] * C_j + [3] * C_i^2 + [4] * C_j^2 + [5] * C_i * C_j
 * 
 *   The SM point is identified and set at (C_i = 0, C_j = 0).
 * 
 * Input:
 *   - ROOT file with TTree named "tree_observables" containing EFT observable values.
 *   - List of Wilson coefficients to scan.
 *   - List of observables to be fitted.
 * 
 * Output:
 *   - Plain text file with fit coefficients *and their errors* for each observable and operator pair.
 *     Format: observable,operator_x,min_x,max_x,operator_y,min_y,max_y,
 *             [1],[1]_err,[2],[2]_err,[3],[3]_err,[4],[4]_err,[5],[5]_err,[0],[0]_err,
 *             sm_val,sm_val_error
 */

// ================================== INCLUDE DEPENDENCIES =====================================

#include <TFile.h>
#include <TTree.h>
#include <TGraph2DErrors.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>

// ================================== FUNCTION DECLARATION =====================================
void fill_graph2d(TGraph2DErrors* g2, const char* observable, const char* input_operator1, 
                  const char* input_operator2, const char* input_file, int& n_points, int& match_order);


// =================================================================================================================
//                                             M A I N    C O D E
// =================================================================================================================

void get_fit2d() {

    const char* input_file = "../results/observables_1500GeV_costheta_02.root";
    std::string output_path = "../jupyter/2d_fits_1500GeV_costheta_02.txt";
    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        std::cerr << "[ERROR] Cannot open output file." << std::endl;
        return;
    }

    // Header with errors included
    outfile << "observable,operator_x,min_x,max_x,operator_y,min_y,max_y,"
            << "[1],[1]_err,[2],[2]_err,[3],[3]_err,[4],[4]_err,[5],[5]_err,[0],[0]_err,"
            << "sm_val,sm_val_error\n";

    // WC available in the sample
    std::vector<std::string> wilson_names = {
        "tGRe", "tGIm", 
        "td1", "tu1", "tj1", "Qu1", "Qd1", "Qj11", "Qj31",
        "td8", "tu8", "tj8", "Qu8", "Qd8", "Qj18", "Qj38"
    };

    // Observables studied
    std::vector<std::string> observable_names = {
        "number_events", 
        "Ckk", "Crr", "Cnn", "Ckr", 
        "D3", "A_plus"
    };

    // Loop over all operator pairs
    for (size_t i = 0; i < wilson_names.size(); ++i) {
        for (size_t j = i + 1; j < wilson_names.size(); ++j) {

            std::string operator1 = wilson_names[i];
            std::string operator2 = wilson_names[j];
            for (const auto& observable : observable_names) {

                TGraph2DErrors* g2 = new TGraph2DErrors();
                int count, match_order;
                fill_graph2d(g2, observable.c_str(), operator1.c_str(), operator2.c_str(),
                             input_file, count, match_order);

                if (count != 13) { // require full 2D grid
                    delete g2;
                    continue;
                }

                std::string x_operator = match_order ? operator1 : operator2;
                std::string y_operator = match_order ? operator2 : operator1;

                double sm_x, sm_y, sm_z;
                g2->GetPoint(0, sm_x, sm_y, sm_z);
                double sm_z_error = g2->GetErrorZ(0);

                // Scan range
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

                // Quadratic 2D fit
                TF2* fit = new TF2("fit", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y",
                                   min_x, max_x, min_y, max_y);
                g2->Fit(fit, "WQ");

                // Extract parameters and errors
                std::vector<double> p(6), e(6);
                for (int k = 0; k < 6; ++k) {
                    p[k] = fit->GetParameter(k);
                    e[k] = fit->GetParError(k);
                }

                // Save to CSV
                outfile << observable << "," << x_operator << "," << min_x << "," << max_x << ","
                        << y_operator << "," << min_y << "," << max_y << ","
                        << p[1] << "," << e[1] << ","
                        << p[2] << "," << e[2] << ","
                        << p[3] << "," << e[3] << ","
                        << p[4] << "," << e[4] << ","
                        << p[5] << "," << e[5] << ","
                        << p[0] << "," << e[0] << ","
                        << sm_z << "," << sm_z_error << "\n";

                std::cout << observable << ": Δ(SM - fit(0,0)) = " << sm_z - p[0] << std::endl;

                delete fit;
                delete g2;
            }
        }
    }

    outfile.close();
    std::cout << "\n[INFO] Fit results written to '" << output_path << "'\n";
}




// =============================== AUXILIARY FUNCTION ===========================================

void fill_graph2d(TGraph2DErrors* g2, const char* observable, const char* input_operator1, 
                  const char* input_operator2, const char* input_file, int& n_points, int& match_order) {

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
    double sm_val, sm_err;
    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") {
            sm_val = obs_val;
            sm_err = obs_err;
            break;
        }
    }

    g2->SetPoint(0, 0.0, 0.0, sm_val);
    g2->SetPointError(0, 0.0, 0.0, sm_err);
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
                g2->SetPoint(index, wilson1, wilson2, obs_val);
                g2->SetPointError(index, 0.0, 0.0, obs_err);
                ++index;
            } else if (match2) {
                g2->SetPoint(index, wilson2, wilson1, obs_val);
                g2->SetPointError(index, 0.0, 0.0, obs_err);
                ++index;
            }
        } else {
            if (match1) {
                g2->SetPoint(index, wilson2, wilson1, obs_val);
                g2->SetPointError(index, 0.0, 0.0, obs_err);
                ++index;
            } else if (match2 || match4) {
                g2->SetPoint(index, wilson1, wilson2, obs_val);
                g2->SetPointError(index, 0.0, 0.0, obs_err);
                ++index;
            }
        }
    }

    n_points = index;
    file->Close();
}
