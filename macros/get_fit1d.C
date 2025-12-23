// =============================================================================================
// SMEFT 1D Quadratic Fit Tool
// =============================================================================================
/*
 * Author: Josue Elizalde
 * Date  : 28.05.2025
 * 
 * Description:
 *   Performs quadratic fits of quantum observables as functions of single Wilson coefficients,
 *   using TGraphErrors. The observable vs coefficient is fitted as:
 * 
 *       O(C_i) = [0] + [1] * C_i + [2] * C_i^2
 *
 *   The SM point is identified at C_i = 0.
 * 
 * Input:
 *   - ROOT file with TTree named "tree_observables" containing EFT observable values.
 * 
 * Output:
 *   - Plain text file with fit results for each observable and operator.
 */

// ==================================  DEPENDENCIES =====================================
// ROOT core I/O and histogramming
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TF1.h>

// ROOT plotting 
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>

// Standard C++ libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>


// ================================== FUNCTION DECLARATIONS =====================================
void fill_graph1d(TGraphErrors* graph, const char* observable, const char* input_operator, const char* input_file);


// =================================================================================================================
//                                             M A I N    C O D E
// =================================================================================================================

void get_fit1d() {

    // Data root file w observables
    const char* input_file = "../results/observables_1500GeV_costheta_02.root";

    // Open output file
    std::string output_path = "../jupyter/1d_fits_1500GeV_costheta_02.txt";
    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        std::cerr << "[ERROR] Cannot open output file: " << output_path << std::endl;
        return;
    }

    // List of Wilson coefficient names (as extracted from branch names)
    std::vector<std::string> wilson_names = {
        "tGRe", "tGIm",
        "td1", "tu1", "tj1", "Qu1", "Qd1", "Qj11", "Qj31",
        "td8", "tu8", "tj8", "Qu8", "Qd8", "Qj18", "Qj38"
    };

    // List of quantum observables to be analyzed 
    std::vector<std::string> observable_names = {
        "number_events", 
        "Ckk", "Crr", "Cnn", "Ckr", 
        "D3", "A_plus"
    };

    // Write header of output file .txt
    outfile << "observable,operator_x,min_x,max_x,"
            << "p1,p1_err,p2,p2_err,p0,p0_err,"
            << "sm_val,sm_val_error\n";

    // FITTING LOOP OVER ALL OPERATOR AND OBSERVABLES
    for (const auto& wilson_operator : wilson_names) {
        std::cout << "[INFO] Operator: " << wilson_operator << std::endl;

        for (const auto& observable : observable_names) {

            // Build graph from tree data
            auto* graph = new TGraphErrors();
            fill_graph1d(graph, observable.c_str(), wilson_operator.c_str(), input_file);

            // Fit with quadratic function
            TF1* fit = new TF1("fit_func", "[0] + [1]*x + [2]*x*x", -4.5, 4.5);
            graph->Fit(fit, "WQ");

            // Extract SM point (always first)
            double sm_x, sm_y;
            graph->GetPoint(0, sm_x, sm_y);
            double sm_y_error = graph->GetErrorY(0);

            // Get min and max of x (Wilson values)
            double min_x = std::numeric_limits<double>::max();
            double max_x = -std::numeric_limits<double>::max();
            for (int ip = 0; ip < graph->GetN(); ++ip) {
                double x, y;
                graph->GetPoint(ip, x, y);
                if (x < min_x) min_x = x;
                if (x > max_x) max_x = x;
            }

            // Fit parameters and errors
            double p0 = fit->GetParameter(0);
            double p1 = fit->GetParameter(1);
            double p2 = fit->GetParameter(2);
            double e0 = fit->GetParError(0);
            double e1 = fit->GetParError(1);
            double e2 = fit->GetParError(2);

            // Save to output file
            outfile << observable << "," << wilson_operator << ",";
            outfile << min_x << "," << max_x << ",";
            outfile << p1 << "," << e1 << ",";
            outfile << p2 << "," << e2 << ",";
            outfile << p0 << "," << e0 << ",";
            outfile << sm_y << "," << sm_y_error << "\n";

            // Optional: print SM deviation
            std::cout << observable << ":  Δ(SM - fit(0)) = " << sm_y - p0 << std::endl;

            delete graph;
            delete fit;
        }
    }


    outfile.close();
    std::cout << "\n[INFO] Fit results written to " << output_path << std::endl;
}

// ================================== AUXILIARY FUNCTION ========================================
// =============================================================================================
// PURPOSE : Fills a TGraphErrors with observable values vs Wilson coefficient (1D scan)
// INPUTS  :
//   - graph: pointer to an empty TGraphErrors object
//   - observable: name of the observable (e.g., "Ckk", "D3", "A_plus")
//   - input_operator: Wilson coefficient name (e.g., "tGRe", "td1", etc.)
//   - input_file: path to the ROOT file with the TTree "tree_observables"
// BEHAVIOR:
//   - Searches the TTree for entries where only the first operator is varied,
//     and the second operator is empty (1D scan).
//   - Extracts the observable value and its uncertainty for each point.
//   - Identifies the Standard Model point ("weight_mc_NOSYS") and inserts it at x = 0.
//
// NOTE:
//   - Assumes all observables have a corresponding "_err" branch for uncertainties.
// =============================================================================================

void fill_graph1d(TGraphErrors* graph, const char* observable, const char* input_operator, const char* input_file) {

    // Open the ROOT file
    TFile* file = TFile::Open(input_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Cannot open file: " << input_file << std::endl;
        return;
    }

    // Retrieve the TTree containing the observables
    TTree* tree = (TTree*)file->Get("tree_observables");
    if (!tree) {
        std::cerr << "[ERROR] TTree 'tree_observables' not found.\n";
        file->Close();
        return;
    }

    // Set up tree branches for metadata
    std::string* weight_name = nullptr;
    std::string* operator1 = nullptr;
    std::string* operator2 = nullptr;
    double wilson1, wilson2;

    tree->SetBranchAddress("weight_name", &weight_name);
    tree->SetBranchAddress("wilson1_value", &wilson1);
    tree->SetBranchAddress("wilson1_operator", &operator1);
    tree->SetBranchAddress("wilson2_value", &wilson2);
    tree->SetBranchAddress("wilson2_operator", &operator2);

    // Prepare branches for observable value and error
    double observable_value = 0.0;
    std::string err_branch = std::string(observable) + "_err";
    double observable_err = 0.0;

    if (tree->SetBranchAddress(observable, &observable_value) != 0 ||
        tree->SetBranchAddress(err_branch.c_str(), &observable_err) != 0) {
        std::cerr << "[ERROR] Observable or error branch not found: " << observable << std::endl;
        file->Close();
        return;
    }

    // Search for the Standard Model point (used as reference at x = 0)
    double sm_val = 0.0;
    double sm_err = 0.0;
    bool found_sm = false;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") {
            sm_val = observable_value;
            sm_err = observable_err;
            found_sm = true;
            break;
        }
    }

    if (!found_sm) {
        std::cerr << "[WARNING] SM point not found for observable '" << observable << "'\n";
        file->Close();
        return;
    }

    // Set SM point at x = 0
    graph->SetPoint(0, 0.0, sm_val);
    graph->SetPointError(0, 0.0, sm_err);

    // Fill graph with EFT points for the selected Wilson operator (1D scan only)
    int index = 1;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Include only 1D scans: first operator matches input and second is empty
        if (*operator1 == input_operator && operator2->empty()) {
            graph->SetPoint(index, wilson1, observable_value);
            graph->SetPointError(index, 0.0, observable_err);
            ++index;
        }
    }

    // Close the file to free memory
    file->Close();
}