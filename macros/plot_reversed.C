// =============================================================================================
// SMEFT Observable Plotter (1D fits over Wilson coefficients)
// =============================================================================================
/*
 * Author: Josue Elizalde
 * Date  : 28.05.2025
 * 
 * Description:
 *   Plots quantum observables as functions of single Wilson coefficients,
 *   extracting data from a ROOT TTree and overlaying polynomial fits.
 * 
 *   Each observable is plotted against the Wilson coefficient using TGraphErrors,
 *   normalized to the SM prediction (O / O_SM), and fitted with a quadratic function:
 * 
 *       O(C) / O_SM = [0] + [1] * C + [2] * C²
 *
 * Inputs:
 *   - A ROOT file with a TTree named "tree_observables", containing:
 *       * observable values and uncertainties,
 *       * Wilson coefficients and their names,
 *       * weight metadata (e.g., weight_mc_NOSYS for the SM point).
 * 
 */


// ================================== INCLUDE DEPENDENCIES =====================================

// ROOT core I/O and histogramming
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

// ROOT plotting and fitting
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TApplication.h>

// Standard C++ libraries
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// ======================================= DECLARATIONS =========================================

void fill_graph2d(TGraphErrors* graph, const char* observable, const char* input_operator, const char* input_observables_file);


// ======================================= MAIN FUNCTION =========================================

void plot_reversed() {

    // Input file containing the SMEFT observable scan results
    const char* input_file = "../results/observables_800GeV_costheta_04.root";

    // List of Wilson coefficients (as labeled in the tree)
    std::vector<std::string> wilson_names = {
        "tGRe", "tGIm",
        "td1", "tu1", "tj1", "Qu1", "Qd1", "Qj11", "Qj31",
        "td8", "tu8", "tj8", "Qu8", "Qd8", "Qj18", "Qj38"
    };

    // List of quantum observables to plot
    std::vector<std::string> observable_names = {
        "number_events", 
        "Ckk", "Crr", "Cnn",
        "Ckr",
        "D3", "A_plus"
    };

    // Marker and color styles for each operator
    int marker_styles[] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
    std::vector<int> colors = {
        kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kViolet+1, kCyan+2, kGray+2,
        kSpring+4, kPink+6, kTeal+2, kAzure+10, kYellow+3, kBlue+3, kGreen+3, kOrange+6
    };

    // Loop over all observables
    for (const auto& observable : observable_names) {

        std::cout << "[INFO] Plotting observable: " << observable << std::endl;

        // Map observable to math-style label
        std::string observable_math;
        if (observable == "Ckk")          observable_math = "C_{kk}";
        else if (observable == "Crr")     observable_math = "C_{rr}";
        else if (observable == "Cnn")     observable_math = "C_{nn}";
        else if (observable == "Ckr")     observable_math = "C_{kr}";
        else if (observable == "Crk")     observable_math = "C_{rk}";
        else if (observable == "D3")      observable_math = "D_{3}";
        else if (observable == "A_plus")  observable_math = "A^{+}";
        else if (observable == "A_minus") observable_math = "A^{-}";
        else if (observable == "number_events") observable_math = "N_{events}";
        else observable_math = observable; // fallback

        // Create a canvas and multigraph
        TCanvas* canvas = new TCanvas(observable.c_str(), observable.c_str(), 1000, 700);
        TMultiGraph* multi_graph = new TMultiGraph();

        

        // Create a legend (2 columns, more compact)
        TLegend* legend = new TLegend(0.70, 0.68, 0.85, 0.88);  // ancho = 0.15
        legend->SetTextSize(0.04);
        legend->SetNColumns(2);
        legend->SetBorderSize(0);

        // Loop over each Wilson coefficient and add its graph
        for (size_t i = 0; i < wilson_names.size(); ++i) {
            const std::string& wilson = wilson_names[i];

            // Create graph for this operator
            TGraphErrors* graph = new TGraphErrors();
            graph->SetMarkerStyle(marker_styles[i % 14]);
            graph->SetMarkerSize(1.2);
            graph->SetMarkerColor(colors[i % colors.size()]);
            graph->SetLineColor(colors[i % colors.size()]);
            graph->SetLineStyle(1);
            graph->SetLineWidth(2);
            graph->SetTitle(wilson.c_str());

            // Fill graph with observable values
            fill_graph2d(graph, observable.c_str(), wilson.c_str(), input_file);

            // Add graph to multigraph and legend
            multi_graph->Add(graph, "P");
            legend->AddEntry(graph, wilson.c_str(), "lp");

            // Perform quadratic fit over Wilson coefficient values
            TF1* fit = new TF1(("fit_" + wilson).c_str(), "[0] + [1]*x + [2]*x*x", -4, 4);
            fit->SetLineColor(colors[i % colors.size()]);
            fit->SetLineStyle(1);
            fit->SetLineWidth(2);

            graph->Fit(fit, "RQ");
            fit->Draw("same");
        }

        // Add titles and axis labels (math style)
        //std::string title = observable_math + "  (m_{#it{t}#it{t}} < 400 GeV )";
        std::string title = observable_math + "  (m_{#it{t}#it{t}} > 800 GeV  |cos#theta|<0.4)";
        multi_graph->SetTitle((title + ";Wilson coefficient value;Observable/SM").c_str());
        multi_graph->Draw("AP");
        // Only draw legend for the first observable
        if (observable == observable_names.front()) {
            legend->Draw();
        }


        // Make axis titles and labels bigger
        multi_graph->GetXaxis()->SetTitleSize(0.047);  // bigger axis title
        multi_graph->GetYaxis()->SetTitleSize(0.047);
        multi_graph->GetXaxis()->SetLabelSize(0.047); // bigger axis numbers
        multi_graph->GetYaxis()->SetLabelSize(0.047);

        // Save plot as PDF
        std::string filename = "../plots/2d/1500_02/zzz" + observable + "_vs_wilsons_400GeV_last.pdf";
        canvas->SaveAs(filename.c_str());
    }

    std::cout << "\n[INFO] ✅ All plots created successfully with quadratic fits.\n";
}


// ===================================== FILLING FUNCTION ========================================

/*
 * Function: fill_graph2d
 * ----------------------
 *   Fills a TGraphErrors with the observable values normalized to SM, as a function
 *   of a single Wilson coefficient (1D scan). Includes the SM point at x = 0.
 *
 * Inputs:
 *   - graph               : pointer to TGraphErrors to be filled
 *   - observable          : name of the observable (e.g., "Ckk", "D3", ...)
 *   - input_operator      : name of the Wilson coefficient to scan (e.g., "tGRe")
 *   - input_observables_file : path to ROOT file with TTree "tree_observables"
 */
void fill_graph2d(TGraphErrors* graph, const char* observable, const char* input_operator, const char* input_observables_file) {

    double ckksm = -0.0562;	
    double crrsm = -0.1249;
    double cnnsm = 	0.1771;


    TFile* file = TFile::Open(input_observables_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << input_observables_file << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("tree_observables");
    if (!tree) {
        std::cerr << "Error: Tree 'tree_observables' not found!" << std::endl;
        file->Close();
        return;
    }

    // Bind branches for SMEFT metadata
    std::string* weight_name = nullptr;
    std::string* operator1 = nullptr;
    std::string* operator2 = nullptr;
    double wilson1, wilson2;

    tree->SetBranchAddress("weight_name", &weight_name);
    tree->SetBranchAddress("wilson1_value", &wilson1);
    tree->SetBranchAddress("wilson1_operator", &operator1);
    tree->SetBranchAddress("wilson2_value", &wilson2);
    tree->SetBranchAddress("wilson2_operator", &operator2);

    // Bind observable value and error branches
    double observable_value = 0.0;
    double observable_err = 0.0;
    std::string err_name = std::string(observable) + "_err";

    if (tree->SetBranchAddress(observable, &observable_value) != 0 ||
        tree->SetBranchAddress(err_name.c_str(), &observable_err) != 0) {
        std::cerr << "Error: Observable or error branch not found for " << observable << std::endl;
        file->Close();
        return;
    }

    // Find SM point (weight_mc_NOSYS)
    double observable_sm = -1;
    double observable_sm_err = -1;
    bool found_sm = false;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") {
            observable_sm = observable_value;
            observable_sm_err = observable_err;
            found_sm = true;
            break;
        }
    }

    if (!found_sm) {
        std::cerr << "[WARNING] SM point not found for observable: " << observable << std::endl;
        file->Close();
        return;
    }

    // Insert SM point at x = 0
    graph->SetPoint(0, 0.0, 1.0);  // SM normalized to 1
    graph->SetPointError(0, 0.0, std::abs( (observable_sm_err / observable_sm)));

    // Fill graph with points where only input_operator is active (1D scan)
    int index = 1;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (*operator1 == input_operator && operator2->empty()) {
            graph->SetPoint(index, wilson1, observable_value / observable_sm);
            graph->SetPointError(index, 0.0, std::abs(observable_err / observable_sm));
            std::cout << weight_name->c_str() << std::endl;
            ++index;
        }
    }

    file->Close();
}
