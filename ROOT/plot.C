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


// ===================================== INCLUDE DEPENDENCIES ===================================

// ROOT core I/O and histogramming
#include <TFile.h>            // For reading ROOT files
#include <TTree.h>            // To handle TTree structure
#include <TH1F.h>             // Included by default (not used)
#include <TCanvas.h>          // To create and save plots

// ROOT advanced graphics
#include <TGraph2DErrors.h>   // Included (not used here)
#include <TLegend.h>          // For plot legends
#include <TApplication.h>     // Needed when using ROOT interactively (not needed in batch)

#include <TGraphErrors.h>     // For 2D graphs with symmetric errors
#include <TMultiGraph.h>      // For combining multiple TGraphs
#include <TF1.h>              // For polynomial fits

// C++ standard library
#include <iostream>           // For standard input/output
#include <tuple>              // For future multiple return values
#include <regex>              // For parsing operator names (not used here)
#include <fstream>            // For I/O (not used here)


// ==================================== FUNCTION DECLARATION ====================================
void fill_graph2d(TGraphErrors* graph, const char* observable, const char* input_operator, const char* input_observables_file);


// ========================================== MAIN CODE =========================================
void plot() {

    // Input ROOT file containing observables and Wilson coefficient values
    const char* input_file = "../results/observables_800GeV_costheta_04.root";
    
    std::string add = "(m_{tt} > 800 GeV  |cos#theta|<0.4)";
    //const char* input_file = "../results/observables_800GeV_costheta_04.root";


    // List of Wilson coefficients to loop over (one plot per coefficient)
    //std::vector<std::string> wilson_names = {"tGRe", "tGIm",
        //"td1", "tu1", "tj1", "Qu1", "Qd1", "Qj11", "Qj31",
        //"td8", "tu8", "tj8", "Qu8", "Qd8", "Qj18", "Qj38"
    //};
    std::vector<std::string> wilson_names = {
        "Qj31", "Qu1", "Qj11" , "td1", "Qj18", "tj1"
    };
 
    // List of quantum observables to plot per Wilson coefficient
    std::vector<std::string> observable_names = {
        "number_events", 
        "Ckk", "Crr", "Cnn",
        "Ckr",
        "D3", "A_plus"
    };

    // Loop over Wilson coefficients
    for (size_t i = 0; i < wilson_names.size(); ++i) {
        std::string wilson_operator = wilson_names[i];
        std::cout << wilson_operator << "\n";

        // Create a canvas and a multigraph for plotting all observables on the same plot
        auto* canvas = new TCanvas(wilson_operator.c_str(), wilson_operator.c_str(), 800, 600);
        auto* multi_graph = new TMultiGraph();

        // Define a set of colors for different observables
        std::vector<int> colors = {
            kRed, kBlue, kGreen+2, kMagenta, kOrange, kViolet+1, kAzure+10, kCyan+2, kGray+2
        };

                // ====== Legend setup ======
         // Create a legend for the plot
        auto* legend = new TLegend(0.35, 0.66, 0.56, 0.90);
        legend->SetNColumns(2);
        legend->SetTextSize(0.045);
        legend->SetBorderSize(2);   // ← pone un borde fino
        bool added_sm_to_legend = false;  // Flag to only add "SM point" once


        // Loop over observables to draw each one on the same canvas
        for (size_t k = 0; k < observable_names.size(); ++k) {
            std::string observable = observable_names[k];

            // Create a graph for the observable
            auto* graph = new TGraphErrors();
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(1.4);
            graph->SetMarkerColor(colors[k % colors.size()]);
            graph->SetLineColor(colors[k % colors.size()]);
            graph->SetTitle(observable.c_str());

            // Fill graph with data from tree (observable normalized to SM)
            fill_graph2d(graph, observable.c_str(), wilson_operator.c_str(), input_file);

            // Scaling
            double scale = 1.0;
            //if (observable == "Ckk")      scale = -0.0562;
            //else if (observable == "Crr") scale = -0.1249;
            //else if (observable == "Cnn") scale =  0.1771;
            //else if (observable == "D3") scale = 0.119425;
            //else if (observable == "A_plus" ) scale = (M_PI / 16.0) * (-1) * (-0.0562 -0.1249);
            //else if (observable == "Ckr") scale = -0.0775786;


            for (int p = 0; p < graph->GetN(); ++p) {
                double x = graph->GetX()[p];
                double y = graph->GetY()[p];
                double ey = graph->GetErrorY(p);

                graph->SetPoint(p, x, y * scale);
                graph->SetPointError(p, 0.0, ey * std::abs(scale));
           }

            // Add to multigraph and legend
            multi_graph->Add(graph, "P");
            multi_graph->Add(graph, "P");

            // ---- legend label in math style (TLatex) ----
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
            else observable_math = observable;

            legend->AddEntry(graph, observable_math.c_str(), "p");
            // ---------------------------------------------

            legend->SetFillStyle(0);        // transparent background

            // Add the SM point as a star marker
            double x0 = graph->GetX()[0];  // x = 0
            double y0 = graph->GetY()[0];  // y = 1

            TGraph* sm_point = new TGraph(1);
            sm_point->SetPoint(0, x0, y0);
            sm_point->SetMarkerStyle(29);     // star
            sm_point->SetMarkerSize(1.5);
            sm_point->SetMarkerColor(kBlack);
            multi_graph->Add(sm_point, "P");

            // Add SM point to legend (only once)
            if (!added_sm_to_legend) {
                legend->AddEntry(sm_point, "SM", "p");  // math-style label for SM
                added_sm_to_legend = true;
            }

            // Perform a 2nd-degree polynomial fit: [0] + [1]x + [2]x²
            TF1* fit = new TF1(("fit_" + observable).c_str(), "[0] + [1]*x + [2]*x*x", -4, 4);
            fit->SetLineColor(colors[k % colors.size()]);
            fit->SetLineWidth(3);
            graph->Fit(fit, "RQ");  // R = use range; Q = quiet
            fit->Draw("same");
        }

        // Set plot title and axis labels
        //std::string add = " (m_{tt} < 400 GeV)";
        //std::string add = observable_math + "  (m_{tt} > 1500 GeV  |cos#theta|<0.2)";
        std::string title = wilson_operator + add + ";Wilson Coefficient Value;Observable/SM";
        multi_graph->SetTitle(title.c_str());

        // Draw all graphs and legend
        multi_graph->Draw("AP*");
        // ====== Axis style customization ======
        multi_graph->GetXaxis()->SetTitleSize(0.05);   // tamaño del título eje X
        multi_graph->GetXaxis()->SetLabelSize(0.05);  // tamaño de números eje X
        multi_graph->GetYaxis()->SetTitleSize(0.05);   // tamaño del título eje Y
        multi_graph->GetYaxis()->SetLabelSize(0.05);  // tamaño de números eje Y

        //multi_graph->GetXaxis()->SetTitleOffset(1.2);  // distancia título-eje X
        //multi_graph->GetYaxis()->SetTitleOffset(1.4);  // distancia título-eje Y

        legend->Draw();

        // Save plot to file
        std::string filename = "../plots/2d/800/wow/" + wilson_operator + "_800GeV_costheta0p4.pdf";
        canvas->SaveAs(filename.c_str());
    }

    std::cout << "\n[INFO] OK";
}


// ====================================== AUXILIARY FUNCTION ======================================
/**
 * Fill a TGraphErrors with normalized observable values and errors from a TTree.
 * 
 * Parameters:
 *   - graph               : pointer to the graph to fill
 *   - observable          : name of observable (e.g., "Ckk")
 *   - input_operator      : name of Wilson coefficient (e.g., "tu1")
 *   - input_observables_file : ROOT file containing the tree with observables
 */
void fill_graph2d(TGraphErrors* graph, const char* observable, const char* input_operator, const char* input_observables_file) {

    // Open ROOT file
    TFile* file = TFile::Open(input_observables_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << input_observables_file << std::endl;
        return;
    }

    // Get TTree
    TTree* tree = (TTree*)file->Get("tree_observables");
    if (!tree) {
        std::cerr << "Error: Tree 'tree_observables' not found!" << std::endl;
        return;
    }

    // Set up branches
    std::string* weight_name = nullptr;
    std::string* operator1 = nullptr;
    std::string* operator2 = nullptr;
    double wilson1, wilson2;
    tree->SetBranchAddress("weight_name", &weight_name);
    tree->SetBranchAddress("wilson1_value", &wilson1);
    tree->SetBranchAddress("wilson1_operator", &operator1);
    tree->SetBranchAddress("wilson2_value", &wilson2);
    tree->SetBranchAddress("wilson2_operator", &operator2);

    // Bind observable and its error
    double observable_value = 0.0;
    if (tree->SetBranchAddress(observable, &observable_value) != 0) {
        std::cerr << "Error: Observable '" << observable << "' not found.\n";
        file->Close();
        return;
    }

    double observable_err = 0.0;
    std::string err_name = std::string(observable) + "_err";
    if (tree->SetBranchAddress(err_name.c_str(), &observable_err) != 0) {
        std::cerr << "Error: Observable error branch '" << err_name << "' not found.\n";
        file->Close();
        return;
    }

    // Find SM reference point
    double observable_sm = 0.0;
    double observable_sm_err = 0.0;
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") {
            observable_sm = observable_value;
            observable_sm_err = observable_err;
            break;
        }
    }

    // Add SM point at origin (x=0, y=1)
    graph->SetPoint(0, 0.0, 1.0);
    graph->SetPointError(0, 0.0, std::abs(observable_sm_err / observable_sm));

    // Fill graph with normalized observable values for 1D scans
    int index = 1;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Only keep points where input_operator is turned on alone
        bool match = *operator1 == input_operator && operator2->empty();
        if (!match) continue;

        // Fill point (x = wilson1, y = observable / SM)
        graph->SetPoint(index, wilson1, observable_value / observable_sm);
        graph->SetPointError(index, 0.0, std::abs(observable_err / observable_sm));
        ++index;
    }

    file->Close();  // Clean up
}
