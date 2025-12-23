// =============================================================================================
// FUNCTION TO EXTRACT OBSERVABLE RATIOS (EFT/SM) OVER 2D WILSON COEFFICIENT GRIDS FOR EACH OBSERVABLE
// =============================================================================================
/*
 * Author: [Josue Elizalde]
 * Date: [28.05.25]
 * 
 * Function: get_points3d
 * Description:
 *   For each pair of Wilson coefficients and each angular observable, this function:
 *     - Loads the EFT/SM ratio values from a ROOT file
 *     - Constructs a 3D graph 
 *     - Fits a quadratic surface in the Wilson coefficient space (TF2 fit)
 *     - Plots the resulting points and fits
 *     - Saves the plots to disk and writes the graphs and fits to a ROOT file
 * 
 * Operators analyzed:
 *   - Full list: {tGRe, tGIm, td1, td8, tu1, tu8, tj1, tj8, Qu1, Qu8, Qd1, Qd8, Qj11, Qj18, Qj31, Qj38}
 *   - Chromomagnetic dipole (tGRe, tGIm)
 *   - Four-fermion operators with color singlet/octet structure e.g. td8, tu1, Qj11, Qj31
 *   
 * Observables analyzed:
 *   - Angular spin correlations and asymmetries: {N, Ckk, Crr, Cnn, Ckr, Crk, D3, A⁺, A⁻}
 *   - Each observable is analyzed as a function of two Wilson coefficients at a time
 *
 * Use case:
 *   - Visualization of multidimensional SMEFT parameter space
 *   - Sensitivity maps and global EFT fit preparation
 */
// INCLUDE DEPENDENCIES
// ROOT core I/O and histogramming
#include <TFile.h>            // File input/output operations
#include <TTree.h>            // Access and manipulate ROOT TTrees
#include <TH1F.h>             // 1D histograms
#include <TCanvas.h>          // Canvas for plotting

// ROOT advanced data analysis
#include <TGraph2DErrors.h>   // 3D graph with asymmetric error bars
#include <TLegend.h>          // Legend for plotting multiple objects
#include <TApplication.h>     // Required when running ROOT interactively

// C++ standard libraries
#include <iostream>           // Standard console I/O
#include <tuple>              // For returning multiple values from a function
#include <regex>              // Regular expressions (used for parsing Wilson names)


// Function Declarations
void function_get_points3d(TGraph2DErrors* g2, const char* observable,  const char* input_operator1, const char* input_operator2, 
    const char* input_observables_file );



// =================================================================================================================
//                                             M A I N    C O D E
// =================================================================================================================

void plot3d() {

    // Input ROOT file containing observables computed per SMEFT weight
    const char* input_file = "../results/observables_800GeV_costheta_04.root";

    // List of Wilson coefficient names (as extracted from branch names)
    //std::vector<std::string> wilson_names = {
    //    "tGRe", "tGIm",
    //    "td1", "tu1", "tj1" , "Qu1", "Qd1", "Qj11", "Qj31",
    //    "td8", "tu8", "tj8" , "Qu8", "Qd8", "Qj18", "Qj38"
    //};

    std::vector<std::string> wilson_names = { "Qj31", "td8", "Qj18", "Qu8"
    };

        // List of quantum observables to be analyzed (original names used in ROOT file)
    std::vector<std::string> observable_names = {
        "number_events",
        "Ckk", "Crr", "Cnn", 
        "D3", "A_plus"
    };

    // Map from internal observable names to LaTeX labels for plots
    std::map<std::string, std::string> observable_labels = {
        {"number_events", "N_{events}"},
        {"Ckk",           "C_{kk}"},
        {"Crr",           "C_{rr}"},
        {"Cnn",           "C_{nn}"},
        {"D3",            "D_{3}"},
        {"A_plus",        "A^{+}"}
    };


    // Output ROOT file to store fits and graphs
    TFile* output_file = new TFile("../results/fits3d_julia_800GeV_costheta_04.txt", "RECREATE");

    // Set plot style (palette and no statistics box)
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);


    // Loop over all unique Wilson coefficient pairs (i < j)
    for (size_t i = 0; i < wilson_names.size(); ++i) {
        for (size_t j = i + 1; j < wilson_names.size(); ++j) {

            std::string operator1 = wilson_names[i];
            std::string operator2 = wilson_names[j];
            
            // Create a canvas divided into subplots (2x2 for 4 observables)
            std::string canvas_name = "c_multi_" + operator1 + "_" + operator2;
            TCanvas* c_multi = new TCanvas(canvas_name.c_str(), "", 1200, 1000);
            int n_obs = observable_names.size();
            int n_cols = std::ceil(std::sqrt(n_obs));
            int n_rows = std::ceil(double(n_obs) / n_cols);
            c_multi->Divide(n_cols, n_rows, 0.01, 0.01);


            // Loop over all observables
            for (size_t k = 0; k < observable_names.size(); ++k) {
                std::string observable = observable_names[k];
                c_multi->cd(k + 1);

                // Load EFT/SM observable ratios as a 2D graph with error bars
                TGraph2DErrors* g2 = new TGraph2DErrors();
                function_get_points3d(g2, observable.c_str(), operator1.c_str(), operator2.c_str(), input_file);

                g2->SetTitle("");
                g2->SetMarkerColor(kRed + 1);
                g2->SetMarkerStyle(20);
                g2->SetMarkerSize(1.2);

                // If no points, annotate and skip
                if (g2->GetN() == 0) {
                    std::cerr << "[WARN] No data for " << observable << " with " << operator1 << ", " << operator2 << std::endl;
                    TLatex* empty = new TLatex(0.3, 0.5, "No data available");
                    empty->SetNDC();
                    empty->SetTextSize(0.055);
                    empty->SetTextColor(kGray+2);
                    empty->Draw();
                    delete g2;
                    continue;
                }

                // Define axis labels and title
                std::string cut_label = "";  //"(m_{tt} > 800 GeV, |cos#theta| < 0.4)";
                //std::string title = observable + "_{EFT}/" + observable + "_{SM} " + cut_label;
                std::string title = observable_labels[observable] + " " + cut_label;
                std::string x_label = operator1;
                std::string y_label = operator2;
                //std::string z_label = observable + "_{EFT}/" + observable + "_{SM}";
                std::string z_label = observable;

                // Compute min/max for axis ranges
                double min_x = std::numeric_limits<double>::max();
                double max_x = -std::numeric_limits<double>::max();
                double min_y = std::numeric_limits<double>::max();
                double max_y = -std::numeric_limits<double>::max();
                double x, y, z;

                for (int ip = 0; ip < g2->GetN(); ++ip) {
                    g2->GetPoint(ip, x, y, z);
                    if (x < min_x) min_x = x;
                    if (x > max_x) max_x = x;
                    if (y < min_y) min_y = y;
                    if (y > max_y) max_y = y;
                }

                double min_z = std::numeric_limits<double>::max();
                double max_z = -std::numeric_limits<double>::max();

                for (int ip = 0; ip < g2->GetN(); ++ip) {
                    g2->GetPoint(ip, x, y, z);
                    if (z < min_z) min_z = z;
                    if (z > max_z) max_z = z;
                }

                // Agrega margen para mejor visualización (opcional)
                double margin_z = 0.1 * (max_z - min_z);
                min_z -= margin_z;
                max_z += margin_z;


                // Add a margin to axis ranges for better plot coverage
                double margin_factor = 0.3;
                double margin_x = margin_factor * (max_x - min_x);
                double margin_y = margin_factor * (max_y - min_y);
                min_x -= margin_x;
                max_x += margin_x;
                min_y -= margin_y;
                max_y += margin_y;

                // Fit function: quadratic surface in x and y
                TF2* fit = new TF2("fit", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y", min_x, max_x, min_y, max_y);
                //fit->SetNpx(50);
                //fit->SetNpy(50);
                //fit->SetLineColor(kBlack);
                //fit->SetFillColorAlpha(kAzure + 1, 0.35);
                g2->Fit(fit, "WQ");  // , Q = quiet

                //PRINT=============================================================================// Supone que ya hiciste: g2->Fit(fit, "WQ");
                std::vector<double> parameters_julia(6);
                parameters_julia[5] = fit->GetParameter(0);
                parameters_julia[0] = fit->GetParameter(1);
                parameters_julia[1] = fit->GetParameter(2);
                parameters_julia[2] = fit->GetParameter(3);
                parameters_julia[3] = fit->GetParameter(4);
                parameters_julia[4] = fit->GetParameter(5);

                std::cout << "parameters_julia = [";
                bool first = true;
                for (double val : parameters_julia) {
                    if (!first) std::cout << ", ";
                    std::cout << val;
                    first = false;
                }
                std::cout << "]" << std::endl;

                std::cout << "\n--- Fit range ---\n";
                std::cout << "x range: [" << min_x << ", " << max_x << "]" << std::endl;
                std::cout << "y range: [" << min_y << ", " << max_y << "]" << std::endl;



                //PRINT============================================================================

                // Dummy frame to control axis range and labels
                TH3F* h_frame = new TH3F("h_frame", "", 
                         2, min_x, max_x, 
                         2, min_y, max_y, 
                         2, min_z, max_z);


                h_frame->GetXaxis()->SetTitle(x_label.c_str());
                h_frame->GetYaxis()->SetTitle(y_label.c_str());
                //h_frame->GetZaxis()->SetTitle(z_label.c_str());

                h_frame->GetXaxis()->SetTitleSize(0.065);
                h_frame->GetYaxis()->SetTitleSize(0.065);
                //h_frame->GetZaxis()->SetTitleSize(0.045);
                h_frame->GetXaxis()->SetTitleOffset(1.4);
                h_frame->GetYaxis()->SetTitleOffset(1.4);
                //h_frame->GetZaxis()->SetTitleOffset(1.4);
                h_frame->GetXaxis()->SetLabelSize(0.048);
                h_frame->GetYaxis()->SetLabelSize(0.0485);
                h_frame->GetZaxis()->SetLabelSize(0.0485);

                // Set 3D viewing angles
                gPad->SetTheta(30);
                gPad->SetPhi(30);
                gPad->Update();

                // Draw axes, fit surface, and data points
                h_frame->Draw("BOX3");
                fit->Draw("SURF2 SAME");
                g2->Draw("E P SAME");

                // Add legend
                TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9);
                legend->AddEntry(g2, "EFT value", "p");
                //legend->AddEntry(fit, "Fit surface", "f");
                legend->SetTextSize(0.045);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->Draw();

                // Additional text annotations
                TLatex* latex = new TLatex();
                latex->SetNDC();
                
                latex->DrawLatex(0.15, 0.70, "#splitline{m_{tt} > 800 GeV}{|cos#theta| < 0.4}");
                latex->SetTextSize(0.04);
                TPaveText* pt = new TPaveText(0.15, 0.85, 0.85, 0.95, "NDC");
                pt->AddText(title.c_str());
                pt->SetFillColor(0);
                pt->SetBorderSize(0);
                pt->SetTextFont(42);
                pt->SetTextSize(0.08);
                pt->Draw();

                // Save graph and fit in the corresponding observable directory
                std::string graph_name = "graph_" + observable + "_" + operator1 + "_" + operator2;
                g2->SetName(graph_name.c_str());
                fit->SetName(("fit_" + observable + "_" + operator1 + "_" + operator2).c_str());

                output_file->cd();
                output_file->mkdir(observable.c_str());
                output_file->cd(observable.c_str());
                g2->Write();
                fit->Write();

                std::string output_name_pdf = "../plots/3d/" + operator1 + "_" + operator2 + ".pdf";
                c_multi->SaveAs(output_name_pdf.c_str());
            }

            // Save canvas after all observables are drawn
            output_file->cd();
            c_multi->Write(canvas_name.c_str());
        }
    }

    // Close the output file
    output_file->Close();
    std::cout << "\n[INFO] All 2D fits and graphs have been written to output_fits_800GeV.root\n";
}
/// ============================================== END OF MAIN CODE ===================================================






// =======================================================================================================
//       AUXILIARY FUNCTION TO EXTRACT EFT/SM RATIOS OVER 2D WILSON COEFFICIENT GRIDS FOR A GIVEN OBSERVABLE
// =======================================================================================================
/*
 * Given a ROOT file containing a TTree with SMEFT-weighted observable values and metadata,
 * this function extracts the EFT/SM ratio for a specified angular observable across a
 * 2D grid of Wilson coefficient values.
 *
 * For each entry in the TTree:
 *   - It compares the current SMEFT operators to the selected (input) pair
 *   - Normalizes the observable value with respect to the SM prediction (weight_mc_NOSYS)
 *   - Fills a TGraph2DErrors with (C1, C2, EFT/SM) and statistical error bars
 *   - Handles both single-operator and mixed-operator SMEFT weights
 *   - Supports reversed operator orderings automatically
 *
 * Inputs:
 *   - g2                    : pointer to TGraph2DErrors object to be filled
 *   - observable            : name of the observable (e.g. "Ckk", "A_plus")
 *   - input_operator1       : first Wilson coefficient to scan (e.g. "tGRe")
 *   - input_operator2       : second Wilson coefficient to scan (e.g. "td8")
 *   - input_observables_file: path to the ROOT file containing the "observables_tree"
 *
 * Output:
 *   - A filled TGraph2DErrors object (g2) with EFT/SM ratios and associated uncertainties
 *     for the selected observable across the 2D grid defined by (operator1, operator2)
 */

void function_get_points3d(TGraph2DErrors* g2, const char* observable, const char* input_operator1, const char* input_operator2, 
    const char* input_observables_file) {

    // Open the input ROOT file containing observable values and metadata
    TFile* file = TFile::Open(input_observables_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << input_observables_file << std::endl;
        return;
    }

    // Access the TTree named "observables_tree"
    TTree* tree = (TTree*)file->Get("tree_observables");
    if (!tree) {
        std::cerr << "Error: Tree 'tree_observables' not found!" << std::endl;
        return;
    }

    // Declare variables and link branch addresses
    std::string* weight_name = nullptr;
    std::string* operator1 = nullptr;
    std::string* operator2 = nullptr;
    double wilson1, wilson2;
    double N, Ckk, Cnn, Crr, Ckr, Crk, A_plus, A_minus, D3;
    tree->SetBranchAddress("weight_name", &weight_name);
    tree->SetBranchAddress("wilson1_value", &wilson1);
    tree->SetBranchAddress("wilson2_value", &wilson2);
    tree->SetBranchAddress("wilson1_operator", &operator1);
    tree->SetBranchAddress("wilson2_operator", &operator2);

    // Dynamically bind selected observable and its associated error branch
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

    // Search for the SM reference value (corresponding to weight_mc_NOSYS)
    double observable_sm = -999;
    double observable_sm_err;
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (*weight_name == "weight_mc_NOSYS") {
            observable_sm = observable_value;
            observable_sm_err = observable_err;
            break;
        }
    }

    if (observable_sm == -999) {
        std::cerr << "[ERROR] No SM observable found in tree!\n";
        file->Close();
        return;
    }

    // Check if the input operator order matches the order in the tree
    bool match_order = false;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (*operator1 == input_operator1 && *operator2 == input_operator2) {
            match_order = true;
            break;
        }
    }

    // Insert the SM reference point: (0, 0, 1.0)
    g2->SetPoint(0, 0.0, 0.0, observable_sm);
    g2->SetPointError(0, 0.0, 0.0, observable_sm_err);
    int contar = 1;

    // Loop over all tree entries and populate the graph based on operator matching
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Skip SM entry (already included)
        if (*weight_name == "weight_mc_NOSYS") continue;

        // Define matching conditions for single and two-operator weights
        bool match1 = (*operator1 == input_operator1 && operator2->empty());
        bool match2 = (*operator1 == input_operator2 && operator2->empty());
        bool match3 = (*operator1 == input_operator1 && *operator2 == input_operator2);
        bool match4 = (*operator1 == input_operator2 && *operator2 == input_operator1);

        // Compute EFT/SM ratio and its relative error
        double ratio = observable_value;
        double err   = observable_err;

        // Assign point to graph, taking operator order into account
        if (match_order) {
            if (match1 || match3) {
                g2->SetPoint(contar, wilson1, wilson2, ratio);
                g2->SetPointError(contar, 0.0, 0.0, err);
                std::cout << *weight_name << ": (" << wilson1 << ", " << wilson2 << ", " << ratio << ")\n";
                contar++;
            }
            if (match2) {
                g2->SetPoint(contar, wilson2, wilson1, ratio);
                g2->SetPointError(contar, 0.0, 0.0, err);
                std::cout << *weight_name << ": (" << wilson2 << ", " << wilson1 << ", " << ratio << ")\n";
                contar++;
            }
        } else {
            if (match1) {
                g2->SetPoint(contar, wilson2, wilson1, ratio);
                g2->SetPointError(contar, 0.0, 0.0, err);
                std::cout << *weight_name << ": (" << wilson2 << ", " << wilson1 << ", " << ratio << ")\n";
                contar++;
            }
            if (match2 || match4) {
                g2->SetPoint(contar, wilson1, wilson2, ratio);
                g2->SetPointError(contar, 0.0, 0.0, err);
                std::cout << *weight_name << ": (" << wilson1 << ", " << wilson2 << ", " << ratio << ")\n";
                contar++;
            }
        }
    }

    std::cout << contar << std::endl;  // Total number of points inserted (debug)

    // Close input file
    file->Close();
}


