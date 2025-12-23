// =============================================================================================
//        FUNCTION TO COMPUTE QUANTUM OBSERVABLES FROM ANGULAR HISTOGRAMS IN A ROOT FILE
// =============================================================================================
/*
 * Author: [Josue Elizalde]
 * Date: [28.05.25]
 * 
 * Description:
 *   This script reads a ROOT file containing parton-level events from simulated top-quark pair production.
 *   For each SMEFT weight branch, it computes several quantum observables based on angular correlations.
 *   These observables characterize spin correlations and asymmetries between decay products of the top quarks.
 *
 *   Computed observables include:
 *     - Number of events N
 *     - Spin correlation coefficients: Ckk, Crr, Cnn, Ckr, Crk
 *     - D3 observable from the angle between final particles: cos(θ_p) 
 *     - Azimuthal asymmetry A⁺ 
 *     - Statistical uncertainties for all observables
 *
 *   Results are written into a new TTree for further SMEFT interpretation or fitting.
 *
 * Use cases:
 *   - EFT sensitivity studies
 *   - Angular analysis in top physics
 *   - Building input for global SMEFT fits
 */

// =============== ROOT Includes ===============
#include <TFile.h>        // File I/O
#include <TTree.h>        // Tree structures
#include <TH1F.h>         // 1D histograms

#include <iostream>       // Console I/O
#include <vector>         // For storing list of branches
#include <string>         // String handling
#include <regex>          // Regex to extract Wilson coefficients

struct Observables {
    double N;         double N_err;         // Number of Events
    double Ckk;       double Ckk_err;       // Helicity axis correlation
    double Crr;       double Crr_err;       // r-axis correlation
    double Cnn;       double Cnn_err;       // Transverse axis correlation
    double Ckr;       double Ckr_err;       // Mixed helicity × r-axis
    double Crk;       double Crk_err;       // Mixed r-axis × helicity
    double D3;        double D3_err;        // Derived from cos(θ_p)
    double A_plus;    double A_plus_err;    // Symmetric azimuthal asymmetry
    double A_minus;   double A_minus_err;   // Antisymmetric azimuthal asymmetry
};


// Function Declarations
// Main function to compute all observables for a given SMEFT weight
Observables function_get_observables(const char* root_file_name, const char* weight_branch_name,
                                     int number_bins, float spin_analyzing_power_a, float spin_analyzing_power_b);

// Extracts names and values of Wilson coefficients from weight name (e.g., "weight_mc_cQj11_p0p3_cQd8_m1p5")
void extract_wilson_values(const std::string& weight_name,
                           std::string& wilson1_operator, double& wilson1_value,
                           std::string& wilson2_operator, double& wilson2_value);



// =============================================================================================
//                                         MAIN CODE 
// =============================================================================================

void get_observables() {
    // Input ROOT file with simulated ttbar events and SMEFT weights
    const char* input_file = "../data/fourthOutput.root";

    // Open the input file safely
    TFile* file = TFile::Open(input_file, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open input file " << input_file << std::endl;
        return;
    }

    // Retrieve the main TTree
    TTree* tree_input = (TTree*)file->Get("truth");
    if (!tree_input) {
        std::cerr << "Error: 'truth' TTree not found in input file." << std::endl;
        file->Close();
        return;
    }

    // Output ROOT file to store computed observables
    const char* output_file = "../results/observables_400GeV.root";
    TFile* outfile = new TFile(output_file, "RECREATE");

    // Create TTree for saving final observables
    TTree* tree_observables = new TTree("tree_observables", "Tree with EFT observables");

    // Define variables to be stored in TTree branches 
    double n_events, n_events_err;   // number of events                   
    tree_observables->Branch("number_events", &n_events);   tree_observables->Branch("number_events_err", &n_events_err); 

    std::string wilson1_operator, wilson2_operator, weight_name;    // wilson_coefficients
    double wilson1_value, wilson2_value;
    tree_observables->Branch("weight_name", &weight_name);
    tree_observables->Branch("wilson1_operator", &wilson1_operator);
    tree_observables->Branch("wilson1_value", &wilson1_value);
    tree_observables->Branch("wilson2_operator", &wilson2_operator);
    tree_observables->Branch("wilson2_value", &wilson2_value);
    
    double Ckk, Ckk_err, Crr, Crr_err, Cnn, Cnn_err, Ckr, Ckr_err, Crk, Crk_err;    // correlations (with asociated errors)
    tree_observables->Branch("Ckk", &Ckk);   tree_observables->Branch("Ckk_err", &Ckk_err);
    tree_observables->Branch("Crr", &Crr);   tree_observables->Branch("Crr_err", &Crr_err);
    tree_observables->Branch("Cnn", &Cnn);   tree_observables->Branch("Cnn_err", &Cnn_err);
    tree_observables->Branch("Ckr", &Ckr);   tree_observables->Branch("Ckr_err", &Ckr_err);
    tree_observables->Branch("Crk", &Crk);   tree_observables->Branch("Crk_err", &Crk_err);
    
    double D3, D3_err, A_plus, A_plus_err, A_minus, A_minus_err; // other observables (with asociated errors)
    tree_observables->Branch("D3", &D3);             tree_observables->Branch("D3_err", &D3_err);
    tree_observables->Branch("A_plus", &A_plus);     tree_observables->Branch("A_plus_err", &A_plus_err);
    tree_observables->Branch("A_minus", &A_minus);   tree_observables->Branch("A_minus_err", &A_minus_err);

    // Loop over branches to find SMEFT weights
    TObjArray* branches = tree_input->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        std::string branch_name = branches->At(i)->GetName();

        // Only keep weight branches
        if (branch_name.find("weight_mc_") != 0) continue;

        // Extract metadata from branch name
        weight_name = branch_name;
        extract_wilson_values(weight_name, wilson1_operator, wilson1_value,
                              wilson2_operator, wilson2_value);

        std::cout << "\n[INFO] Processing " << weight_name
                  << " with: " << wilson1_operator << " = " << wilson1_value
                  << ", " << wilson2_operator << " = " << wilson2_value << std::endl;

        // Compute observables from angular histograms
        Observables o = function_get_observables(input_file, branch_name.c_str(), 22, 1.0, -1.0);

        // Assign results
        n_events = o.N;  n_events_err = o.N_err;   
        Ckk = o.Ckk;     Ckk_err = o.Ckk_err;
        Crr = o.Crr;     Crr_err = o.Crr_err;
        Cnn = o.Cnn;     Cnn_err = o.Cnn_err;
        Ckr = o.Ckr;     Ckr_err = o.Ckr_err;
        Crk = o.Crk;     Crk_err = o.Crk_err;
        D3  = o.D3;      D3_err  = o.D3_err;
        A_plus = o.A_plus;   A_plus_err = o.A_plus_err;
        A_minus = o.A_minus; A_minus_err = o.A_minus_err;

        std::cout << "  → Ckk = " << Ckk << " ± " << Ckk_err
                  << ", N = " <<  n_events  << std::endl;

        // Save to output TTree
        tree_observables->Fill();
    }

    // Finalize output
    outfile->cd();
    tree_observables->Write();
    outfile->Close();

    std::cout << "\n[INFO] Saved observables to: " << output_file << std::endl;
}




/// =================================================END OF MAIN CODE=====================================================

// =======================================================================================================
//        AUXILIARY FUNCTION TO EXTRACT QUANTUM OBSERVABLES FROM SMEFT-WEIGHTED EVENTS IN A ROOT TTREE
// =======================================================================================================
/*
 * Given a ROOT file containing parton-level truth information with SMEFT weight branches,
 * this function computes several angular observables and spin correlations relevant for
 * top-quark pair production analyses.
 *
 * Observables include:
 *   - Spin correlation coefficients (C_ij) along helicity (k), r, and transverse (n) axes
 *   - Azimuthal asymmetries A⁺ and A⁻
 *   - The D3 observable derived from cos(θₚ)
 *
 * Inputs:
 *   - root_file_name          : name of the input ROOT file
 *   - weight_branch_name      : branch containing the SMEFT weight
 *   - number_bins             : binning for angular distributions
 *   - spin_analyzing_power_a  : spin-analyzing power for particle a (e.g. lepton)
 *   - spin_analyzing_power_b  : spin-analyzing power for particle b (e.g. lepton or b-jet)
 *
 * Output:
 *   - A struct containing the computed observables and their statistical uncertainties
 */

Observables function_get_observables(const char* root_file_name, const char* weight_branch_name, 
    int number_bins, float spin_analyzing_power_a, float spin_analyzing_power_b) {

    // Open input ROOT file
    TFile *file = TFile::Open(root_file_name, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << root_file_name << std::endl;
        return {};
    }

    gROOT->cd(); // Important: create histograms in memory, not inside input file!

    // Retrieve TTree containing parton-level information
    TTree *tree = (TTree*)file->Get("truth");
    if (!tree) { 
        std::cerr << "Error: TTree 'truth' not found in file\n"; 
        return {}; 
    }

    // Verify if the TTree contains branch with the specified weight
    if (tree->GetBranch(weight_branch_name) == nullptr) {
        std::cerr << "Error: Branch " << weight_branch_name << " not found in TTree 'truth'\n";
        return {};
    }

    // Get message for succesfull opening and to get the barnch 
    std::cout << "Successfully opened tree and branch: " << weight_branch_name << std::endl;

    // Define kinematic variables in different spin quantization axes: (p = particle, m = antiparticle)
    float cos_theta_helicity_p, cos_theta_raxis_p, cos_theta_transverse_p;  // helicity (k), r-axis (r), transverse (n)
    tree->SetBranchAddress("spin_SDM_cos_theta_helicity_p_NOSYS", &cos_theta_helicity_p);
    tree->SetBranchAddress("spin_SDM_cos_theta_transverse_p_NOSYS", &cos_theta_transverse_p);
    tree->SetBranchAddress("spin_SDM_cos_theta_raxis_p_NOSYS", &cos_theta_raxis_p);

    float cos_theta_helicity_m, cos_theta_raxis_m, cos_theta_transverse_m; // helicity (k), r-axis (r), transverse (n)
    tree->SetBranchAddress("spin_SDM_cos_theta_helicity_m_NOSYS", &cos_theta_helicity_m);
    tree->SetBranchAddress("spin_SDM_cos_theta_transverse_m_NOSYS", &cos_theta_transverse_m);
    tree->SetBranchAddress("spin_SDM_cos_theta_raxis_m_NOSYS", &cos_theta_raxis_m);

    //Define kinematic variables
    float cos_theta_top;    // Define Cos theta top
    tree->SetBranchAddress("spin_SDM_cos_theta_top_NOSYS", &cos_theta_top);
    float cos_theta_p;      // Angular observable cosθab (angle between final particles)  
    tree->SetBranchAddress("spin_SDM_theta_p_NOSYS", &cos_theta_p);  
    float m_ttbar;          // Invariant mass of the ttbar system (after final state radiation FSR)
    tree->SetBranchAddress("Ttbar_MC_ttbar_afterFSR_m", &m_ttbar); 
    float weight;           // Event weight (depends on selected SMEFT contribution or MC variation)
    tree->SetBranchAddress(weight_branch_name, &weight);
    
    // Define Histograms for angular correlations in helicity base: cos(theta_i^p) * cos(theta_j^m), i,j = k, r, n
    int n_bins = number_bins;
    TH1F *hist_kk = new TH1F(Form("hist_kk_%s", weight_branch_name), Form("kk-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);
    TH1F *hist_kr = new TH1F(Form("hist_kr_%s", weight_branch_name), Form("kr-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);
    TH1F *hist_kn = new TH1F(Form("hist_kn_%s", weight_branch_name), Form("kn-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);

    TH1F *hist_rr = new TH1F(Form("hist_rr_%s", weight_branch_name), Form("rr-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);
    TH1F *hist_rk = new TH1F(Form("hist_rk_%s", weight_branch_name), Form("rk-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);
    TH1F *hist_rn = new TH1F(Form("hist_rn_%s", weight_branch_name), Form("rn-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);

    TH1F *hist_nn = new TH1F(Form("hist_nn_%s", weight_branch_name), Form("nn-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);
    TH1F *hist_nk = new TH1F(Form("hist_nk_%s", weight_branch_name), Form("nk-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);
    TH1F *hist_nr = new TH1F(Form("hist_nr_%s", weight_branch_name), Form("nr-axis (%s)", weight_branch_name), n_bins, -1.15, 1.15);

    hist_kk->GetXaxis()->SetTitle("cos(#theta_{k}^{+})cos(#theta_{k}^{-})");
    hist_kr->GetXaxis()->SetTitle("cos(#theta_{k}^{+})cos(#theta_{r}^{-})");
    hist_kn->GetXaxis()->SetTitle("cos(#theta_{k}^{+})cos(#theta_{n}^{-})");

    hist_rr->GetXaxis()->SetTitle("cos(#theta_{r}^{+})cos(#theta_{r}^{-})");
    hist_rk->GetXaxis()->SetTitle("cos(#theta_{r}^{+})cos(#theta_{k}^{-})");
    hist_rn->GetXaxis()->SetTitle("cos(#theta_{r}^{+})cos(#theta_{n}^{-})");

    hist_nn->GetXaxis()->SetTitle("cos(#theta_{n}^{+})cos(#theta_{n}^{-})");
    hist_nk->GetXaxis()->SetTitle("cos(#theta_{n}^{+})cos(#theta_{k}^{-})");
    hist_nr->GetXaxis()->SetTitle("cos(#theta_{n}^{+})cos(#theta_{r}^{-})");

    // Define Histogram for cos(theta_p)
    TH1F *hist_cos_theta_p = new TH1F(Form("hist_cos_theta_p_%s", weight_branch_name), Form("cos(#theta_{p}) (%s)", weight_branch_name), n_bins, -1.15, 1.15);  
    hist_cos_theta_p->GetXaxis()->SetTitle("cos(#theta_{p})");

    // Define Histogram for cos_theta_top
    TH1F *hist_cos_theta_top = new TH1F(Form("hist_cos_theta_top_%s", weight_branch_name), Form("cos(#theta_{top}) (%s)", weight_branch_name), n_bins, -1.15, 1.15);  
    hist_cos_theta_top->GetXaxis()->SetTitle("cos(#theta_{top})");

    const double m_ttbar_cut = 400e3;       // Minimum m_ttbar cut (400-800 GeV)
    const double abs_cos_theta_cut = 0.2;   // minimum angle for top
    double n_events = 0.0;                  // number of events
    double n_events_err = 0.0; 
    // Loop over events and fill histograms
    Long64_t n_entries = tree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);

        //Apply cuts (threshold or boosted regions)
        //if (m_ttbar < m_ttbar_cut || fabs(cos_theta_top) > abs_cos_theta_cut ) continue;  // top mass and angle cut =========== CUT ========
        if (m_ttbar > m_ttbar_cut) continue;

        // Calculate number of events
        n_events += weight;
        n_events_err += weight * weight;

        // Fill histograms
        hist_kk->Fill(cos_theta_helicity_p * cos_theta_helicity_m,   weight);
        hist_kr->Fill(cos_theta_helicity_p * cos_theta_raxis_m,      weight);
        hist_kn->Fill(cos_theta_helicity_p * cos_theta_transverse_m, weight);

        hist_rr->Fill(cos_theta_raxis_p * cos_theta_raxis_m,      weight);
        hist_rk->Fill(cos_theta_raxis_p * cos_theta_helicity_m,   weight);
        hist_rn->Fill(cos_theta_raxis_p * cos_theta_transverse_m, weight);

        hist_nn->Fill(cos_theta_transverse_p * cos_theta_transverse_m, weight);
        hist_nk->Fill(cos_theta_transverse_p * cos_theta_helicity_m,   weight);
        hist_nr->Fill(cos_theta_transverse_p * cos_theta_raxis_m,      weight);

        hist_cos_theta_p->Fill(cos_theta_p, weight);

        hist_cos_theta_top->Fill(cos_theta_top, weight);
    }

    // Normalize histograms
    hist_kk->Scale(1.0 / hist_kk->Integral());
    hist_rr->Scale(1.0 / hist_rr->Integral());
    hist_nn->Scale(1.0 / hist_nn->Integral());

    hist_cos_theta_p->Scale(1.0 / hist_cos_theta_p->Integral());
    hist_cos_theta_top->Scale(1.0 / hist_cos_theta_top->Integral());

    hist_kr->Scale(1.0 / hist_kr->Integral());
    hist_kn->Scale(1.0 / hist_kn->Integral());

    hist_rk->Scale(1.0 / hist_rk->Integral());
    hist_rn->Scale(1.0 / hist_rn->Integral());

    hist_nk->Scale(1.0 / hist_nk->Integral());
    hist_nr->Scale(1.0 / hist_nr->Integral());

    // Calculate obsevables
    double Ckk = (9.0 / (spin_analyzing_power_a * spin_analyzing_power_b)) * hist_kk->GetMean();
    double Crr = (9.0 / (spin_analyzing_power_a * spin_analyzing_power_b)) * hist_rr->GetMean();
    double Cnn = (9.0 / (spin_analyzing_power_a * spin_analyzing_power_b)) * hist_nn->GetMean();
    double Ckr = (9.0 / (spin_analyzing_power_a * spin_analyzing_power_b)) * hist_kr->GetMean();
    double Crk = (9.0 / (spin_analyzing_power_a * spin_analyzing_power_b)) * hist_rk->GetMean();
    double D3 = 3.0 * hist_cos_theta_p->GetMean();


    // asymmetry: from Improved tests of entanglement and Bell inequalities with LHC tops- Aguilar-Casas
    double A_plus = (M_PI / 16.0) * (spin_analyzing_power_a * spin_analyzing_power_b) * (Ckk + Crr);

    // Calculate statistical uncertainties
    double Ckk_err = (9.0 / std::abs((spin_analyzing_power_a * spin_analyzing_power_b))) * hist_kk->GetMeanError();
    double Crr_err = (9.0 / std::abs((spin_analyzing_power_a * spin_analyzing_power_b))) * hist_rr->GetMeanError();
    double Cnn_err = (9.0 / std::abs((spin_analyzing_power_a * spin_analyzing_power_b))) * hist_nn->GetMeanError();
    double Ckr_err = (9.0 / std::abs((spin_analyzing_power_a * spin_analyzing_power_b))) * hist_kr->GetMeanError();
    double Crk_err = (9.0 / std::abs((spin_analyzing_power_a * spin_analyzing_power_b))) * hist_rk->GetMeanError();
    double D3_err  = 3.0 * hist_cos_theta_p->GetMeanError();

    double A_plus_err  = (M_PI / 16.0) * std::abs((spin_analyzing_power_a * spin_analyzing_power_b)) * std::sqrt(Crr_err * Crr_err + Ckk_err * Ckk_err);

    // Close File
    file->Close();
    delete file;

    // Return values with struct Observables format
    Observables result;
    result.N          = n_events;
    result.N_err      = sqrt(n_events_err);
    result.Ckk        = Ckk;
    result.Ckk_err    = Ckk_err;
    result.Crr        = Crr;
    result.Crr_err    = Crr_err;
    result.Cnn        = Cnn;
    result.Cnn_err    = Cnn_err;
    result.Ckr        = Ckr;
    result.Ckr_err    = Ckr_err;
    result.Crk        = Crk;
    result.Crk_err    = Crk_err;
    result.D3         = D3;
    result.D3_err     = D3_err;
    result.A_plus     = A_plus;
    result.A_plus_err = A_plus_err;
    return result;

}


// =======================================================================================================================
//                           FUNCTION TO EXTRACT WILSON COEFFICIENT VALUES AND NAMES
// =======================================================================================================================
void extract_wilson_values(const std::string& weight_name,
                           std::string& wilson1_operator, double& wilson1_value,
                           std::string& wilson2_operator, double& wilson2_value) {
    // Initialize output values
    wilson1_operator = "";
    wilson2_operator = "";
    wilson1_value = 0.0;
    wilson2_value = 0.0;

    // Regular expression to match patterns like: cOperator_m1p0 or cOperator_p0p3
    std::regex pattern("(c[[:alnum:]]+)_([mp]\\d+p\\d+)");
    std::sregex_iterator iter(weight_name.begin(), weight_name.end(), pattern);
    std::sregex_iterator end;

    int count = 0;
    for (; iter != end && count < 2; ++iter, ++count) {
        // Extract operator name, e.g., "cQj31"
        std::string op = (*iter)[1].str();

        // Remove the leading 'c' to get a cleaner name like "Qj31"
        if (!op.empty() && op[0] == 'c') {
            op = op.substr(1);
        }

        // Extract encoded value string, e.g., "m1p5" or "p0p3"
        std::string val = (*iter)[2].str();

        // Convert value string to double
        // Example: "m1p5" → -1.5, "p0p3" → 0.3
        double value = 0.0;
        std::string sign = val.substr(0, 1);                // 'm' or 'p'
        std::string val_no_sign = val.substr(1);            // e.g., "1p5"
        std::replace(val_no_sign.begin(), val_no_sign.end(), 'p', '.'); // e.g., "1.5"
        value = std::stod(val_no_sign);                     // convert to double
        if (sign == "m") value *= -1;                       // apply sign

        // Assign values to output variables
        if (count == 0) {
            wilson1_operator = op;
            wilson1_value = value;
        } else if (count == 1) {
            wilson2_operator = op;
            wilson2_value = value;
        }
    }
}




