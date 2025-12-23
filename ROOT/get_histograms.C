// =============================================================================================
//         FUNCTION TO EXTRACT AND SAVE HISTOGRAMS OF ANGULAR OBSERVABLES FROM A ROOT FILE
// =============================================================================================
/*
 * Author: [Josue Elizalde]
 * Date: [28.05.25]
 * 
 * Description:
 *   This code processes a ROOT file containing parton-level information
 *   from simulated top-quark pair events. For each SMEFT weight branch found
 *   in the TTree, it computes histograms of angular observables and
 *   spin correlation, applying event selection cuts.
 *   Histograms are normalized and saved to a separate ROOT file for further analysis.
 *
 *   Histograms extracted:
 *     - Spin correlation distributions: products of cos(theta_top) * cos(theta_antitop) projected along the helicity, r-axis, and transverse bases
 *     - Angular distribution between final-state particles: cos(theta_p)
 *     - Angular distribution of the top quark direction: cos(theta_top)
 *
 *   Target use:
 *     - SMEFT sensitivity studies
 *     - Comparison of angular observables across different operator contributions and regions
 * 
 *
 */
// ROOT headers for file handling, histograms, and plotting
#include <TFile.h>        // For TFile (reading/writing ROOT files)
#include <TTree.h>        // For TTree (event data structure)
#include <TH1F.h>         // For TH1F (1D histograms)
#include <TCanvas.h>      // For TCanvas (optional: plotting purposes)
#include <TLegend.h>      // For legends in plots (optional)

#include <iostream>       // For std::cout and std::cerr
#include <TApplication.h> // Required only if using graphical display (optional)


// Histograms needed to extract values for the spin density matrix
struct Histograms {
    TH1F* hist_kk;             // cos(θ_k⁺) × cos(θ_k⁻)
    TH1F* hist_rr;             // cos(θ_r⁺) × cos(θ_r⁻)
    TH1F* hist_nn;             // cos(θ_n⁺) × cos(θ_n⁻)

    TH1F* hist_kr;             // cos(θ_k⁺) × cos(θ_r⁻)
    TH1F* hist_kn;             // cos(θ_k⁺) × cos(θ_n⁻)
    TH1F* hist_rk;             // cos(θ_r⁺) × cos(θ_k⁻)
    TH1F* hist_rn;             // cos(θ_r⁺) × cos(θ_n⁻)
    TH1F* hist_nk;             // cos(θ_n⁺) × cos(θ_k⁻)
    TH1F* hist_nr;             // cos(θ_n⁺) × cos(θ_r⁻)

    TH1F* hist_cos_theta_p;    // Distribution cosθ_ab (angle between final particles)  
    TH1F* hist_cos_theta_top;  // Distribution of cos(θ_top)
};

// Function Declarations
Histograms extract_histograms(const char* root_file_name, const char* weight_branch_name, int number_bins);


// =================================================================================================================
//                                             M A I N    C O D E
// =================================================================================================================
void get_histograms() {

    // Path to the input ROOT file
    const char* input_file = "../data/lastData.root";
    std::cout << "Opening file: " << input_file << std::endl;

    // Load file and retrieve the TTree named "truth"
    TFile* file = TFile::Open(input_file, "READ");
    TTree* tree = (TTree*)file->Get("truth");

    // Extract list of all weight branches (i.e., those containing "weight_mc")
    std::vector<std::string> weight_list;
    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        std::string branch_name = branches->At(i)->GetName();
        if (branch_name.find("weight_mc") != std::string::npos) {
            weight_list.push_back(branch_name);
            std::cout << branch_name << "\n";
        }
    }

    // Output ROOT file to store histograms
    std::string filename = "../histograms/histograms_800GeV.root"; // Customize filename based on cuts
    TFile* output_file = new TFile(filename.c_str(), "RECREATE");
    if (!output_file->IsOpen()) {
        std::cerr << "Error: Could not create output ROOT file!" << std::endl;
        return;
    }

    // Loop over each SMEFT weight branch and compute histograms
    for (const auto& weight_name : weight_list) {
        std::cout << "Processing " << weight_name << std::endl;
        Histograms h = extract_histograms(input_file, weight_name.c_str(), 22);

        // Write all histograms to the output file
        output_file->cd();
        h.hist_kk->Write(); h.hist_rr->Write(); h.hist_nn->Write();
        h.hist_kr->Write(); h.hist_kn->Write(); h.hist_rk->Write();
        h.hist_rn->Write(); h.hist_nk->Write(); h.hist_nr->Write();
        h.hist_cos_theta_p->Write(); h.hist_cos_theta_top->Write();

        // Free dynamically allocated memory
        delete h.hist_kk; delete h.hist_rr; delete h.hist_nn;
        delete h.hist_kr; delete h.hist_kn; delete h.hist_rk;
        delete h.hist_rn; delete h.hist_nk; delete h.hist_nr;
        delete h.hist_cos_theta_p; delete h.hist_cos_theta_top;
    }

    std::cout << filename << " created!" << std::endl;
    output_file->Close();
    delete output_file;
}


// =============================== AUXILIARY FUNCTION ========================================
/*
 * Extract and return a set of angular histograms for a specific SMEFT weight branch.
 * Applies selection cuts and computes normalized angular correlation distributions.
 */

Histograms extract_histograms(const char* root_file_name, const char* weight_branch_name, int number_bins) {

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

    std::cout << "Successfully opened tree and branch: " << weight_branch_name << std::endl;

    // Define kinematic variables in different spin quantization axes: (p = particle, m = antiparticle)
    // helicity (k), r-axis (r), transverse (n)
    float cos_theta_helicity_p, cos_theta_raxis_p, cos_theta_transverse_p; 
    tree->SetBranchAddress("spin_SDM_cos_theta_helicity_p_NOSYS", &cos_theta_helicity_p);
    tree->SetBranchAddress("spin_SDM_cos_theta_transverse_p_NOSYS", &cos_theta_transverse_p);
    tree->SetBranchAddress("spin_SDM_cos_theta_raxis_p_NOSYS", &cos_theta_raxis_p);

    float cos_theta_helicity_m, cos_theta_raxis_m, cos_theta_transverse_m;
    tree->SetBranchAddress("spin_SDM_cos_theta_helicity_m_NOSYS", &cos_theta_helicity_m);
    tree->SetBranchAddress("spin_SDM_cos_theta_transverse_m_NOSYS", &cos_theta_transverse_m);
    tree->SetBranchAddress("spin_SDM_cos_theta_raxis_m_NOSYS", &cos_theta_raxis_m);

    //More kinematic variables
    // Define Cos theta top
    float cos_theta_top;
    tree->SetBranchAddress("spin_SDM_cos_theta_top_NOSYS", &cos_theta_top);

    // Angular observable cosθab (angle between final particles)  
    float cos_theta_p;
    tree->SetBranchAddress("spin_SDM_theta_p_NOSYS", &cos_theta_p);  
    
    // Invariant mass of the ttbar system (after final state radiation FSR)
    float m_ttbar;
    tree->SetBranchAddress("Ttbar_MC_ttbar_afterFSR_m", &m_ttbar); 

    // Event weight (depends on selected SMEFT contribution or MC variation)
    float weight;
    tree->SetBranchAddress(weight_branch_name, &weight);
    
    // Define Histograms for angular correlations in helicity base: cos(theta_i^p) * cos(theta_j^m), i,j = k, r, n
    int n_bins = number_bins;
    // Naming convention: hist_[axis1][axis2]_[weight]
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

    // Event selection cuts:
    const double m_ttbar_cut = 800e3; // Selects boosted region
    const double abs_cos_theta_cut = 0.4; // Top quark angular cut
    
    // Loop over events and fill histograms
    Long64_t n_entries = tree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);
        
        // Apply cuts (threshold or boosted regions)  
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                  
        //if (m_ttbar < m_ttbar_cut || fabs(cos_theta_top) > abs_cos_theta_cut ) continue;  //  m_ttbar and top angle cut  
        if (m_ttbar < m_ttbar_cut) continue;  //  m_ttbar and top angle cut  
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Fill histograms with angular correlations weighted by SMEFT contribution
        hist_kk->Fill(cos_theta_helicity_p * cos_theta_helicity_m,   weight); // Example: helicity × helicity → cos(θ^+_k) × cos(θ^-_k)
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


    // Close File
    file->Close();
    delete file;

    return {hist_kk, hist_rr, hist_nn, hist_kr, hist_kn,
        hist_rk, hist_rn, hist_nk, hist_nr, hist_cos_theta_p, hist_cos_theta_top};

}



