import ROOT
import pandas as pd
import numpy as np
import xgboost as xgb
import os

# Use English for comments in the code.

# --- 1. Configuration ---
MODEL_PATH = "bs_tautau_xgb_model.bin"
TREE_NAME = "tree"

# List of input and output files
# Format: (input_file, output_file)
FILES_TO_PROCESS = [
    ("final_bstautau_vtx.root",  "applied_final_bstautau_vtx.root"),
    ("final_ttsemileptonic_vtx.root",  "applied_final_ttsemileptonic_vtx.root"),

    # Add more files as needed
]

FEATURES = [
    "alpha1", "alpha2", "x1", "x2", "m_vis", "m_exact",
    "pair_pt", "pair_eta", "n_extra_trks", "iso_ratio",
    "tau1_vtxProb", "tau2_vtxProb", "tau1_lxySig", "tau2_lxySig",
#    "j_ParTRawTauhtaue", "j_ParTRawTauhtaumu", 
#    "j_ParTRawTauhtauh", "j_ParTRawSingletau",
    "tau1_m_rho", "tau2_m_rho",
    "tau1_deltaPV_dr", "tau2_deltaPV_dr",
    "tau1_mass", "tau2_mass",
    "tau1_pt", "tau2_pt",
    "tau1_eta", "tau2_eta",
    "j_deepflavB"
]


#FEATURES = [
##    "alpha1", "alpha2", "x1", "x2", "m_vis", "m_exact",
#    "pair_pt", "pair_eta", "n_extra_trks", "iso_ratio",
#    "tau1_vtxProb", "tau2_vtxProb", "tau1_lxySig", "tau2_lxySig",
##    "j_ParTRawTauhtaue", "j_ParTRawTauhtaumu", 
##    "j_ParTRawTauhtauh", "j_ParTRawSingletau",
#    "tau1_m_rho", "tau2_m_rho",
#    "tau1_deltaPV_dr", "tau2_deltaPV_dr",
#    "tau1_mass", "tau2_mass",
#    "tau1_pt", "tau2_pt",
#    "tau1_eta", "tau2_eta",
#    "j_deepflavB"
#]



# v3: full variables
#FEATURES = [
#    "alpha1", "alpha2", "x1", "x2", "m_vis", "m_exact",
#    "pair_pt", "pair_eta", "n_extra_trks", "iso_ratio",
#    "tau1_vtxProb", "tau2_vtxProb", "tau1_lxySig", "tau2_lxySig",
#    "j_ParTRawTauhtaue", "j_ParTRawTauhtaumu", 
#    "j_ParTRawTauhtauh", "j_ParTRawSingletau",
#    "tau1_m_rho", "tau2_m_rho",
##    "tau1_deltaPV_dx", "tau2_deltaPV_dx",
##    "tau1_deltaPV_dy", "tau2_deltaPV_dy",
##    "tau1_deltaPV_dz", "tau2_deltaPV_dz",
#    "tau1_mass", "tau2_mass",
#    "tau1_pt", "tau2_pt",
#    "tau1_eta", "tau2_eta",
#    "j_deepflavB"
#]


# v1
#FEATURES = [
#    "pair_pt", "pair_eta", "n_extra_trks", "iso_ratio",
#    "tau1_vtxProb", "tau2_vtxProb", "tau1_lxySig", "tau2_lxySig",
##    "j_ParTRawTauhtaue", "j_ParTRawTauhtaumu", 
##    "j_ParTRawTauhtauh", "j_ParTRawSingletau",
#    "tau1_m_rho", "tau2_m_rho",
##    "tau1_deltaPV_dx", "tau2_deltaPV_dx",
##    "tau1_deltaPV_dy", "tau2_deltaPV_dy",
##    "tau1_deltaPV_dz", "tau2_deltaPV_dz",
#    "tau1_mass", "tau2_mass",
#    "tau1_pt", "tau2_pt",
#    "tau1_eta", "tau2_eta",
#    "j_deepflavB"
#]
#
## v2: v1 + adding di-tau score
#FEATURES = [
#    "pair_pt", "pair_eta", "n_extra_trks", "iso_ratio",
#    "tau1_vtxProb", "tau2_vtxProb", "tau1_lxySig", "tau2_lxySig",
#    "j_ParTRawTauhtaue", "j_ParTRawTauhtaumu", 
#    "j_ParTRawTauhtauh", "j_ParTRawSingletau",
#    "tau1_m_rho", "tau2_m_rho",
##    "tau1_deltaPV_dx", "tau2_deltaPV_dx",
##    "tau1_deltaPV_dy", "tau2_deltaPV_dy",
##    "tau1_deltaPV_dz", "tau2_deltaPV_dz",
#    "tau1_mass", "tau2_mass",
#    "tau1_pt", "tau2_pt",
#    "tau1_eta", "tau2_eta",
#    "j_deepflavB"
#]






def apply_mva():
    # Load the trained model
    print("Loading model: {}".format(MODEL_PATH))
    bst = xgb.Booster()
    bst.load_model(MODEL_PATH)

    for in_file, out_file in FILES_TO_PROCESS:
        if not os.path.exists(in_file):
            print("Warning: {} not found. Skipping...".format(in_file))
            continue

        print("Processing: {} -> {}".format(in_file, out_file))
        
        # Open original file and get tree
        f_in = ROOT.TFile.Open(in_file)
        
        hist = f_in.Get("h_genEventSumw")

        tree = f_in.Get(TREE_NAME)
        
        # 1. Convert Tree to Dataframe for fast inference
        # In Python 2.7/PyROOT, we iterate and fill a dict
        data_dict = {f: [] for f in FEATURES}
        for event in tree:
            for f in FEATURES:
                data_dict[f].append(getattr(event, f))
        
        df = pd.DataFrame(data_dict)
        
        # 2. Inference
        dmatrix = xgb.DMatrix(df[FEATURES])
        scores = bst.predict(dmatrix) # This returns the probability [0, 1]

        # 3. Write to New ROOT File
        f_out = ROOT.TFile(out_file, "RECREATE")
        new_tree = tree.CloneTree(0) # Clone structure without entries
        
        # Create a buffer for the new branch
        mva_score = np.array([0], dtype=np.float32)
        new_tree.Branch("mva_score", mva_score, "mva_score/F")
        
        print("Writing scores to new tree...")
        num_entries = tree.GetEntries()
        for i in range(num_entries):
            tree.GetEntry(i)
            # Fill the buffer with the calculated score
            mva_score[0] = scores[i]
            new_tree.Fill()
            
            if i % 10000 == 0:
                print("Progress: {}/{}".format(i, num_entries))

        new_tree.Write()
        hist.Write()
        f_out.Close()
        f_in.Close()
        print("Successfully created: {}".format(out_file))

if __name__ == "__main__":
    apply_mva()
