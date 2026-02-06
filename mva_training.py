import ROOT
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import os
from xgboost import plot_importance

# Use English for comments in the code.

def load_data(input_file, tree_name, var_list, n_events=-1):
    """Convert ROOT TTree to Pandas DataFrame with filtering and weights.
       n_events: Number of events to load. -1 means all events.
    """
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        print("Error: Could not open {0}".format(input_file))
        return pd.DataFrame()
        
    tree = f.Get(tree_name)
    if not tree:
        print("Error: Tree {0} not found in {1}".format(tree_name, input_file))
        f.Close()
        return pd.DataFrame()

    # Include weight components and pair_pt for filtering
    needed_vars = list(set(var_list + ['is_true_signal', 'pair_pt', 'L1PreFiringWeight_Nom', 'genWeight', 'puWeight']))
    
    data_dict = {var: [] for var in needed_vars}
    
    count = 0
    for event in tree:
        # Stop if we reached the requested number of events
        if n_events > 0 and count >= n_events:
            break

        # Filter: Only use events where reconstruction was successful
        if event.pair_pt == -1:
            continue
            
        for var in needed_vars:
            data_dict[var].append(getattr(event, var))
        
        count += 1
    
    f.Close()
    df = pd.DataFrame(data_dict)
    
    if not df.empty:
        # Calculate event weight: PreFiring * GenWeight * Pileup
        df['event_weight'] = df['L1PreFiringWeight_Nom'] * df['genWeight'] * df['puWeight']
    
    print("Loaded {0} events from {1}".format(len(df), input_file))
    return df

# Features configuration
FEATURES = [
    "alpha1", "alpha2", "x1", "x2", "m_vis", "m_exact",
    "pair_pt", "pair_eta", "n_extra_trks", "iso_ratio",
    "tau1_vtxProb", "tau2_vtxProb", "tau1_lxySig", "tau2_lxySig",
    "tau1_m_rho", "tau2_m_rho",
    "tau1_deltaPV_dr", "tau2_deltaPV_dr",
    "tau1_mass", "tau2_mass",
    "tau1_pt", "tau2_pt",
    "tau1_eta", "tau2_eta",
    "j_deepflavB"
]

def train(sig_file, bkg_file, n_sig=-1, n_bkg=-1):
    if not os.path.exists('plots'):
        os.makedirs('plots')

    # 1. Load Data
    print("Loading Signal data...")
    df_sig = load_data(sig_file, "tree", FEATURES, n_events=n_sig)
    print("Loading Background data...")
    df_bkg = load_data(bkg_file, "tree", FEATURES, n_events=n_bkg)
    
    # Combine datasets
    df = pd.concat([df_sig, df_bkg], ignore_index=True)
    
    if df.empty:
        print("Error: DataFrame is empty. Check input files.")
        return

    X = df[FEATURES]
    y = df['is_true_signal']
    weights = df['event_weight']
    
    # Split into Training and Testing (80% / 20%)
    X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
        X, y, weights, test_size=0.2, random_state=42
    )

    # 2. Configure XGBoost with scale_pos_weight
    # Ratio of Backgrounds to Signals to balance the classes
    n_pos = np.sum(y_train == 1)
    n_neg = np.sum(y_train == 0)
    spw = float(n_neg) / n_pos if n_pos > 0 else 1.0
    print("Calculated scale_pos_weight: {0:.2f}".format(spw))

    model = xgb.XGBClassifier(
        max_depth=5,
        learning_rate=0.1,
        n_estimators=500,
        subsample=0.8,
        colsample_bytree=0.8,
        objective='binary:logistic',
        scale_pos_weight=spw,
        n_jobs=-1,
        random_state=42
    )
    
    # 3. Training
    print("Training XGBoost model with {0} training events...".format(len(X_train)))
    model.fit(
        X_train, y_train,
        sample_weight=w_train,
        eval_set=[(X_test, y_test)],
        sample_weight_eval_set=[w_test],
        eval_metric='logloss',
        early_stopping_rounds=20,
        verbose=True
    )
    
    # 4. Evaluation and Overtraining Check
    y_pred_test = model.predict_proba(X_test)[:, 1]
    y_pred_train = model.predict_proba(X_train)[:, 1]
    
    fpr, tpr, _ = roc_curve(y_test, y_pred_test, sample_weight=w_test)
    roc_auc = auc(fpr, tpr)
    print("Weighted ROC AUC Score: {0:.4f}".format(roc_auc))

    # --- Plot: Overtraining Check ---
    plt.figure(figsize=(8, 6))
    bins = np.linspace(0, 1, 40)
    
    for label, color, name in [(1, 'red', 'Signal'), (0, 'blue', 'Background')]:
        # Test (Points)
        mask_test = (y_test == label)
        hist_test, _ = np.histogram(y_pred_test[mask_test], bins=bins, weights=w_test[mask_test], density=True)
        raw_counts, _ = np.histogram(y_pred_test[mask_test], bins=bins)
        scale = hist_test / np.where(raw_counts > 0, raw_counts, 1)
        err = np.sqrt(raw_counts) * scale
        
        center = (bins[:-1] + bins[1:]) / 2
        plt.errorbar(center, hist_test, yerr=err, fmt='o', c=color, label='{0} (Test)'.format(name))
        
        # Train (Steps)
        mask_train = (y_train == label)
        plt.hist(y_pred_train[mask_train], bins=bins, weights=w_train[mask_train], 
                 histtype='step', fill=False, color=color, alpha=0.3, density=True, 
                 label='{0} (Train)'.format(name), linewidth=2)

    plt.xlabel('XGBoost Score')
    plt.ylabel('Arbitrary Units (Normalized)')
    plt.yscale('log')
    plt.title('Overtraining Check')
    plt.legend(loc='upper center')
    plt.grid(alpha=0.3)
    plt.savefig('plots/overtraining_check.png')
    plt.close()
    print("Saved: plots/overtraining_check.png")

    # --- Plot: ROC Curves ---
    plt.figure(figsize=(8, 8))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = {0:.2f})'.format(roc_auc))
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.savefig('plots/roc_curve.png')
    plt.close()

    # --- Feature Importance ---
    plt.figure(figsize=(10, 8))
    xgb.plot_importance(model, importance_type='gain', max_num_features=15)
    plt.title("Feature Importance (Gain)")
    plt.tight_layout()
    plt.savefig('plots/feature_importance_gain.png')
    plt.close()

    plt.figure(figsize=(10, 8))
    plot_importance(
        model, 
        importance_type='weight', # Specify 'weight' to rank by frequency
        title='Feature Importance (Weight)',
        xlabel='Weight'
    )
    plt.tight_layout()
    plt.savefig('plots/feature_importance_weight.png')
    plt.close() 

    # 5. Save Model
    model.save_model("bs_tautau_xgb_model.bin")
    print("Model saved as bs_tautau_xgb_model.bin")

if __name__ == "__main__":
    SIGNAL_FILE = "signalOnl.root"
#    SIGNAL_FILE = "test_bstautau_vtx.root"
    BKG_FILE = "test_bstautau_vtx.root"
    
    # -1 means use all events
    train(SIGNAL_FILE, BKG_FILE, n_sig=-1, n_bkg=-1)

