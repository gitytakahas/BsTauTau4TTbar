import ROOT
import numpy as np
import sys
import os
from array import array
import glob
import math
import argparse

def deltaR(eta1, phi1, eta2, phi2):
    dphi = phi1 - phi2
    while dphi >  math.pi: dphi -= 2*math.pi
    while dphi < -math.pi: dphi += 2*math.pi
    return math.sqrt((eta1 - eta2)**2 + dphi**2)


# --- 1. Configuration ---
TAU_CONFIG = [
    ("BsTau_pt", 'F'), ("BsTau_eta", 'F'), ("BsTau_phi", 'F'),
    ("BsTau_trkIdx1", 'I'), ("BsTau_trkIdx2", 'I'), ("BsTau_trkIdx3", 'I'),
    ("BsTau_pion1_charge", 'I'), ("BsTau_pion2_charge", 'I'), ("BsTau_pion3_charge", 'I'),
    ("BsTau_pion1_pt", 'F'), ("BsTau_pion2_pt", 'F'), ("BsTau_pion3_pt", 'F'),
    ("BsTau_pion1_eta", 'F'), ("BsTau_pion2_eta", 'F'), ("BsTau_pion3_eta", 'F'),
    ("BsTau_pion1_phi", 'F'), ("BsTau_pion2_phi", 'F'), ("BsTau_pion3_phi", 'F'),
    ("BsTau_vtxProb", 'F'), ("BsTau_flightLen", 'F'), ("BsTau_lxySig", 'F'), 
    ("BsTau_nExtra", 'I'), ("BsTau_mass", 'F'),
    ("BsTau_genMatchId", 'I'), 
    ("BsTau_isFromBsTau", 'I'),
    ("BsTau_charge", 'I'),
    ("BsTau_fitMass", 'F'),
    ("BsTau_fitPt", 'F'),
    ("BsTau_maxDoca", 'F'),
    ("BsTau_minDoca", 'F'),
    ("BsTau_pvips", 'F'),
    ("BsTau_deltaChi2", 'F'),
]


J_BRANCH_LIST = [
    "j_ParTRawB", "j_ParTRawC", "j_ParTRawOther", 
    "j_ParTRawSingletau", "j_ParTRawTauhtaue", "j_ParTRawTauhtauh", "j_ParTRawTauhtaumu",
    "j_deepflavB",
    "BsTau_jetNCharged", "BsTau_jetNNeutral", "BsTau_deepFlavBB",
    "BsTau_jetEnergyFraction", "BsTau_jetPt", "BsTau_jetEta", "BsTau_jetPhi", "BsTau_jetMass"
]

def get_genEventSumw(file_path):
    """Retrieve the sum of genEventSumw from the Runs tree in a ROOT file."""
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie(): return 1.0
    runs_tree = f.Get("Runs")
    if not runs_tree:
        f.Close()
        return 1.0
    sumw = 0
    for entry in runs_tree:
        sumw += entry.genEventSumw
    f.Close()
    return float(sumw)

def calculate_reconstruction(p_vis, e_vis, sv, pv):
    dx, dy, dz = sv['x']-pv['x'], sv['y']-pv['y'], sv['z']-pv['z']
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    if dist < 1e-5: return None
    nx, ny, nz = dx/dist, dy/dist, dz/dist
    px, py, pz = p_vis['px'], p_vis['py'], p_vis['pz']
    p_vis_vec = np.array([px, py, pz])
    p_vis_mag = np.linalg.norm(p_vis_vec)
    p_vis_n = px*nx + py*ny + pz*nz
    
    m_tau, m_vis2 = 1.77686, max(0, e_vis**2 - p_vis_mag**2)
    denom = 2.0 * (e_vis - p_vis_n)
    if abs(denom) < 1e-8: return None
    
    p_nu = (m_tau**2 - m_vis2) / denom
    if p_nu < 0: return None
    
    p_full_mag = p_vis_n + p_nu
    x = p_vis_n / p_full_mag if p_full_mag > 0 else -1
    alpha = np.arccos(np.clip(p_vis_n / p_vis_mag, -1.0, 1.0)) if p_vis_mag > 0 else 0
    
    return {
        'x': x, 
        'alpha': alpha, 
        'p_full': p_vis_vec + (p_nu * np.array([nx, ny, nz])), 
        'e_full': e_vis + p_nu
    }

def analyze_and_save(input_path, output_file):

# --- 1. Identify Input Files ---
    input_files = []
    
    if os.path.isdir(input_path):
        input_files = glob.glob(os.path.join(input_path, "*.root"))
    else:
        input_files = glob.glob(input_path)

    if not input_files:
        print("Error: No ROOT files found for input: {0}".format(input_path))
        return

    # --- 2. Setup TChain ---
    chain = ROOT.TChain("Events")
    for f in input_files:
        print("Adding to chain: {0}".format(f))
        chain.Add(f)
    
    total_entries = chain.GetEntries()

    # --- 3. SumW calculation (from all matched files) ---
    total_sumw = 0
    for f_path in input_files:
        total_sumw += get_genEventSumw(f_path)
    
    print("\nTotal Files: {0}, Total Entries: {1}, Total SumW: {2}".format(
        len(input_files), total_entries, total_sumw))



    f_out = ROOT.TFile(output_file, "RECREATE")
    new_tree = ROOT.TTree("tree", "BsTauTau_Full_Reconstruction")


#    if not os.path.exists(input_file):
#        print("Error: {0} not found.".format(input_file))
#        return
#
#    total_sumw = get_genEventSumw(input_file)
#    print("Opening input: {0} (Total SumW: {1})".format(input_file, total_sumw))
#        
#    f_in = ROOT.TFile.Open(input_file)
#    tree = f_in.Get("Events")
#    total_entries = tree.GetEntries()
#
#    f_out = ROOT.TFile(output_file, "RECREATE")
#    new_tree = ROOT.TTree("tree", "BsTauTau_Full_Reconstruction")

    # --- Setup Buffers ---
    m_vis, m_exact = array('f', [0.0]), array('f', [0.0])
    x1, x2, alpha1, alpha2 = array('f', [0.0]), array('f', [0.0]), array('f', [0.0]), array('f', [0.0])
    pair_pt, pair_eta, pair_phi = array('f', [0.0]), array('f', [0.0]), array('f', [0.0])
    n_extra_trks = array('i', [0])
    iso_ratio, dr_taus = array('f', [0.0]), array('f', [0.0])
    tau1_m_rho, tau2_m_rho = array('f', [0.0]), array('f', [0.0])
    is_true_signal = array('i', [0])
    puWeight, L1PreFiringWeight_Nom, genWeight = array('f', [1.0]), array('f', [1.0]), array('f', [1.0])
    mu1_pt, mu1_eta = array('f',[-1.0]), array('f',[-1.0])
    nj, nbj = array('i',[-1]), array('i',[-1])
    PuppiMET_pt = array('f',[-1.0])
    

    new_tree.Branch("m_vis", m_vis, "m_vis/F")
    new_tree.Branch("m_exact", m_exact, "m_exact/F")
    new_tree.Branch("x1", x1, "x1/F"); new_tree.Branch("x2", x2, "x2/F")
    new_tree.Branch("alpha1", alpha1, "alpha1/F"); new_tree.Branch("alpha2", alpha2, "alpha2/F")
    new_tree.Branch("pair_pt", pair_pt, "pair_pt/F")
    new_tree.Branch("pair_eta", pair_eta, "pair_eta/F")
    new_tree.Branch("pair_phi", pair_phi, "pair_phi/F")
    new_tree.Branch("n_extra_trks", n_extra_trks, "n_extra_trks/I")
    new_tree.Branch("iso_ratio", iso_ratio, "iso_ratio/F")
    new_tree.Branch("dr_taus", dr_taus, "dr_taus/F")
    new_tree.Branch("tau1_m_rho", tau1_m_rho, "tau1_m_rho/F")
    new_tree.Branch("tau2_m_rho", tau2_m_rho, "tau2_m_rho/F")
    new_tree.Branch("is_true_signal", is_true_signal, "is_true_signal/I")
    new_tree.Branch("puWeight", puWeight, "puWeight/F")
    new_tree.Branch("L1PreFiringWeight_Nom", L1PreFiringWeight_Nom, "L1PreFiringWeight_Nom/F")
    new_tree.Branch("genWeight", genWeight, "genWeight/F")
    new_tree.Branch("mu1_pt", mu1_pt, "mu1_pt/F")
    new_tree.Branch("mu1_eta", mu1_eta, "mu1_eta/F")
    new_tree.Branch("nj", nj, "nj/I")
    new_tree.Branch("nbj", nbj, "nbj/I")
    new_tree.Branch("PuppiMET_pt", PuppiMET_pt, "PuppiMET_pt/F")

    outputs = {}
    for prefix in ["tau1_", "tau2_"]:
        v = "{0}deltaPV_dr".format(prefix)
        outputs[v] = array('f', [0.0])
        new_tree.Branch(v, outputs[v], "{0}/F".format(v))
        
        for name, type_char in TAU_CONFIG:
            s = name.replace("BsTau_", "")
            key = prefix + s
            if type_char == 'I':
                outputs[key] = array('i', [0])
                new_tree.Branch(key, outputs[key], "{0}/I".format(key))
            else:
                outputs[key] = array('f', [0.0])
                new_tree.Branch(key, outputs[key], "{0}/F".format(key))

    for bname in J_BRANCH_LIST:
        outputs[bname] = array('f', [0.0])
        new_tree.Branch(bname, outputs[bname], "{0}/F".format(bname))

    # --- Processing ---
    for i in range(total_entries):
#        tree.GetEntry(i)
        chain.GetEntry(i)
        if i % 10000 == 0:
            sys.stdout.write("\rProgress: {0:.1f}%".format(100.0 * i / total_entries))
            sys.stdout.flush()

        if len(chain.BsTau_pt) < 2: continue

        
        # Group tau candidates by jet using original rounding logic
        taus_by_jet = {}
        for it in range(len(chain.BsTau_pt)):
            if chain.BsTau_jetPt[it] < 0: continue
            # Use jet kinematics as key to identify tau candidates belonging to the same jet
            key = (round(chain.BsTau_jetPt[it], 3), round(chain.BsTau_jetEta[it], 3))
            if key not in taus_by_jet: taus_by_jet[key] = []
            taus_by_jet[key].append(it)

        # Pair building within each jet
        for key, indices in taus_by_jet.items():
#	    print("check", key, indices, len(indices))

            if len(indices) < 2: continue
            indices_sorted = sorted(indices, key=lambda idx: chain.BsTau_vtxProb[idx], reverse=True)
            
            # Sort by vertex probability (descending)           
            
            i1, i2, found_pair = -1, -1, False
            for idx_a in range(len(indices_sorted)):
                for idx_b in range(idx_a + 1, len(indices_sorted)):
                    ia, ib = indices_sorted[idx_a], indices_sorted[idx_b]
                    
                    # Check for track overlap using trkIdx
                    t1 = {chain.BsTau_trkIdx1[ia], chain.BsTau_trkIdx2[ia], chain.BsTau_trkIdx3[ia]}
                    t2 = {chain.BsTau_trkIdx1[ib], chain.BsTau_trkIdx2[ib], chain.BsTau_trkIdx3[ib]}
#                    print("This is check:", t1, t2, chain.BsTau_charge[ia], chain.BsTau_charge[ib])
 
                    # OS charge requirement and ensure no shared tracks
                    if t1.isdisjoint(t2) and int(chain.BsTau_charge[ia] + chain.BsTau_charge[ib]) == 0:
                        i1, i2, found_pair = ia, ib, True
#                        print('pair found!')
                        break
                if found_pair: break
            
            if not found_pair: continue
	
#	    print("pair found !!!!!!!!!!!!!!!!")

            # Match jet branch index for the selected tau pair
            matched_j_idx = -1
            for j_idx in range(len(chain.j_pt)):

                dr = deltaR(chain.j_eta[j_idx], chain.j_phi[j_idx], chain.BsTau_jetEta[i1], chain.BsTau_jetPhi[i1])
                if dr < 0.1:
#                if abs(chain.j_pt[j_idx] - chain.BsTau_jetPt[i1]) < 0.01 and abs(chain.j_eta[j_idx] - chain.BsTau_jetEta[i1]) < 0.01:
                    matched_j_idx = j_idx; break
            if matched_j_idx == -1: 
                print('no matching jet found ...')
                continue

            # Kinematic reconstruction for tau1 and tau2
            def get_kin_full(it):
                pt, eta, phi, m = chain.BsTau_pt[it], chain.BsTau_eta[it], chain.BsTau_phi[it], chain.BsTau_mass[it]
                px, py, pz = pt*np.cos(phi), pt*np.sin(phi), pt*np.sinh(eta)
                e = np.sqrt(px**2 + py**2 + pz**2 + m**2)
                
                # Best m_rho candidate from 3 pions
                m_rho_best, min_diff = -1.0, 999.0
                pions = []
                for idx in [1, 2, 3]:
                    p = ROOT.TLorentzVector()
                    p.SetPtEtaPhiM(getattr(chain, "BsTau_pion{0}_pt".format(idx))[it], 
                                  getattr(chain, "BsTau_pion{0}_eta".format(idx))[it], 
                                  getattr(chain, "BsTau_pion{0}_phi".format(idx))[it], 0.13957)
                    pions.append({'p': p, 'q': getattr(chain, "BsTau_pion{0}_charge".format(idx))[it]})
                
                for a in range(3):
                    for b in range(a+1, 3):
                        if pions[a]['q'] + pions[b]['q'] == 0:
                            m_rho_cand = (pions[a]['p'] + pions[b]['p']).M()
                            if abs(m_rho_cand - 0.770) < min_diff:
                                min_diff, m_rho_best = abs(m_rho_cand - 0.770), m_rho_cand
                
                res = calculate_reconstruction({'px':px, 'py':py, 'pz':pz}, e, 
                                              {'x':chain.BsTau_svX[it], 'y':chain.BsTau_svY[it], 'z':chain.BsTau_svZ[it]}, 
                                              {'x':chain.BsTau_pvX[it], 'y':chain.BsTau_pvY[it], 'z':chain.BsTau_pvZ[it]})
                return res, px, py, pz, e, eta, phi, m_rho_best


#            print("this is good event")
            res1, px1, py1, pz1, e1, eta1, phi1, rho1 = get_kin_full(i1)
            res2, px2, py2, pz2, e2, eta2, phi2, rho2 = get_kin_full(i2)
            if not res1 or not res2: continue

#            print("this is good event2")
            # Fill event weights
            genWeight[0] = float(getattr(chain, "genWeight", 1.0))
            puWeight[0] = float(getattr(chain, "puWeight", 1.0))
            L1PreFiringWeight_Nom[0] = float(getattr(chain, "L1PreFiringWeight_Nom", 1.0))
            mu1_pt[0] = float(getattr(chain, "mu1_pt", -1.0))
            mu1_eta[0] = float(getattr(chain, "mu1_eta", -1.0))
            nj[0] = int(getattr(chain, "nj", -1.0))
            nbj[0] = int(getattr(chain, "nbj", -1.0))
            PuppiMET_pt[0] = float(getattr(chain, "PuppiMET_pt", -1.0))

            # Visible and Rho mass
            tau1_m_rho[0], tau2_m_rho[0] = float(rho1), float(rho2)
            x1[0], x2[0] = float(res1['x']), float(res2['x'])
            alpha1[0], alpha2[0] = float(res1['alpha']), float(res2['alpha'])
            m_vis[0] = float(np.sqrt(max(0, (e1+e2)**2 - ((px1+px2)**2 + (py1+py2)**2 + (pz1+pz2)**2))))
            
            # Full mass reconstruction (collinear approximation)
            if 0 < res1['x'] < 1 and 0 < res2['x'] < 1:
                p_sum = res1['p_full'] + res2['p_full']
                e_sum = res1['e_full'] + res2['e_full']
                m_exact[0] = float(np.sqrt(max(0, e_sum**2 - np.sum(p_sum**2))))
                
                v_pair = ROOT.TLorentzVector()
                v_pair.SetPxPyPzE(p_sum[0], p_sum[1], p_sum[2], e_sum)
                pair_pt[0], pair_eta[0], pair_phi[0] = float(v_pair.Pt()), float(v_pair.Eta()), float(v_pair.Phi())
            else: 
                m_exact[0] = -1.0; pair_pt[0] = -1.0; pair_eta[0] = 0.0; pair_phi[0] = 0.0

            # Isolation and geometry
            selected_trks = {chain.BsTau_trkIdx1[i1], chain.BsTau_trkIdx2[i1], chain.BsTau_trkIdx3[i1], 
                             chain.BsTau_trkIdx1[i2], chain.BsTau_trkIdx2[i2], chain.BsTau_trkIdx3[i2]}
            # Count total tracks in jet minus the 6 tracks belonging to the tau pair
            all_tau_trks = set()
            for idx in indices:
                all_tau_trks.update([chain.BsTau_trkIdx1[idx], chain.BsTau_trkIdx2[idx], chain.BsTau_trkIdx3[idx]])
            n_extra_trks[0] = int(len(all_tau_trks - selected_trks))
            
            #            iso_ratio[0] = float(np.sqrt((px1+px2)**2 + (py1+py2)**2) / chain.j_pt[matched_j_idx])
#            print('jet1, jet2=', chain.BsTau_jetPt[i1], chain.BsTau_jetPt[i2])
            iso_ratio[0] = float(np.sqrt((px1+px2)**2 + (py1+py2)**2) / chain.BsTau_jetPt[i1])

            dr_taus[0] = float(np.sqrt((eta1 - eta2)**2 + (ROOT.TVector2.Phi_mpi_pi(phi1 - phi2))**2))
            is_true_signal[0] = int(1 if (chain.BsTau_isFromBsTau[i1] and chain.BsTau_isFromBsTau[i2]) else 0)

            # Store individual tau data
            for p, idx in [("tau1_", i1), ("tau2_", i2)]:
                dx, dy, dz = chain.BsTau_pvX[idx]-chain.BsTau_svX[idx], chain.BsTau_pvY[idx]-chain.BsTau_svY[idx], chain.BsTau_pvZ[idx]-chain.BsTau_svZ[idx]
                outputs[p+"deltaPV_dr"][0] = float(np.sqrt(dx**2 + dy**2 + dz**2))
                
                for name, type_char in TAU_CONFIG:
                    s = name.replace("BsTau_", "")
                    val = getattr(chain, name)[idx]
                    outputs[p+s][0] = int(val) if type_char == 'I' else float(val)

            # Fill jet-related branches
            for b in J_BRANCH_LIST:
                val = -999.0
                try:
                    data = getattr(chain, b)
                    target_idx = i1 if b.startswith("BsTau_jet") else matched_j_idx
                    if target_idx < len(data): val = data[target_idx]
                except Exception: pass
                outputs[b][0] = float(val if val is not None else -999.0)

            new_tree.Fill()

    # Finalize files
    h_sumw = ROOT.TH1D("h_genEventSumw", "Total GenEventSumw", 1, 0, 1)
    h_sumw.SetBinContent(1, total_sumw)
    h_sumw.Write()
    f_out.Write()
    f_out.Close()
    print("\nProcessing completed successfully.")

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Process ROOT files.")
    parser.add_argument("-i", "--input", required=True, 
                        help="Input path: can be a directory, a single file, or a wildcard (e.g. 'data/*.root')")
    parser.add_argument("-o", "--output", required=True, 
                        help="Output file absolute path")

  
    args = parser.parse_args()
    
    analyze_and_save(args.input, args.output)
    
#    analyze_and_save("all_Skim.root", "sig_reco.root")


