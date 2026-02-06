import ROOT
import os
import math

# Use English for comments in the code.

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

try:
    from officialStyle import officialStyle
except ImportError:
    def officialStyle(style):
        pass

# --- ROOT Global Settings ---
ROOT.gROOT.SetBatch(True)
officialStyle(ROOT.gStyle)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

isMVA = True
prefix = ""
if isMVA: prefix = "applied_"

produce = True

# Define the variable to be used for combine fit
#FINAL_VARIABLE = "m_exact"
FINAL_VARIABLE = "mva_score"
BIN_NAME = "bin1"  # Name of the directory in ROOT file and bin in datacard

base = "pair_pt!=-1"

if isMVA:
#    base += "&&mva_score > 0.7"
#    base += "&& mva_score > 0.1"
    base += "&&1"

weight = "L1PreFiringWeight_Nom*genWeight*puWeight"
lumi = 59.7 # 2018 only

# --- 1. Configuration ---
# Format: (file_path, cut_string, weight_string, label, color, style, cross_section, process_name)
DATA_SAMPLES = [
#    (prefix + "test_bstautau.root",      base + " && is_true_signal==1", weight, "Sig (Sig File add)", ROOT.kRed+1,  1, 830*0.16, "signal"),
#    (prefix + "test_bstautau_vtx.root",      base + " && is_true_signal==1", weight, "Sig (Sig File add 2)", ROOT.kBlue+1,  1, 830*0.16, "signal2"),
#    (prefix + "final_bstautau_vtx.root",      base + " && is_true_signal == 1", weight, "Sig (Match)", ROOT.kRed+1,  1, 830*0.16, "signal"),
#    (prefix + "final_bstautau_vtx.root",      base + " && is_true_signal != 1", weight, "Sig (No-Match)", ROOT.kBlue+1, 1, 830*0.16, "sigfile_truth"),
    (prefix + "final_bstautau_vtx.root",      base + " && 1", weight, "Sig File", ROOT.kBlue+1, 1, 830*0.16, "signal"),
    (prefix + "final_ttsemileptonic_vtx.root", base + " && is_true_signal != 1", weight, "Bkg (Bkg File)", ROOT.kGreen+2, 2, 366.29, "bkg_ttbar"),
]

# --- 2. Variables ---
# Format: (name, nbins, xmin, xmax, xtitle, isLog)
VAR_LIST = [
    ("m_vis", 50, 0, 10, "m_{vis} [GeV]", False),
    ("m_exact", 30, 3.5, 8, "m_{exact} [GeV]", False),
    ("x1", 50, 0, 1.2, "x_{1}", False),
    ("x2", 50, 0, 1.2, "x_{2}", False),
    ("alpha1", 50, 0, 0.5, "#alpha_{1} [rad]", False),
    ("alpha2", 50, 0, 0.5, "#alpha_{2} [rad]", False),
    ("pair_pt", 50, 0, 150, "B_{s} p_{T} [GeV]", False),
    ("pair_eta", 50, -3, 3, "B_{s} #eta", False),
    ("pair_phi", 50, -math.pi, math.pi, "B_{s} #eta", False),
    ("n_extra_trks", 11, -0.5, 10.5, "Number of extra tracks", False),
    ("iso_ratio", 50, 0, 1.2, "Isolation ratio", False),
    ("dr_taus", 50, 0, 0.8, "#DeltaR(#tau_{1}, #tau_{2})", False),
    ("tau1_m_rho", 50, 0, 2.0, "#tau_{1} m_{#rho} [GeV]", False),
    ("tau2_m_rho", 50, 0, 2.0, "#tau_{2} m_{#rho} [GeV]", False),
    ("tau1_vtxProb", 50, 0, 1.0, "#tau_{1} vtxProb", False),
    ("tau2_vtxProb", 50, 0, 1.0, "#tau_{2} vtxProb", False),
    ("tau1_lxySig", 50, 0, 50, "#tau_{1} L_{xy} sig", False),
    ("tau2_lxySig", 50, 0, 50, "#tau_{2} L_{xy} sig", False),
    ("j_ParTRawTauhtaue", 50, 0, 1.0, "ParT: e", False),
    ("j_ParTRawTauhtaumu", 50, 0, 1.0, "ParT: #mu", False),
    ("j_ParTRawTauhtauh", 50, 0, 1.0, "ParT: h", False),
    ("j_ParTRawSingletau", 50, 0, 1.0, "ParT: SingleTau", False),
    ("tau1_deltaPV_dr", 50, -2.0, 2.0, "Tau1 PV: Delta(r)", False),
    ("tau2_deltaPV_dr", 50, -2.0, 2.0, "Tau2 PV: Delta(r)", False),
    ("tau1_mass", 50, 0.4, 1.75, "Tau1 mass", False),
    ("tau2_mass", 50, 0.4, 1.75, "Tau2 mass", False),
    ("tau1_pt", 50, 0, 60, "Tau1 pt", False),
    ("tau1_eta", 50, -2.4, 2.4, "Tau1 eta", False),
    ("tau1_phi", 50, -math.pi, math.pi, "Tau1 phi", False),
    ("tau2_pt", 50, 0, 60, "Tau2 pt", False),
    ("tau2_eta", 50, -2.4, 2.4, "Tau2 eta", False),
    ("tau2_phi", 50, -math.pi, math.pi, "Tau2 phi", False),
    ("BsTau_jetNCharged", 30, 0, 30, "NCharged", False),
    ("BsTau_jetNNeutral", 30, 0, 30, "NNeutral", False),
    ("j_deepflavB", 30, 0, 1., "deep flav. B", False),
]

if isMVA:
    # Added mva_score with isLog = True
    VAR_LIST.append(("mva_score", 50, 0, 1.0, "XGBS score", True))

def write_datacard(filename, shapes_file, processes, rates):
    """Generates a CMS Combine datacard."""
    n_bkg = len(processes) - 1
    with open(filename, 'w') as f:
        f.write("imax 1  number of channels\n")
        f.write("jmax {}  number of backgrounds\n".format(n_bkg))
        f.write("kmax * number of nuisance parameters\n")
        f.write("-" * 60 + "\n")
        f.write("shapes * * {} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n".format(os.path.basename(shapes_file)))
        f.write("-" * 60 + "\n")
        f.write("bin {}\n".format(BIN_NAME))
        f.write("observation -1\n")
        f.write("-" * 60 + "\n")
        
        f.write("{:<15} {:<15}".format("bin", BIN_NAME))
        for _ in range(n_bkg): f.write(" {:<15}".format(BIN_NAME))
        f.write("\n")

        f.write("{:<15}".format("process"))
        for p in processes: f.write(" {:<15}".format(p))
        f.write("\n")

        f.write("{:<15}".format("process"))
        for i in range(len(processes)): f.write(" {:<15}".format(i))
        f.write("\n")

        f.write("{:<15}".format("rate"))
        for r in rates: f.write(" {:<15.4f}".format(r))
        f.write("\n")
        f.write("-" * 60 + "\n")
        
        f.write("{:<12} {:<6} {:<10}".format("lumi", "lnN", "1.025"))
        for _ in range(n_bkg): f.write(" {:<10}".format("-"))
        f.write("\n")

def run_plotting():
    ensureDir('plots')
    ensureDir('combine')

    for var_name, nbins, xmin, xmax, xtitle, isLog in VAR_LIST:
        f_combine = None
        combine_dir = None
        combine_processes = []
        combine_rates = []
        h_data_obs = None
        
        if produce and var_name == FINAL_VARIABLE:
            shapes_path = "combine/shapes.root"
            f_combine = ROOT.TFile(shapes_path, "RECREATE")
            combine_dir = f_combine.mkdir(BIN_NAME)

        canvas = ROOT.TCanvas("c_{}".format(var_name), "", 800, 700)
        canvas.SetLeftMargin(0.15)
        # Apply log scale if option is True
        if isLog:
            canvas.SetLogy(1)

        leg = ROOT.TLegend(0.45, 0.70, 0.90, 0.88)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.03)

        hists_for_canvas = []
        max_y = 0

        for i, (fpath, cut, weight, label, color, style, crosssection, proc_name) in enumerate(DATA_SAMPLES):
            f_in = ROOT.TFile.Open(fpath)
            if not f_in or f_in.IsZombie():
                print("Warning: Could not open {}. Skipping...".format(fpath))
                continue
            
            h_norm = f_in.Get("h_genEventSumw")
            genEventSumw = h_norm.GetBinContent(1)

            print('genEventSumw = ', genEventSumw)

            tree = f_in.Get("tree")
            h_temp = ROOT.TH1F("h_temp_{}_{}".format(var_name, i), "", nbins, xmin, xmax)
            tree.Project(h_temp.GetName(), var_name, '('+cut+')*'+weight)

            h_temp.GetXaxis().SetTitle(xtitle)
            # Apply scale factor
            scale = crosssection*lumi*1000./genEventSumw
            h_temp.Scale(scale)

            # Combine processing
            if f_combine:
                combine_dir.cd()
                h_comb = h_temp.Clone(proc_name)
                h_comb.Write()
                combine_processes.append(proc_name)
                combine_rates.append(h_temp.Integral())

                # Create data_obs (sum of all processes)
                if h_data_obs is None:
                    h_data_obs = h_temp.Clone("data_obs")
                    h_data_obs.SetDirectory(0)
                else:
                    h_data_obs.Add(h_temp)

            # Plotting preparation (copy for normalization)
            h_plot = h_temp.Clone(h_temp.GetName() + "_plot")
            h_plot.SetDirectory(0)
            h_plot.SetLineColor(color)
            h_plot.SetLineWidth(3)
            h_plot.SetLineStyle(style)

            raw_entries = h_temp.GetSumOfWeights()
            if h_plot.Integral() > 0:
                h_plot.Scale(1.0 / h_plot.Integral())
            
            leg.AddEntry(h_plot, "{} ({:.0f})".format(label, raw_entries), "l")
            if h_plot.GetMaximum() > max_y: max_y = h_plot.GetMaximum()
            hists_for_canvas.append(h_plot)
            
            f_in.Close()

        # Save data_obs if generated
        if f_combine and h_data_obs:
            combine_dir.cd()
            h_data_obs.Write()
            f_combine.Close()
            write_datacard("combine/datacard.txt", "shapes.root", combine_processes, combine_rates)
            print(">>> Generated combine/shapes.root (with data_obs) and combine/datacard.txt")

        # Drawing
        canvas.cd()
        for i, h in enumerate(hists_for_canvas):
            if isLog:
                h.SetMaximum(max_y * 10.0) # More margin for log scale
                h.SetMinimum(1e-4)        # Avoid 0 in log scale
            else:
                h.SetMaximum(max_y * 1.4)
                h.SetMinimum(0)

            h.Draw("HIST" if i == 0 else "HIST SAME")
        leg.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.15, 0.92, "#bf{CMS} #it{Internal}")
        canvas.Print("plots/{}_multi_comp.png".format(var_name))

if __name__ == "__main__":
    run_plotting()

