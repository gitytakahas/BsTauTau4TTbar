import ROOT

def get_genEventSumw(file_path):
    """Retrieve the sum of genEventSumw from the Runs tree in a ROOT file."""

    print(file_path)

    f = ROOT.TFile.Open(file_path)
    runs_tree = f.Get("Runs")
    if not runs_tree:
        raise RuntimeError("No Runs tree found in file", file_path)

    sumw = 0
    sumcnt = 0
    for entry in runs_tree:
        sumw += entry.genEventSumw
        sumcnt += entry.genEventCount

    print(">> The sumw = ",sumw, '. event count = ', sumcnt)
    f.Close()
    return sumw

if __name__ == "__main__":
    get_genEventSumw("all_Skim_bstautau.root")
    get_genEventSumw("all_Skim_ttsemeleptonic.root")

    
