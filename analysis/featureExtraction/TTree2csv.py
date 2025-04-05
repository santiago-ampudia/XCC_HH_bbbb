import ROOT

# open the file
f = ROOT.TFile.Open("kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBttHHbbbbESpreadDurham1034BSplitSampleN.root")

# get the tree
t = f.Get("TreeBtt")

# create a csv file
csv_file = open("pyTorchAnalysis/Btt.csv", "w")

# define header names (in the desired order) excluding the skipped branches
header_names = [
    "aplanarity",
    "invMassB1",
    "invMassB2",
    "minJetM",
    "sphericity",
    "cosThetaB1",
    "cosThetaB2",
    "cosThetaB3",
    "cosThetaB4",
    "sumPt",
    "jetB1Pt",
    "jetB2Pt",
    "jetB3Pt",
    "jetB4Pt",
    "jetB1M",
    "jetB2M",
    "jetB3M",
    "jetB4M",
    "jetNObjects",
    "minJetNObjects",
    "invMassB1AntiKt",
    "invMassB2AntiKt",
    "nJetsAntiKt",
    "invMassB11Best",
    "invMassB21Best",
    "invMassB12Best",
    "invMassB22Best",
    "invMassB13Best",
    "invMassB23Best",
    "invMassB14Best",
    "invMassB24Best",
    "invMassB15Best",
    "invMassB25Best",
    "invMassB16Best",
    "invMassB26Best",
    "invMassB17Best",
    "invMassB27Best",
    "invMassB18Best",
    "invMassB28Best",
    "distanceZ1MinChiSquaredZZMass",
    "distanceZ2MinChiSquaredZZMass",
    "exclYmerge12",
    "exclYmerge23",
    "exclYmerge34",
    "exclYmerge45",
    "exclYmerge56",
    "invMassZZ1",
    "invMassZZ2",
    "thrust",
    "boostB1",
    "boostB2",
    "boostB3",
    "boostB4",
    "boostSystem",
    "missingET",
    "invMass4Jets",
    "deltaRJetPairs",
    "pTAssymetry",
    "jetPairMassDelta",
    "invMassB1FitBest",
    "invMassB2FitBest",
    "chi2ndfBest",
    "jetPairMassDeltaFit"
]

# write the header row to the CSV file
csv_file.write(", ".join(header_names) + "\n")

# loop over all events
for i in range(t.GetEntries()):
    t.GetEntry(i)
    # retrieve the leaf values in the order defined by header_names
    leaves = []
    for name in header_names:
        leaves.append(t.GetLeaf(name).GetValue())
    # write line to csv file
    csv_file.write(", ".join(map(str, leaves)) + "\n")

# close the csv file
csv_file.close()
