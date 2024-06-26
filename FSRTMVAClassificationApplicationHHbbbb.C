 #include <cstdlib>
 #include <vector>
 #include <iostream>
 #include <map>
 #include <string>
 #include <map>
 #include <algorithm>
 #include <random>
 
 #include "TFile.h"
 #include "TTree.h"
 #include "TString.h"
 #include "TSystem.h"
 #include "TROOT.h"
 #include "TRandom.h"
 #include "TStopwatch.h"
 
 #include "TMVA/Tools.h"
 #include "TMVA/Reader.h"
 #include "TMVA/MethodCuts.h"
 
 
 using namespace TMVA;
 
 void FSRTMVAClassificationApplicationHHbbbbHelper( TString myMethodList, int topology, int nBack, vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, int& sizeTree, int nbin, TH1F& histOutput, TH1F& histOutputW, double weight, string inputMethod, TString NNVars, string rtdCut)
 {
	if(nBack == 0) cout<<"Helper 0"<<endl;
 	if(nBack == 1) cout<<"Helper 1"<<endl;
 	if(nBack == 2) cout<<"Helper 2"<<endl;
 	if(nBack == 3) cout<<"Helper 3"<<endl;
 	
 	//---------------------------------------------------------------
    // This loads the library
    TMVA::Tools::Instance();
 
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
 
    Use[inputMethod] = 1;
    /*// Cut optimisation
    Use["Cuts"]            = 0;
    Use["CutsD"]           = 0;
    Use["CutsPCA"]         = 0;
    Use["CutsGA"]          = 0;
    Use["CutsSA"]          = 0;
    //
    // 1-dimensional likelihood ("naive Bayes estimator")
    Use["Likelihood"]      = 0;
    Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
    Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
    Use["LikelihoodKDE"]   = 0;
    Use["LikelihoodMIX"]   = 0;
    //
    // Mutidimensional likelihood and Nearest-Neighbour methods
    Use["PDERS"]           = 0;
    Use["PDERSD"]          = 0;
    Use["PDERSPCA"]        = 0;
    Use["PDEFoam"]         = 0;
    Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
    Use["KNN"]             = 0; // k-nearest neighbour method
    //
    // Linear Discriminant Analysis
    Use["LD"]              = 0; // Linear Discriminant identical to Fisher
    Use["Fisher"]          = 0;
    Use["FisherG"]         = 0;
    Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
    Use["HMatrix"]         = 0;
    //
    // Function Discriminant analysis
    Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
    Use["FDA_SA"]          = 0;
    Use["FDA_MC"]          = 0;
    Use["FDA_MT"]          = 0;
    Use["FDA_GAMT"]        = 0;
    Use["FDA_MCMT"]        = 0;
    //
    // Neural Networks (all are feed-forward Multilayer Perceptrons)
    Use["MLP"]             = 0; // Recommended ANN
    Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
    Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
    Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
    Use["TMlpANN"]         = 0; // ROOT's own ANN
    Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
    Use["DNN_GPU"] = 0;         // Multi-core accelerated DNN.
    //
    // Support Vector Machine
    Use["SVM"]             = 0;
    //
    // Boosted Decision Trees
    Use["BDT"]             = 0; // uses Adaptive Boost
    Use["BDTG"]            = 0; // uses Gradient Boost
    Use["BDTB"]            = 0; // uses Bagging
    Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
    Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
    //
    // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
    Use["RuleFit"]         = 0;
    // ---------------------------------------------------------------
    Use["Plugin"]          = 0;
    Use["Category"]        = 0;
    Use["SVM_Gauss"]       = 0;
    Use["SVM_Poly"]        = 0;
    Use["SVM_Lin"]         = 0;*/
 
    std::cout << std::endl;
    std::cout << "==> Start TMVAClassificationApplication" << std::endl;
 
    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
       for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
 
       std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
       for (UInt_t i=0; i<mlist.size(); i++) {
          std::string regMethod(mlist[i]);
 
          if (Use.find(regMethod) == Use.end()) {
             std::cout << "Method \"" << regMethod
                       << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
             for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
                std::cout << it->first << " ";
             }
             std::cout << std::endl;
             return;
          }
          Use[regMethod] = 1;
       }
    }
 
    // --------------------------------------------------------------------------------------------------
 
    // Create the Reader object
 
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Float_t aplanarity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, invMassB1, invMassB2, jetB1M, jetB2M, jetB3M, jetB4M, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, minJetM, sphericity, sumPt, nParticles, minConstSize, jetNObjects, minJetNObjects;
    reader->AddVariable( "aplanarity", &aplanarity );
    reader->AddVariable( "cosThetaB1", &cosThetaB1 );
    reader->AddVariable( "cosThetaB2", &cosThetaB2 );
    reader->AddVariable( "cosThetaB3", &cosThetaB3 );
    reader->AddVariable( "cosThetaB4", &cosThetaB4 );
    reader->AddVariable( "invMassB1", &invMassB1 );
    reader->AddVariable( "invMassB2", &invMassB2 );
    reader->AddVariable( "jetB1M", &jetB1M );
    reader->AddVariable( "jetB2M", &jetB2M );
    reader->AddVariable( "jetB3M", &jetB3M );
    reader->AddVariable( "jetB4M", &jetB4M );
    reader->AddVariable( "jetB1Pt", &jetB1Pt );
    reader->AddVariable( "jetB2Pt", &jetB2Pt );
    reader->AddVariable( "jetB3Pt", &jetB3Pt );
    reader->AddVariable( "jetB4Pt", &jetB4Pt );
    reader->AddVariable( "minJetM", &minJetM );
    reader->AddVariable( "sphericity", &sphericity );
    reader->AddVariable( "sumPt", &sumPt );
    //reader->AddVariable( "nParticles", &nParticles );
    //reader->AddVariable( "minConstSize", &minConstSize );
    //reader->AddVariable( "jetNObjects", &jetNObjects );
    //reader->AddVariable( "minJetNObjects", &minJetNObjects );
 
 
    // Prepare input tree (this must be replaced by your data source)
    // in this example, there is a toy tree with signal and one with background events
    // we'll later on use only the "signal" events for the test in this example.
    //
    string inputSText = "outputTreeSHHbbbbESpreadDurham"+rtdCut+"Test.root";
    TFile* inputS = new TFile(inputSText.c_str());
    string inputBqqText = "outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
    TFile* inputBqq = new TFile(inputBqqText.c_str());
    string inputBttbarText = "outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Test.root";
    TFile* inputBtt = new TFile(inputBttbarText.c_str());
    string inputBZZText = "outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Test.root";
    TFile* inputBZZ = new TFile(inputBZZText.c_str());
    string inputBWWText = "outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Test.root";
    TFile* inputBWW = new TFile(inputBWWText.c_str());
    
    std::cout << "--- TMVAClassificationApplication   : Using input file for signal: " << inputS->GetName() << std::endl;
 
    // Event loop
 
    // Prepare the event tree
    // - Here the variable names have to corresponds to your tree
    // - You can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
    
    TTree* theTree;
    
    if(topology == 1) theTree = (TTree*)inputS->Get("TreeSTest");
    else if(topology == 2) theTree = (TTree*)inputBqq->Get("TreeBqqTest");
    else if(topology == 3) theTree = (TTree*)inputBtt->Get("TreeBttTest");
    else if(topology == 4) theTree = (TTree*)inputBZZ->Get("TreeBZZTest");
    else if(topology == 5) theTree = (TTree*)inputBWW->Get("TreeBWWTest");
    
    theTree->SetBranchAddress( "aplanarity", &aplanarity );
    theTree->SetBranchAddress( "cosThetaB1", &cosThetaB1 );
    theTree->SetBranchAddress( "cosThetaB2", &cosThetaB2 );
    theTree->SetBranchAddress( "cosThetaB3", &cosThetaB3 );
    theTree->SetBranchAddress( "cosThetaB4", &cosThetaB4 );
    theTree->SetBranchAddress( "invMassB1", &invMassB1 );
    theTree->SetBranchAddress( "invMassB2", &invMassB2 );
    theTree->SetBranchAddress( "jetB1M", &jetB1M );
    theTree->SetBranchAddress( "jetB2M", &jetB2M );
    theTree->SetBranchAddress( "jetB3M", &jetB3M );
    theTree->SetBranchAddress( "jetB4M", &jetB4M );
    theTree->SetBranchAddress( "jetB1Pt", &jetB1Pt );
    theTree->SetBranchAddress( "jetB2Pt", &jetB2Pt );
    theTree->SetBranchAddress( "jetB3Pt", &jetB3Pt );
    theTree->SetBranchAddress( "jetB4Pt", &jetB4Pt );
    theTree->SetBranchAddress( "minJetM", &minJetM );
    theTree->SetBranchAddress( "sphericity", &sphericity );
    theTree->SetBranchAddress( "sumPt", &sumPt );
    //theTree->SetBranchAddress( "nParticles", &nParticles );
    //theTree->SetBranchAddress( "minConstSize", &minConstSize );
    //theTree->SetBranchAddress( "jetNObjects", &jetNObjects );
    //theTree->SetBranchAddress( "minJetNObjects", &minJetNObjects );
    
    
    
    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;
 
    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
 
    if(topology == 1) std::cout << "--- Processing: " << theTree->GetEntries() << " HH events" << std::endl;
    if(topology == 2) std::cout << "--- Processing: " << theTree->GetEntries() << " qq events" << std::endl;
    if(topology == 3) std::cout << "--- Processing: " << theTree->GetEntries() << " ttbar events" << std::endl;
    if(topology == 4) std::cout << "--- Processing: " << theTree->GetEntries() << " ZZ events" << std::endl;
    if(topology == 5) std::cout << "--- Processing: " << theTree->GetEntries() << " WW events" << std::endl;
    
    TStopwatch sw;
    sw.Start();
   
	   if(nBack == 0) cout<<"Starting application of BDTqq!"<<endl;
	   if(nBack == 1) cout<<"Starting application of BDTttbar!"<<endl;
	   if(nBack == 2) cout<<"Starting application of BDTZZ!"<<endl;
	   if(nBack == 3) cout<<"Starting application of BDTWW!"<<endl;
	   
	    // Book the MVA methods
	    TString dir;
	    if(nBack == 0) dir    = "datasetqq" + rtdCut + NNVars + "/weights/";
	    else if(nBack == 1) dir    = "datasetttbar" + rtdCut + NNVars + "/weights/";
	    else if(nBack == 2) dir    = "datasetZZ" + rtdCut + NNVars + "/weights/";
	    else if(nBack == 3) dir    = "datasetWW" + rtdCut + NNVars + "/weights/";
	    TString prefix = "TMVAClassification";
	    // Book method(s)
	    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	       if (it->second) {
		  TString methodName = TString(it->first) + TString(" method");
		  TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
		  reader->BookMVA( methodName, weightfile );
	       }
	    }
    
	    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	 
	       //if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
	 
	       theTree->GetEntry(ievt);
	       
	       TString methodName = inputMethod + TString(" method");
	       
	       ///Printing for qualitative feeling
	       if(topology == 1 && ievt == 0) cout<<"POS 0: "<<reader->EvaluateMVA( methodName )<<endl;
	       if(topology == 1 && ievt == 70) cout<<"POS 70: "<<reader->EvaluateMVA( methodName )<<endl;
	       
	       ////Filling out arrays for plane of NN outputs
	       if(nBack == 0) BDTqqOutput.push_back(reader->EvaluateMVA(methodName)); 
	       if(nBack == 1) BDTttbarOutput.push_back(reader->EvaluateMVA(methodName));
	       if(nBack == 2) BDTZZOutput.push_back(reader->EvaluateMVA(methodName));
	       if(nBack == 3) BDTWWOutput.push_back(reader->EvaluateMVA(methodName));
	       
	       /////Filling out histograms 
	       histOutput.Fill(reader->EvaluateMVA(methodName));
	       histOutputW.Fill(reader->EvaluateMVA(methodName), weight);
	       
	       sizeTree = ievt+1;
	       
	    }
    
    
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();
 
    // Get efficiency for cuts classifier
    if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                 << " (for a required signal efficiency of " << effS << ")" << std::endl;
 
    if (Use["CutsGA"]) {
 
       // test: retrieve cuts for particular signal efficiency
       // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
       TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;
 
       if (mcuts) {
          std::vector<Double_t> cutsMin;
          std::vector<Double_t> cutsMax;
          mcuts->GetCuts( 0.7, cutsMin, cutsMax );
          std::cout << "--- -------------------------------------------------------------" << std::endl;
          std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
          for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
             std::cout << "... Cut: "
                       << cutsMin[ivar]
                       << " < \""
                       << mcuts->GetInputVar(ivar)
                       << "\" <= "
                       << cutsMax[ivar] << std::endl;
          }
          std::cout << "--- -------------------------------------------------------------" << std::endl;
       }
    }
 
    // Write histograms
 
    TFile *target  = new TFile( "TMVApp.root","RECREATE" );
    
    target->Close();
 
    //std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
 
    delete reader;
 
    std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
    
    /*
    if(topology == 1)
    {
	    if(nBack == 0)
	    {
	    	histBdtqqG->SetLineColor(kBlue);
	    	TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper1 = new TCanvas();
	   	cFSRTMVAClassificationApplicationHHbbbbHelper1->cd();
	    	histBdtqqG->Draw("HIST");
	    	histOutput.Add(histBdtqqG);
	    }
	    
	    if(nBack == 1)
	    {
		histBdtttbar->SetLineColor(kBlue);
		TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper2 = new TCanvas();
	    	cFSRTMVAClassificationApplicationHHbbbbHelper2->cd();
	    	histBdtttbar->Draw("HIST");
	    	histOutput.Add(histBdtttbarG);
	    }
    }	
    if(topology == 2)
    {
	    if(nBack == 0)
	    {
	    	histBdtqqG->SetLineColor(kRed+2);
	    	TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper3 = new TCanvas();
	   	cFSRTMVAClassificationApplicationHHbbbbHelper3->cd();
	    	histBdtqqG->Draw("HIST");
	    	histOutput.Add(histBdtqqG);
	    }
	    
	    if(nBack == 1)
	    {
		histBdtttbar->SetLineColor(kRed+2);
		TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper4 = new TCanvas();
	    	cFSRTMVAClassificationApplicationHHbbbbHelper4->cd();
	    	histBdtttbar->Draw("HIST");
	    	histOutput.Add(histBdtttbarG);
	    }
    }	
    if(topology == 3)
    {
	    if(nBack == 0)
	    {
	    	histBdtqq->SetLineColor(kPink+6);
	    	TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper5 = new TCanvas();
	   	cFSRTMVAClassificationApplicationHHbbbbHelper5->cd();
	    	histBdtqq->Draw("HIST");
	    	histOutput.Add(histBdtqqG);
	    }
	    
	    if(nBack == 1)
	    {
		histBdtttbar->SetLineColor(kPink+6);
		TCanvas *cFSRTMVAClassificationApplicationHHbbbbHelper6 = new TCanvas();
	    	cFSRTMVAClassificationApplicationHHbbbbHelper6->cd();
	    	histBdtttbar->Draw("HIST");
	    	histOutput.Add(histBdtttbarG);
	    }
    }
    */	
 	
 }
 
 ////Function that generates the 2D hist with the outputs of the BDTs.
 void BDTOutputPlane(double bottomHistLimit, double topHistLimit, vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, int sizeHH, int sizeqq, int sizettbar, TH2F& histBDTsOutputHH)
 {
 	TCanvas *cBDTOutputPlane1 = new TCanvas();
	cBDTOutputPlane1->cd();
 	histBDTsOutputHH.Draw("col");
 	for(int i=0; i<BDTqqOutput.size(); i++)
 	{
		double radius = 0.005;
		TEllipse* r = new TEllipse(BDTqqOutput[i], BDTttbarOutput[i], radius, radius); 
		if(i<sizeHH) r->SetLineColor(kBlue);	
		else if(i<(sizeHH+sizeqq)) r->SetLineColor(kRed+2);
		else if(i<(sizeHH+sizeqq+sizettbar)) r->SetLineColor(kPink+6);
		r->SetLineStyle(2);
              	r->Draw("same");
 	}
 }
 
 ////Function that generates a histogram with the overlap of two given histograms.
 void BDTOutputOverlap(double bottomHistLimit, double topHistLimit, int nbin, TH1F& histBDTHH, Int_t colorHH, TH1F& histBDTBack, Int_t colorBack, string title)
 {
 	TH1F *histBdtOverlap     = new TH1F( "MVA_BDT_Overlap", title.c_str(), nbin, bottomHistLimit-.05, topHistLimit+.05);
 	histBdtOverlap->SetFillColor(kViolet);
 	
 	TH1F *histBDTHHClone = (TH1F*)histBDTHH.Clone("histBDTHHClone");
 	TH1F *histBDTBackClone = (TH1F*)histBDTBack.Clone("histBDTBackClone");
 	histBDTHHClone->SetFillColor(colorHH);
 	histBDTBackClone->SetFillColor(colorBack);
 	
 	
 	double contEventsHH=0, contEventsBack=0;
 	
 	for (int i=1; i<=nbin; i++) 
  	{
     		double content1 = histBDTHHClone->GetBinContent(i);
     		double content2 = histBDTBackClone->GetBinContent(i);
	     	if(content1!=0 && content2!=0)
	    	{
	     		double overlappingContent = min(content1, content2);
	     		histBdtOverlap->SetBinContent(i, overlappingContent);
	     	}
          	contEventsHH+=content1;
          	contEventsBack+=content2;
   	}
   	
   	/*TCanvas *cBDTOutputOverlap1 = new TCanvas();
   	cBDTOutputOverlap1->cd();
   	histBDTHH->Draw("HIST");
   	histBDTBack->Draw("HIST same");
   	histBdtOverlap->Draw("HIST same");*/
   	
   	
   	double maxHH = 	histBDTHHClone->GetMaximum();
   	double maxBack = histBDTBackClone->GetMaximum();
   	
   	//cout<<"maxHH: "<<maxHH<<endl<<"maxBack: "<<maxBack<<endl;
   	
   	TCanvas *cBDTOutputOverlap1 = new TCanvas();
   	cBDTOutputOverlap1->cd();
   	if(maxHH>maxBack)
   	{
   		histBDTHHClone->SetTitle(title.c_str());
   		histBDTHHClone->Draw("HIST");
   		histBDTBackClone->Draw("HIST same");
   		histBdtOverlap->Draw("HIST same");
   	}
   	else
   	{
   		histBDTBackClone->SetTitle(title.c_str());
   		histBDTBackClone->Draw("HIST");
   		histBDTHHClone->Draw("HIST same");
   		histBdtOverlap->Draw("HIST same");
   	} 
   	
   	
   	return;
 }
 
 ////Function that, given two vectors, sorts the first one by increasing order and the entries of the second one to follow the first one. For example, if I have vector1 2,3,1,4,5 and vector2 7,2,9,3,1, I need vector1 to be 1,2,3,4,5 and vector2 to be 9,7,2,3,1. 
 void sortVectorPair(vector<double>& vector1, vector<double>& vector2)
 {
 	    std::vector<std::pair<double, double>> pairedVectors(vector1.size());
	    for (size_t i = 0; i < vector1.size(); ++i) {
		pairedVectors[i] = std::make_pair(vector1[i], vector2[i]);
	    }
	    // Step 2: Sort the vector of pairs based on the first element of the pair
	    std::sort(pairedVectors.begin(), pairedVectors.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
		return a.first < b.first;
	    });
	    // Step 3: Unpack the sorted pairs back into two separate vectors
	    for (size_t i = 0; i < pairedVectors.size(); ++i) {
		vector1[i] = pairedVectors[i].first;
		vector2[i] = pairedVectors[i].second;
	    }
	    return;
 }
 
 
 ////Function that, given five vectors, sorts the first one by increasing order and rearranges the entries of the second, third... ones to follow the first one. For example, if I have vector1 3,1,2, vector2 9,7,8, and vector3 4,2,3, I need vector1 to be 1,2,3 and vector2 to be 7,8,9 and vector3 to be 2,3,4.
void sortVectorQuintuple(std::vector<double>& vector1, std::vector<double>& vector2, std::vector<double>& vector3, std::vector<double>& vector4, std::vector<double>& vector5)
{
    // Check that all vectors are the same size
    if (vector1.size() != vector2.size() || vector1.size() != vector3.size() || vector1.size() != vector4.size() || vector1.size() != vector5.size()) {
        throw std::invalid_argument("All vectors must have the same size.");
    }

    // Step 1: Pair all vectors into a vector of tuples
    std::vector<std::tuple<double, double, double, double, double>> pairedVectors(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        pairedVectors[i] = std::make_tuple(vector1[i], vector2[i], vector3[i], vector4[i], vector5[i]);
    }

    // Step 2: Sort the vector of tuples based on the first element of the tuple
    std::sort(pairedVectors.begin(), pairedVectors.end(), 
        [](const std::tuple<double, double, double, double, double>& a, const std::tuple<double, double, double, double, double>& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

    // Step 3: Unpack the sorted tuples back into the four separate vectors
    for (size_t i = 0; i < pairedVectors.size(); ++i) {
        vector1[i] = std::get<0>(pairedVectors[i]);
        vector2[i] = std::get<1>(pairedVectors[i]);
        vector3[i] = std::get<2>(pairedVectors[i]);
        vector4[i] = std::get<3>(pairedVectors[i]);
        vector5[i] = std::get<4>(pairedVectors[i]);
    }

    return;
}
 
 
///Function that plots the significance as a function of the cuts on BDTqq and BDTttbar; it then returns the maximum significance.
 double findMaxSignificanceAndROC(double bottomHistLimit, double topHistLimit, double nbin, TH1F& histBDTHH, TH1F& histBDTBack, double& n, double& b, double& optCut)
 {
 	 double low=bottomHistLimit, high=topHistLimit;
 	 double cut;
 	 Int_t binLow = histBDTHH.FindBin(low);
 	 Int_t binHigh = histBDTHH.FindBin(high);
 	 double totalIntegral = histBDTHH.Integral(binLow, binHigh);
 	 double totalIntegralBack = histBDTBack.Integral(binLow, binHigh);
 	 double contEvents=totalIntegral;
	 double significance=0.0, maxSignificance=0.0, signalRight=0.0, BDTEffi=0.0;
         //cout<<endl<<endl<<"Total integral for histBdt: "<<totalIntegral<<endl<<"Right integral for histBdt: "<<rightIntegral<<endl;
         
         TH2F *histROC = new TH2F("ROC", "Signal and Background Acceptance", nbin, 0.0, 1.0, nbin, 0, 1.0);
  	 TH2F *histROCRej = new TH2F("ROCRej", "Signal Acceptance and Background Rejection", nbin, 0.5, 1.0, nbin, 0, 10000.0);
  	 histROC->GetXaxis()->SetTitle("Signal accepted");
   	 histROC->GetYaxis()->SetTitle("Background rejected (1/background accepted)");
   
  	 TH2F *histSignificance     = new TH2F( "hist_significance", "Significance s/sqrt(s+b)", nbin, bottomHistLimit-.1, topHistLimit+.1, nbin, 0, 10.0);
  	 histSignificance->GetXaxis()->SetTitle("Cut value for BDT");
   	 histSignificance->GetYaxis()->SetTitle("Significance");
   
   	 for (double i=bottomHistLimit; i<=topHistLimit; i+=0.0001) 
   	 {
   		cut=i;
   		Int_t binCut = histBDTHH.FindBin(cut);
   		double rightIntegral = histBDTHH.Integral(binCut, binHigh); 
   		double rightIntegralBack = histBDTBack.Integral(binCut, binHigh);
	   	/*cout<<"Cut: "<<i<<endl;
	   	cout<<"total integral: "<<totalIntegral<<endl;
	   	cout<<"total integral back: "<<totalIntegralBack<<endl;
	   	cout<<"Right integral: "<<rightIntegral<<endl;
	   	cout<<"Right integral back: "<<rightIntegralBack<<endl<<endl;*/
	   	histROC->Fill(rightIntegral/totalIntegral, rightIntegralBack/totalIntegralBack);
	   	histROCRej->Fill(rightIntegral/totalIntegral, 1/(rightIntegralBack/totalIntegralBack));
	   	significance=rightIntegral/(sqrt(rightIntegral+rightIntegralBack));
	   	histSignificance->Fill(cut, significance);
	   	if(significance>maxSignificance) 
	   	{
	   		maxSignificance=significance;
	   		optCut=cut;
	   		n=rightIntegral+rightIntegralBack;
	   		b=rightIntegralBack;
	   		signalRight=rightIntegral;
	   	}
	   	//if((rightIntegral/totalIntegral)>.599 && (rightIntegral/totalIntegral)<.601) cout<<"For "<<(rightIntegral/totalIntegral)<<" signal acceptance, background acceptance: "<<(rightIntegralBack/totalIntegralBack)<<endl<<endl<<endl;
	  }
	  
	  TCanvas *cfindMaxSignificanceAndROC1 = new TCanvas();
	  TCanvas *cfindMaxSignificanceAndROC2 = new TCanvas();
	  TCanvas *cfindMaxSignificanceAndROC3 = new TCanvas();
	  
	  cfindMaxSignificanceAndROC1->cd();
	  histROC->Draw("BOX");
	  
	  //gStyle->SetOptStat(0);
	  cfindMaxSignificanceAndROC2->cd();
	  cfindMaxSignificanceAndROC2->SetLogy();
	  histROCRej->Draw("BOX");
	  
	  cfindMaxSignificanceAndROC3->cd();
	  histSignificance->Draw("HIST");
	  
	  return maxSignificance;
 }
 
 /////Function that generates ROC curve without using histograms and integrals, as well as calculates the max. significance
 void findSignificance(double bottomHistLimit, double topHistLimit, double nbin, vector<double> NNOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW, int nBack, TH1F& histBDTHH, TH1F& histBDTBack, double& defCut, double& maxSignificance, double weightHH, double weightBack, double targetFractionHH, TH2F& histROC, TH2F& histROCRej, TH2F& histSignificance, double& totalRemaining, double& HHRemaining, double& backRemaining)
 {
 	defCut=bottomHistLimit;
	double significance; 	
  	 
	double size = NNOutput.size();
	
	//cout<<endl<<endl<<endl<<endl<<endl<<"NNOutputSize: "<<size<<endl<<endl<<endl<<endl;
	
	vector<double> NNOutputHH(NNOutput.begin(), NNOutput.begin()+sizeHH);
	double bottomLimitBack, topLimitBack;
	if(nBack == 0) 
	{
		bottomLimitBack = sizeHH;
		topLimitBack = sizeHH+sizeqq;
	}
	else if(nBack == 1)
	{
		bottomLimitBack = sizeHH+sizeqq;
		topLimitBack = sizeHH+sizeqq+sizettbar;
	}
	else if(nBack == 2)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar;
		topLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ;
	}
	else if(nBack == 3)
	{
		bottomLimitBack = sizeHH+sizeqq+sizettbar+sizeZZ;
		topLimitBack = size;
	}
	//cout<<"bottomLimitBack: "<<bottomLimitBack<<"      topLimitBack: "<<topLimitBack<<endl<<endl<<endl<<endl;
	vector<double> NNOutputBack(NNOutput.begin()+bottomLimitBack, NNOutput.end()-(size-topLimitBack));
	//cout<<"No crashea en la creacion del vector de back ouputs"<<endl<<endl<<endl;
	sort(NNOutputHH.begin(), NNOutputHH.end());
	sort(NNOutputBack.begin(), NNOutputBack.end());
	
	int sizeBack = NNOutputBack.size();
	
	cout<<"sizeBack: "<<sizeBack<<endl<<endl<<endl;
	
	//cout<<"Size: "<<size<<endl<<"SizeHH: "<<NNOutputHH.size()<<endl<<"SizeBack: "<<NNOutputBack.size()<<endl<<endl;
	/*for(int i=0; i<100; i++) cout<<NNOutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<100; i++) cout<<NNOutputBack[i]<<", ";
	cout<<endl;*/
	
	int indexHH=0, indexBack=0;
	double elementNHH=NNOutputHH[indexHH],elementN1HH, elementNBack=NNOutputBack[indexBack], elementN1Back;
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack=0, eventsCutBack=0;
	bool flagFraction=false;
	//for(double cut=bottomHistLimit; cut<=topHistLimit; cut+=0.000001)
	for(double cut=bottomHistLimit; cut<=topHistLimit; cut+=0.0001)
	{
		
		
		while(cut>elementNHH && indexHH<sizeHH)
		{
			indexHH++;
			if(indexHH<sizeHH) elementNHH = NNOutputHH[indexHH];
		}
		while(cut>elementNBack && indexBack<sizeBack)
		{
			indexBack++;
			if(indexBack<sizeBack) elementNBack = NNOutputBack[indexBack];
		}
		
		eventsRemainingHH = sizeHH-indexHH;
		eventsRemainingBack = sizeBack-indexBack;
		
		double fractionHH = eventsRemainingHH/sizeHH;
		double fractionBack = eventsRemainingBack/sizeBack;
		histROCRej.Fill(fractionHH, 1/(fractionBack));
		histROC.Fill(eventsRemainingHH/sizeHH, eventsRemainingBack/sizeBack);
		
		significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBack*weightBack)));
		/*cout<<"CUT: "<<cut<<endl<<endl;
		cout<<"IndexHH: "<<indexHH<<"    elem: "<<elementNHH<<endl<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl<<endl;
		cout<<"IndexBack: "<<indexBack<<"    elem: "<<elementNBack<<endl<<"Back remaining: "<<eventsRemainingBack<<"    back cut: "<<sizeBack-eventsRemainingBack<<endl;*/
		//cout<<"Significance: "<<significance<<"         maxSignificnace: "<<maxSignificance<<endl<<endl;
		//cout<<"--------"<<endl;
	   	histSignificance.Fill(cut, significance);
	   	if(significance>maxSignificance) 
	   	{
	   		maxSignificance=significance;
	   		defCut=cut;
	   		backRemaining=eventsRemainingBack*weightBack;
	   		HHRemaining=eventsRemainingHH*weightHH;
	   		totalRemaining=eventsRemainingHH+eventsRemainingBack;
	   	}
		
		/*if(fractionHH < targetFractionHH+0.001 && fractionHH > targetFractionHH-0.001 && !flagFraction)
		{
			flagFraction=true;
			cout<<endl<<"(WITH ARRAYS) For "<<fractionHH*100<<"% of signal ("<<eventsRemainingHH<<" out of "<<sizeHH<<" events), background: "<<fractionBack*100<<"% ("<<eventsRemainingBack<<" out of "<<sizeBack<<" events)"<<endl;
			double low=bottomHistLimit, high=topHistLimit;
 	 		Int_t binLow = histBDTHH.FindBin(cut);
 	 		Int_t binHigh = histBDTHH.FindBin(high);
 	 		double eventsRemainingHHIntegral = histBDTHH.Integral(binLow, binHigh);
 	 		double eventsRemainingBackIntegral = histBDTBack.Integral(binLow, binHigh);
 	 		double fractionHHIntegral = eventsRemainingHHIntegral/sizeHH;
			double fractionBackIntegral = eventsRemainingBackIntegral/sizeBack;
 	 		cout<<endl<<"(WITH INTEGRALS) For "<<fractionHHIntegral*100<<"% of signal ("<<eventsRemainingHHIntegral<<" out of "<<sizeHH<<" events), background: "<<fractionBackIntegral*100<<"% ("<<eventsRemainingBackIntegral<<" out of "<<sizeBack<<" events)"<<endl;
 	 		defCut = cut;
 	 		cout<<"Cut for this %: "<<cut<<endl<<endl;
 	 		
 	 		cout<<"Significance for this: "<<significance<<endl<<endl;
		}*/
	}
	cout<<"HHRemainingUnweighted: "<<HHRemaining/weightHH<<endl<<"backRemainingUnweighted: "<<backRemaining/weightBack<<endl; 
	cout<<"HHRemaining: "<<HHRemaining<<endl<<"backRemaining: "<<backRemaining<<endl; 
	
	TCanvas *cfindSignificance1 = new TCanvas();
	TCanvas *cfindSignificance2 = new TCanvas();
	TCanvas *cfindSignificance3 = new TCanvas();
	
	cfindSignificance1->cd();
	histROC.Draw("BOX");
	  
	//gStyle->SetOptStat(0);
	cfindSignificance2->cd();
	cfindSignificance2->SetLogy();
	histROCRej.Draw("BOX");
	
	cfindSignificance3->cd();
	histSignificance.Draw("BOX");
	
 }
 
 
 ///Function that calculates the significance of the data after applying multiple NNs.
 void findSignificanceCombined(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NN1Output, vector<double>& NN2Output, int sizeHH, int sizeBack1, int sizeBack2, int topology1, int topology2, double& defCutNN1, double& defCutNN2, double& maxSignificanceCombined, double weightHH, double weightBack1, double weightBack2, TH2F& histROCCombined, TH2F& histROCRejCombined, TH2F& histSignificanceCombined, double& totalRemainingCombined, double& HHRemainingCombined, double& backRemainingCombined)
{
	defCutNN1=bottomHistLimit;
	defCutNN2=bottomHistLimit;
	double significance;
	
	vector<double> NN1OutputHH(NN1Output.begin(), NN1Output.begin()+sizeHH);
	vector<double> NN2OutputHH(NN2Output.begin(), NN2Output.begin()+sizeHH);
	vector<double> NN1OutputBack1(NN1Output.begin()+sizeHH, NN1Output.begin()+sizeHH+sizeBack1);
	vector<double> NN2OutputBack1(NN2Output.begin()+sizeHH, NN2Output.begin()+sizeHH+sizeBack1);
	vector<double> NN1OutputBack2(NN1Output.begin()+sizeHH+sizeBack1, NN1Output.begin()+sizeHH+sizeBack1+sizeBack2);
	vector<double> NN2OutputBack2(NN2Output.begin()+sizeHH+sizeBack1, NN2Output.begin()+sizeHH+sizeBack1+sizeBack2);
	
	
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<"    NN2OutputHH size: "<<NN2OutputHH.size()<<endl;
	cout<<"NN1OutputBack1 size: "<<NN1OutputBack1.size()<<"    NN2OutputBack1 size: "<<NN2OutputBack1.size()<<endl;
	cout<<"NN1OutputBack2 size: "<<NN1OutputBack2.size()<<"    NN2OutputBack2 size: "<<NN2OutputBack2.size()<<endl<<endl<<endl;
	
	for(int i=0; i<6; i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<6; i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN1OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN2OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN1OutputBack2[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN2OutputBack2[i]<<", ";
	cout<<endl;*/
	
	sortVectorPair(NN1OutputHH, NN2OutputHH);
	sortVectorPair(NN1OutputBack1, NN2OutputBack1);
	sortVectorPair(NN1OutputBack2, NN2OutputBack2);
	
	/*cout<<"Now sorted: "<<endl;
	for(int i=0; i<6; i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<6; i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN1OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<4; i++) cout<<NN2OutputBack1[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN1OutputBack2[i]<<", ";
	cout<<endl;
	for(int i=0; i<3; i++) cout<<NN2OutputBack2[i]<<", ";
	cout<<endl;*/
	
	int index1HH=0, index2HH=0, index1Back1=0, index2Back1=0, index1Back2=0, index2Back2=0;
	double elementN1HH=NN1OutputHH[index1HH], elementN2HH=NN2OutputHH[index2HH], elementN1Back1=NN1OutputBack1[index1Back1], elementN2Back1=NN2OutputBack1[index2Back1], elementN1Back2=NN1OutputBack2[index1Back2], elementN2Back2=NN2OutputBack2[index2Back2];
	double eventsRemainingHH=0, eventsCutHH=0, eventsRemainingBack1=0, eventsCutBack1=0, eventsRemainingBack2=0, eventsCutBack2=0, eventsRemainingBack=0, eventsCutBack=0;
	double sizeBack=sizeBack1+sizeBack2;
	
	for(double cut1=bottomHistLimit; cut1<=topHistLimit; cut1+=0.001)
	{
		while(cut1>elementN1HH && index1HH<sizeHH)
		{
			index1HH++;
			if(index1HH<sizeHH) elementN1HH = NN1OutputHH[index1HH];
		}
		while(cut1>elementN1Back1 && index1Back1<sizeBack1)
		{
			index1Back1++;
			if(index1Back1<sizeBack1) elementN1Back1 = NN1OutputBack1[index1Back1];
		}
		while(cut1>elementN1Back2 && index1Back2<sizeBack2)
		{
			index1Back2++;
			if(index1Back2<sizeBack2) elementN1Back2 = NN1OutputBack2[index1Back2];
		}
		
		vector<double> NNOutputSurvivorsHH(NN2OutputHH.begin()+index1HH, NN2OutputHH.end());
		vector<double> NNOutputSurvivorsBack1(NN2OutputBack1.begin()+index1Back1, NN2OutputBack1.end());
		vector<double> NNOutputSurvivorsBack2(NN2OutputBack2.begin()+index1Back2, NN2OutputBack2.end());
		
		sort(NNOutputSurvivorsHH.begin(), NNOutputSurvivorsHH.end());
		sort(NNOutputSurvivorsBack1.begin(), NNOutputSurvivorsBack1.end());
		sort(NNOutputSurvivorsBack2.begin(), NNOutputSurvivorsBack2.end());
		
		index2HH=0;
		index2Back1=0;
		index2Back2=0;
		if(NNOutputSurvivorsHH.size()>0) elementN2HH=NNOutputSurvivorsHH[index2HH];
		if(NNOutputSurvivorsBack1.size()>0) elementN2Back1=NNOutputSurvivorsBack1[index2Back1];
		if(NNOutputSurvivorsBack2.size()>0) elementN2Back2=NNOutputSurvivorsBack2[index2Back2];
		
		for(double cut2=bottomHistLimit; cut2<=topHistLimit; cut2+=0.001)
		{
			while(cut2>elementN2HH && index2HH<NNOutputSurvivorsHH.size())
			{
				index2HH++;
				if(index2HH<NNOutputSurvivorsHH.size()) elementN2HH = NNOutputSurvivorsHH[index2HH];
			}
			while(cut2>elementN2Back1 && index2Back1<NNOutputSurvivorsBack1.size())
			{
				index2Back1++;
				if(index2Back1<NNOutputSurvivorsBack1.size()) elementN2Back1 = NNOutputSurvivorsBack1[index2Back1];
			}
			while(cut2>elementN2Back2 && index2Back2<NNOutputSurvivorsBack2.size())
			{
				index2Back2++;
				if(index2Back2<NNOutputSurvivorsBack2.size()) elementN2Back2 = NNOutputSurvivorsBack2[index2Back2];
			}
			
			eventsRemainingHH = (sizeHH-index1HH)-(NNOutputSurvivorsHH.size()-(NNOutputSurvivorsHH.size()-index2HH));
			eventsRemainingBack1 = (sizeBack1-index1Back1)-(NNOutputSurvivorsBack1.size()-(NNOutputSurvivorsBack1.size()-index2Back1));
			eventsRemainingBack2 = (sizeBack2-index1Back2)-(NNOutputSurvivorsBack2.size()-(NNOutputSurvivorsBack2.size()-index2Back2));
			
			
			/*weightHH=1;
			weightBack1=1;
			weightBack2=1;*/
			
			eventsRemainingBack = eventsRemainingBack1*weightBack1 + eventsRemainingBack2*weightBack2;
			
			significance=(eventsRemainingHH*weightHH)/(sqrt((eventsRemainingHH*weightHH)+(eventsRemainingBack1*weightBack1+eventsRemainingBack2*weightBack2)));
			
			if(significance>maxSignificanceCombined) 
		   	{
		   		maxSignificanceCombined=significance;
		   		defCutNN1=cut1;
		   		defCutNN2=cut2;
		   		totalRemainingCombined=eventsRemainingHH*weightHH+eventsRemainingBack;
		   		backRemainingCombined=eventsRemainingBack;
		   		HHRemainingCombined=eventsRemainingHH*weightHH;
		   	}
		   	
		   	
		   	/*cout<<"CUT1: "<<cut1<<endl;
			cout<<"Index1HH: "<<index1HH<<"    elementN1HH: "<<elementN1HH<<endl;
			cout<<"Index1Back1: "<<index1Back1<<"    elementN1Back1: "<<elementN1Back1<<endl;
			cout<<"Index1Back2: "<<index1Back2<<"    elementN1Back2: "<<elementN1Back2<<endl;
			cout<<"CUT2: "<<cut2<<endl;
			cout<<"Index2HH: "<<index2HH<<"    elementN2HH: "<<elementN2HH<<endl;
			cout<<"Index2Back1: "<<index2Back1<<"    elementN2Back1: "<<elementN2Back1<<endl;
			cout<<"Index2Back2: "<<index2Back2<<"    elementN2Back2: "<<elementN2Back2<<endl;
			cout<<"Signal remaining: "<<eventsRemainingHH<<"    signal cut: "<<sizeHH-eventsRemainingHH<<endl;
			cout<<"Back1 remaining: "<<eventsRemainingBack1<<"    back1 cut: "<<sizeBack1-eventsRemainingBack1<<endl;
			cout<<"Back2 remaining: "<<eventsRemainingBack2<<"    back2 cut: "<<sizeBack2-eventsRemainingBack2<<endl;
			for(int i=0; i<NNOutputSurvivorsHH.size(); i++) cout<<NNOutputSurvivorsHH[i]<<", ";
			cout<<endl;
			for(int i=0; i<NNOutputSurvivorsBack1.size(); i++) cout<<NNOutputSurvivorsBack1[i]<<", ";
			cout<<endl;
			for(int i=0; i<NNOutputSurvivorsBack2.size(); i++) cout<<NNOutputSurvivorsBack2[i]<<", ";
			cout<<endl;
			cout<<"Significance: "<<significance<<"         maxSignificance: "<<maxSignificanceCombined<<endl<<endl;
			cout<<"--------"<<endl<<endl;*/
		
		
		}
	}
	
	/*cout<<"maxSignificanceCombined: "<<maxSignificanceCombined<<endl;
	cout<<"For cutNN1: "<<defCutNN1<<"    and cutNN2: "<<defCutNN2<<endl;
	cout<<"totalRemainingCombined: "<<totalRemainingCombined<<endl;
	cout<<"HHRemainingCombined: "<<HHRemainingCombined<<endl;
	cout<<"backRemainingCombined: "<<backRemainingCombined<<endl;*/

}

void findSignificanceCutsCombined(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NN1Output, vector<double>& NN2Output, vector<double>& NN3Output, vector<double>& NN4Output, int sizeHH, int sizeBack1, int sizeBack2, int sizeBack3, int sizeBack4, double defCutNN1, double defCutNN2, double defCutNN3, double defCutNN4, double& maxSignificanceCombined, double weightHH, double weightBack1, double weightBack2, double weightBack3, double weightBack4, TH2F& histROCCombined, TH2F& histROCRejCombined, TH2F& histSignificanceCombined, double& totalRemainingCombined, double& HHRemainingCombined, double& backRemainingCombined)
{
	vector<double> topologyTracker;
	for(int i=0; i<NN1Output.size(); i++)
	{
		if(i<sizeHH) topologyTracker.push_back(0);
		else if(i<sizeHH+sizeBack1) topologyTracker.push_back(1); 
		else if(i<sizeHH+sizeBack1+sizeBack2) topologyTracker.push_back(2);
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3) topologyTracker.push_back(3);  
		else if(i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) topologyTracker.push_back(4);  
	}
	
	/*cout<<endl<<"topologyTracker: ";
	for(int i=0;i<topologyTracker.size();i++) cout<<topologyTracker[i]<<", ";
	
	cout<<endl<<endl<<endl<<"NN1Output size: "<<NN1Output.size()<<"    NN2Output size: "<<NN2Output.size()<<"    NN3Output size: "<<NN3Output.size()<<"    NN4Output size: "<<NN3Output.size()<<endl<<endl;
	
	cout<<"NN1Output: ";
	for(int i=0;i<NN1Output.size();i++) cout<<NN1Output[i]<<", ";
	cout<<endl<<"NN2Output: ";
	for(int i=0;i<NN2Output.size();i++) cout<<NN2Output[i]<<", ";
	cout<<endl<<"NN3Output: ";
	for(int i=0;i<NN3Output.size();i++) cout<<NN3Output[i]<<", ";
	cout<<endl<<"NN4Output: ";
	for(int i=0;i<NN4Output.size();i++) cout<<NN4Output[i]<<", ";
	cout<<endl<<"topologyTracker: ";
	for(int i=0;i<topologyTracker.size();i++) cout<<topologyTracker[i]<<", ";*/
	
	vector<double> NN1OutputHH(NN1Output.begin(), NN1Output.begin()+sizeHH);
	vector<double> NN1OutputBack(NN1Output.begin()+sizeHH, NN1Output.end());
	vector<double> NN2OutputHH(NN2Output.begin(), NN2Output.begin()+sizeHH);
	vector<double> NN2OutputBack(NN2Output.begin()+sizeHH, NN2Output.end());
	vector<double> NN3OutputHH(NN3Output.begin(), NN3Output.begin()+sizeHH);
	vector<double> NN3OutputBack(NN3Output.begin()+sizeHH, NN3Output.end());
	vector<double> NN4OutputHH(NN4Output.begin(), NN4Output.begin()+sizeHH);
	vector<double> NN4OutputBack(NN4Output.begin()+sizeHH, NN4Output.end());
	vector<double> topologyTrackerHH(topologyTracker.begin(), topologyTracker.begin()+sizeHH);
	vector<double> topologyTrackerBack(topologyTracker.begin()+sizeHH, topologyTracker.end());
	
	/*cout<<endl<<endl<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<"topologyTrackerHH size: "<<topologyTrackerHH.size()<<endl<<"topologyTrackerBack size: "<<topologyTrackerBack.size()<<endl<<endl;
	
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	sortVectorQuintuple(NN1OutputHH, NN2OutputHH, NN3OutputHH, NN4OutputHH, topologyTrackerHH);
	sortVectorQuintuple(NN1OutputBack, NN2OutputBack, NN3OutputBack, NN4OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after qq sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	//cout<<"HH starting: "<<topologyTrackerHH.size()<<endl;
	
	double event;
	
	int indexNN1OutputHH = 0;
	event = NN1OutputHH[indexNN1OutputHH];
	while(event<defCutNN1 && indexNN1OutputHH<NN1OutputHH.size())
	{
		indexNN1OutputHH++;
		event = NN1OutputHH[indexNN1OutputHH];
	}
	
	int indexNN1OutputBack = 0;
	event = NN1OutputBack[indexNN1OutputBack];
	while(event<defCutNN1 && indexNN1OutputBack<NN1OutputBack.size())
	{
		indexNN1OutputBack++;
		event = NN1OutputBack[indexNN1OutputBack];
	}
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN1OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN1OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN1OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN1OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN1OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN1OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN1OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN1OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN1OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN1OutputBack);
	

	/*cout<<endl<<endl<<"After qq cutting: "<<endl;
	cout<<"indexNN1OutputHH: "<<indexNN1OutputHH<<"    indexNN1OutputBack: "<<indexNN1OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;
	
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	
	
	sortVectorQuintuple(NN2OutputHH, NN3OutputHH, NN4OutputHH, NN1OutputHH, topologyTrackerHH);
	sortVectorQuintuple(NN2OutputBack, NN3OutputBack, NN4OutputBack, NN1OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after ttbar sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	int indexNN2OutputHH = 0;
	event = NN2OutputHH[indexNN2OutputHH];
	while(event<defCutNN2 && indexNN2OutputHH<NN2OutputHH.size())
	{
		indexNN2OutputHH++;
		event = NN2OutputHH[indexNN2OutputHH];
	}
	
	int indexNN2OutputBack = 0;
	event = NN2OutputBack[indexNN2OutputBack];
	while(event<defCutNN2 && indexNN2OutputBack<NN2OutputBack.size())
	{
		indexNN2OutputBack++;
		event = NN2OutputBack[indexNN2OutputBack];
	}
	
	//cout<<endl<<endl<<"indexNN2OutputHH: "<<indexNN2OutputHH<<"    indexNN2OutputBack: "<<indexNN2OutputBack<<endl;
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN2OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN2OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN2OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN2OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN2OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN2OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN2OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN2OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN2OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN2OutputBack);
	
	/*cout<<endl<<endl<<"After ttbar cutting: "<<endl;
	cout<<"indexNN2OutputHH: "<<indexNN2OutputHH<<"    indexNN2OutputBack: "<<indexNN2OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	
	//cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;
	
	
	sortVectorQuintuple(NN3OutputHH, NN4OutputHH, NN1OutputHH, NN2OutputHH, topologyTrackerHH);
	sortVectorQuintuple(NN3OutputBack, NN4OutputBack, NN1OutputBack, NN2OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after ZZ sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	
	
	int indexNN3OutputHH = 0;
	event = NN3OutputHH[indexNN3OutputHH];
	while(event<defCutNN3 && indexNN3OutputHH<NN3OutputHH.size())
	{
		indexNN3OutputHH++;
		event = NN3OutputHH[indexNN3OutputHH];
	}
	
	int indexNN3OutputBack = 0;
	event = NN3OutputBack[indexNN3OutputBack];
	while(event<defCutNN3 && indexNN3OutputBack<NN3OutputBack.size())
	{
		indexNN3OutputBack++;
		event = NN3OutputBack[indexNN3OutputBack];
	}
	
	//cout<<endl<<endl<<"indexNN3OutputHH: "<<indexNN3OutputHH<<"    indexNN3OutputBack: "<<indexNN3OutputBack<<endl;
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN3OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN3OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN3OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN3OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN3OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN3OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN3OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN3OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN3OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN3OutputBack);
	
	/*cout<<endl<<endl<<"After ZZ cutting: "<<endl;
	cout<<"indexNN3OutputHH: "<<indexNN3OutputHH<<"    indexNN3OutputBack: "<<indexNN3OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;*/
	
	sortVectorQuintuple(NN4OutputHH, NN1OutputHH, NN2OutputHH, NN3OutputHH, topologyTrackerHH);
	sortVectorQuintuple(NN4OutputBack, NN1OutputBack, NN2OutputBack, NN3OutputBack, topologyTrackerBack);
	
	/*cout<<endl<<endl<<"Now after WW sorting: ";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	
	
	int indexNN4OutputHH = 0;
	event = NN4OutputHH[indexNN4OutputHH];
	while(event<defCutNN4 && indexNN4OutputHH<NN4OutputHH.size())
	{
		indexNN4OutputHH++;
		event = NN4OutputHH[indexNN4OutputHH];
	}
	
	int indexNN4OutputBack = 0;
	event = NN4OutputBack[indexNN4OutputBack];
	while(event<defCutNN4 && indexNN4OutputBack<NN4OutputBack.size())
	{
		indexNN4OutputBack++;
		event = NN4OutputBack[indexNN4OutputBack];
	}
	
	//cout<<endl<<endl<<"indexNN4OutputHH: "<<indexNN4OutputHH<<"    indexNN4OutputBack: "<<indexNN4OutputBack<<endl;
	
	NN1OutputHH.erase(NN1OutputHH.begin(), NN1OutputHH.begin()+indexNN4OutputHH);
	NN1OutputBack.erase(NN1OutputBack.begin(), NN1OutputBack.begin()+indexNN4OutputBack);
	NN2OutputHH.erase(NN2OutputHH.begin(), NN2OutputHH.begin()+indexNN4OutputHH);
	NN2OutputBack.erase(NN2OutputBack.begin(), NN2OutputBack.begin()+indexNN4OutputBack);
	NN3OutputHH.erase(NN3OutputHH.begin(), NN3OutputHH.begin()+indexNN4OutputHH);
	NN3OutputBack.erase(NN3OutputBack.begin(), NN3OutputBack.begin()+indexNN4OutputBack);
	NN4OutputHH.erase(NN4OutputHH.begin(), NN4OutputHH.begin()+indexNN4OutputHH);
	NN4OutputBack.erase(NN4OutputBack.begin(), NN4OutputBack.begin()+indexNN4OutputBack);
	topologyTrackerHH.erase(topologyTrackerHH.begin(), topologyTrackerHH.begin()+indexNN4OutputHH);
	topologyTrackerBack.erase(topologyTrackerBack.begin(), topologyTrackerBack.begin()+indexNN4OutputBack);
	
	/*cout<<endl<<endl<<"After WW cutting: "<<endl;
	cout<<"indexNN4OutputHH: "<<indexNN4OutputHH<<"    indexNN4OutputBack: "<<indexNN4OutputBack<<endl;
	cout<<"HH remaining: "<<topologyTrackerHH.size()<<endl;*/
	
	/*cout<<"NN1OutputHH size: "<<NN1OutputHH.size()<<endl<<"NN1OutputBack size: "<<NN1OutputBack.size()<<endl<<"NN2OutputHH size: "<<NN2OutputHH.size()<<endl<<"NN2OutputBack size: "<<NN2OutputBack.size()<<endl<<"NN3OutputHH size: "<<NN3OutputHH.size()<<endl<<"NN3OutputBack size: "<<NN3OutputBack.size()<<endl<<"NN4OutputHH size: "<<NN4OutputHH.size()<<endl<<"NN4OutputBack size: "<<NN4OutputBack.size()<<endl<<endl;*/
	
	/*cout<<endl<<endl<<"Final vectors";
	cout<<endl<<"NN1OutputHH: ";
	for(int i=0;i<NN1OutputHH.size();i++) cout<<NN1OutputHH[i]<<", ";
	cout<<endl<<"NN1OutputBack: ";
	for(int i=0;i<NN1OutputBack.size();i++) cout<<NN1OutputBack[i]<<", ";
	cout<<endl<<"NN2OutputHH: ";
	for(int i=0;i<NN2OutputHH.size();i++) cout<<NN2OutputHH[i]<<", ";
	cout<<endl<<"NN2OutputBack: ";
	for(int i=0;i<NN2OutputBack.size();i++) cout<<NN2OutputBack[i]<<", ";
	cout<<endl<<"NN3OutputHH: ";
	for(int i=0;i<NN3OutputHH.size();i++) cout<<NN3OutputHH[i]<<", ";
	cout<<endl<<"NN3OutputBack: ";
	for(int i=0;i<NN3OutputBack.size();i++) cout<<NN3OutputBack[i]<<", ";
	cout<<endl<<"NN4OutputHH: ";
	for(int i=0;i<NN4OutputHH.size();i++) cout<<NN4OutputHH[i]<<", ";
	cout<<endl<<"NN4OutputBack: ";
	for(int i=0;i<NN4OutputBack.size();i++) cout<<NN4OutputBack[i]<<", ";
	cout<<endl<<"topologyTrackerHH: ";
	for(int i=0;i<topologyTrackerHH.size();i++) cout<<topologyTrackerHH[i]<<", ";
	cout<<endl<<"topologyTrackerBack: ";
	for(int i=0;i<topologyTrackerBack.size();i++) cout<<topologyTrackerBack[i]<<", ";*/
	
	int contBack1=0, contBack2=0, contBack3=0, contBack4=0;
	for(int i=0; i<topologyTrackerBack.size(); i++)
	{
		if(topologyTrackerBack[i]==1) contBack1++;
		else if(topologyTrackerBack[i]==2) contBack2++;
		else if(topologyTrackerBack[i]==3) contBack3++;
		else if(topologyTrackerBack[i]==4) contBack4++;
	}
	
	HHRemainingCombined = NN1OutputHH.size()*weightHH;
	backRemainingCombined = contBack1*weightBack1+contBack2*weightBack2+contBack3*weightBack3+contBack4*weightBack4;
	int HHRemainingCombinedUnweighted = NN1OutputHH.size();
	int backRemainingCombinedUnweighted = contBack1+contBack2+contBack3+contBack4;
	
	maxSignificanceCombined = HHRemainingCombined/sqrt(HHRemainingCombined+backRemainingCombined);
	
	cout<<endl<<endl<<endl<<"remainingHH unweighted: "<<HHRemainingCombinedUnweighted<<endl<<"remainingBackUnweighted: "<<backRemainingCombinedUnweighted<<endl<<endl;
	cout<<"remainingHH: "<<HHRemainingCombined<<endl<<"remainingBack: "<<backRemainingCombined<<endl<<endl<<endl;
	
	
	return;
}

void findSignificanceCutsCombinedCheck(double bottomHistLimit, double topHistLimit, double nbin, vector<double>& NN1Output, vector<double>& NN2Output, vector<double>& NN3Output, vector<double>& NN4Output, int sizeHH, int sizeBack1, int sizeBack2, int sizeBack3, int sizeBack4, double defCutNN1, double defCutNN2, double defCutNN3, double defCutNN4, double& maxSignificanceCombined, double weightHH, double weightBack1, double weightBack2, double weightBack3, double weightBack4, TH2F& histROCCombined, TH2F& histROCRejCombined, TH2F& histSignificanceCombined, double& totalRemainingCombined, double& HHRemainingCombined, double& backRemainingCombined)
{
	vector<double> NN1OutputCopy(NN1Output.begin(), NN1Output.end());
	vector<double> NN2OutputCopy(NN2Output.begin(), NN2Output.end());
	vector<double> NN3OutputCopy(NN3Output.begin(), NN3Output.end());
	vector<double> NN4OutputCopy(NN4Output.begin(), NN4Output.end());
	
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]<defCutNN1 && NN1OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
		}
	}
	for(int i=0; i<NN2OutputCopy.size(); i++)
	{
		if(NN2OutputCopy[i]<defCutNN2 && NN2OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
		}
	}
	for(int i=0; i<NN3OutputCopy.size(); i++)
	{
		if(NN3OutputCopy[i]<defCutNN3 && NN3OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
		}
	}
	for(int i=0; i<NN4OutputCopy.size(); i++)
	{
		if(NN4OutputCopy[i]<defCutNN4 && NN4OutputCopy[i]!=-999)
		{
			NN1OutputCopy[i]=-999;
			NN2OutputCopy[i]=-999;
			NN3OutputCopy[i]=-999;
			NN4OutputCopy[i]=-999;
		}
	}
	
	int contHH=0, contBack1=0, contBack2=0, contBack3=0, contBack4=0;
	for(int i=0; i<NN1OutputCopy.size(); i++)
	{
		if(NN1OutputCopy[i]!=-999 && i<sizeHH) contHH++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1) contBack1++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2) contBack2++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3) contBack3++;
		else if(NN1OutputCopy[i]!=-999 && i<sizeHH+sizeBack1+sizeBack2+sizeBack3+sizeBack4) contBack4++;
	}
	
	HHRemainingCombined = contHH*weightHH;
	backRemainingCombined = contBack1*weightBack1+contBack2*weightBack2+contBack3*weightBack3+contBack4*weightBack4;
	int HHRemainingCombinedUnweighted = contHH;
	int backRemainingCombinedUnweighted = contBack1+contBack2+contBack3+contBack4;
	
	maxSignificanceCombined = HHRemainingCombined/sqrt(HHRemainingCombined+backRemainingCombined);
	
	cout<<endl<<endl<<"From check: "<<endl<<"remainingHH unweighted: "<<HHRemainingCombinedUnweighted<<endl<<"remainingBackUnweighted: "<<backRemainingCombinedUnweighted<<endl<<endl;
	cout<<"remainingHH: "<<HHRemainingCombined<<endl<<"remainingBack: "<<backRemainingCombined<<endl<<endl<<endl;
}

 
 void splitHistogram(TH1F& hist, double bottomHistLimit, double topHistLimit, double cut, string title)
 {
 	TH1F *histRemaining = new TH1F("histRemaining", title.c_str(), hist.GetNbinsX(), bottomHistLimit, topHistLimit);
 	histRemaining->SetFillColor(hist.GetFillColor());
 	
 	
 	int binLow = hist.FindBin(cut);
 	int binHigh = hist.FindBin(topHistLimit);
 	for(int i=binLow; i<=binHigh; i++)
 	{
 	 	double binContent = hist.GetBinContent(i);
        	histRemaining->SetBinContent(i, binContent);
 	} 
 	
 	TH1F *histRemainingZoommed = (TH1F*)histRemaining->Clone("histRemainingZoommed");
	histRemainingZoommed->GetYaxis()->SetRangeUser(0, histRemaining->GetMaximum()*1.05);
 	
 	TLine* linex = new TLine(cut,0.0,cut,(hist.GetMaximum())*1.05);
	linex->SetLineColor(kOrange);
	linex->SetLineWidth(2);
	linex->SetLineStyle(7);
 	
 	TCanvas *csplitHistogram1 = new TCanvas();
 	TCanvas *csplitHistogram2 = new TCanvas();
 	TCanvas *csplitHistogram3 = new TCanvas();
 	
 	csplitHistogram1->cd();
 	hist.Draw("HIST");
	linex->Draw("same");
	
	histRemaining->GetYaxis()->SetRangeUser(0, hist.GetMaximum()*1.05);
	csplitHistogram2->cd();
 	histRemaining->Draw("HIST");
	linex->Draw("same");
	
	csplitHistogram3->cd();
 	histRemainingZoommed->Draw("HIST");
 	TLine* linex2 = new TLine(cut,0.0,cut,(histRemainingZoommed->GetMaximum())*1.05);
	linex2->SetLineColor(kOrange);
	linex2->SetLineWidth(2);
	linex2->SetLineStyle(7);
	linex2->Draw("same");
	binLow = hist.FindBin(bottomHistLimit);
	int entries = histRemainingZoommed->Integral(binLow, binHigh);
	string entriesText = "Entries: " + to_string(entries);
	TPaveText *pave = new TPaveText(0.3, 100, 0.5, 150, "");	
        pave->AddText(entriesText.c_str());
        pave->SetBorderSize(0);  // optionally set border size
        pave->SetFillColor(0);   // set fill color (0 is transparent)
        pave->SetTextFont(1);   // optional: set the font
        pave->SetTextSize(0.03); // optional: set the text size
        pave->Draw("same");
	
	
 }
 
void pairMassCut(int topology, double targetFractionHH, double& effi)
{
    	TFile* input = new TFile("outputTreesHHbbbbESpreadDurham.root");
    
 	std::cout << "--- pairMassCut   : Using input file: " << input->GetName() << std::endl;
 
    	TTree* theTree;
    
    	if(topology == 1) theTree = (TTree*)input->Get("TreeS");
    	else if(topology == 2) theTree = (TTree*)input->Get("TreeBqq");
    	else if(topology == 3) theTree = (TTree*)input->Get("TreeBtt");
    
    	Float_t aplanarity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, invMassB1, invMassB2, jetB1M, jetB2M, jetB3M, jetB4M, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, minJetM, sphericity, sumPt;
    	theTree->SetBranchAddress( "aplanarity", &aplanarity );
    	theTree->SetBranchAddress( "cosThetaB1", &cosThetaB1 );
    	theTree->SetBranchAddress( "cosThetaB2", &cosThetaB2 );
    	theTree->SetBranchAddress( "cosThetaB3", &cosThetaB3 );
    	theTree->SetBranchAddress( "cosThetaB4", &cosThetaB4 );
    	theTree->SetBranchAddress( "invMassB1", &invMassB1 );
    	theTree->SetBranchAddress( "invMassB2", &invMassB2 );
    	theTree->SetBranchAddress( "jetB1M", &jetB1M );
    	theTree->SetBranchAddress( "jetB2M", &jetB2M );
    	theTree->SetBranchAddress( "jetB3M", &jetB3M );
    	theTree->SetBranchAddress( "jetB4M", &jetB4M );
    	theTree->SetBranchAddress( "jetB1Pt", &jetB1Pt );
    	theTree->SetBranchAddress( "jetB2Pt", &jetB2Pt );
    	theTree->SetBranchAddress( "jetB3Pt", &jetB3Pt );
    	theTree->SetBranchAddress( "jetB4Pt", &jetB4Pt );
    	theTree->SetBranchAddress( "minJetM", &minJetM );
    	theTree->SetBranchAddress( "sphericity", &sphericity );
    	theTree->SetBranchAddress( "sumPt", &sumPt );
    	
    	double nEntries = theTree->GetEntries();
    	double invMassB1Remaining=0, invMassB2Remaining=0, eventsRemaining=0;
    	for (Long64_t entry = 0; entry < nEntries; entry++) 
    	{
		theTree->GetEntry(entry);
		//std::cout << "Entry " << entry << ": invMassB = " << invMassB << std::endl;
		
		if(50<invMassB1 && invMassB1<180 && 50<invMassB2 && invMassB2<180)
		{ 
			invMassB1Remaining++;
			invMassB2Remaining++;
			eventsRemaining++;
		}
		
		/*if(73<invMassB2 && invMassB2<140)
		{ 
			invMassB2Remaining++;
		}*/
    	}
    	
    	cout<<"invMassB1Remaining: "<<eventsRemaining<<"      fraction: "<<eventsRemaining/nEntries<<endl<<"invMassB2Remaining: "<<invMassB2Remaining<<"      fraction: "<<invMassB2Remaining/nEntries<<endl;
    	
    	effi = eventsRemaining/nEntries;
    	
    	
    	return;
}

void addPointTo2DHistogram(TH2F& hist, double x, double y, int color, double radiusX, double radiusY)
{
	TEllipse* point = new TEllipse(x,y,radiusX,radiusY); 
        point->SetLineStyle(1);
        point->SetFillColor(color);
        point->SetLineColor(color);
   	point->SetLineWidth(1);
   	
   	TCanvas *caddPointTo2DHistogram1 = new TCanvas();
   	
   	//gStyle->SetOptStat(0);
   	caddPointTo2DHistogram1->cd();
	caddPointTo2DHistogram1->SetLogy();
   	hist.Draw("BOX");
	point->Draw("same");
}

void addPointToCanvas(double x, double y, int color, double radiusX, double radiusY)
{
	TEllipse* point = new TEllipse(x,y,radiusX,radiusY); 
        point->SetLineStyle(1);
        point->SetFillColor(color);
        point->SetLineColor(color);
   	point->SetLineWidth(1);
	point->Draw("same");
}

void add2DHistToFile(TH2F& hist, string fileName, string histName)
{
	int color = hist.GetLineColor();
	TH2F *histClone = (TH2F*)hist.Clone(histName.c_str());
	histClone->SetName(histName.c_str());
	histClone->SetLineColor(color);
	TFile *file = TFile::Open(fileName.c_str(), "UPDATE");
	histClone->Write();
	file->Close();
	return;
}

void draw2DHistsFromFile(string fileName, string canvasScale, bool point, double x, double y, int pointColor, double radiusX, double radiusY)
{
	TFile *file = TFile::Open(fileName.c_str(), "READ");
	//gStyle->SetOptStat(0);
	TCanvas *cdraw2DHistsFromFile1 = new TCanvas();
	cdraw2DHistsFromFile1->cd();
	if(canvasScale == "log") cdraw2DHistsFromFile1->SetLogy();
	TIter next(file->GetListOfKeys());
    	TKey *key;
    	bool isFirstHistogram = true;
    	while ((key = (TKey*)next()))
    	{
    		if (strstr(key->GetClassName(), "TH2F"))
    		{
    			TH2F *hist = (TH2F*)key->ReadObj();
    			string str(hist->GetName());
    			if(str == "ROCRejBDT") hist->SetFillColor(kRed);
    			else if(str == "ROCRejBDTInvMass") hist->SetFillColor(kRed-9);
    			else if(str == "ROCRejBDTG") hist->SetFillColor(kBlue);
    			else if(str == "ROCRejBDTGInvMass") hist->SetFillColor(kAzure+5);
    			else if(str == "ROCRejMLP") hist->SetFillColor(kViolet);
    			else if(str == "ROCRejMLPInvMass") hist->SetFillColor(kMagenta-4);
    			/*hist->SetFillColor(kRed);
    			hist->SetLineColor(kRed);*/
    			hist->Draw("BOX same");
    			/*hist->Draw(isFirstHistogram ? "" : "BOX");
    			isFirstHistogram = false;
    			int color = hist->GetFillColor();
    			cout<<endl<<endl<<endl<<"Hist name: "<<hist->GetName()<<"       Color: "<<color<<endl<<endl<<endl;*/
    		}
    	}
    	if(point == true) addPointToCanvas(x, y, pointColor, radiusX, radiusY);
}


void findCrossSectionHHbbbb(double totalRemaining, double HHRemaining, double backRemaining, double luminosity, double& crossSection, double& errorTop, double& errorBottom, int nbin)
{ 
	double xMin = -2000.0;   // Minimum x value
   	double xMax = 10000.0;  // Maximum x value
	
	double branchingRatioHbb=0.5824; ////////////at Higgs mass of 125.00 GeV
	double totalBR = branchingRatioHbb * branchingRatioHbb; 
	double totalEffi=HHRemaining/(1462612*0.001225*totalBR);
	//cout<<"HHRemaining: "<<HHRemaining<<"    (1462612*0.001225): "<<1462612*0.001225<<"    totalEffi: "<<totalEffi<<endl;
	
	double s, b=backRemaining, n=totalRemaining; 
	cout<<"b: "<<b<<"    n: "<<n<<"    n-b: "<<n-b<<endl;
	double Lb, Lsb;
	double chiSquared, minChiSquared, minX;
	int contDebug=0;
	
	////Loop to find minChiSquared
	cout<<"Empieza minChiSquared"<<endl;
	for(double x=xMin; x<=xMax; x+=.01)
	{
		//if(contDebug<5) cout<<"Entra"<<endl<<"x: "<<x<<"    totalEffi: "<<totalEffi<<"    totalBR: "<<totalBR<<endl;
		s=x*totalEffi*totalBR;
		//if(contDebug<5) cout<<"s: "<<s<<"    b: "<<b<<"    n: "<<n<<endl;
      		//chiSquared=-2*TMath::Log((pow(TMath::E(), -(s+b)) * pow((s+b), n))/(pow(TMath::E(), -(b)) * pow((b), n)));
      		chiSquared = -2*(-s+n*TMath::Log((s+b)/b));
      		//if(contDebug<5) cout<<"Lsb: "<<Lsb<<"    Lb: "<<Lb<<"    ChiSquared: "<<chiSquared<<endl<<endl;
	      	if(chiSquared<minChiSquared)
	      	{
	      		minChiSquared=chiSquared;
	      		minX=x;
	      	}
	      	contDebug++;
	}
	TH2F *histChiSquared = new TH2F("hChiSquared", "ChiSquared as a function of cross-section * lumi", nbin, xMin, xMax, nbin, minChiSquared*1.1, 20);
	///Loop to plot chiSquared
	cout<<"Termina minChiSquared e inicia plot chiSquared"<<endl;
	for(double x=xMin; x<=xMax; x+=.01)
	{
		s=x*totalEffi*totalBR;
      		//chiSquared=-2*TMath::Log((pow(TMath::E(), -(s+b)) * pow((s+b), n))/(pow(TMath::E(), -(b)) * pow((b), n)));
      		chiSquared = -2*(-s+n*TMath::Log((s+b)/b));
      		if (chiSquared != 1000000000) histChiSquared->Fill(x, chiSquared);
	}
	cout<<"Termina plot chiSquared"<<endl;
		
	TCanvas *cfindCrossSectionHHbbbb1 = new TCanvas();
	cfindCrossSectionHHbbbb1->cd();
	histChiSquared->Draw();
	
	crossSection = minX / luminosity;
	
	cout<<"Signal efficiency: "<<totalEffi<<endl<<"Branching ratio: "<<totalBR<<endl;
	cout<<"min ChiSquared: "<<minChiSquared<<"     for x: "<<minX<<endl<<"Cross section: "<<crossSection<<endl;
	
	////To find uncertainty
	TLine* lineError = new TLine(xMin,minChiSquared+1,xMax,minChiSquared+1);
	lineError->SetLineColor(kRed);
	lineError->SetLineWidth(2);
	lineError->SetLineStyle(1);
	lineError->Draw("same");
	double xErrorLeft, errorLeft, distanceMinLeft=100000, xErrorRight, errorRight, distanceMinRight=100000; 
	for(double x=xMin; x<=xMax; x+=0.1)
	{
		s=x*totalEffi*totalBR;
      		//chiSquared=-2*TMath::Log((pow(TMath::E(), -(s+b)) * pow((s+b), n))/(pow(TMath::E(), -(b)) * pow((b), n)));
      		chiSquared = -2*(-s+n*TMath::Log((s+b)/b));
	      	if(chiSquared<minChiSquared+1.1 && chiSquared>minChiSquared+0.9)
	      	{
	      		if(x<minX)
	      		{
	      			if(abs(chiSquared-(minChiSquared+1))<distanceMinLeft)
	      			{
	      				distanceMinLeft=abs(chiSquared-(minChiSquared+1));
	      				xErrorLeft=x;
	      			}
	      		}
	      		else
	      		{
	      			if(abs(chiSquared-(minChiSquared+1))<distanceMinRight)
	      			{
	      				distanceMinRight=abs(chiSquared-(minChiSquared+1));
	      				xErrorRight=x;
	      			}
	      		}
	      	}
	      	contDebug++;
	}
	errorLeft=abs((xErrorLeft-minX)/(luminosity*crossSection));	
	errorRight=abs((xErrorRight-minX)/(luminosity*crossSection));	
	cout<<"errorLeft: "<<errorLeft<<"    for x: "<<xErrorLeft<<"    with distance: "<<distanceMinLeft<<endl;
	cout<<"errorRight: "<<errorRight<<"    for x: "<<xErrorRight<<"    with distance: "<<distanceMinRight<<endl<<endl;
	////////
	
}

/////Function that generates the files to train NN to find optimal cut for NNs outputs
void generateFilesOutputNN(vector<double>& BDTqqOutput, vector<double>& BDTttbarOutput, vector<double>& BDTZZOutput, vector<double>& BDTWWOutput, int sizeHH, int sizeqq, int sizettbar, int sizeZZ, int sizeWW)
{
	TFile *outputTreeSNNTrain = new TFile("outputTreeSNNTrain.root", "recreate");
 	TTree TreeSNNTrain("TreeSNNTrain","a simple Tree with simple variables (Train)");
  	TFile *outputTreeSNNTest = new TFile("outputTreeSNNTest.root", "recreate");
  	TTree TreeSNNTest("TreeSNNTest","a simple Tree with simple variables (Test)");
  	TFile *outputTreeBqqNNTrain = new TFile("outputTreeBqqNNTrain.root", "recreate");
 	TTree TreeBqqNNTrain("TreeBqqNNTrain","a bqqimple Tree with simple variables (Train)");
  	TFile *outputTreeBqqNNTest = new TFile("outputTreeBqqNNTest.root", "recreate");
  	TTree TreeBqqNNTest("TreeBqqNNTest","a bqqimple Tree with simple variables (Test)");
  	TFile *outputTreeBttbarNNTrain = new TFile("outputTreeBttbarNNTrain.root", "recreate");
 	TTree TreeBttbarNNTrain("TreeBttbarNNTrain","a bttbarimple Tree with simple variables (Train)");
  	TFile *outputTreeBttbarNNTest = new TFile("outputTreeBttbarNNTest.root", "recreate");
  	TTree TreeBttbarNNTest("TreeBttbarNNTest","a bttbarimple Tree with simple variables (Test)");
  	TFile *outputTreeBZZNNTrain = new TFile("outputTreeBZZNNTrain.root", "recreate");
 	TTree TreeBZZNNTrain("TreeBZZNNTrain","a bZZimple Tree with simple variables (Train)");
  	TFile *outputTreeBZZNNTest = new TFile("outputTreeBZZNNTest.root", "recreate");
  	TTree TreeBZZNNTest("TreeBZZNNTest","a bZZimple Tree with simple variables (Test)");
  	TFile *outputTreeBWWNNTrain = new TFile("outputTreeBWWNNTrain.root", "recreate");
 	TTree TreeBWWNNTrain("TreeBWWNNTrain","a bWWimple Tree with simple variables (Train)");
  	TFile *outputTreeBWWNNTest = new TFile("outputTreeBWWNNTest.root", "recreate");
  	TTree TreeBWWNNTest("TreeBWWNNTest","a bWWimple Tree with simple variables (Test)");
  	
  	float NN1Output, NN2Output, NN3Output, NN4Output;
  	TreeSNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeSNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeSNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeSNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeBqqNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBqqNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBqqNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBqqNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeBttbarNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBttbarNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBttbarNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBttbarNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeBZZNNTrain.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBZZNNTrain.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBZZNNTrain.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBZZNNTrain.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	
  	TreeSNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeSNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeSNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeSNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeBqqNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBqqNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBqqNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBqqNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeBttbarNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBttbarNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBttbarNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBttbarNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	TreeBZZNNTest.Branch("NN1Output",&NN1Output,"NN1Output/F");
  	TreeBZZNNTest.Branch("NN2Output",&NN2Output,"NN2Output/F");
  	TreeBZZNNTest.Branch("NN3Output",&NN3Output,"NN3Output/F");
  	TreeBZZNNTest.Branch("NN4Output",&NN4Output,"NN4Output/F");
  	
  	for(int i=0; i<BDTqqOutput.size(); i++)
  	{
  		NN1Output = BDTqqOutput[i];
  		NN2Output = BDTttbarOutput[i];
  		NN3Output = BDTZZOutput[i];
  		NN4Output = BDTWWOutput[i];
  		
  		random_device rd;  // Will be used to obtain a seed for the random number engine
		mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
		uniform_int_distribution<> distrib(1, 4); // Uniform distribution between 1 and 4
		if(distrib(gen) != 1)
		{ 
			if(i<sizeHH) TreeSNNTrain.Fill();
			else if(i<sizeHH+sizeqq) TreeBqqNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar) TreeBttbarNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ) TreeBZZNNTrain.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW) TreeBWWNNTrain.Fill();
		}
		else
		{ 
			if(i<sizeHH) TreeSNNTest.Fill();
			else if(i<sizeHH+sizeqq) TreeBqqNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar) TreeBttbarNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ) TreeBZZNNTest.Fill();
			else if(i<sizeHH+sizeqq+sizettbar+sizeZZ+sizeWW) TreeBWWNNTest.Fill();
		}
						
  	}
  	
  	outputTreeSNNTrain->cd();
   	TreeSNNTrain.Write();
   	outputTreeBqqNNTrain->cd();
   	TreeBqqNNTrain.Write();
   	outputTreeBttbarNNTrain->cd();
   	TreeBttbarNNTrain.Write();
   	outputTreeBZZNNTrain->cd();
   	TreeBZZNNTrain.Write();
   	outputTreeBWWNNTrain->cd();
   	TreeBWWNNTrain.Write();
   	
   	outputTreeSNNTest->cd();
   	TreeSNNTest.Write();
   	outputTreeBqqNNTest->cd();
   	TreeBqqNNTest.Write();
   	outputTreeBttbarNNTest->cd();
   	TreeBttbarNNTest.Write();
   	outputTreeBZZNNTest->cd();
   	TreeBZZNNTest.Write();
   	outputTreeBWWNNTest->cd();
   	TreeBWWNNTest.Write();
   	
}


 ///Function that finds the limits of the NN output
 void findHistLimits(string method, double& bottomHistLimit, double& topHistLimit, int& methodColor)
 {
 	if(method == "Likelihood")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "LikelihoodD")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 0.9999;
 	}
 	if(method == "LikelihoodPCA")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "LikelihoodKDE")
 	{
 		bottomHistLimit = -0.00001;
 		topHistLimit = 0.99999;
 	}
 	if(method == "LikelihoodMIX")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "PDERS")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "PDERSD")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "PDERSPCA")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "KNN")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	if(method == "HMatrix")
 	{
 		bottomHistLimit = -0.95;
 		topHistLimit = 1.55;
 	}
 	if(method == "Fisher")
 	{
 		bottomHistLimit = -4.0;
 		topHistLimit = 4.0;
 	}
 	else if(method == "FisherG")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "BoostedFisher")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "LD")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "MLP")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 		methodColor = kViolet;
 	}
 	else if(method == "MLPBFGS")
 	{
 		bottomHistLimit = -1.25;
 		topHistLimit = 1.5;
 	}
 	else if(method == "MLPBNN")
 	{
 		bottomHistLimit = -1.25;
 		topHistLimit = 1.5;
 	}
 	else if(method == "CFMlpANN")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "TMlpANN")
 	{
 		bottomHistLimit = -1.3;
 		topHistLimit = 1.3;
 	}
 	else if(method == "DNN_GPU")
 	{
 		bottomHistLimit = -0.1;
 		topHistLimit = 1.1;
 	}
 	else if(method == "DNN_CPU")
 	{
 		bottomHistLimit = -0.1;
 		topHistLimit = 1.1;
 	}
 	else if(method == "BDT")
 	{
 		bottomHistLimit = -0.8;
 		topHistLimit = 0.8;
 		methodColor = kRed;
 	}
 	else if(method == "BDTG")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 		methodColor = kBlue;
 	}
 	else if(method == "BDTB")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "BDTD")
 	{
 		bottomHistLimit = -0.8;
 		topHistLimit = 0.8;
 	}
 	else if(method == "BDTF")
 	{
 		bottomHistLimit = -1.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "RuleFit")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "SVM_Gauss")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "SVM_Poly")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "SVM_Lin")
 	{
 		bottomHistLimit = 0.0;
 		topHistLimit = 1.0;
 	}
 	else if(method == "FDA_MT")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 3.0;
 	}
 	else if(method == "FDA_GA")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 3.0;
 	}
 	else if(method == "Category")
 	{
 		bottomHistLimit = -2.0;
 		topHistLimit = 2.0;
 	}
 	else if(method == "Plugin")
 	{
 		bottomHistLimit = -0.8;
 		topHistLimit = 0.8;
 	}
 	else cout<<"ERROOOOOOOOOOOOOOOOOOR UNKNOWN METHOD"<<endl;
 	
 	return;
 }
 
 void FSRTMVAClassificationApplicationHHbbbb()
 {
 	cout<<"Aquì empieza main()"<<endl;
 	gStyle->SetOptStat(0);
 	
 	string method = "BDTG";
 	string variables = "All1V34BSplit";
 	string rtdCut = "10";
 	
 	string inputTrainHHText="outputTreeSHHbbbbESpreadDurham"+rtdCut+"Train.root";
 	string inputTestHHText="outputTreeSHHbbbbESpreadDurham"+rtdCut+"Test.root";
 	TFile* inputTrainHH = new TFile(inputTrainHHText.c_str());
    	TTree* theTreeTrainHH = (TTree*)inputTrainHH->Get("TreeSTrain");
 	TFile* inputTestHH = new TFile(inputTestHHText.c_str());
    	TTree* theTreeTestHH = (TTree*)inputTestHH->Get("TreeSTest");
    	double weightFactorHH = (theTreeTestHH->GetEntries()+theTreeTrainHH->GetEntries());
    	weightFactorHH = weightFactorHH/(theTreeTestHH->GetEntries());
    	string inputTrainqqText="outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Train.root";
 	string inputTestqqText="outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
    	TFile* inputTrainqq = new TFile(inputTrainqqText.c_str());
    	TTree* theTreeTrainqq = (TTree*)inputTrainqq->Get("TreeBqqTrain");
 	TFile* inputTestqq = new TFile(inputTestqqText.c_str());
    	TTree* theTreeTestqq = (TTree*)inputTestqq->Get("TreeBqqTest");
    	double weightFactorqq = (theTreeTestqq->GetEntries()+theTreeTrainqq->GetEntries());
    	weightFactorqq = weightFactorqq/(theTreeTestqq->GetEntries());
    	string inputTrainttbarText="outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Train.root";
 	string inputTestttbarText="outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Test.root";
    	TFile* inputTrainttbar = new TFile(inputTrainttbarText.c_str());
    	TTree* theTreeTrainttbar = (TTree*)inputTrainttbar->Get("TreeBttTrain");
 	TFile* inputTestttbar = new TFile(inputTestttbarText.c_str());
    	TTree* theTreeTestttbar = (TTree*)inputTestttbar->Get("TreeBttTest");
    	double weightFactorttbar = (theTreeTestttbar->GetEntries()+theTreeTrainttbar->GetEntries());
    	weightFactorttbar = weightFactorttbar/(theTreeTestttbar->GetEntries());
    	string inputTrainZZText="outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Train.root";
 	string inputTestZZText="outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Test.root";
    	TFile* inputTrainZZ = new TFile(inputTrainZZText.c_str());
    	TTree* theTreeTrainZZ = (TTree*)inputTrainZZ->Get("TreeBZZTrain");
 	TFile* inputTestZZ = new TFile(inputTestZZText.c_str());
    	TTree* theTreeTestZZ = (TTree*)inputTestZZ->Get("TreeBZZTest");
    	double weightFactorZZ = (theTreeTestZZ->GetEntries()+theTreeTrainZZ->GetEntries());
    	weightFactorZZ = weightFactorZZ/(theTreeTestZZ->GetEntries());
    	string inputTrainWWText="outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Train.root";
 	string inputTestWWText="outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Test.root";
    	TFile* inputTrainWW = new TFile(inputTrainWWText.c_str());
    	TTree* theTreeTrainWW = (TTree*)inputTrainWW->Get("TreeBWWTrain");
 	TFile* inputTestWW = new TFile(inputTestWWText.c_str());
    	TTree* theTreeTestWW = (TTree*)inputTestWW->Get("TreeBWWTest");
    	double weightFactorWW = (theTreeTestWW->GetEntries()+theTreeTrainWW->GetEntries());
    	weightFactorWW = weightFactorWW/(theTreeTestWW->GetEntries());
    	cout<<"WeightFactorHH: "<<weightFactorHH<<endl;
    	cout<<"WeightFactorqq: "<<weightFactorqq<<endl;
    	cout<<"WeightFactorttbar: "<<weightFactorttbar<<endl;
    	cout<<"WeightFactorZZ: "<<weightFactorZZ<<endl;
    	cout<<"WeightFactorWW: "<<weightFactorWW<<endl;
    	
 	
 	double bottomHistLimit=0.0, topHistLimit=0.0;
 	int methodColor;
 	findHistLimits(method, bottomHistLimit, topHistLimit, methodColor);
 	
 	int sizeHH=0, sizeqq=0, sizettbar=0, sizeZZ=0, sizeWW=0;
 	int nbin=10000, nbinNN=100;
 	double weightHH=0.001225*weightFactorHH, weightqq=0.0349*weightFactorqq, weightttbar=0.503*weightFactorttbar, weightZZ=0.8167*weightFactorZZ, weightWW=0.5149*weightFactorWW;
 	vector<double> BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput;
 	
 	
 	TH1F *histBdtqqHH     = new TH1F( "MVA_BDTqqHH",           "MVA_BDTqq on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtqqqq     = new TH1F( "MVA_BDTqqqq",           "MVA_BDTqq on qq",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarHH     = new TH1F( "MVA_BDTttbarHH",           "MVA_BDTttbar on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarttbar     = new TH1F( "MVA_BDTttbarttbar",           "MVA_BDTttbar on ttbar",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtZZHH     = new TH1F( "MVA_BDTZZHH",           "MVA_BDTZZ on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtZZZZ     = new TH1F( "MVA_BDTZZZZ",           "MVA_BDTZZ on ZZ",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtWWHH     = new TH1F( "MVA_BDTWWHH",           "MVA_BDTWW on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtWWWW     = new TH1F( "MVA_BDTWWWW",           "MVA_BDTWW on WW",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtqqHHW    = new TH1F( "MVA_BDTqqHHW",           "MVA_BDTqq on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtqqqqW     = new TH1F( "MVA_BDTqqqqW",           "MVA_BDTqq on qq",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarHHW     = new TH1F( "MVA_BDTttbarHHW",           "MVA_BDTttbar on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtttbarttbarW     = new TH1F( "MVA_BDTttbarttbarW",           "MVA_BDTttbar on ttbar",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtZZHHW     = new TH1F( "MVA_BDTZZHHW",           "MVA_BDTZZ on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtZZZZW     = new TH1F( "MVA_BDTZZZZW",           "MVA_BDTZZ on ZZ",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histBdtWWHHW     = new TH1F( "MVA_BDTWWHHW",           "MVA_BDTWW on HH",           nbinNN, bottomHistLimit-.05, topHistLimit+.05);
 	TH1F *histBdtWWWWW     = new TH1F( "MVA_BDTWWWWW",           "MVA_BDTWW on WW",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	TH1F *histDummy     = new TH1F( "histDummy",           "histDummy",           nbinNN, bottomHistLimit-.05, topHistLimit+.05 );
 	
 	
 	histBdtqqHH->SetLineColor(kBlue);
 	histBdtqqqq->SetLineColor(kRed+2);
 	histBdtttbarHH->SetLineColor(kBlue);
 	histBdtttbarttbar->SetLineColor(kPink+6);
 	histBdtZZHH->SetLineColor(kBlue);
 	histBdtZZZZ->SetLineColor(kGreen+1);
 	histBdtWWHH->SetLineColor(kBlue);
 	histBdtWWWW->SetLineColor(kYellow-2);
 	histBdtqqHHW->SetLineColor(kBlue);
 	histBdtqqqqW->SetLineColor(kRed+2);
 	histBdtttbarHHW->SetLineColor(kBlue);
 	histBdtttbarttbarW->SetLineColor(kPink+6);
 	histBdtZZHHW->SetLineColor(kBlue);
 	histBdtZZZZW->SetLineColor(kGreen+1);
 	histBdtWWHHW->SetLineColor(kBlue);
 	histBdtWWWWW->SetLineColor(kYellow-2);
 
 	histBdtqqHH->SetFillColor(kBlue);
 	histBdtqqqq->SetFillColor(kRed+2);
 	histBdtttbarHH->SetFillColor(kBlue);
 	histBdtttbarttbar->SetFillColor(kPink+6);
 	histBdtZZHH->SetFillColor(kBlue);
 	histBdtZZZZ->SetFillColor(kGreen+1);
 	histBdtWWHH->SetFillColor(kBlue);
 	histBdtWWWW->SetFillColor(kYellow-2);
 	histBdtqqHHW->SetFillColor(kBlue);
 	histBdtqqqqW->SetFillColor(kRed+2);
 	histBdtttbarHHW->SetFillColor(kBlue);
 	histBdtttbarttbarW->SetFillColor(kPink+6);
 	histBdtZZHHW->SetFillColor(kBlue);
 	histBdtZZZZW->SetFillColor(kGreen+1);
 	histBdtWWHHW->SetFillColor(kBlue);
 	histBdtWWWWW->SetFillColor(kYellow-2);
 	
 	double targetFractionHH=0.8;
 	double effiHH=0, effiqq=0, effittbar=0, effiZZ=0, effiWW=0;
 	//pairMassCut(1, targetFractionHH, effiHH);
 	//pairMassCut(2, targetFractionHH, effiqq);
 	//pairMassCut(3, targetFractionHH, effittbar);
 	
 	/////For arguments: empty string, topology (1: HH, 2: qq, 3:ttbar, 4: ZZ, 5: WW), nBack (0: qq, 1: ttbar, 2: ZZ, 3: WW), vector for BDTqq output, vector for BDTttbar output, var for size of tree, nbin for hists, histogram for BDT output, weighted histogram for BDT output, weight, method, NNVars (e.g. All if using all variables, PairMasses if only using the pair masses, etc)
 	//////////
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeHH, nbinNN, *histBdtqqHH, *histBdtqqHHW, weightHH, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeHH, nbinNN, *histBdtttbarHH, *histBdtttbarHHW, weightHH, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeHH, nbinNN, *histBdtZZHH, *histBdtZZHHW, weightHH, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 1, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeHH, nbinNN, *histBdtWWHH, *histBdtWWHHW, weightHH, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeqq, nbinNN, *histBdtqqqq, *histBdtqqqqW, weightqq, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 2, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeqq, nbinNN, *histDummy, *histDummy, weightqq, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizettbar, nbinNN, *histBdtttbarttbar, *histBdtttbarttbarW, weightttbar, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 3, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizettbar, nbinNN, *histDummy, *histDummy, weightttbar, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightZZ, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightZZ, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeZZ, nbinNN, *histBdtZZZZ, *histBdtZZZZW, weightZZ, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 4, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeZZ, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 0, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 1, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 2, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeWW, nbinNN, *histDummy, *histDummy, weightWW, method, variables, rtdCut);
 	FSRTMVAClassificationApplicationHHbbbbHelper("", 5, 3, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeWW, nbinNN, *histBdtWWWW, *histBdtWWWWW, weightWW, method, variables, rtdCut);
 	//////////////
 	
 	cout<<"BDTqqOutput size: "<<BDTqqOutput.size()<<endl;
 	cout<<"BDTttbarOutput size: "<<BDTttbarOutput.size()<<endl;
 	cout<<"BDTZZOutput size: "<<BDTZZOutput.size()<<endl;
 	cout<<"BDTWWOutput size: "<<BDTWWOutput.size()<<endl;
 	
 	TH2F *histBDTsOutputHH = new TH2F("btsoutputHH", "NNqq vs NNttbar", 100, bottomHistLimit-0.1, topHistLimit+0.1, 100, bottomHistLimit-0.1, topHistLimit+0.1);
 	histBDTsOutputHH->GetXaxis()->SetTitle("NNqq output");
   	histBDTsOutputHH->GetYaxis()->SetTitle("NNttbar output");
 	//BDTOutputPlane(bottomHistLimit, topHistLimit, BDTqqOutput, BDTttbarOutput, sizeHH, sizeqq, sizettbar, *histBDTsOutputHH);
 	
 	///////////////////////
 	BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqHH, kBlue, *histBdtqqqq, kRed+2, method + "qq on HH and qq (unweighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtqqHHW, kBlue, *histBdtqqqqW, kRed+2, method + "qq on HH and qq (weighted)");
 	BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtttbarHH, kBlue, *histBdtttbarttbar, kPink+6, method + "ttbar on HH and ttbar (unweighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtttbarHHW, kBlue, *histBdtttbarttbarW, kPink+6, method + "ttbar on HH and ttbar (weighted)");
 	BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtZZHH, kBlue, *histBdtZZZZ, kGreen+1, method + "ZZ on HH and ZZ (unweighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtZZHHW, kBlue, *histBdtZZZZW, kGreen+1, method + "ZZ on HH and ZZ (weighted)");
 	BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtWWHH, kBlue, *histBdtWWWW, kYellow-2, method + "WW on HH and WW (unweighted)");
 	//BDTOutputOverlap(bottomHistLimit, topHistLimit, nbinNN, *histBdtWWHHW, kBlue, *histBdtWWWWW, kYellow-2, method + "WW on HH and WW (weighted)");
 	/////////////////////////
 	
 	//double n, b, optCut;
 	//double maxSignificance=0;
 	//maxSignificance = findMaxSignificanceAndROC(bottomHistLimit, topHistLimit, nbin, *histBdtqqHHW, *histBdtqqqqW, n, b, optCut);
 	
 	
 	double maxSignificanceqq=0, defCutqq, totalRemainingqq=0, HHRemainingqq=0, backRemainingqq=0;
 	double maxSignificancettbar=0, defCutttbar, totalRemainingttbar=0, HHRemainingttbar=0, backRemainingttbar=0;
 	double maxSignificanceZZ=0, defCutZZ, totalRemainingZZ=0, HHRemainingZZ=0, backRemainingZZ=0;
 	double maxSignificanceWW=0, defCutWW, totalRemainingWW=0, HHRemainingWW=0, backRemainingWW=0;
 	double maxSignificanceCombined=0, defCutqqCombined, defCutttbarCombined=0, totalRemainingCombined=0, HHRemainingCombined=0, backRemainingCombined=0;
 	
 	TH2F *histROCqq = new TH2F("hROCqq", "Signal and Background Acceptance (NNqq)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCttbar = new TH2F("hROCttbar", "Signal and Background Acceptance (NNttbar)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCZZ = new TH2F("hROCZZ", "Signal and Background Acceptance (NNZZ)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCWW = new TH2F("hROCWW", "Signal and Background Acceptance (NNWW)", nbin, 0.0, 1.0, nbin, 0, 1.0);
 	TH2F *histROCqqttbar = new TH2F("hROCqqttbar", "Signal and Background Acceptance (NNqq and NNttbar)", nbin, 0.0, 1.0, nbin, 0, 1.0);
  	TH2F *histROCRejqq = new TH2F("hROCRejqq", "Signal Acceptance and Background Rejection (NNqq)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejttbar = new TH2F("hROCRejttbar", "Signal Acceptance and Background Rejection (NNttbar)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejZZ = new TH2F("hROCRejZZ", "Signal Acceptance and Background Rejection (NNZZ)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejWW = new TH2F("hROCRejWW", "Signal Acceptance and Background Rejection (NNWW)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	TH2F *histROCRejqqttbar = new TH2F("hROCRejqqttbar", "Signal Acceptance and Background Rejection (NNqq and NNttbar)", nbin, 0.5, 1.0, nbin, 0, 100000.0);
  	
  	TH2F *histSignificanceqq     = new TH2F( "hist_significanceqq", "Significance s/sqrt(s+b) (NNqq)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceqq->GetXaxis()->SetTitle("Cut value for NNqq");
   	histSignificanceqq->GetYaxis()->SetTitle("Significance for NNqq");
   	TH2F *histSignificancettbar     = new TH2F( "hist_significancettbar", "Significance s/sqrt(s+b) (NNttbar)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificancettbar->GetXaxis()->SetTitle("Cut value for NNttbar)");
   	histSignificancettbar->GetYaxis()->SetTitle("Significance for NNttbar");
   	TH2F *histSignificanceZZ     = new TH2F( "hist_significanceZZ", "Significance s/sqrt(s+b) (NNZZ)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceZZ->GetXaxis()->SetTitle("Cut value for NNZZ");
   	histSignificanceZZ->GetYaxis()->SetTitle("Significance for NNZZ");
   	TH2F *histSignificanceWW     = new TH2F( "hist_significanceWW", "Significance s/sqrt(s+b) (NNWW)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceWW->GetXaxis()->SetTitle("Cut value for NNWW");
   	histSignificanceWW->GetYaxis()->SetTitle("Significance for NNWW");
   	TH2F *histSignificanceqqttbar     = new TH2F( "hist_significanceqqttbar", "Significance s/sqrt(s+b) (NNqq and NNttbar)", nbin, bottomHistLimit-.05, topHistLimit+.05, nbin, 0, 20.0);
  	histSignificanceqqttbar->GetXaxis()->SetTitle("Cut value for NNqq");
   	histSignificanceqqttbar->GetYaxis()->SetTitle("Significance for NNqq"); //////////CHECAR COMO DIMENSIONALIZAR
   	
   	///////////////////
   	cout<<endl<<"Significance for "<<method<<"qq applied only to HH and qq: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, 0, *histBdtqqHH, *histBdtqqqq, defCutqq, maxSignificanceqq, weightHH, weightqq, targetFractionHH, *histROCqq, *histROCRejqq, *histSignificanceqq, totalRemainingqq, HHRemainingqq, backRemainingqq);
 	cout<<"maxSignificance (qq): "<<maxSignificanceqq<<endl<<"Cut in NNqq for max significance (qq): "<<defCutqq<<endl;
 	
 	cout<<endl<<"Significance for "<<method<<"ttbar applied only to HH and ttbar: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTttbarOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, 1, *histBdtttbarHH, *histBdtttbarttbar, defCutttbar, maxSignificancettbar, weightHH, weightttbar, targetFractionHH, *histROCttbar, *histROCRejttbar, *histSignificancettbar, totalRemainingttbar, HHRemainingttbar, backRemainingttbar);
 	cout<<"maxSignificance (ttbar): "<<maxSignificancettbar<<endl<<"Cut in NNttbar for max significance (ttbar): "<<defCutttbar<<endl;
 	
 	cout<<endl<<"Significance for "<<method<<"ZZ applied only to HH and ZZ: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTZZOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, 2, *histBdtZZHH, *histBdtZZZZ, defCutZZ, maxSignificanceZZ, weightHH, weightZZ, targetFractionHH, *histROCZZ, *histROCRejZZ, *histSignificanceZZ, totalRemainingZZ, HHRemainingZZ, backRemainingZZ);
 	cout<<"maxSignificance (ZZ): "<<maxSignificanceZZ<<endl<<"Cut in NNZZ for max significance (ZZ): "<<defCutZZ<<endl;
 	
 	cout<<endl<<"Significance for "<<method<<"WW applied only to HH and WW: "<<endl;
 	findSignificance(bottomHistLimit, topHistLimit, nbin, BDTWWOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, 3, *histBdtWWHH, *histBdtWWWW, defCutWW, maxSignificanceWW, weightHH, weightWW, targetFractionHH, *histROCWW, *histROCRejWW, *histSignificanceWW, totalRemainingWW, HHRemainingWW, backRemainingWW);
 	cout<<"maxSignificance (WW): "<<maxSignificanceWW<<endl<<"Cut in NNWW for max significance (WW): "<<defCutWW<<endl;
 	
 	///////////////////////
 	
 	//cout<<endl<<"(using crazy for loops) Significance for "<<method<<"qq and "<<method<<"ttbar combined applied to HH, qq and, ttbar: "<<endl;
 	//findSignificanceCombined(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, sizeHH, sizeqq, sizettbar, 0, 1, defCutqqCombined, defCutttbarCombined, maxSignificanceCombined, weightHH, weightqq, weightttbar, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
 	//cout<<"maxSignificance (combined): "<<maxSignificanceCombined<<endl<<"Cut in NNqq for max significance (combined): "<<defCutqqCombined<<endl<<"Cut in NNttbar for max significance (combined): "<<defCutttbarCombined<<endl;
 	
 	//////////////
 	defCutqq=0.636;
 	defCutttbar=0.802;
 	cout<<endl<<"Significance for "<<method<<"qq, "<<method<<"ttbar, "<<method<<"ZZ, and "<<method<<"WW combined applied to HH, qq, ttbar, ZZ, and WW (trained and applied independently): "<<endl;
 	findSignificanceCutsCombined(bottomHistLimit, topHistLimit, nbin, BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW, defCutqq, defCutttbar, defCutZZ, defCutWW, maxSignificanceCombined, weightHH, weightqq, weightttbar, weightZZ, weightWW, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
 	cout<<"maxSignificance (combined): "<<maxSignificanceCombined<<endl;
 	//////////////
 	
 	/*//////////DEBUGGING PURPOSES
 	vector<double> testNN1Output = {0.01, 0.61, 0.11, 0.81, 0.91, 0.98, 0.03, 0.933, 0.13, 0.23, 0.98, 0.85, 0.95, 0.92, 0.92};
 	vector<double> testNN2Output = {0.02, 0.03, 0.12, 0.82, 0.92, 0.98, 0.94, 0.95, 0.84, 0.94, 0.06, 0.066, 0.96, 0.93, 0.91};
 	vector<double> testNN3Output = {0.08, 0.16, 0.65, 0.89, 0.96, 0.98, 0.95, 0.96, 0.93, 0.82, 0.95, 0.84, 0.97, 0.91, 0.904};
 	//findSignificanceCombined(bottomHistLimit, topHistLimit, nbin, testNN1Output, testNN2Output, 6, 4, 3, 0, 1, defCutqqCombined, defCutttbarCombined, maxSignificanceCombined, weightHH, weightqq, weightttbar, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined); 
 	findSignificanceCutsCombined(bottomHistLimit, topHistLimit, nbin, testNN1Output, testNN2Output, testNN3Output,  6, 4, 3, 2, 0.6, 0.7, 0.8, maxSignificanceCombined, weightHH, weightqq, weightttbar, weightZZ, *histROCqqttbar, *histROCRejqqttbar, *histSignificanceqqttbar, totalRemainingCombined, HHRemainingCombined, backRemainingCombined);
 	/////////DEBUGGING PURPOSES*/
 	
 	
 	//add2DHistToFile(*histROCRej, "ROCHists.root", "ROCRej"+method);
 	
 	//splitHistogram(*histBdtqqqq, bottomHistLimit, topHistLimit, defCut, "BDTqq on qq remaining after cut");
 	
 	//addPointTo2DHistogram(*histROCRej, effiHH, 1/effiqq, kRed, .005, 1);
 	
 	///For arguments: (name of file with hists, scale (log or linear), bool for whether a point will be graffed, x coord, y coord, point color, radius x, radius y 
 	//draw2DHistsFromFile("ROCHists.root", "log", true, effiHH, 1/effiqq, kGreen, .005, 1);
 	
 	
 	
 	
 	
 	
 	double luminosity = 4900, crossSectionqq, errorTopqq, errorBottomqq;
 	double crossSectionttbar, errorTopttbar, errorBottomttbar;
 	double crossSectionZZ, errorTopZZ, errorBottomZZ;
 	double crossSectionWW, errorTopWW, errorBottomWW;
 	double crossSectionCombined, errorTopCombined, errorBottomCombined;
 	
 	cout<<endl<<endl;
 	
 	/*cout<<"Cross-section and error for the output of "<<method<<"qq applied only to HH and qq: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingqq+backRemainingqq, HHRemainingqq, backRemainingqq, luminosity, crossSectionqq, errorTopqq, errorBottomqq, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"ttbar applied only to HH and ttbar: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingttbar+backRemainingttbar, HHRemainingttbar, backRemainingttbar, luminosity, crossSectionttbar, errorTopttbar, errorBottomttbar, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"ZZ applied only to HH and ZZ: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingZZ+backRemainingZZ, HHRemainingZZ, backRemainingZZ, luminosity, crossSectionttbar, errorTopZZ, errorBottomZZ, nbin);
 	cout<<"Cross-section and error for the output of "<<method<<"WW applied only to HH and WW: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingWW+backRemainingWW, HHRemainingWW, backRemainingWW, luminosity, crossSectionttbar, errorTopWW, errorBottomWW, nbin);*/
 	cout<<"Cross-section and error for the output of "<<method<<" combined applied to HH, qq, ttbar, ZZ, and WW: "<<endl;
 	findCrossSectionHHbbbb(HHRemainingCombined+backRemainingCombined, HHRemainingCombined, backRemainingCombined, luminosity, crossSectionCombined, errorTopCombined, errorBottomCombined, nbin);
 	
 	
 	
 	
 	
 	
 	
 	/*TCanvas *cMain2 = new TCanvas();
 	TCanvas *cMain3 = new TCanvas();
 	TCanvas *cMain4 = new TCanvas();
 	TCanvas *cMain5 = new TCanvas();
 	
 	cMain2->cd();
 	histBdtqqHH->Draw("HIST");
 	
 	cMain3->cd();
 	histBdtqqqq->Draw("HIST");
 	
 	cMain4->cd();
 	histBdtttbarHH->Draw("HIST");
 	
 	cMain5->cd();
 	histBdtttbarttbar->Draw("HIST");*/
 	
 	
 	cout<<endl<<endl<<endl<<"General info of input files: "<<endl;
 	cout<<"w HH size: "<<sizeHH*weightHH<<endl<<"uw HH size: "<<sizeHH<<endl<<endl;
 	cout<<"w qq size: "<<sizeqq*weightqq<<endl<<"uw qq size: "<<sizeqq<<endl<<endl;
 	cout<<"w ttbar size: "<<sizettbar*weightttbar<<endl<<"uw ttbar size: "<<sizettbar<<endl<<endl;
 	cout<<"w ZZ size: "<<sizeZZ*weightZZ<<endl<<"uw ZZ size: "<<sizeZZ<<endl<<endl;
 	cout<<"w WW size: "<<sizeWW*weightWW<<endl<<"uw WW size: "<<sizeWW<<endl<<endl;
 	
 	
 	cout<<"BDTqqOutput size: "<<BDTqqOutput.size()<<endl;
 	cout<<"BDTttbarOutput size: "<<BDTttbarOutput.size()<<endl;
 	cout<<"BDTZZOutput size: "<<BDTZZOutput.size()<<endl;
 	cout<<"BDTWWOutput size: "<<BDTWWOutput.size()<<endl;
 	
 	
 	///////Filling out tree for output NN
 	//generateFilesOutputNN(BDTqqOutput, BDTttbarOutput, BDTZZOutput, BDTWWOutput, sizeHH, sizeqq, sizettbar, sizeZZ, sizeWW);
  	///////Filling out tree for output NN
 	
 	
 	
 	return;
 }
 
 /*int main( int argc, char** argv )
 {
    cout<<"Aquì empieza main()"<<endl;
    TString methodList;
    for (int i=1; i<argc; i++) {
       TString regMethod(argv[i]);
       if(regMethod=="-b" || regMethod=="--batch") continue;
       if (!methodList.IsNull()) methodList += TString(",");
       methodList += regMethod;
    }
    FSRTMVAClassificationApplicationHHbbbb(methodList);
    return 0;
 }*/
