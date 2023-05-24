#include "Anfangsverteilung.h"
#include "MersenneTwister.h"

#include <random>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstdlib>
# include <cstdio>
# include <vector>
# include <fstream>
#include <string>
#include <math.h>
#include <algorithm>  
#include <boost/filesystem.hpp>

using std::vector;
using std::cout;
using std::cin;


// Initiate random number generator
time_t seconds ;
MTRand mtrand1(time(NULL));
std::random_device rd;
std::mt19937 gen(rd());

// Define a function to generate a random number for the Gillespie Algorithm (Defined at the bottom)
double exprand( double lambda );

// Define a function to print vectors
void print(vector <double> const &a);

// Define a function to calculate total number of cells with ecDNA
int numCellsWithEcdna(int numCellsWithEcdnaA, int numCellsWithEcdnaB);
int numCellsWithEcdna(std::vector <std::vector <double> > v);

// Define a function to count number of nonzero entries in a vector
int getNonZeroSize (vector <double> v);

// Create a directory for simulation outputs 
void createOutputDir (std::string outputFolder);

// Define a function for ecDNA barcoding/labelling
vector < vector <double> > barcodeAll(vector <double> v);

// Define a function to evolve unlabelled ecDNA
void ecDNAEvolveWithoutLabels(int NumCells, int NumNeutral, int amplify, double fitness, int initialcopies, int runs, std::string outputFolder);

// Define a function to evolve labelled ecDNA
void ecDNAEvolveWithLabels(int NumCells, int NumNeutral, int amplify, double fitness, int initialcopies, int runs, std::string outputFolder);

int main(int argc, char* argv[])
{
    if (argc < 8) { // We expect 8 arguments: the program name, followed by 7 args
        std::cerr << "Usage: " << argv[0] << " NumCells(int) NumNeutral(int)"<<
        " amplify(int) fitness(double) initialcopies(int) runs(int) outputFolder(str)" << std::endl;
        return 1;
    }
    else {
    	// convert to correct types
    	int x1 = atoi(argv[1]);
    	int x2 = atoi(argv[2]);
    	int x3 = atoi(argv[3]);
    	double x4 = atof(argv[4]);
    	int x5 = atoi(argv[5]);
    	int x6 = atoi(argv[6]);
    	
    	cout << "Running simulations with NumCells="<<x1<<", NumNeutral="<<x2<<", amplify="<<x3<<", fitness="<<x4<<
    	", initialcopies="<<x5<<", runs="<<x6<<", outputFolder="<<argv[7]<<"\n";					

    	createOutputDir(argv[7]);
        ecDNAEvolveWithLabels(x1, x2, x3, x4, x5, x6, argv[7]);
    	return 0;
    }
}     

// Auxiliary functions
double exprand( double lambda)
{
    double y = mtrand1.randExc();
    return (-log(1.-y)/lambda);
}

int numCellsWithEcdna(int numCellsWithEcdnaA, int numCellsWithEcdnaB)
{
    if (numCellsWithEcdnaA != numCellsWithEcdnaB) return -1;
    return numCellsWithEcdnaA;
}

int numCellsWithEcdna(std::vector <std::vector <double> > v)
{
    // calculate number of cells with each ecDNA type
    std::vector <double> numCellsEcdnaType (v.size(), 0);
    for (int i=0; i<v.size(); i++){
        numCellsEcdnaType[i] = getNonZeroSize(v[i]);
    }
    // get maximum of vector numCellsEcdnaType; this is the number of cells with
    // any type of ecDNA
    std::vector<double>::iterator numCellsWithEcdna = std::max_element(numCellsEcdnaType.begin(), numCellsEcdnaType.end());
    return *numCellsWithEcdna;
}

int getNonZeroSize (std::vector<double> v)
{
    int nonZeroSize = 0;
    for (unsigned i=0; i<v.size(); i++) {
      if (v[i]!=0.0) nonZeroSize++;
  }
  return nonZeroSize;
}

void createOutputDir (std::string outputFolder)
{

    // std::string outputFolder = "../exps/";
    const char* path = outputFolder.c_str();
    boost::filesystem::path dir(path);
    if(boost::filesystem::create_directory(dir))
    {
        std::cerr<< "Directory Created: "<<outputFolder<<std::endl;
    }
    
}

vector < vector <double> > barcodeAll (vector <double> v)
{
    // This function assigns a unique label to each ecDNA in the current population of cells. 
   
    int numCells = v.size();
    // cout<<"\nTotal number of cells = "<<N<<"\n";

    // compute total number of labels
    int numLabels = std::accumulate(v.begin(), v.end(), decltype(v)::value_type(0));
    // cout<<"Total number of labels = "<<T<<"\n";

    // initialize vector of vectors for labelled data
    vector< vector<double> > vLabelled (numLabels, vector<double> (numCells,0));
    // for(int k=0; k<numLabels; k++) {
        // cout<<"v["<<k<<"] = ";
        // for(int i=0; i<numCells; i++){
        //     cout<<vLabelled[k][i]<<" ";
        // }
        // cout<<"\n";
    // }

    // do labelling for i=0
    for(int k=0; k<numLabels+1; k++){
        if(k < v[0]){
            vLabelled[k][0] = 1;
        }
    }
    // do labelling for i>0
    for(int i=1; i<numCells; i++){
        int sum2 = 0;
        for(int j=0; j<i+1; j++){
           sum2 += v[j];
          }
          int sum1 = sum2 - v[i];
        for(int k=0; k<numLabels+1; k++){
           if(k>sum1 && k<=sum2){
                vLabelled[k-1][i] = 1;
            }
        }
     }
    return vLabelled;
    // print labelled vectors
    // for(int k=0; k<T; k++) {
    //     for(int i=0; i<N; i++){
    //         cout<<sL[k][i]<<" ";
    //     }
    //     cout<<"\n";
   // }
}

void ecDNAEvolveWithoutLabels(int NumCells, int NumNeutral, int amplify, double fitness, int initialcopies, int runs, std::string outputFolder)
{
   // Initiate a bunch of vectors to store cell states throughout the simulation
   vector <double> State (1,initialcopies); // this initializes a vector of size 1, with value=initialcopies  
   vector < vector <double> > FinalOutput (runs ,vector <double> (NumCells+1,0)); // final number of ecDNA in each cell
   vector < vector <double> > Neutral (runs ,vector <double> (NumCells,0));
   double pFrac; // fraction of cells positive for ecDNA
   double nFrac; // fraction of cells negative for ecDNA (i.e. neutral cells)
   std::fstream dataFracs; // file to store cell fractions
   double totalN; // total number of cells

   std::string fractionsBaseFileName = "cellFractions_";
   std::string fractionsFileName = outputFolder + fractionsBaseFileName + std::to_string((int)fitness) + "_" + std::to_string(initialcopies) + ".txt";
   dataFracs.open (outputFolder+fractionsFileName ,std::ios::out);
   //vector < vector <double> > cellFractions (runs, NumCells+1, aFrac, bFrac, nFrac); // how to best store all this stuff???? can vector be more than 1-dimensional?
   
   int count1 = 0;  // Dummy variable to count number of simulation repeats
    
    while (count1 < runs)
    {
	int count = 0; // Dummy Variable to count until the number of maximal cells
               
        // Reset some of the vectors of the simulation
        State.resize(1);
        State.at(0) = initialcopies;
        NumNeutral = 0;
        

        while (count < NumCells)
        {
            // cout << count1 << " " << count << "\n"; // output current state of the simulation
            
            // Compute fractions of cells with each type of ecDNA
            
            totalN = count+1; 
            pFrac = getNonZeroSize(State)/totalN;
            nFrac = NumNeutral/totalN; 

            // Write fractions to file
            dataFracs << count1 << " " << totalN << " " << pFrac << " " << nFrac << "\n";

            // Define a bunch of dummy variables to store intermediate ecDNA copy number and cell fitness            
            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            
            double s1 = 0;
            double s2 = 0;
            double s3 = 0;

            
            // Run the Gillespie algorithm to decide which cell should divide next 
            // (here it is only 2 possible divisions, cells with or without ecDNA
            // The Gillespie algorithm can be more complicated and can contain cell death 
            // negative selection or copy number dependent selection etc.
            
            a1 = exprand( NumNeutral ); // Random number of neutral cell population
            a2 = fitness * State.size();
            a3 = exprand( a2 );         // Random number of ecDNA cell population
            
            Neutral.at(count1).at(count) = NumNeutral ;
            
            // If a1 < a3 a neutral cell is picked for proliferation. 
            // In contrast if a3 >= a1 a cell with ecDNA is picked for proliferation
            if (a1 <= a3)
            {
                NumNeutral++;
            }
            else
            {
                if (State.size() > 0)
                {  
            	    s1 = mtrand1.randInt(State.size()-1); // Pick a random cell with ecDNA to proliferate
                    double c4 = amplify*State.at(s1);  // double the ecDNA copies in the mother cell
                    // Random binomial trial to distribute the ecDNA copies into daughter cells
                    std::binomial_distribution<> d(c4, 0.5); 
                    s2 = d(gen);
                    
                    // Below is just a few statements to make sure to count all possible cases of daughter cells correctly
                    
		    s3 = State.at(s1);
                }
                else
                {
                    s1 = 0;
                    double c4 = amplify*State.at(s1);
                    std::binomial_distribution<> d(c4 ,0.5);
                    s2 = d(gen);
                    
                    s3 = State.at(s1);
                }
                
                if (s2 == 0)
                {
                    // s2 is the output of the binomial: s2=0 means one daughter cell is neutral.
                    NumNeutral++;
                    State.at(s1) = amplify*s3;
                }
                else
                { 
                    State.at(s1)=s2;
                    if ((amplify*s3) == s2)
                    {
                        // Then the other daughter cell is neutral
                        NumNeutral++;
                    }
                    else
                    {
                    	// Here neither daughter cells are neutral.
                        double c3 = (amplify*s3) - s2;
                        State.push_back(c3);
                    }
                }
            }
            count++;
        }
        for (int i=0; i<State.size(); i++)
        {
            FinalOutput.at(count1).at(i)= State.at(i);
        }
        for (int i=State.size(); i<(NumCells); i++)
        {
            FinalOutput.at(count1).at(i)= 0;
        }
        count1++;
    }


    dataFracs.close();
    
    std::fstream datei ;
    
    // This gives the ecDNA copy number of each type for each cell at the end of the simulation (all measures can be constructed from here)
     
    std::string baseFileName = "Summary_";
    std::string fileName = outputFolder + baseFileName + std::to_string((int)fitness) + "_" + std::to_string(initialcopies) + ".txt";
    datei.open (fileName ,std::ios::out);

    for ( int i=0 ; i<runs ; i++)
    { 
        if (i!=0)
        { 
            datei << "\n";
        }
        for ( int j=0 ; j< NumCells ; j++)
        {
            datei << FinalOutput.at(i).at(j) << " " ;
        }
    }
    
    datei.close() ;
}

void ecDNAEvolveWithLabels(int NumCells, int NumNeutral, int amplify, double fitness, int initialcopies, int runs, std::string outputFolder)
{
   // Initiate a bunch of vectors to store cell states throughout the simulation
   vector <double> State (1,initialcopies); // this initializes a vector of size 1, with value=initialcopies
   // create a variable to store the final number of ecDNA copies of each type in each cell. note that here we assume that all ecDNA are assigned unique labels
   vector <vector <vector <double> > > FinalOutput (runs , vector < vector<double>> (initialcopies, vector <double> (NumCells+1,0))); 
   vector < vector <double> > Neutral (runs ,vector <double> (NumCells,0)); // check if size must be NumCells+1
   
   
   int count1 = 0;  // Dummy variable to count number of simulation repeats
    
    while (count1 < runs)
    {
	int count = 0; // Dummy Variable to count until the number of maximal cells
               
        // Reset some of the vectors of the simulation
        State.resize(1);
        State.at(0) = initialcopies;
        // Label all ecDNA
        vector <vector <double> > stateLabelled = barcodeAll(State);
        int numLabels = stateLabelled.size();
        // compute number of cells with ecDNA
        int numEcdnaCells = numCellsWithEcdna(stateLabelled);
        NumNeutral = 0;
        

        while (count < NumCells)
        {
            cout << count1 << " " << count << "\n"; // output current state of the simulation

            // Define a bunch of dummy variables to store intermediate ecDNA copy number and cell fitness            
            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            
            double s1 = 0;
            vector<double> s2 (numLabels,0);
            vector<double> s3 (numLabels,0);

            
            // Run the Gillespie algorithm to decide which cell should divide next 
            // (here it is only 2 possible divisions, cells with or without ecDNA
            // The Gillespie algorithm can be more complicated and can contain cell death 
            // negative selection or copy number dependent selection etc.
            
            a1 = exprand( NumNeutral ); // Random number of neutral cell population
            a2 = fitness * numEcdnaCells; // fitness times total number of cells with ecDNA
            a3 = exprand( a2 );         // Random number of ecDNA cell population
            
            Neutral.at(count1).at(count) = NumNeutral ;
            // If a1 < a3 a neutral cell is picked for proliferation. 
            // In contrast if a3 >= a1 a cell with ecDNA is picked for proliferation
            if (a1 <= a3)
            {
                NumNeutral++;
            }
            else
            {
                if (numEcdnaCells > 0) // check if there are any cells with ecDNA
                {  
                    // Pick a random cell with ecDNA to proliferate (note! -1)
            	    s1 = mtrand1.randInt(numEcdnaCells - 1); 
                }
                else
                {
                    s1 = 0;
                }

                vector <double> c4 (numLabels, 0);
                for (int k=0; k<numLabels; k++)
                {
                    c4[k] = amplify*stateLabelled[k].at(s1); // double ecDNA copies of each type in the mother cell
                    // Random binomial trials to distribute the ecDNA copies into daughter cells
                    std::binomial_distribution<> d(c4[k], 0.5); // does this give independent distributions each time?
                    s2[k] = d(gen);
                    s3[k] = stateLabelled[k].at(s1); // number of ecDNA of each type in cell s1
                }
   
                // Below is just a few statements to make sure to count all possible cases of daughter cells correctly

                // check if all s2's are zero; this would mean that one daughter cell is neutral
		if (std::all_of(s2.begin(), s2.end(), [](int i) { return i==0; }) == 0)
                {
                    NumNeutral++;
                    for (int k=0; k<numLabels; k++) stateLabelled[k].at(s1) = amplify*s3[k];
                }
                
                else
                {   
		    for (int k=0; k<numLabels; k++) stateLabelled[k].at(s1) = s2[k];
                    // check if the other daughter cell is neutral. note that in the following condition we assume
                    // that the amplification factor is the same for all ecDNA types
                    if (amplify*std::accumulate(s3.begin(), s3.end(), 0) == std::accumulate(s2.begin(), s2.end(), 0))
                    {
                        NumNeutral++;
                    }
                    else // neither daughter cell is neutral
                    {
                        vector <double> c3 (numLabels, 0);
                        for (int k=0; k<numLabels; k++) 
                        {
                            c3[k] = amplify*s3[k] - s2[k];
                            stateLabelled[k].push_back(c3[k]);
                        }
                    }
                }
            }
            count++;
        }
        for (int i=0; i<numEcdnaCells; i++)
        {
            for (int k=0; k<numLabels; k++)
            {
                FinalOutput.at(count1).at(k).at(i)= stateLabelled.at(k).at(i);
            }
        }
        for (int i=numEcdnaCells; i<NumCells; i++)
        {
            for (int k=0; k<numLabels; k++)
            {
                FinalOutput.at(count1).at(k).at(i)= 0;
            }
        }
        count1++;
    }
    
    std::fstream datei ;
    
    // This gives the ecDNA copy number of each type for each cell at the end of the simulation (all measures can be constructed from here)
     
    std::string baseFileName = "Summary_";
    std::string fileName = outputFolder + baseFileName + std::to_string((int)fitness) + "_" + std::to_string(initialcopies) + ".txt";
    datei.open (fileName ,std::ios::out);

    for ( int i=0 ; i<runs ; i++)
    { 
        if (i!=0)
        { 
            datei << "\n";
        }
        for ( int j=0 ; j< NumCells ; j++)
        {   
            for (int k=0; k<initialcopies; k++)
            {
                datei << FinalOutput.at(i).at(k).at(j)<< " " ;
            }
        }
    }
    
    datei.close() ;
}

