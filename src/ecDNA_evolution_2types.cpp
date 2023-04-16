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
#include <math.h>
#include <algorithm>  

using std::vector;
using std::cout;
using std::cin;

int NumCells = 10;//100000;    // Maximum number of cells
int NumNeutral = 0;      // Initial number of cells with no ecDNA
int amplify = 2;         // factor of ecDNA amplification upon cell division X -> amplify * X
double fitness = 1;//3;      // relative fitness of cells with ecDNA ( fitness =1 corresponds to neutral dynamics)
int initialcopies = 10;//45; // Initial copies of ecDNA in the first founder cell
int initialcopies_a = 3; // Initial copies of type 'a' ecDNA in founder cell (set to nonzero)
int initialcopies_b = initialcopies - initialcopies_a; // Initial copies of type 'b' ecDNA in founder cell
int runs = 1;//20;         // Number of simulation repeats
// Initiate a bunch of vectors to store cell states throughout the simulation
vector <double> State_a (1,initialcopies_a); // this initializes a vector of size 1, with value=initialcopies for type 'a'
vector <double> State_b (1,initialcopies_b); // this initializes a vector of size 1, with value=initialcopies for type 'b'
vector < vector <double> > FinalOutput_a (runs ,vector <double> (NumCells+1,0)); // final number of ecDNA of type 'a' for each cell
vector < vector <double> > FinalOutput_b (runs ,vector <double> (NumCells+1,0));
vector < vector <double> > Neutral (runs ,vector <double> (NumCells,0));
// vector <double> Rate (NumCells ,0); // vector of zeros of size NumCells
// vector <double> Rate1 (NumCells ,0);

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

int main()
{
    int count1 = 0;  // Dummy variable to count number of simulation repeats
    
    while (count1 < runs)
    {
	int count = 0; // Dummy Variable to count until the number of maximal cells
               
        // Reset some of the vectors of the simulation
        State_a.resize(1);
        State_a.at(0) = initialcopies_a;
        State_b.resize(1);
        State_b.at(0) = initialcopies_b;
        NumNeutral = 0;

        while (count < NumCells)
        {
            cout << count1 << " " << count << "\n"; // output current state of the simulation
            
            // Define a bunch of dummy variables to store intermediate ecDNA copy number and cell fitness            
            double a1 = 0;
            double a2 = 0;
            double a3 = 0;
            
            double s1 = 0;
            double s2_a = 0;
            double s2_b = 0;
            double s3_a = 0;
            double s3_b = 0;
            
            // Run the Gillespie algorithm to decide which cell should divide next 
            // (here it is only 2 possible divisions, cells with or without ecDNA
            // The Gillespie algorithm can be more complicated and can contain cell death 
            // negative selection or copy number dependent selection etc.
            
            a1 = exprand( NumNeutral ); // Random number of neutral cell population
            a2 = fitness * numCellsWithEcdna(State_a.size(), State_b.size()); // fitness times total number of ecDNA
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
                if (numCellsWithEcdna(State_a.size(), State_b.size()) > 0)
                {  
                    // Pick a random cell with ecDNA to proliferate (note! -1)
            	    s1 = mtrand1.randInt(numCellsWithEcdna(State_a.size(), State_b.size()) - 1); 
                }
                else
                {
                    s1 = 0;
                }
                double c4_a = amplify*State_a.at(s1);  // double the ecDNA 'a' copies in the mother cell
                double c4_b = amplify*State_b.at(s1);  // double the ecDNA 'b' copies in the mother cell
                
                // Random binomial trials to distribute the ecDNA copies into daughter cells
                std::binomial_distribution<> d_a(c4_a, 0.5); 
                s2_a = d_a(gen);
                
                std::binomial_distribution<> d_b(c4_b, 0.5); 
                s2_b = d_b(gen);
                    
                // Below is just a few statements to make sure to count all possible cases of daughter cells correctly
                // if (s2 == 2)
                // {
                //     Rate1.at(count)++;
                // }
                
                // number of ecDNA of types a and b in cell s1
		s3_a = State_a.at(s1);
                s3_b = State_b.at(s1);
                
                if (s2_a == 0 && s2_b == 0)
                {
                    // s2 is the output of the binomial: s2=0 means one daughter cell is neutral.
                    NumNeutral++;
                    State_a.at(s1) = amplify*s3_a;
                    State_b.at(s1) = amplify*s3_b;
                    // Rate.at(count)++;
                }
                
                else
                { 
                    State_a.at(s1)=s2_a;
                    State_b.at(s1)=s2_b;
                    if ((amplify*s3_a + amplify*s3_b) == s2_a + s2_b)
                    {
                        // Then the other daughter cell is neutral
                        NumNeutral++;
                        // Rate.at(count)++ ;
                    }
                    else
                    {
                    	// Here neither daughter cells are neutral
                        double c3_a = (amplify*s3_a) - s2_a;
                        double c3_b = (amplify*s3_b) - s2_b;
                        State_a.push_back(c3_a);
                        State_b.push_back(c3_b);
                    }
                }
            }
            count++;
        }
        for (int i=0; i<numCellsWithEcdna(State_a.size(), State_b.size()); i++)
        {
            FinalOutput_a.at(count1).at(i)= State_a.at(i);
            FinalOutput_b.at(count1).at(i)= State_b.at(i);
        }
        for (int i=numCellsWithEcdna(State_a.size(), State_b.size()); i<(NumCells); i++)
        {
            FinalOutput_a.at(count1).at(i)= 0;
            FinalOutput_b.at(count1).at(i)= 0;
        }
        count1++;
    }

    std::fstream datei_a ;
    std::fstream datei_b ;
    /*datei.open ("NonNeutral().txt" ,std::ios::out);
    for ( int i=0 ; i < State.size() ; i++)
    {		datei << State.at(i) << " " ;}
    datei.close() ;*/
    
    /*datei.open ("Rate.txt" ,std::ios::out);
    for ( int i=0 ; i < NumCells ; i++)
    {		datei << Rate.at(i) << "\n " ;}
    datei.close() ;*/
    
    /* datei.open ("Rate1.txt" ,std::ios::out);
    for ( int i=0 ; i < NumCells ; i++)
    {		datei << Rate1.at(i) << "\n " ;}
    datei.close() ;*/

    /*datei.open ("Neutral2.txt" ,std::ios::out);
    for ( int i=0 ; i<runs ; i++)
    { if (i!=0)
    { datei << "\n" ;}
        for ( int j=0 ; j< NumCells ; j++)
            datei << Neutral.at(i).at(j) << " " ;}
    datei.close() ;*/
    
    // This gives the ecDNA copy number of each type for each cell at the end of the simulation (all measures can be constructed from here)
    datei_a.open ("NonNeutralSummary_a.txt" ,std::ios::out);
    datei_b.open ("NonNeutralSummary_b.txt" ,std::ios::out);
    for ( int i=0 ; i<runs ; i++)
    { 
        if (i!=0)
        { 
            datei_a << "\n";
            datei_b << "\n";
        }
        for ( int j=0 ; j< NumCells ; j++)
        {
            datei_a << FinalOutput_a.at(i).at(j) << " " ;
            datei_b << FinalOutput_b.at(i).at(j) << " " ;
        }
    }
    
    datei_a.close() ;
    datei_b.close() ;
}     

// Auxiliary functions
double exprand( double lambda)
{
    double y = mtrand1.randExc();
    return (-log(1.-y)/lambda);
}

void print(vector <double> const &a) 
{
   cout << "The vector elements are : ";

   for(int i=0; i < a.size(); i++)
   {cout << a.at(i) << ' ';}
   cout << std::endl;
}

int numCellsWithEcdna(int numCellsWithEcdnaA, int numCellsWithEcdnaB)
{
    if (numCellsWithEcdnaA != numCellsWithEcdnaB) return -1;
    return numCellsWithEcdnaA;
}
