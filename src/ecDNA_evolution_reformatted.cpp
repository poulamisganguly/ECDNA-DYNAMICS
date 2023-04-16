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
int initialcopies_a = 3; // Initial copies of type a ecDNA 
int runs = 1;//20;         // Number of simulation repeats
// Initiate a bunch of vectors to store cell states throughout the simulation
vector <int> Labels; // this is a label vector containing the labels of all ecDNA in the cell
vector <double> State (1,initialcopies); // this initializes a vector of size 1, with value=initialcopies
vector < vector <double> > FinalOutput (runs ,vector <double> (NumCells+1,0));
vector < vector <double> > Neutral (runs ,vector <double> (NumCells,0));
vector <double> Rate (NumCells ,0);
vector <double> Rate1 (NumCells ,0);

// Initiate random number generator
time_t seconds ;
MTRand mtrand1(time(NULL));
std::random_device rd;
std::mt19937 gen(rd());

// Define a function to generate a random number for the Gillespie Algorithm (Defined at the bottom)
double exprand( double lambda );

// Define a function to print vectors
void print(vector <double> const &a);


int main()
{
    int count1 = 0;  // Dummy variable to count number of simulation repeats
    
    // generate as many labels as initial ecDNA copies
    // for (int n=0; n<initialcopies; n++)
    // {
    //	Labels.push_back(n);
    // }
    //cout << Labels.size() << std::endl;
    
    // State.insert( State.end(), Labels.begin(), Labels.end() ); // concatenate Labels with State vector
    //cout << "State size = " << State.size() << std::endl;
    //print(State);
    
    while (count1 < runs)
    {
	int count = 0; // Dummy Variable to count until the number of maximal cells
               
        // Reset some of the vectors of the simulation
        State.resize(1);
        State.at(0) = initialcopies;
        NumNeutral = 0;

        while (count < NumCells)
        {
            cout << count1 << " " << count << "\n"; // output current state of the simulation
            
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
                    if (s2 == 2)
                    {
                        Rate1.at(count)++;
                    }
		    s3 = State.at(s1);
                }
                else
                {
                    s1 = 0;
                    double c4 = amplify*State.at(s1);
                    std::binomial_distribution<> d(c4 ,0.5);
                    s2 = d(gen);
                    if (s2 == 2)
                    {
                        Rate1.at(count)++;
                    }
                    s3 = State.at(s1);
                }
                
                if (s2 == 0)
                {
                    // s2 is the output of the binomial: s2=0 means one daughter cell is neutral.
                    NumNeutral++;
                    State.at(s1) = amplify*s3;
                    Rate.at(count)++;
                }
                else
                { 
                    State.at(s1)=s2;
                    if ((amplify*s3) == s2)
                    {
                        // Then the other daughter cell is neutral
                        NumNeutral++;
                        Rate.at(count)++ ;
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

    std::fstream datei ;
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
    
    // This gives the ecDNA copy number for each cell at the end of the simulation (all measures can be constructed from here)
    datei.open ("NonNeutralSummaryPC3.txt" ,std::ios::out);
    for ( int i=0 ; i<runs ; i++)
    { 
        if (i!=0)
        { 
            datei << "\n";
        }
        for ( int j=0 ; j< NumCells ; j++)
            datei << FinalOutput.at(i).at(j) << " " ;
    }
    datei.close() ;
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

