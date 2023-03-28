

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

using std::vector;
using std::cout;
using std::cin;


int Samplesize = 1000000;


vector < vector <double> > Sample (Samplesize ,vector <double> (2,0));




time_t seconds ;
MTRand mtrand1(time(NULL));

std::random_device rd;
std::mt19937 gen(rd());



int main()
{
   
    for (int i = 0; i < Samplesize ; i++)
    {
        
        /*std::binomial_distribution<> d(100 ,0.9);*/
        
        double n = mtrand1.randInt(70)+2;
        
        std::binomial_distribution<> d(n ,0.5);
        
        double s = d(gen);
        Sample.at(i).at(0) = s;
        Sample.at(i).at(1) = n-s;}

    
     std::fstream datei ;
   datei.open ("Sample_0.5_COLO.txt" ,std::ios::out);
      for ( int i=0 ; i<Samplesize ; i++)
      { if (i!=0)
      { datei << "\n" ;}
          for ( int j=0 ; j< 2 ; j++)
              datei << Sample.at(i).at(j) << " " ;}
      datei.close() ;
    
}     











