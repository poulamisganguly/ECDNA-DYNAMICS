

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
# include <math.h>
# include <algorithm>  

using namespace std;

class Ecdna {
  public:
    int label; // label before or after barcoding
    double fitness; // fitness of each ecDNA type 
    
    // constructor for class ecDNA
    Ecdna (int a, double b) : label{a}, fitness{b} {}
};

class Cell {
  public:
    bool type; // with or without ecDNA
    double fitness, amp; // s parameter and amplification factor before division
    int copyNum; //number of ecDNA copies 
    std::vector<Ecdna> ecdna; // vector of ecDNA 
    
    // constructor for class Cell
    Cell (bool a, double b, int n, double m) : type{a}, fitness{b}, copyNum{n}, amp{m} {}
    
    // assign ecDNA with label=0 and fitness=1
    void initialize_cell () {
        ecdna.assign (copyNum, Ecdna(0,1.0));
        }
    
    // output all ecDNA in cell in format (label, fitness)
    void output_cell_state () {
        for (Ecdna i: ecdna) {std::cout << i.label << ", " << i.fitness << "\n";}
        
    }
        
    // amplify ecDNA in the cell and do random segregation
    void divide () {
        
        
        }
    
};

class Population
{
int number; // number of cells

};


int main () {


return 0;
}


