/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Headers, libraries, and data structures for VHCI

*/

//Make including safe
#ifndef VCI_HEADERS
#define VCI_HEADERS

//Header Files
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <string>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <sys/stat.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <Eigen/StdList>
#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <boost/functional/hash.hpp>
//Set namespaces for common libraries
using namespace Eigen;
using namespace std;
using namespace Spectra;
using namespace boost;


//Global exact constants
const double pi = 4*atan(1); //Pi
const double rt2pi = sqrt(2*pi); //Needed for Gaussian broadening

//Global measured constants (NIST, CODATA 2010)
const double cs = 2.99792458e-10; //Speed of light (cm)
const double k = 0.69503476; //Boltzmann constant (cm^-1)

//Global derived constants
const double h = 2*pi; //Planck constant (cm^-1)

//Timers
int StartTime = 0; //Time the calculation starts
int EndTime = 0; //Time the calculation ends

//Custom data structures
struct HOFunc
{
    //Data structure for harmonic oscillator basis functions
    double Freq; //Frequency
    int Quanta; //Number of quanta in the mode
    // condition for equality
    bool operator == (const HOFunc &other) const {
        return (Freq == other.Freq) && (Quanta == other.Quanta);
    }
    // condition for inequality
    bool operator != (const HOFunc &other) const {
        return (Freq != other.Freq) || (Quanta != other.Quanta);
    }
};

struct WaveFunction
{
    //Data structure for storing a VCI wavefunction
    int M; //Number of modes
    vector<HOFunc> Modes; // WaveFunctions
    bool operator == (const WaveFunction &other) const {
        bool same_wfn = 1;
        if(other.M != M){
            return 0;
        }
        for(int n=0; n<other.M; n++){
            if(other.Modes[n] != Modes[n]){
                return 0;
            }
        }
        return same_wfn;
    }
    // condition for inequality
    bool operator != (const WaveFunction &other) const {
        bool diff_wfn = 0;
        if(other.M != M){
            return 1;
        }
        for(int n=0; n<other.M; n++){
            if(other.Modes[n] != Modes[n]){
                return 1;
            }
        }
        return diff_wfn;
    }

};

struct FConst
{
    //Data structure for anharmonic force constants
    //Note: fc should include the permutation term
    double fc; //Value of the force constant
    vector<int> fcpow; //Modes and powers for the force constant
    vector<int> ShortModes; // Shortened list of modes only including those affected
    vector<int> ModePowers; // Power for each affected mode
};

struct WfnHasher
{
    size_t operator () (const WaveFunction& key) const 
    {
        // function to generate unique hash for a WaveFunction using Boost
        size_t seed = 0;
        for(int n=0; n<key.M; n++){
            hash_combine(seed, hash_value(key.Modes[n].Quanta));
            //hash_combine(seed, key.Modes[n].Freq);
        }
        return seed;
    }
};


//Global variables
int Ncpus = 0; //Number of CPUs for the calculations
int NEig; //Number of eigenstates to include in HB optimization
vector<WaveFunction> BasisSet; //Full basis set
vector<FConst> AnharmFC; //List of force constants
vector<FConst> AnharmHB; //List of force constants for HB procedure with singles and doubles terms
vector<FConst> CubicFC; //List of cubic force constants
vector<FConst> QuarticFC; //List of quartic force constants
vector<FConst> QuinticFC; //List of quintic force constants
vector<FConst> SexticFC; //List of sextic force constants

double HCI_Eps; // Variational CI energy parameter epsilon_1
bool perturb; // Should I use PT2 correction?
bool restart; // Am I restarting from a checkpoint file?
double PT2_Eps; // PT2 basis energy parameter epsilon_2

typedef unordered_set<WaveFunction, WfnHasher> HashedStates;
typedef SparseMatrix<double, 0, ptrdiff_t> SpMat;
typedef Triplet<double, ptrdiff_t> Trip;
//Function declarations
void AnharmHam(MatrixXd&);

double AnharmPot(int,int,const FConst&);

void AnnihilationLO(double&,int&);

bool CheckFile(const string&);

void IsRestart(const string&);

void WriteCheckpoint(VectorXd&,MatrixXd&,fstream&);

void ReadCheckpoint(VectorXd&,MatrixXd&,fstream&);

void CreationLO(double&,int&);

double Fact(int);

int FindMaxThreads();

void PrintFancyTitle();

void PrintFreqs(VectorXd&,MatrixXd&,fstream&);

void ReadCIArgs(int,char*,fstream&,fstream&,string&);

void ReadCIInput(fstream&);

void ScaleFC();

bool ScreenState(int,int,const vector<int>&,const FConst&);

void DenseDiagonalize(MatrixXd&,MatrixXd&,VectorXd&);

void ArnoldiDiagonalize(MatrixXd&,MatrixXd&,VectorXd&);

void SparseDiagonalize(SpMat&,MatrixXd&,VectorXd&);

void ZerothHam(MatrixXd&);

void ZerothHamSparse(vector<Trip>&);

void AnharmHamSparse(vector<Trip>&);

void MakeHamSparse(SpMat&);

bool sortByFC(const FConst&, const FConst&);

void HeatBath_Sort_FC();

void AddStatesHB(HashedStates&,HashedStates&,int,double,double);

void Perform_HCI(MatrixXd&, VectorXd&, const string&);

void QDiffVec(int, int, int&, int&, vector<int>&);

void DoPT2(MatrixXd&, VectorXd&);

void AddASCI(vector<WaveFunction>&, MatrixXd&, VectorXd&);

void AddConnected(HashedStates&, HashedStates& NewStates, int);

void Perform_ASCI(MatrixXd&, VectorXd&, const string&);
//Function definitions 
#include "Core_functions.cpp"
#include "Input_Reader.cpp"
#include "Ham.cpp"
#include "HB.cpp"
#include "ASCI.cpp"
#include "PT2.cpp"

#endif
