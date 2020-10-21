
//Main header
#include "VCI_headers.h"

int main(int argc, char* argv[])
{
    //Misc. initialization
    StartTime = (unsigned)time(0); //Time the program starts
    cout.precision(12);
    cout << fixed;
    //End of section
    //Initialize local variables
    fstream vcidata,freqfile; //File streams for the input and output files
    string cpname;
    MatrixXd CIVec; //CI eigenvectors (wavefunction)
    VectorXd CIFreq; //CI eigenvalues (frequencies)
    CIVec = MatrixXd(1,1); CIVec.setZero();
    CIFreq = VectorXd(1); CIFreq.setZero();
    double Ezpe = 0; //CI zero-point energy
    double RunTime; //Run time for the calculations
    string TimeUnits; //Seconds, minutes, or hours for RunTime
    //End of section
    
    
    //Print title and compile date
    PrintFancyTitle();
    cout << "Last modification: ";
    cout << __TIME__ << " on ";
    cout << __DATE__ << '\n';
    cout << '\n';
    cout.flush();
    //End of section
    
    //Gather input and check for errors
    cout << "Reading input..." << '\n';
    ReadCIArgs(argc,argv,vcidata,freqfile,cpname); //Read arguments
    ReadCIInput(vcidata); //Read input files
    //End of section
    if(restart){ // Load checkpoint file if it exists
        cout << '\n' << "Loading checkpoint file...";
        cout.flush();
        fstream checkfile;
        checkfile.open(cpname,ios_base::in);
        ReadCheckpoint(CIFreq,CIVec,checkfile);
        checkfile.close();
        cout << "done." << '\n';
        cout << "  Basis functions: " << BasisSet.size() << '\n' << endl;
    }
    if(HCI_Eps == 0.){ // One-shot regular VCI procedure
        if(!restart){
            SpMat HSparse(BasisSet.size(),BasisSet.size());
            cout << "Constructing the Vibrational CI Hamiltonian...";
            cout.flush();
            MakeHamSparse(HSparse);
            if(NEig==1){
                cout << "done." << "\n" << "\n" << "Solving for the ground state..." << "\n" << endl;
            }else{
                cout << "done." << "\n" << "\n" << "Solving for the lowest " 
                    << NEig << " eigenstates..." << "\n" << endl;
            }
            SparseDiagonalize(HSparse,CIVec,CIFreq); //Diagonalization wrapper
            if(BasisSet.size()>100000){ 
                cout << "Writing checkpoint file...";
                cout.flush();
                fstream checkfile;
                checkfile.open(cpname,ios_base::out | ios_base::trunc);
                if (!checkfile.good() || cpname=="none"){
                    cout << '\n' << "Error: checkpoint file could not be found or no name was given." << '\n' << endl;
                }else{
                    WriteCheckpoint(CIFreq,CIVec,checkfile);
                    cout << "done." << '\n' << endl;
                }
                checkfile.close();
            }
        }
    }
    else{ // Iterative HCI procedure
        cout << "Beginning Heat Bath CI procedure..." << '\n';
        Perform_HCI(CIVec,CIFreq,cpname);
    }
    if(perturb){
        EndTime = (unsigned)time(0); //Time the calculation stops
        RunTime = (double)(EndTime-StartTime); //Total run time
        if (RunTime >= 3600)
        {
            //Switch to hours
            RunTime /= 3600;
            TimeUnits = "hours";
        }
        else if (RunTime >= 60)
        {
            //Switch to minutes
            RunTime /= 60;
            TimeUnits = "minutes";
        }
        else
        {
            //Stick with seconds
            TimeUnits = "seconds";
        }
        cout.precision(2); //Truncate time
        cout << " Variational stage took " << RunTime;
        cout << " " << TimeUnits;
        cout << '\n' << '\n';
        
        cout << "Performing Epstein-Nesbet perturbation theory..." << endl;
        if(PT2_Eps>0){    
            cout << " Using heat bath sorting of perturbative states with cutoff: " << scientific << PT2_Eps << fixed << endl;
        }
        else{
            cout << " Including all connected states in PT2 calculation." << endl;
        }
        DoPT2(CIVec,CIFreq);
    }
    
        //Print freqs (do before removing ZPE)
    cout.precision(10); //Replace settings
//    sort(CIFreq.data(),CIFreq.data()+CIFreq.size());//Sort eigenvalues in ascending order
    PrintFreqs(CIFreq,CIVec,freqfile); //Print eigenvalues to file
    
    //Calculate and remove ZPE
    Ezpe = CIFreq.minCoeff(); //Find ZPE
    if (Ezpe < 0)
    {
        //Print error
        cout << "Error: The zero-point energy is negative.";
        cout << " Something is wrong.";
        cout << '\n';
        cout << "Try reducing the number of quanta or removing";
        cout << " large anharmonic";
        cout << '\n';
        cout << "force constants.";
        cout << '\n' << '\n';
        cout.flush();
        //Quit
        exit(0);
    }
    
    //Print results
    cout << "Printing the results...";
    cout << '\n' << '\n';
    cout << "Results:" << '\n';
    cout << "  Zero-point energy: ";
    cout.precision(2); //Truncate energy
    cout << Ezpe << " (1/cm)"; //Print energy
    cout.precision(12); //Replace settings
    cout << '\n';
    cout.flush(); //Print progress
    EndTime = (unsigned)time(0); //Time the calculation stops
    RunTime = (double)(EndTime-StartTime); //Total run time
    if (RunTime >= 3600)
    {
        //Switch to hours
        RunTime /= 3600;
        TimeUnits = "hours";
    }
    else if (RunTime >= 60)
    {
        //Switch to minutes
        RunTime /= 60;
        TimeUnits = "minutes";
    }
    else
    {
        //Stick with seconds
        TimeUnits = "seconds";
    }
    cout.precision(2); //Truncate time
    cout << "  Run time: " << RunTime;
    cout.precision(12); //Replace settings
    cout << " " << TimeUnits;
    cout << '\n' << '\n';
    cout << "Done.";
    cout << '\n' << '\n';
    cout.flush(); //Print progress
    //End of section
    
    //Empty memory and delete junk files
    vcidata.close();
    freqfile.close();
    //End of section
    
    //Quit
    return 0;
};
