/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Basic routines for VHCI

*/

//Core utility functions
void PrintFancyTitle()
{
    cout << '\n';
    cout << "######################################";
    cout << "#######################################";
    cout << '\n';
    cout << "#                                       ";
    cout << "                                    #";
    cout << '\n';
    cout << "#              ";
    cout << "Vibrational Heat-Bath Configuration Interaction";
    cout << "              #";
    cout << '\n';
    cout << "#                         ";
    cout << "By: Jonathan H. Fetherolf";
    cout << "                         #";
    cout << '\n';
    cout << "#                                       ";
    cout << "                                    #";
    cout << '\n';
    cout << "#                      ";
    cout << "Based on LOVCI by Eric G. Kratz";
    cout << "                      #";
    cout << '\n';
    cout << "#                                       ";
    cout << "                                    #";
    cout << '\n';
    cout << "######################################";
    cout << "#######################################";
    cout << '\n' << '\n';
    cout.flush();
    return;
};

bool CheckFile(const string& file)
{
    //Checks if a file exists
    struct stat buffer;
    if (stat(file.c_str(),&buffer) != -1)
    {
        return 1;
    }
    return 0;
};

int FindMaxThreads()
{
    //Function to count the number of allowed threads
    int ct = 0; //Generic counter
    #pragma omp parallel reduction(+:ct)
    ct += 1; //Add one for each thread
    #pragma omp barrier
    //Return total count
    return ct;
};

double Fact(int n)
{
    //Calculate a factorial
    double val = 1;
    while (n > 0)
    {
        val *= n;
        n -= 1;
    }
    return val;
};

void IsRestart(const string& cpname){
    fstream checkfile;
    checkfile.open(cpname,ios_base::in);
    if ( checkfile.peek() == std::ifstream::traits_type::eof() ){
        restart=0; 
    }else{restart=1;}
    checkfile.close();
}

void WriteCheckpoint( VectorXd& Freqs, MatrixXd& Vecs, fstream& cpfile ){
    // Function for creating checkpoint file that stores basis and eigenvectors
    // Save basic dimensions of system
    cpfile << BasisSet.size() << " " << BasisSet[0].M << " " << NEig << '\n';
    // Save fundamental freqs and eigenvalues
    for(unsigned int j=0; j<BasisSet[0].M; j++){
        cpfile << fixed << setprecision(17) << BasisSet[0].Modes[j].Freq << " ";
    }
    for(unsigned int j=0; j<Freqs.size(); j++){
        cpfile << fixed << setprecision(17) << Freqs(j) << " ";
    }
    cpfile << '\n';
    for(unsigned int i=0; i<BasisSet.size(); i++){
        //Save basis states
        for(unsigned int j=0; j<BasisSet[i].M; j++){
            cpfile << BasisSet[i].Modes[j].Quanta << " ";
        }
        for(unsigned int j=0; j<Vecs.cols(); j++){
            cpfile << fixed << setprecision(17) << Vecs(i,j) << " ";
        }
        cpfile << '\n';
    }
    cpfile.flush();
}

void ReadCheckpoint( VectorXd& Freqs, MatrixXd& Vecs, fstream& cpfile ){
    // Function for reading checkpoint file upon restart
    // Read fundamental freqs and eigenvalues
    int NBasis;
    cpfile >> NBasis;
    int param;
    cpfile >> param;
    if(param!=BasisSet[0].M){
        cout << "Error: Checkpoint file contains wrong number of fundamental modes!" << endl;
        exit(0);
    }
    cpfile >> param;
    if(param!=NEig){
        cout << "Error: Checkpoint file contains wrong number of eigenstates!" << endl;
        exit(0);
    }
    double freq;
    for(unsigned int i=0; i<BasisSet[0].M; i++){
        cpfile >> freq;
        if(freq!=BasisSet[0].Modes[i].Freq){
            cout << "Error: Checkpoint file contains wrong fundamental frequencies!" << endl;
            exit(0);
        }
    }
    Freqs = VectorXd(NEig);
    for(unsigned int i=0; i<NEig; i++){
        cpfile >> freq;
        Freqs(i) = freq;

    }
    vector<WaveFunction> CPBasis(NBasis);
    WaveFunction wftemp = BasisSet[0];
    int qtemp;
    double citemp;
    Vecs = MatrixXd(NBasis,NEig);
    for(unsigned int i=0; i<NBasis; i++){
        //Reads the number of quanta and CI coefficents for each basis state
        for(unsigned int j=0; j<BasisSet[0].M; j++){
            cpfile >> qtemp;
            wftemp.Modes[j].Quanta = qtemp;
        }
        CPBasis[i] = wftemp;
        for(unsigned int j=0; j<NEig; j++){
            cpfile >> citemp;
            Vecs(i,j) = citemp;
        }
    }
    BasisSet = CPBasis;
}
