/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Functions for building and diagonalizing the Hamiltonian matrix

*/

inline void CreationLO(double& ci, int& ni)
{
    //Creation ladder operator
    ci *= sqrt(ni+1); //Update coefficient
    ni += 1; //Update state
    return;
};

inline void AnnihilationLO(double& ci, int& ni)
{
    //Annihilation ladder operator
    ci *= sqrt(ni); //Update coefficient
    ni -= 1; //Update state
    //Check for impossible states
    if (ni < 0)
    {
        ni = -1; //Used later to remove the state
        ci = 0; //Delete state
    }
    return;
};

double AnharmPot(int n, int m, const FConst& fc)
{
    //Calculate anharmonic matrix elements for <m|H|n>
    double Vnm = 0;
    // Initialize new states
    vector<vector<int> > NewStates(fc.ShortModes.size());
    vector<vector<double> > StateCoeffs(fc.ShortModes.size());
    //Create new states
    for (unsigned int i=0;i<fc.ShortModes.size();i++)
    {
        //Apply operators
        NewStates[i].push_back(BasisSet[n].Modes[fc.ShortModes[i]].Quanta);
        StateCoeffs[i].push_back(1.0);
        for (int j=0;j<fc.ModePowers[i];j++)
        {
            //Create new state for mode i
            vector<int> stateupdate;
            vector<double> coeffupdate;
            for (unsigned int k=0;k<NewStates[i].size();k++)
            {
                int quant;
                double coeff;
                //Creation
                quant = NewStates[i][k];
                coeff = StateCoeffs[i][k];
                CreationLO(coeff,quant);
                stateupdate.push_back(quant);
                coeffupdate.push_back(coeff);
                //Annihilation
                quant = NewStates[i][k];
                coeff = StateCoeffs[i][k];
                AnnihilationLO(coeff,quant);
                if (quant >= 0)
                {
                    stateupdate.push_back(quant);
                    coeffupdate.push_back(coeff);
                }
            }
            //Save states
            NewStates[i] = stateupdate;
            StateCoeffs[i] = coeffupdate; // Accounting for permutations
        }
    }
    //Sum energies
    vector<double> CoeffSum;
    for (unsigned int i=0;i<fc.ShortModes.size();i++)
    {
        CoeffSum.push_back(0.0);
        for (unsigned int j=0;j<NewStates[i].size();j++)
        {
            int quantn = NewStates[i][j];
            int quantm = BasisSet[m].Modes[fc.ShortModes[i]].Quanta;
            if (quantn == quantm)
            {
                CoeffSum[i] += StateCoeffs[i][j];
            }
        }
    }
    //Scale by the force constant
    Vnm = fc.fc;
    //Combine coeffcients
    for (unsigned int i=0;i<CoeffSum.size();i++)
    {
        Vnm *= CoeffSum[i];
    }
    return Vnm;
};

//Hamiltonian operators
void ZerothHam(MatrixXd& H)
{
    //Calculate the harmonic Hamiltonian matrix elements
    //Harmonic matrix elements
    #pragma omp parallel for
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        //Loop over all modes
        double Ei = 0.; //Hii matrix element
        for (int j=0;j<BasisSet[i].M;j++)
        {
            //Calculate partial energies
            double Ej = 0.5;
            Ej += BasisSet[i].Modes[j].Quanta;
            Ej *= BasisSet[i].Modes[j].Freq;
            //Update matrix element
            Ei += Ej;
        }
        //Update Hamiltonian
        H(i,i) += Ei;
    }
    return;
};

void AnharmHam(MatrixXd& H)
{
    //Add anharmonic terms to the Hamiltonian
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }
    #pragma omp parallel for
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        vector<int> qdiffvec(BasisSet[0].M,0);
        for (unsigned int j=i;j<BasisSet.size();j++) // starts with j=i to exploit Hermiticity (Hij = Hji)
        {
            double Vij = 0;
            int qdiff = 0; // total number of quanta difference between states
            int mchange = 0; // number of modes with nonzero change in quanta
            QDiffVec(i,j,qdiff,mchange,qdiffvec);
            if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0){ 
                // States cannot differ by more than fcmax quanta 
                for (unsigned int k=0;k<QuarticFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(i,j,QuarticFC[k]);
                    }
                }
                for (unsigned int k=0;k<SexticFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,SexticFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(i,j,SexticFC[k]);
                    }
                }
                H(i,j) += Vij;
                if(i!=j){
                    H(j,i) += Vij; // Exploiting Hermiticity (Hij=Hji)
                }
            }
            if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1){ 
                // fcmax-1 assumes max order is even (4th or 6th)
                // States cannot differ by more than fcmax quanta 
                for (unsigned int k=0;k<CubicFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,CubicFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(i,j,CubicFC[k]);
                    }
                }
                for (unsigned int k=0;k<QuinticFC.size();k++)
                {
                    if (ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k]))
                    {// Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        Vij += AnharmPot(i,j,QuinticFC[k]);
                    }
                }
                H(i,j) += Vij;
                if(i!=j){
                    H(j,i) += Vij; // Exploiting Hermiticity (Hij=Hji)
                }
            }
        }
    }
    return;
};

//Hamiltonian operators
void ZerothHamSparse(vector<Trip>& HTrip)
{
    //Calculate the harmonic Hamiltonian matrix elements in sparse form of Eigen::Triplet
    //Harmonic matrix elements
    #pragma omp parallel for 
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        //Loop over all modes
        double Ei = 0.; //Hii matrix element
        for (int j=0;j<BasisSet[i].M;j++)
        {
            //Calculate partial energies
            double Ej = 0.5;
            Ej += BasisSet[i].Modes[j].Quanta;
            Ej *= BasisSet[i].Modes[j].Freq;
            //Update matrix element
            Ei += Ej;
        }
        //Update Hamiltonian
        #pragma omp critical
        HTrip.push_back(Trip(i,i,Ei/2.)); // Dividing by 2 so I can do H=H+H*
    }
    return;
};

void AnharmHamSparse(vector<Trip>& HTrip)
{
    //Add anharmonic terms to the Sparse Hamiltonian in sparse form of Eigen::Triplet
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }
    #pragma omp parallel for
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
        vector<int> qdiffvec(BasisSet[0].M,0);
        for (unsigned int j=i;j<BasisSet.size();j++) // starts with j=i to exploit Hermiticity (Hij = Hji)
        {
            double Vij = 0;
            int mchange = 0; // number of modes with nonzero change in quanta
            int qdiff = 0; // total number of quanta difference between states
            QDiffVec(i,j,qdiff,mchange,qdiffvec);
            if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0){ 
                // States cannot differ by more than fcmax quanta
                double W = 0.;
                double V = 0.;
                for (unsigned int k=0;k<QuarticFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        double val = AnharmPot(i,j,QuarticFC[k]);
                        Vij += val; 
                        if(abs(val) > abs(W)){
                            W = val;
                            V = QuarticFC[k].fc;
                        }
                    }
                }
                for (unsigned int k=0;k<SexticFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,SexticFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        double val = AnharmPot(i,j,SexticFC[k]); 
                        Vij += val;
                        if(abs(val) > abs(W)){
                            W = val;
                            V = SexticFC[k].fc;
                        }
                    }
                }
                #pragma omp critical
                if(i==j){
                    HTrip.push_back(Trip(i,j,Vij/2.)); // Dividing by 2 so I can do H=H+H*

                }else{
                    HTrip.push_back(Trip(i,j,Vij)); // Exploiting Hermiticity (Hij=Hji)
                }
            }
            if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1){
                // fcmax-1 assumes max order is even (4th or 6th)
                // States cannot differ by more than fcmax quanta 
                double W = 0.;
                double V = 0.;
                for (unsigned int k=0;k<CubicFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,CubicFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        
                        double val = AnharmPot(i,j,CubicFC[k]);
                        Vij += val; 
                        if(abs(val) > abs(W)){
                            W = val;
                            V = CubicFC[k].fc;
                        }
                        //Vij += AnharmPot(i,j,CubicFC[k]);
                    }
                }
                for (unsigned int k=0;k<QuinticFC.size();k++){
                    if (ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k])){
                        // Screen force constants that cannot connect basis states i and j
                        //Add anharmonic matrix elements
                        double val = AnharmPot(i,j,QuinticFC[k]);
                        Vij += val; 
                        if(abs(val) > abs(W)){
                            W = val;
                            V = QuinticFC[k].fc;
                        }
                        //Vij += AnharmPot(i,j,QuinticFC[k]);
                    }
                }
                #pragma omp critical
                if(i==j){
                    HTrip.push_back(Trip(i,j,Vij/2.)); // Dividing by 2 so I can do H=H+H*

                }else{
                    HTrip.push_back(Trip(i,j,Vij)); // Exploiting Hermiticity (Hij=Hji)
                }
            }
        }
    }
    return;
};

inline void MakeHamSparse(SpMat& HSp){
    //Build the sparse CI Hamiltonian
    vector< Trip > HTrip;
    ZerothHamSparse(HTrip);
    AnharmHamSparse(HTrip);
    HSp.setFromTriplets(HTrip.begin(),HTrip.end());
    HSp.makeCompressed();
    HTrip = vector< Trip >(); // Free memory
    SpMat HSpT = HSp.transpose();
    HSpT.makeCompressed();
    HSp += HSpT; // Complete symmetric matrix
    HSpT = SpMat(1,1); // Free memory
    HSp.makeCompressed();
    cout << "The Hamiltonian is " << fixed << setprecision(2) << 
        100*(1.-(double)HSp.nonZeros()/(double)HSp.size()) << "% sparse." << endl;
    return; 
}

//Utility functions
inline void SparseDiagonalize(SpMat& H, MatrixXd& Psi, VectorXd& E)
{
    typedef SparseSymMatProd<double,Eigen::Lower,0,ptrdiff_t> SparseMVProd;
    SparseMVProd op(H);
    // Construct eigen solver object, requesting the largest three eigenvalues
    int NCV = 0;
//    int NState = 0;
//    if(NEig > BasisSet.size()){
//        NState = BasisSet.size()-1;
//    }else{
//        NState = NEig;
//    }
    NCV = max(2*NEig+1,20); // Default from Scipy's Lanczos/Arnoldi implementation
    SymEigsSolver< double, SMALLEST_ALGE, SparseMVProd > eigs(&op, NEig, NCV);
    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(1000,1e-10,SMALLEST_ALGE);
    if(eigs.info() == SUCCESSFUL){
        E = eigs.eigenvalues().real();
        Psi = eigs.eigenvectors().real();
    }else{
        cout << "Error: Eigenvalues did not converge." << endl; exit(0);}
    return;
};

inline void DenseDiagonalize(MatrixXd& H, MatrixXd& Psi, VectorXd& E){
    //Wrapper for the Eigen diagonalization
    SelfAdjointEigenSolver<MatrixXd> SE; //Schrodinger equation
    SE.compute(H); //Diagonalize the matrix
    E = SE.eigenvalues().real(); //Extract frequencies
    Psi = SE.eigenvectors().real(); //Extract CI vectors
    return;
};

inline void QDiffVec(int n, int m, int& qtot, int& mchange, vector<int>& DiffVec)
{
    //Function that calculates ifference in quanta between basis states 
    // Pre-computing this saves time on inner loop over AnharmFC
    //qtot is the total number of changed quanta, mchange is number modes with changed quanta, DiffVec is a vector with number of changed quanta per mode
    for (unsigned int i=0;i<BasisSet[0].M;i++)
    {
        int qdiff = 0; //Number of quanta between states
        qdiff += BasisSet[n].Modes[i].Quanta;
        qdiff -= BasisSet[m].Modes[i].Quanta;
        DiffVec[i] = abs(qdiff);
//        DiffVec.push_back(abs(qdiff));
        qtot += abs(qdiff);
        if(qdiff!=0){
            mchange += 1;
        }
    }
    return;
};

inline bool ScreenState(int qdiff, int mchange, const vector<int>& QDiffVec, const FConst& fc)
{
    //Function for ignoring states that have no overlap
    //qdiff is the number of changed quanta, mchange is number modes with changed quanta, QDiffVec is a vector with number of changed quanta per mode, fc is force constant
    bool keepstate = 1; //Assume the state is good
    if(qdiff > fc.fcpow.size() || 
            mchange > fc.ShortModes.size()
            || qdiff%2 != fc.fcpow.size()%2
            ){
        return 0;
    }
    //Check based on force constant powers (check that raising and lowering results in quanta match)
    for (unsigned int i=0;i<fc.ShortModes.size();i++)
    {
        if ( QDiffVec[fc.ShortModes[i]] > fc.ModePowers[i] || 
                QDiffVec[fc.ShortModes[i]] % 2 != fc.ModePowers[i] % 2){
            //Impossible for the states to overlap if mode power is too small or wrong even/odd parity
            //Skip the rest of the checks
            return 0;
        }
    }
    //Check overlap of all other modes (modes not involved in FC)
    for (int i=0;i<QDiffVec.size();i++)
    {
        bool cont = 1; //Continue the check
        for (unsigned int j=0;j<fc.ShortModes.size();j++)
        {
            if (fc.ShortModes[j] == i)
            {
                //Ignore this mode since it is in the FC
                cont = 0;
                break; // No need to continue checking mode against fc.ShortModes
            }
        }
        if (cont)
        {
            if ( QDiffVec[i] != 0)
            {
                //Remove state due to zero overlap in mode i
                return 0; // No need to check if the other modes match 
            }
        }
    }
    //Return decision
    return keepstate;
};

void PrintFreqs(VectorXd& Freqs, MatrixXd& Vecs, fstream& outfile)
{
    //Function to print the first NEig eigenvalues
    //And their modes assignments
    int NFreqs = NEig; // set maximum number of printed eigenvalues
    if (Freqs.size() < NEig){ NFreqs = Freqs.size();} // set to number of eigenvalues if < 100
    outfile << fixed << setprecision(10) << Freqs(0) << "   "; // Print ZPE
    double maxcoeff = max(abs(Vecs.col(0).maxCoeff()),abs(Vecs.col(0).minCoeff()));
    for( unsigned int m=0; m<Vecs.rows(); m++ ){ // Looping over elements of eigenvector
        if(abs(Vecs(m,0))>=maxcoeff/2){ // Output the mode identities of the largest element(s)
            outfile << fixed << setprecision(2) << abs(Vecs(m,0)) << " (";
            int nmode = 0;
            for( unsigned int l=0; l<BasisSet[m].M; l++){
                if(BasisSet[m].Modes[l].Quanta>0){
                    if(nmode>0){
                    outfile << "+";
                    }
                    outfile << BasisSet[m].Modes[l].Quanta << "w" << l;
                    nmode += 1;
                }
            }
            outfile << ")";
        }
    }
    outfile << '\n';
    for( unsigned int n=1; n<NFreqs; ++n){
        outfile << fixed << setprecision(12) << Freqs(n)-Freqs(0) << "   "; // Print eigenvalue minus ZPE
        double maxcoeff = max(abs(Vecs.col(n).maxCoeff()),abs(Vecs.col(n).minCoeff()));
        int state = 0;
        for( unsigned int m=0; m<Vecs.rows(); m++ ){
            if(abs(Vecs(m,n))>=maxcoeff/2){
                if(state>0){
                    outfile << "  ";
                }
                outfile << fixed << setprecision(2) << abs(Vecs(m,n)) << " (";
                state += 1;
                int nmode = 0;
                for( unsigned int l=0; l<BasisSet[m].M; l++){
                    if(BasisSet[m].Modes[l].Quanta>0){
                        if(nmode>0){
                        outfile << "+";
                        }
                        outfile << BasisSet[m].Modes[l].Quanta << "w" << l;
                        nmode += 1;
                    }
                }
                outfile << ")";
            }
        }
        outfile << '\n';
    }
    //Print and return
    outfile.flush();
    return;
};
