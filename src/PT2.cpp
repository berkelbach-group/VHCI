/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

 Implementation of Epstein-Nesbet PT2 with Heat Bath sorting

*/

void DoPT2(MatrixXd& Evecs, VectorXd& Evals){
    if(HCI_Eps==0){
        HeatBath_Sort_FC();
    }
    cout << " Finding connected states..." << endl;
    int N_opt;
    if(NEig > BasisSet.size()){ // If we don't have enough states to optimize for yet
        N_opt = BasisSet.size();
    }else{ // If number of states exceeds the number selected to optimize
        N_opt = NEig;
    }
    vector<double> Cmax; // Vector of maximum coefficients over excited states
        // MatrixXd absEvecs = abs(Evecs);
    for(unsigned int i=0; i<Evecs.rows(); i++){
        double Cij = Evecs(i,0);
        for(unsigned int j=1; j<N_opt; j++){
            if(abs(Evecs(i,j))>abs(Cij)){
                Cij = Evecs(i,j);
            }
        }
        Cmax.push_back(Cij);
    }
    HashedStates HashedBasisInit; // hashed unordered_set containing BasisSet to check for duplicates
    HashedStates HashedPTBasis; // hashed unordered_set of new states that only allows unique states to be inserted
    for(const WaveFunction& wfn : BasisSet){
        HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
    }
    for(unsigned int n=0; n<Cmax.size(); n++){ // Loop over max CI coefficients and add basis states
        AddStatesHB(HashedBasisInit,HashedPTBasis,n,Cmax[n],PT2_Eps);
    }
    int PTBasisSize = HashedPTBasis.size();
    cout << " Perturbative space contains " << PTBasisSize << " states." << endl;
    for(const WaveFunction& wfn : HashedPTBasis){ // Add NewStates to the BasisSet
        BasisSet.push_back(wfn);
    }
    HashedPTBasis = HashedStates(); // Free from memory
    HashedBasisInit = HashedStates(); // Free from memory
    //// PT Basis has been constructed; now we calculate matrix elements ////
    if(NEig==1){
        cout << " Calculating the 2nd-order perturbative correction on the ground state energy" << endl;
    }else{
        cout << " Calculating the 2nd-order perturbative correction for on first " << N_opt << " eigenvalues." << endl;
    }
    vector<double> DeltaE(N_opt,0.);  // Vector will contain the PT correction for each eigenvalue
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }  
    #pragma omp parallel for
    for(unsigned int a=Evecs.rows(); a<Evecs.rows()+PTBasisSize;a++){
        vector<double> SumHaiCi(N_opt,0.);
        vector<int> qdiffvec(BasisSet[0].M,0);
        for(unsigned int i=0; i<Evecs.rows(); i++){
            double Hai = 0;
            int mchange = 0; // number of modes with nonzero change in quanta
            int qdiff = 0; // total number of quanta difference 
            QDiffVec(a,i,qdiff,mchange,qdiffvec);
            if(qdiff <= fcmax && mchange <= fcmax && qdiff%2==0){ 
                // States cannot differ by more than fcmax quanta
                for (unsigned int k=0;k<QuarticFC.size();k++){
                    if ( ScreenState(qdiff,mchange,qdiffvec,QuarticFC[k]) ){
                        //Screen force constants for connection                      
                        //Add anharmonic matrix elements
                        Hai += AnharmPot(a,i,QuarticFC[k]);
                    }
                }
                for (unsigned int k=0;k<SexticFC.size();k++){
                    if ( ScreenState(qdiff,mchange,qdiffvec,SexticFC[k]) ){
                        //Screen force constants for connection                      
                        //Add anharmonic matrix elements
                        Hai += AnharmPot(a,i,SexticFC[k]);
                    }
                }
                for(int n=0; n<N_opt; n++){
               //     #pragma omp critical
                    SumHaiCi[n] += Hai*Evecs(i,n); // C_i Hai for each eigenvalue of interest
                }
            }
            if(qdiff <= fcmax-1 && mchange <= fcmax-1 && qdiff%2==1){ 
                // fcmax-1 assumes 4th or 6th max order 
                // States cannot differ by more than fcmax quanta
                for (unsigned int k=0;k<CubicFC.size();k++){
                    if ( ScreenState(qdiff,mchange,qdiffvec,CubicFC[k]) ){
                        //Screen force constants for connection                      
                        //Add anharmonic matrix elements
                        Hai += AnharmPot(a,i,CubicFC[k]);
                    }
                }
                for (unsigned int k=0;k<QuinticFC.size();k++){
                    if ( ScreenState(qdiff,mchange,qdiffvec,QuinticFC[k]) ){
                        //Screen force constants for connection                      
                        //Add anharmonic matrix elements
                        Hai += AnharmPot(a,i,QuinticFC[k]);
                    }
                }
                for(int n=0; n<N_opt; n++){
               //     #pragma omp critical
                    SumHaiCi[n] += Hai*Evecs(i,n); // C_i Hai for each eigenvalue of interest
                }
            }
        }
        double Ea = 0.; //Hii matrix element
        for (unsigned int j=0;j<BasisSet[a].M;j++){
          //Calculate partial energies
          double Ej = 0.5;
          Ej += BasisSet[a].Modes[j].Quanta;
          Ej *= BasisSet[a].Modes[j].Freq;
          //Update matrix element
          Ea += Ej;
        }
        vector<int> zerodiffvec(BasisSet[0].M,0);
        int qdiff=0;
        int mchange=0;
        for (unsigned int k=0;k<QuarticFC.size();k++){ // Only even-ordered fc can affect this
            if ( ScreenState(qdiff,mchange,zerodiffvec,QuarticFC[k]) ){
            // Screen force constants that cannot connect basis states a and a
            //Add anharmonic matrix elements
                Ea += AnharmPot(a,a,QuarticFC[k]);
            }
        }
        for (unsigned int k=0;k<SexticFC.size();k++){ // Only even-ordered fc can affect this
            if ( ScreenState(qdiff,mchange,zerodiffvec,SexticFC[k]) ){
            // Screen force constants that cannot connect basis states a and a
            //Add anharmonic matrix elements
                Ea += AnharmPot(a,a,SexticFC[k]);
            }
        }
        for(int n=0; n<N_opt; n++){
            #pragma omp atomic // Will cause floating point error if blindly done in parallel
            DeltaE[n] += pow(SumHaiCi[n],2)/(Evals(n)-Ea);
        }
    }
    for(int n=0; n<N_opt; n++){
        Evals(n) += DeltaE[n];
    }
    cout << " New ZPE after PT2 correction is: " << Evals(0) << endl;
}
