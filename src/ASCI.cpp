/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Functions for generating the product basis via Adaptive Sampling CI algorithm

*/


void AddASCI(vector<WaveFunction>& NewBasis, MatrixXd& Evecs, VectorXd& Evals){
    int Ncore = 20000;
    int Ntarget = 20000;
    vector<WaveFunction> OldBasis = BasisSet; // Preserve the real basis
    cout << " Finding connected states..." << endl;
    vector<double> Cmax; // Vector of maximum coefficients over excited states
   // MatrixXd absEvecs = abs(Evecs);
    VectorXd maxVal = (Evecs.cwiseAbs()).rowwise().maxCoeff();
    for(unsigned int i=0; i<Evecs.rows(); i++){
        Cmax.push_back(maxVal[i]); 
    }
    maxVal.resize(0); // clear from memory
    int N_opt;
    if(NEig > BasisSet.size()){ // If we don't have enough states to optimize for yet
        N_opt = BasisSet.size();
    }else{ // If number of states exceeds the number selected to optimize
        N_opt = NEig;
    }
    if(Ncore > BasisSet.size()){ // Basis is small enough that all states are core states 
        Ncore = BasisSet.size();
    }
	vector<int> maxidx(Cmax.size()); // Contains indices, will be sorted by CI coefficient size
    int x=0;
	std::iota(maxidx.begin(),maxidx.end(),x++); //Initializing indices
	sort( maxidx.begin(),maxidx.end(), [&](int i,int j){return abs(Cmax[i])>abs(Cmax[j]);} );
    HashedStates HashedConnected; // hashed unordered_set of new states that only allows unique states to be inserted
    HashedStates HashedBasisInit; // hashed unordered_set of new states that only allows unique states to be inserted
    
    for(const WaveFunction& wfn : BasisSet){
        HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
    }
    for(unsigned int m=0; m<Ncore; m++){
        AddConnected(HashedBasisInit,HashedConnected, maxidx[m]); // Adds states connected to Ncore states with largest CI coefficients
    }
    vector<int>().swap(maxidx); // clear from memory
    int ConnectedSize = HashedConnected.size();
    cout << " Performing ASCI addition over " << ConnectedSize << " new connected states." << endl;
    for(const WaveFunction& wfn : HashedConnected){ // Add connected states to basis set 
        BasisSet.push_back(wfn);
    }
    HashedConnected = HashedStates(); // Free from memory
    HashedBasisInit = HashedStates(); // Free from memory
    swap(HashedConnected,HashedBasisInit); // Free from memory
    //// PT Basis has been constructed; now we calculate matrix elements ////
    int fcmax=0;
    for (unsigned int k=0;k<AnharmFC.size();k++){ //Order of maximum anharmonic term
        if(AnharmFC[k].fcpow.size()>fcmax){
            fcmax = AnharmFC[k].fcpow.size();
        }
    }
    //#pragma omp parallel for
    
    vector<double> newCmax; // Vector of maximum estimated coefficients over excited states
    for(unsigned int a=0; a<Evecs.rows()+ConnectedSize;a++){
        vector<double> SumHaiCi(N_opt,0.); // Estimated coefficient for basis state a in eigenvector for excited state n
        vector<int> qdiffvec(BasisSet[0].M,0);
        for(unsigned int i=0; i<Evecs.rows(); i++){
            if(i==a){continue;}
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
        //#pragma omp atomic // Will cause floating point error if blindly done in parallel
        double Aia = SumHaiCi[0]/(Evals(0)-Ea);
        for(int n=1; n<N_opt; n++){
            double Aia2 = SumHaiCi[n]/(Evals(n)-Ea);
            if(abs(Aia2)>abs(Aia)){
                Aia = Aia2;
            }
        }
        newCmax.push_back(Aia);
    }

    vector<int> newmaxidx(newCmax.size()); // Contains indices, will be sorted by CI coefficient size
    x=0;
	std::iota(newmaxidx.begin(),newmaxidx.end(),x++); //Initializing indices
	sort( newmaxidx.begin(),newmaxidx.end(), [&](int i,int j){return abs(newCmax[i])>abs(newCmax[j]);} );
    for(unsigned int i=0; i<Ntarget; i++){
        NewBasis.push_back(BasisSet[newmaxidx[i]]); // Add Ntarget states with largest wfn coefficients
    }
    BasisSet = OldBasis;
}

void AddConnected(HashedStates& BasisInit, HashedStates& NewStates, int n){ // Find all connected basis states
    for(unsigned int i=0; i<AnharmFC.size(); ++i){ // Loop over sorted force constants
        if(AnharmFC[i].ShortModes.size()==1){ // Stupid way to enumerate new states from F_ij
            for( int a=-(AnharmFC[i].ModePowers[0]); a<(AnharmFC[i].ModePowers[0]+1); a+=2 ){
                if( a != 0){// Skip if no change
                    WaveFunction tmp = BasisSet[n];
                    tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta += a;
                    int ntot = 0;
                    for(unsigned int p=0; p<tmp.M; p++){
                        ntot += tmp.Modes[p].Quanta;
                    }
                    if( tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta >=0 &&
                            BasisInit.count(tmp) == 0 &&
                            ntot < 11
                            ){ //make sure a|0> = 0 and tmp does not exist in original basis
                        NewStates.insert(tmp); // add new state to set
                    }
                }
            }
        }
        if(AnharmFC[i].ShortModes.size()==2){ // Stupid way to enumerate new states from F_ij
            for( int a=-(AnharmFC[i].ModePowers[0]); a<(AnharmFC[i].ModePowers[0]+1); a+=2 ){
                for( int b=-(AnharmFC[i].ModePowers[1]); b<(AnharmFC[i].ModePowers[1]+1); b+=2 ){
                    if(abs(a)+abs(b) != 0){// Skip if no change
                        WaveFunction tmp = BasisSet[n];
                        tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta += a;
                        tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta += b;
                        int ntot = 0;
                        for(unsigned int p=0; p<tmp.M; p++){
                            ntot += tmp.Modes[p].Quanta;
                        }
                        if( tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta >=0 &&
                               tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta >=0 && 
                               BasisInit.count(tmp) == 0 &&
                               ntot < 11
                               ){ //make sure a|0> = 0
                            NewStates.insert(tmp); // add new state to set
                        }
                    }
                }
            }
        }
        if(AnharmFC[i].ShortModes.size()==3){ // Stupid way to enumerate new states from F_ijk
            for( int a=-(AnharmFC[i].ModePowers[0]); a<(AnharmFC[i].ModePowers[0]+1); a+=2 ){
                for( int b=-(AnharmFC[i].ModePowers[1]); b<(AnharmFC[i].ModePowers[1]+1); b+=2 ){
                    for( int c=-(AnharmFC[i].ModePowers[2]); c<(AnharmFC[i].ModePowers[2]+1); c+=2 ){
                        if(abs(a)+abs(b)+abs(c) != 0){// Skip if no change
                            WaveFunction tmp = BasisSet[n];
                            tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta += a;
                            tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta += b;
                            tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta += c;
                            int ntot = 0;
                            for(unsigned int p=0; p<tmp.M; p++){
                                ntot += tmp.Modes[p].Quanta;
                            }
                            if( tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta >=0 &&
                                   tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta >=0  &&
                                   tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta >=0 &&
                                   BasisInit.count(tmp) == 0 &&
                                   ntot < 11
                                   ){ //make sure a|0> = 0
                                NewStates.insert(tmp); // add new state
                            }
                        }      
                    }
                }
            }
        }
        if(AnharmFC[i].ShortModes.size()==4){ // Stupid way to enumerate new states from F_ijkl
            for( int a=-(AnharmFC[i].ModePowers[0]); a<(AnharmFC[i].ModePowers[0]+1); a+=2 ){
                for( int b=-(AnharmFC[i].ModePowers[1]); b<(AnharmFC[i].ModePowers[1]+1); b+=2 ){
                    for( int c=-(AnharmFC[i].ModePowers[2]); c<(AnharmFC[i].ModePowers[2]+1); c+=2 ){
                        for( int d=-(AnharmFC[i].ModePowers[3]); d<(AnharmFC[i].ModePowers[3]+1); d+=2 ){
                            if(abs(a)+abs(b)+abs(c)+abs(d) != 0){ // Skip if no change
                                WaveFunction tmp = BasisSet[n];
                                tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta += a;
                                tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta += b;
                                tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta += c;
                                tmp.Modes[AnharmFC[i].ShortModes[3]].Quanta += d;
                                int ntot = 0;
                                for(unsigned int p=0; p<tmp.M; p++){
                                    ntot += tmp.Modes[p].Quanta;
                                }
                                if( tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta >=0 &&
                                       tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta >=0  &&
                                       tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta >=0 &&
                                       tmp.Modes[AnharmFC[i].ShortModes[3]].Quanta >=0 &&
                                       BasisInit.count(tmp) == 0 &&
                                       ntot < 11
                                       ){ //make sure a|0> = 0
                                    NewStates.insert(tmp); // add new state
                                }
                            }
                        }
                    }
                }
            }
        }   
        if(AnharmFC[i].ShortModes.size()==5){ // Stupid way to enumerate new states from F_ijklm 
            for( int a=-(AnharmFC[i].ModePowers[0]); a<(AnharmFC[i].ModePowers[0]+1); a+=2 ){
                for( int b=-(AnharmFC[i].ModePowers[1]); b<(AnharmFC[i].ModePowers[1]+1); b+=2 ){
                    for( int c=-(AnharmFC[i].ModePowers[2]); c<(AnharmFC[i].ModePowers[2]+1); c+=2 ){
                        for( int d=-(AnharmFC[i].ModePowers[3]); d<(AnharmFC[i].ModePowers[3]+1); d+=2 ){
                            for( int f=-(AnharmFC[i].ModePowers[4]); f<(AnharmFC[i].ModePowers[4]+1); f+=2 ){
                                
                                if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f) != 0){ // Skip if no change
                                    WaveFunction tmp = BasisSet[n];
                                    tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta += a;
                                    tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta += b;
                                    tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta += c;
                                    tmp.Modes[AnharmFC[i].ShortModes[3]].Quanta += d;
                                    tmp.Modes[AnharmFC[i].ShortModes[4]].Quanta += f;
                                    int ntot = 0;
                                    for(unsigned int p=0; p<tmp.M; p++){
                                        ntot += tmp.Modes[p].Quanta;
                                    }
                                    if( tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta >=0 &&
                                           tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta >=0  &&
                                           tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta >=0 &&
                                           tmp.Modes[AnharmFC[i].ShortModes[3]].Quanta >=0 &&
                                           tmp.Modes[AnharmFC[i].ShortModes[4]].Quanta >=0 && 
                                           BasisInit.count(tmp) == 0 &&
                                           ntot < 11
                                           ){ //make sure a|0> = 0
                                        NewStates.insert(tmp); // add new state
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(AnharmFC[i].ShortModes.size()==6){ // Stupid way to enumerate new states from F_ijklmn 
            for( int a=-(AnharmFC[i].ModePowers[0]); a<(AnharmFC[i].ModePowers[0]+1); a+=2 ){
                for( int b=-(AnharmFC[i].ModePowers[1]); b<(AnharmFC[i].ModePowers[1]+1); b+=2 ){
                    for( int c=-(AnharmFC[i].ModePowers[2]); c<(AnharmFC[i].ModePowers[2]+1); c+=2 ){
                        for( int d=-(AnharmFC[i].ModePowers[3]); d<(AnharmFC[i].ModePowers[3]+1); d+=2 ){
                            for( int f=-(AnharmFC[i].ModePowers[4]); f<(AnharmFC[i].ModePowers[4]+1); f+=2 ){
                                for( int g=-(AnharmFC[i].ModePowers[5]); g<(AnharmFC[i].ModePowers[5]+1); g+=2 ){
                                    if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f)+abs(g) != 0){ // Skip if no change
                                        WaveFunction tmp = BasisSet[n];
                                        tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta += a;
                                        tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta += b;
                                        tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta += c;
                                        tmp.Modes[AnharmFC[i].ShortModes[3]].Quanta += d;
                                        tmp.Modes[AnharmFC[i].ShortModes[4]].Quanta += f;
                                        tmp.Modes[AnharmFC[i].ShortModes[5]].Quanta += g;
                                        int ntot = 0;
                                        for(unsigned int p=0; p<tmp.M; p++){
                                            ntot += tmp.Modes[p].Quanta;
                                        }
                                        if( tmp.Modes[AnharmFC[i].ShortModes[0]].Quanta >=0 &&
                                               tmp.Modes[AnharmFC[i].ShortModes[1]].Quanta >=0 &&
                                               tmp.Modes[AnharmFC[i].ShortModes[2]].Quanta >=0 &&
                                               tmp.Modes[AnharmFC[i].ShortModes[3]].Quanta >=0 &&
                                               tmp.Modes[AnharmFC[i].ShortModes[4]].Quanta >=0 &&
                                               tmp.Modes[AnharmFC[i].ShortModes[5]].Quanta >=0 &&
                                               BasisInit.count(tmp) == 0 &&
                                               ntot < 11
                                               ){ //make sure a|0> = 0
                                            NewStates.insert(tmp); // add new state
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(AnharmFC[i].fcpow.size() == 0 || AnharmFC[i].fcpow.size() > 6 ){
            //Print error message
            cout << "Error: HCI only works with force constants up to 6th order." << endl;
            cout.flush();
            exit(0);
        }
    }
}

void Perform_ASCI(MatrixXd& Evecs, VectorXd& Evals, const string& cpname){
    bool is_converged = 0;
    int iter_count = 1;
    if(!restart){
        // Initialize calculation from beginning
        cout << " Initial basis size is: " << BasisSet.size() << "\n";
        MatrixXd VCIHam(BasisSet.size(),BasisSet.size());
        cout << "got this far!" << endl;
        VCIHam.setZero();
        cout << " Constructing Hamiltonian...";
        ZerothHam(VCIHam);
        AnharmHam(VCIHam);
        cout << " Diagonalizing..." << endl;
        DenseDiagonalize(VCIHam,Evecs,Evals); //Diagonalization wrapper
        cout.precision(2);
        cout << " ZPE is: " << Evals(0) << "\n" << endl;
        //Free memory
        VCIHam = MatrixXd(1,1);
    }
    while(!is_converged){
        cout << " Performing ASCI iteration: " << iter_count << endl;
        vector<WaveFunction> NewStates; // New states being added via CIPSI algorithm
        AddASCI(NewStates, Evecs, Evals);
        BasisSet = NewStates; // REPLACE Basis set with new one rather than combine
        cout << " Basis size after iteration " << iter_count << " is: " << BasisSet.size() << endl;
        ++iter_count;
        cout << " Constructing Hamiltonian...";
        cout.flush();
        SpMat HSparse(BasisSet.size(),BasisSet.size());
        MakeHamSparse(HSparse);
        if(NEig==1){
            cout << " Finding the ground state..." << endl;
        }else{
            cout << " Finding the first " << NEig << " eigenstates..." << endl;
        }
        vector<double> old_evals;
        for(unsigned int i=0; i<Evals.size(); i++){
            old_evals.push_back(Evals(i));
        }
        SparseDiagonalize(HSparse,Evecs,Evals); //Diagonalization wrapper
        cout.precision(2);
        cout << " ZPE is: " << Evals(0) << endl << endl;
        cout.flush();
        if(BasisSet.size()>100000){ 
            cout << " Writing checkpoint file...";
            cout.flush();
            fstream checkfile;
            checkfile.open(cpname,ios_base::out | ios_base::trunc);
            if (!checkfile.good() || cpname=="none"){
                cout << '\n' << " Error: checkpoint file could not be found or no name was given." << '\n' << endl;
            }else{
                WriteCheckpoint(Evals,Evecs,checkfile);
                cout << "done." << '\n' << endl;
            }
            checkfile.close();
        }
        bool converge_test = true;
        for(unsigned int i=0; i<Evals.size(); i++){  // Check if any eigenvalues changed by more than 0.5
            double diff = abs(Evals(i)-old_evals[i]);
            if(diff>0.5){
                converge_test=false;
                break;
            }
        }
        if(converge_test){
            is_converged = true;
        }
    }
    return;
}
