/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Functions for generating the heat bath product basis

*/


bool sortByFC(const FConst &lhs, const FConst &rhs) { return abs(lhs.fc) > abs(rhs.fc); } 
// Sort force constants from large to small magnitude

void HeatBath_Sort_FC(){
    sort(AnharmHB.begin(), AnharmHB.end(), sortByFC); // Sort force constants from large to small magnitude
}

void AddStatesHB(HashedStates& BasisInit, HashedStates& NewStates, int n, double Cn, double eps){ // Expand basis via Heat Bath algorithm
    for(unsigned int i=0; i<AnharmHB.size(); ++i){ // Loop over sorted force constants
        if(abs(Cn*AnharmHB[i].fc) >= eps){ // States connected by fc will be added if |fc*Cn| >= eps
            if(AnharmHB[i].ShortModes.size()==1){ // Stupid way to enumerate new states from F_ij
                for( int a=-(AnharmHB[i].ModePowers[0]); a<(AnharmHB[i].ModePowers[0]+1); a+=2 ){
                    if( a != 0){// Skip if no change
                        WaveFunction tmp = BasisSet[n];
                        tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta += a;
                        if( tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta >=0 &&
							    BasisInit.count(tmp) == 0
                                ){ //make sure a|0> = 0 and tmp does not exist in original basis
                            NewStates.insert(tmp); // add new state to set
                        }
                    }
                }
            }
            if(AnharmHB[i].ShortModes.size()==2){ // Stupid way to enumerate new states from F_ij
                for( int a=-(AnharmHB[i].ModePowers[0]); a<(AnharmHB[i].ModePowers[0]+1); a+=2 ){
                    for( int b=-(AnharmHB[i].ModePowers[1]); b<(AnharmHB[i].ModePowers[1]+1); b+=2 ){
                        if(abs(a)+abs(b) != 0){// Skip if no change
                            WaveFunction tmp = BasisSet[n];
                            tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta += a;
                            tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta += b;
                            if( tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta >=0 &&
                                   tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta >=0 && 
								   BasisInit.count(tmp) == 0
                                   ){ //make sure a|0> = 0
                                NewStates.insert(tmp); // add new state to set
                            }
                        }
                    }
                }
            }
            if(AnharmHB[i].ShortModes.size()==3){ // Stupid way to enumerate new states from F_ijk
                for( int a=-(AnharmHB[i].ModePowers[0]); a<(AnharmHB[i].ModePowers[0]+1); a+=2 ){
                    for( int b=-(AnharmHB[i].ModePowers[1]); b<(AnharmHB[i].ModePowers[1]+1); b+=2 ){
                        for( int c=-(AnharmHB[i].ModePowers[2]); c<(AnharmHB[i].ModePowers[2]+1); c+=2 ){
                            if(abs(a)+abs(b)+abs(c) != 0){// Skip if no change
                                WaveFunction tmp = BasisSet[n];
                                tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta += a;
                                tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta += b;
                                tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta += c;
                                if( tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta >=0 &&
                                       tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta >=0  &&
                                       tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta >=0 &&
                                       BasisInit.count(tmp) == 0
                                       ){ //make sure a|0> = 0
                                    NewStates.insert(tmp); // add new state
                                }
                            }      
                        }
                    }
                }
            }
            if(AnharmHB[i].ShortModes.size()==4){ // Stupid way to enumerate new states from F_ijkl
                for( int a=-(AnharmHB[i].ModePowers[0]); a<(AnharmHB[i].ModePowers[0]+1); a+=2 ){
                    for( int b=-(AnharmHB[i].ModePowers[1]); b<(AnharmHB[i].ModePowers[1]+1); b+=2 ){
                        for( int c=-(AnharmHB[i].ModePowers[2]); c<(AnharmHB[i].ModePowers[2]+1); c+=2 ){
                            for( int d=-(AnharmHB[i].ModePowers[3]); d<(AnharmHB[i].ModePowers[3]+1); d+=2 ){
                                if(abs(a)+abs(b)+abs(c)+abs(d) != 0){ // Skip if no change
                                    WaveFunction tmp = BasisSet[n];
                                    tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta += a;
                                    tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta += b;
                                    tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta += c;
                                    tmp.Modes[AnharmHB[i].ShortModes[3]].Quanta += d;
                                    if( tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta >=0 &&
                                           tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta >=0  &&
                                           tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta >=0 &&
                                           tmp.Modes[AnharmHB[i].ShortModes[3]].Quanta >=0 &&
                                           BasisInit.count(tmp) == 0
                                           ){ //make sure a|0> = 0
                                        NewStates.insert(tmp); // add new state
                                    }
                                }
                            }
                        }
                    }
                }
            }   
            if(AnharmHB[i].ShortModes.size()==5){ // Stupid way to enumerate new states from F_ijklm 
                for( int a=-(AnharmHB[i].ModePowers[0]); a<(AnharmHB[i].ModePowers[0]+1); a+=2 ){
                    for( int b=-(AnharmHB[i].ModePowers[1]); b<(AnharmHB[i].ModePowers[1]+1); b+=2 ){
                        for( int c=-(AnharmHB[i].ModePowers[2]); c<(AnharmHB[i].ModePowers[2]+1); c+=2 ){
                            for( int d=-(AnharmHB[i].ModePowers[3]); d<(AnharmHB[i].ModePowers[3]+1); d+=2 ){
                                for( int f=-(AnharmHB[i].ModePowers[4]); f<(AnharmHB[i].ModePowers[4]+1); f+=2 ){
                                    if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f) != 0){ // Skip if no change
                                        WaveFunction tmp = BasisSet[n];
                                        tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta += a;
                                        tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta += b;
                                        tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta += c;
                                        tmp.Modes[AnharmHB[i].ShortModes[3]].Quanta += d;
                                        tmp.Modes[AnharmHB[i].ShortModes[4]].Quanta += f;
                                        if( tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta >=0  &&
                                               tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].ShortModes[3]].Quanta >=0 &&
                                               tmp.Modes[AnharmHB[i].ShortModes[4]].Quanta >=0 && 
                                               BasisInit.count(tmp) == 0
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
            if(AnharmHB[i].ShortModes.size()==6){ // Stupid way to enumerate new states from F_ijklmn 
                for( int a=-(AnharmHB[i].ModePowers[0]); a<(AnharmHB[i].ModePowers[0]+1); a+=2 ){
                    for( int b=-(AnharmHB[i].ModePowers[1]); b<(AnharmHB[i].ModePowers[1]+1); b+=2 ){
                        for( int c=-(AnharmHB[i].ModePowers[2]); c<(AnharmHB[i].ModePowers[2]+1); c+=2 ){
                            for( int d=-(AnharmHB[i].ModePowers[3]); d<(AnharmHB[i].ModePowers[3]+1); d+=2 ){
                                for( int f=-(AnharmHB[i].ModePowers[4]); f<(AnharmHB[i].ModePowers[4]+1); f+=2 ){
                                    for( int g=-(AnharmHB[i].ModePowers[5]); g<(AnharmHB[i].ModePowers[5]+1); g+=2 ){
                                        if(abs(a)+abs(b)+abs(c)+abs(d)+abs(f)+abs(g) != 0){ // Skip if no change
                                            WaveFunction tmp = BasisSet[n];
                                            tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta += a;
                                            tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta += b;
                                            tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta += c;
                                            tmp.Modes[AnharmHB[i].ShortModes[3]].Quanta += d;
                                            tmp.Modes[AnharmHB[i].ShortModes[4]].Quanta += f;
                                            tmp.Modes[AnharmHB[i].ShortModes[5]].Quanta += g;
                                            if( tmp.Modes[AnharmHB[i].ShortModes[0]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].ShortModes[1]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].ShortModes[2]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].ShortModes[3]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].ShortModes[4]].Quanta >=0 &&
                                                   tmp.Modes[AnharmHB[i].ShortModes[5]].Quanta >=0 &&
												   BasisInit.count(tmp) == 0
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
            if(AnharmHB[i].fcpow.size() == 0 || AnharmHB[i].fcpow.size() > 6 ){
                //Print error message
                cout << "Error: HCI only works with force constants up to 6th order." << endl;
                cout.flush();
                exit(0);
            }
        }
        else{break;}// break loop if you've reached element < eps, since all future elements will be smaller (HB sorting)
    }
}

void Perform_HCI(MatrixXd& Evecs, VectorXd& Evals, const string& cpname){
    HeatBath_Sort_FC();
    bool is_converged = 0;
    int iter_count = 1;
    if(!restart){
        // Initialize calculation from beginning
        cout << " Initial basis size is: " << BasisSet.size() << "\n";
        MatrixXd VCIHam(BasisSet.size(),BasisSet.size());
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
        cout << " Performing HCI iteration: " << iter_count << endl;
        vector<WaveFunction> NewStates; // New states being added via HB algorithm
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
        //VectorXd maxVal = (Evecs.cwiseAbs()).rowwise().maxCoeff();
        HashedStates HashedBasisInit; // hashed unordered_set containing BasisSet to check for duplicates
        HashedStates HashedNewStates; // hashed unordered_set of new states that only allows unique states to be inserted
        for( WaveFunction& wfn : BasisSet){
            HashedBasisInit.insert(wfn); // Populate hashed unordered_set with initial basis states
        }
        for(unsigned int n=0; n<Cmax.size(); ++n){ // Loop over max CI coefficients and add basis states
            AddStatesHB(HashedBasisInit, HashedNewStates,n,Cmax[n],HCI_Eps);
        }
        if((double)HashedNewStates.size()/(double)BasisSet.size() <= 0.01){
            is_converged = 1; // Is converged in basis only grows by 1%
            cout << " Basis converged after " << iter_count-1 << " iterations and contains " << BasisSet.size() << " product states." << endl;
            break;
        }else{
            for(const WaveFunction& wfn : HashedNewStates){ // Add NewStates to the BasisSet
                BasisSet.push_back(wfn);
            }
        }
        HashedBasisInit = HashedStates(); // Free from memory
        HashedNewStates = HashedStates(); // Free from memory
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
    }
    return;
}
