/*

#############################################################################
#                                                                           #
#              Vibrational Heat-Bath Configuration Interaction              # 
#                         By: Jonathan H. Fetherolf                         #
#                                                                           #
#                      Based on LOVCI by Eric G. Kratz                      #
#                                                                           #
#############################################################################

Parses the input command and file and initializes the calculation

*/

//Input reader
void ReadCIArgs(int argc, char* argv[], fstream& vcidata, fstream& freqfile, string& cpfilename)
{
    //Function to read the command line arguments
    cpfilename = "none";
    bool DoQuit = 0; //Exit if an error is detected
    string dummy; //Generic string
    stringstream call; //Stream for system calls and reading/writing files
    //Read command line arguments
    if (argc == 1)
    {
        //Escape if there are no arguments
        cout << '\n';
        cout << "Missing arguments...";
        cout << '\n' << '\n';
        cout << "Usage: vhci -n Ncpus -i Modes.inp -o Freqs.dat -c Checkpoint.txt";
        cout << '\n' << '\n';
        cout << "Use -h or --help for detailed instructions.";
        cout << '\n' << '\n';
        cout.flush();
        exit(0);
    }
    if ((argc % 2) != 1)
    {
        dummy = string(argv[1]);
        if ((dummy != "-h") and (dummy != "--help"))
        {
            //Escape if there are missing arguments
            cout << '\n';
            cout << "Wrong number of arguments...";
            cout << '\n' << '\n';
            cout << "Usage: vhci -n Ncpus -i Modes.inp -o Freqs.dat -c Checkpoint.txt";
            cout << '\n' << '\n';
            cout << "Use -h or --help for detailed instructions.";
            cout << '\n' << '\n';
            cout.flush();
            exit(0);
        }
    }
    for (int i=0;i<argc;i++)
    {
        //Read file names and CPUs
        dummy = string(argv[i]);
        if ((dummy == "-h") or (dummy == "--help"))
        {
            //Print helpful information and exit
            cout << '\n';
            cout << "Usage: vhci -n Ncpus -i Modes.inp -o Freqs.txt";
            cout << '\n' << '\n';
            cout << "Command line arguments:";
            cout << '\n' << '\n';
            cout << "  -n    Number of CPUs to use for the calculation.";
            cout << '\n' << '\n';
            cout << "  -i    Input file.";
            cout << '\n' << '\n';
            cout << "  -o    Output file for the first NEig VCI frequencies.";
            cout << '\n' << '\n';
            cout << "  -c    Checkpoint file for HCI iterations on large systems.";
            cout << '\n' << '\n';
            cout.flush();
            exit(0);
        }
        if (dummy == "-n")
        {
            Ncpus = atoi(argv[i+1]);
        }
        if (dummy == "-i")
        {
            vcidata.open(argv[i+1],ios_base::in);
        }
        if (dummy == "-o")
        {
            freqfile.open(argv[i+1],ios_base::out);
        }
        if (dummy == "-c"){
            cpfilename = argv[i+1];
        }
    }
    //Check for argument errors
    if (Ncpus < 1)
    {
        //Checks the number of threads and continue
        cout << " Warning: Calculations cannot run with ";
        cout << Ncpus << " CPUs.";
        cout << '\n';
        cout << "  Do you know how computers work?";
        cout << " Ncpus set to 1";
        cout << '\n';
        Ncpus = 1;
        cout.flush(); //Print warning
    }
    if (Ncpus > FindMaxThreads())
    {
        cout << " Warning: Too many threads requested.";
        cout << '\n';
        cout << "  Requested: " << Ncpus << ",";
        cout << " Available: " << FindMaxThreads();
        cout << '\n';
        cout << "  Ncpus set to " << FindMaxThreads();
        cout << '\n';
        Ncpus = FindMaxThreads();
        cout.flush(); //Print warning
    }
    if (!vcidata.good())
    {
        //Check input file
        cout << " Error: Could not open input file.";
        cout << '\n';
        DoQuit = 1;
    }
    if (!freqfile.good())
    {
        //Check output file
        cout << " Error: Could not output file.";
        cout << '\n';
        DoQuit = 1;
    }

    if (DoQuit)
    {
        //Quit if there is an error
        cout << '\n';
        cout.flush();
        exit(0);
    }
    else
    {
        //Sarcastically continue
        cout << '\n';
        cout << "No fatal errors detected.";
        cout << '\n';
        cout << " And there was much rejoicing. Yay...";
        cout << '\n' << '\n';
        cout.flush();
    }
    if(cpfilename=="none"){
        cout << "No checkpoint filename specified. Checkpoints will not be saved." << '\n' << '\n'; 
    }
    else{
        cout << "Locating checkpoint file..." << '\n';
        IsRestart(cpfilename);
        if(restart==1){
            cout << "Checkpoint file found. Calculation will be resumed." << '\n' << '\n';
        }else{
            cout << "No checkpoint file found. Calculation will start from the beginning." << '\n' << '\n';
        }
        cout.flush();
    }
    //Set threads
    omp_set_num_threads(Ncpus);
    Eigen::setNbThreads(Ncpus);
    return;
};

void ReadCIInput(fstream& vcidata)
{
    //Function to read the input files
    string dummy; //Generic sting
    //Count basis functions and read modes
    int Nmodes = 0; //Number of different modes
    int Ntot = 0; //Maximum quanta in a single product state
    vector<HOFunc> BasisCount; //Temp. storage of modes
    int Nfc = 0; //Number of force constants
    //Heat bath parameters 
    getline(vcidata,dummy); //Junk (String with comment for input file)
    vcidata >> dummy; //Junk
    vcidata >> HCI_Eps; // Variational HCI cutoff energy
    vcidata >> dummy;
    vcidata >> NEig; // Number of eigenstates to include in Heat Bath optimization
    // Set PT2 parameters
    vcidata >> dummy;
    vcidata >> perturb;
    vcidata >> dummy;
    vcidata >> PT2_Eps;
    if (perturb != 0 && perturb != 1)
    {
        //Print an error message
        cout << "Error: perturb must be either 0 or 1 ";
        cout << '\n' << '\n';
        cout.flush(); //Print message
        //Quit
        exit(0);
    }
    if(perturb==0){
        PT2_Eps=0;
    }
    //Set maximum simultaneous excitations 
    vcidata >> dummy; //Junk
    vcidata >> Ntot; //Maximum quanta in a single product state
    //Read active modes
    vcidata >> dummy; //Clear junk
    vcidata >> Nmodes; //Read modes
    for (int i=0;i<Nmodes;i++)
    {
        //Actual vibrational modes
        HOFunc tmp;
        int modeid;
        vcidata >> modeid; //Read mode ID
        //Check mode order
        if (modeid != i)
        {
            //Print an error message
            cout << "Error: Expected mode " << i;
            cout << " but read data for mode " << modeid;
            cout << '\n' << '\n';
            cout.flush(); //Print message
            //Quit
            exit(0);
        }
        vcidata >> tmp.Freq; //Frequency
        vcidata >> tmp.Quanta; //Max number of quanta
        BasisCount.push_back(tmp);
    }
    //Read anharmonic force constants
    vcidata >> dummy; //Clear junk
    vcidata >> Nfc; //Read number of anharmonic force constants
    for (int i=0;i<Nfc;i++)
    {
        //Save force constant data
        FConst tmp;
        int fcpower = 0;
        vcidata >> fcpower;
        for (int j=0;j<fcpower;j++)
        {
            int modej = 0;
            vcidata >> modej;
            tmp.fcpow.push_back(modej);
        }
        vcidata >> tmp.fc; //Read force constant value
        AnharmFC.push_back(tmp);
    }

    //Sort the modes and powers into ShortModes and ModePowers
    for(unsigned int n=0; n<AnharmFC.size(); n++){
        for (unsigned int i=0;i<AnharmFC[n].fcpow.size();i++){
            if (i == 0){
                //Initialize the vector
                AnharmFC[n].ShortModes.push_back(AnharmFC[n].fcpow[i]);
                AnharmFC[n].ModePowers.push_back(1);
            }
            else{
                //Update powers for i > 0
                bool foundmode = 0;
                for (unsigned int j=0;j<AnharmFC[n].ShortModes.size();j++){
                    if (AnharmFC[n].fcpow[i] == AnharmFC[n].ShortModes[j]){
                        AnharmFC[n].ModePowers[j] += 1;
                        foundmode = 1;
                    }
                }   
                if (!foundmode){
                    AnharmFC[n].ShortModes.push_back(AnharmFC[n].fcpow[i]);
                    AnharmFC[n].ModePowers.push_back(1);
                }
            }
        }
    }
    ScaleFC();
    for(unsigned int n=0; n<AnharmFC.size(); n++){
        switch(AnharmFC[n].fcpow.size()){
            case 3:
                CubicFC.push_back(AnharmFC[n]);
                break;
            case 4:
                QuarticFC.push_back(AnharmFC[n]);
                break;
            case 5:
                QuinticFC.push_back(AnharmFC[n]);
                break;
            case 6:
                SexticFC.push_back(AnharmFC[n]);
                break;
            default:
                cout << "Error: Force constants are only supported up to 6th order!" << endl;
                exit(0);
        }
    }
    

    
    //Create product basis
    if(restart){
        WaveFunction wftemp;
        for(unsigned int i=0; i<BasisCount.size(); i++){
            // Make first basis state
            HOFunc hotemp;
            hotemp.Freq = BasisCount[i].Freq;
            hotemp.Quanta = 0;
            wftemp.Modes.push_back(hotemp);
        }
        wftemp.M = wftemp.Modes.size();
        BasisSet.push_back(wftemp);
    }
    else{
        int Nexc = 0;
        if(BasisCount[0].Quanta+1 > Ntot){
            Nexc = Ntot+1;
        }else{
            Nexc = BasisCount[0].Quanta+1; 
        }
        for (int i=0;i<Nexc;i++)
        {
            WaveFunction temp;
            HOFunc tmp;
            tmp.Freq = BasisCount[0].Freq;
            tmp.Quanta = i;
            temp.Modes.push_back(tmp);
            BasisSet.push_back(temp);
        }
        for (unsigned int i=1;i<BasisCount.size();i++)
        {
            vector<WaveFunction> NewBasis;
            if(BasisCount[i].Quanta+1 > Ntot){
                Nexc = Ntot+1;
            }else{
                Nexc = BasisCount[i].Quanta+1; 
            }
            for (int j=0;j<Nexc;j++)
            {
                vector<WaveFunction> BasisCopy = BasisSet;
                for (unsigned int k=0;k<BasisSet.size();k++)
                {
                    HOFunc tmp;
                    //Copy data
                    tmp.Freq = BasisCount[i].Freq;
                    //Set quanta
                    tmp.Quanta = j;
                    //Add to arrays
                    BasisCopy[k].Modes.push_back(tmp);
                }
                for (unsigned int n=0;n<BasisCopy.size();n++)
                {
                    int Modes_tot = BasisCopy[n].Modes.size();
                    int Quanta_tot = 0; //Count total quanta of each product state
                    for(int j=0; j<Modes_tot; ++j){
                        Quanta_tot += BasisCopy[n].Modes[j].Quanta;
                    }
                    if(Quanta_tot<=Ntot){ //Limit total quanta of product state to Ntot
                        NewBasis.push_back(BasisCopy[n]);
                    }
                }
            }
            BasisSet = NewBasis;
        }
        //Correct array lengths
        for (unsigned int i=0;i<BasisSet.size();i++)
        {
            BasisSet[i].M = BasisSet[i].Modes.size();
        }
    }

    AnharmHB = AnharmFC; // AnharmHB contains all the unmodified coupling terms
    
    // Add single excitations
    
    for(unsigned int m=0; m<BasisSet[0].M; m++){ // Add single excitation terms to AnharmHB
        FConst temp;
        temp.fcpow.push_back(m);
        temp.ShortModes.push_back(m);
        temp.ModePowers.push_back(1);
        temp.fc = 0.;
        for(unsigned int n=0; n<CubicFC.size(); n++){ // Add single excitation terms from CubicFC
            if(CubicFC[n].ShortModes.size()==1 && CubicFC[n].ShortModes[0]==m){ // V_mmm
                temp.fc += 3*CubicFC[n].fc;
            }
            if(CubicFC[n].ShortModes.size()==2){
                for( unsigned int i=0; i<2; i++){ // loop over ShortModes
                    if(CubicFC[n].ShortModes[i]==m && CubicFC[n].ModePowers[i]==1){ // V_mnn or V_mmn
                        temp.fc += 2*CubicFC[n].fc;
                        break;
                    }
                }
            }
        }
        if(abs(temp.fc) > 1e-8){
            AnharmHB.push_back(temp);
        }
    }
    
    // Add double excitations

    for(unsigned int m=0; m<BasisSet[0].M; m++){ // Add single excitation terms to AnharmHB
        // m=n case
        FConst temp;
        temp.fcpow.push_back(m);
        temp.fcpow.push_back(m); // twice because m=n
        temp.ShortModes.push_back(m);
        temp.ModePowers.push_back(2); // power = 2 because m=n
        temp.fc = 0.;
        for(unsigned int n=0; n<QuarticFC.size(); n++){ // Add single excitation terms from CubicFC
            if(QuarticFC[n].ShortModes.size()==1 && QuarticFC[n].ShortModes[0]==m){ // V_mmmm
                temp.fc += 4*QuarticFC[n].fc;
            }
            if(QuarticFC[n].ShortModes.size()==2){
                for( unsigned int i=0; i<2; i++){ // loop over ShortModes
                    if(QuarticFC[n].ShortModes[i]==m && QuarticFC[n].ModePowers[i]==2){ // V_mmnn
                        temp.fc += 2*QuarticFC[n].fc;
                        break;
                    }
                }
            }
        }
        if(abs(temp.fc) > 1e-8){
            AnharmHB.push_back(temp);
        }
    }   
    
    for(unsigned int m=0; m<BasisSet[0].M-1; m++){ // Add single excitation terms to AnharmHB
        // m!=n case
        for(unsigned int n=m+1; n<BasisSet[0].M; n++){ //ShortModes and fcpow should be sorted in ascending order
            FConst temp;
            temp.fcpow.push_back(m);
            temp.fcpow.push_back(n);
            temp.ShortModes.push_back(m);
            temp.ShortModes.push_back(n);
            temp.ModePowers.push_back(1);
            temp.ModePowers.push_back(1);
            temp.fc = 0.;
            for(unsigned int i=0; i<QuarticFC.size(); i++){ // Add single excitation terms from CubicFC
                vector<int>::iterator findpow3 = find(QuarticFC[i].ModePowers.begin(), QuarticFC[i].ModePowers.end(), 3);
                vector<int>::iterator findm = find(QuarticFC[i].ShortModes.begin(), QuarticFC[i].ShortModes.end(), m);
                vector<int>::iterator findn = find(QuarticFC[i].ShortModes.begin(), QuarticFC[i].ShortModes.end(), n);
                int mindex = distance(QuarticFC[i].ShortModes.begin(),findm);
                int nindex = distance(QuarticFC[i].ShortModes.begin(),findn);
                if(QuarticFC[i].ShortModes.size()==2 
                        && findpow3 != QuarticFC[i].ModePowers.end()
                        && findm != QuarticFC[i].ShortModes.end()
                        && findn != QuarticFC[i].ShortModes.end() ){ 
                    // V_mnnn and V_mmmn
                    temp.fc += 3*QuarticFC[i].fc;
                }
                if(QuarticFC[i].ShortModes.size()==3 
                        && findm != QuarticFC[i].ShortModes.end()
                        && findn != QuarticFC[i].ShortModes.end() 
                        && QuarticFC[i].ModePowers[mindex] == 1
                        && QuarticFC[i].ModePowers[nindex] == 1
                        ){ 
                    // V_mnll
                    temp.fc += 2*QuarticFC[i].fc;
                }
            }
            if(abs(temp.fc) > 1e-8){
                AnharmHB.push_back(temp);
            }
        }
    } 
    

    //Print settings
    cout << "General settings:" << '\n';
    cout << "  CPU threads: " << Ncpus << '\n' << '\n'; 
    cout << "  Normal modes: " << BasisCount.size() << '\n';
    cout << "  Max quanta per state: " << Ntot << '\n';
    cout << "  Basis functions: " << BasisSet.size() << '\n';
    cout << "  Force constants: ";
    if (Nfc == 0)
    {
        cout << "N/A" << '\n' << '\n';
    }
    else
    {
        cout << AnharmFC.size() << '\n' << '\n';
    }
    if(HCI_Eps != 0.){
        cout.precision(3);
        cout << "  Using HCI with energy cutoff: " << HCI_Eps << '\n';
        if(NEig == 1){
            cout << "  Optimizing just the ground states (ZPE)." << '\n';
        }
        else{
            cout << "  Optimizing the first " << NEig << " eigenstates." << '\n';
        }
    }
    cout << '\n';
    cout.precision(12); //Replace settings
    cout.flush();
    return;
};

void ScaleFC()
{
    //Adjust force constants for permutations and sqrt(2) terms
    for (unsigned int i=0;i<AnharmFC.size();i++)
    {
        //Modify all force constants
        double fcscale;
        //Add sqrt(2) terms
        fcscale = 1.;
        fcscale = ((double)AnharmFC[i].fcpow.size());
        fcscale = pow(2,fcscale);
        fcscale = 1/sqrt(fcscale);
        //Calculate permutations
        for (unsigned int j=0;j<AnharmFC[i].ModePowers.size();j++)
        {
            fcscale /= Fact(AnharmFC[i].ModePowers[j]);
        }
        //Scale FC
        AnharmFC[i].fc *= fcscale;
    }
    return;
};
