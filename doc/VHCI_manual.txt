[//]: # (Mixture of GitHub markdown and HTML. HTML is needed for formatting.)

***
<div align=center> <h2>
Vibrational Heat-Bath Configuration Interaction
</h2> </div>

<div align=center> <h4> By: Jonathan H. Fetherolf </h4> </div>
<div align=center> <h4> Based on LOVCI by Eric G. Kratz </h4> </div>

***

### Introduction

VHCI is a program to calculate anharmonic vibrational spectra from ab initio potential energy surfaces.  VHCI uses the heat bath selected CI criterion to build an optimized CI basis.  It is a based on code by Eric G. Kratz available at https://github.com/kratman/VibCI.

### Installation

Currently, the binary is not included in the repository. The Makefile can be
used to generate VHCI executable. Only a handful of packages are
required to compile the code. An approximate list of packages is given below.
```
VHCI binary: OpenMP, Eigen3, Spectra, Boost
```

To install VHCI, clone the git repository
```
user:$ mkdir VHCI
user:$ git clone https://github.com/berkelbach-group/VHCI.git ./VHCI/
```

or unpack the zipped source code
```
user:$ mkdir VHCI
user:$ cd VHCI/
user:$ unzip VHCI-master.zip
user:$ mv VHCI-master/* .
user:$ rmdir VHCI-master
```

The source code can be compiled with the Makefile provided with VHCI.
On Ubuntu boxes, the Makefile should function without modifications. However,
with other operating systems it may be necessary to change the path to the required packages.
```
Default: -I/usr/include/eigen3/ -I/usr/include/boost/ -I/usr/include/Spectra
```

The Makefile can produce both the binary.
```
user:$ make install
```

### Performing VHCI calculations

Only a single input file is required to run the VHCI program. All keywords
in the input file must be given in the correct order.

Example input can be found in the tests and doc directories.
```
user:$ vhci -n Ncpus -i input.inp -o freqs.dat -c checkpoint.check

 -n  =>    Number of CPUs used for the calculations
 -i  =>    File name for the VHCI input
 -o  =>    File name for the calculated vibrational frequencies
 -c  =>    File name for the checkpoint file (optional)
```
The VHCI program will output the vibrational eigenvalues of the system specified in 
the input file to the specified output file. It will also output the mode 
identities for the largest element(s) of the CI eigenvectors.

If a checkpoint filename is specified, when the VHCI program executes it will check 
if that file already exists.  If it does exist, it will read in the data and restart 
where the calculation was left off.  On each iteration (each time the VCI 
Hamiltonian is built and diagonalized) it will save the basis, eigenstates and 
eigenvalues.  By default, the checkpoint will only be saved if the basis contains 
more than 10^5 states.  If no checkpoint file is needed, this argument can be 
left blank.

#### Input file format 

Comment: [comment] <br>
HCI_Eps: [eps1] <br>
Nstate: [Nstate] <br>
PT2: [bool] <br>
PT2_Eps: [eps2] <br>
Max_Quanta: [Nmax] <br>
Modes: [Nm] <br>
 [id] [freq] [quanta] <br>
 ...
Force_constants: [Nfc] <br>
 [power] [modes] [value] <br>
 ...

#### Keywords

Comment: This line will not be read in by the VHCI program

HCI_Eps: The variational cutoff parameter for heat bath basis selection, 
also referred to as eps1.  A smaller value of eps1 results in more states 
being added to the variational basis.  If eps1 is set to 0, the program 
will perform conventional (noniterative) CI.

Nstate: The number of states to be computed, starting with the ground 
state.  The program will compute and output the lowest Nstate eigenvalues, 
and the heat bath algorithm will optimize the basis for the first Nstate 
states. 

PT2: A boolean (0 or 1) indicating whether a 2nd-order perturbation 
correction should be performed.

PT2_Eps: The perturbative heat bath cutoff parameter, also referred to as 
eps2.  A smaller value of eps2 results in more states being added to the 
perturbative basis.  If eps2 is set to 0, the program will perform a full 
(unscreened) PT2 correction.

Max_quanta:  Maximum total number of quanta per product state, also 
referred to as Ntot.  This determines how many excitations are allowed per 
basis state on the first iteration of heat bath VCI, or how many 
configurations are included in conventional VCI.

Modes: The keyword needs the total number of active vibrational modes (Nm),
followed by the properties of each mode. Every mode is given an integer ID 
(from 0 to Nm-1), a frequency (cm^-1), the maximum number of quanta allowed 
on that mode in the basis set.

Force_constants: The (Nfc) anharmonic force constants are defined by the
total power of the modes (quadratic=2,cubic=3,quartic=4,etc), the list of
integer IDs for the modes, and the value of the force constant. The list of
modes can be given in any order and the force constants can include any
positive power, but heat bath CI is only implemented up to 6th-order. 
The values of the force constants (in cm^-1) are the same as those 
given in the output of common QM packages (e.g. Gaussian, GAMESS).
