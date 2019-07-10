# rates_calcom
Software to estimate inclusive and hadron rates for calibration and commissioning of CLAS12


******************************************************

Producing lund files with Inclusive Electron Events

	  Raffaella De Vita, Silvia Pisano

******************************************************


Aim of the software is to produce events in the LUND format containing inclusive electrons in an extended kinematics, ranging from elastic to DIS regime.

The procedure consists in two steps:

1) First, it makes use of the Misak code (in Fortran), that calculates the inclusive cross-section at given values of (theta, E) of the outgoing electron.

It produces a table that maps the cross-section as a function of (theta, E). It is then used in the next step by the generator to produce inclusive events accordingly.

2) The generator itself, that reads the table produced at the previous point and that generate inclusive events.

------------------------------------------------------------------------------------------------

STEP 1: PRODUCE A TABLE WITH THE VALUE OF THE CROSS-SECTION IN A GIVEN (THETA, E) RANGE

------------------------------------------------------------------------------------------------

Misak software is in the directory

cross_section_calculation

It can be set through the file init_incl.dat, that contains the following text:


INCIDENT   ELECTRON    6.000   0.000   0.000   0.000   0.000   0.000      
TARGET     H2          1.000   1.000   0.000   0.100   0.000   0.000
RAD_EFFECT YES         5.000   0.020   0.010   0.010   0.050   0.000  

SWELLING   V1          0.000   0.000   0.000   0.000   0.000   0.000  
EMC        NES         0.000   0.000   0.000   0.000   0.000   0.000  

ELEC_SPECT             0.000   0.000   0.000   0.000   0.000   0.000
Ee` -RANGE YES         0.500   6.000   0.005   0.000   0.000   0.000
THe -RANGE             5.000  35.000   0.500  15.000   0.000   0.000
Q0  -RANGE NES         0.020   0.100   0.001   0.000   0.000   0.000
W   -RANGE NO          0.900   0.910   0.025   0.000   0.000   0.000
X   -RANGE NES         0.300   1.890   0.050   0.000   0.000   0.000

The relevant elements are the INCIDENT ELECTRON ENERGY (6 GeV in the example) and the ranges in (theta, E): 

Ee` -RANGE YES         0.500   6.000   0.005   0.000   0.000   0.000
THe -RANGE             5.000  35.000   0.500  15.000   0.000   0.000

to be set accordingly to the configuration to be analyzed.

The code can be run through the command

./inclusive

and produces as output the file

fort.26

Its columns represent

Theta  E'  Sigma/dE/dOmega   Sigma_radiated/dE/dOmega   Q2   W Sigma/dW/dQ2  Sigma_radiated/dW/dQ2


Tables for the standard 6 and 11 GeV kinematics have been already produced. They are stored in the directory

cross_section_data/

File for 6 and 11 GeV are, respectively:

data_5_35_0.5_6.0_6.0GeV.dat

data_5_35_0.5_10.9_11.0GeV.dat


---------------------------------------------------------------

STEP 2: GENERATE INCLUSIVE ELECTRON EVENTS IN THE LUND FORMAT

---------------------------------------------------------------

In order to generate a set of inclusive events stored in the LUND format, run the script

generate_lund_events.C

as

root -b -q generate_lund_events.C

It contains the code

{

  gSystem->CompileMacro("DisFunctions.C", "kf");
  gSystem->CompileMacro("GenFunctions.C", "kf");
  gSystem->CompileMacro("InclusiveRateAnalysis.C");
  t=new InclusiveRateAnalysis();
  /* 11 GeV beam energy */
  t->SetKinematicLimits( 5, 35, 0.5, 0.5, 10.9, 0.005, 75674, 11.0 );
  t->SetFile("cross_section_data/data_5_35_0.5_10.9_11.0GeV.dat");
  /* --------- */
  /* 6 GeV beam energy */
  // t->SetKinematicLimits( 5, 35, 0.5, 0.5, 6.0, 0.005, 46872, 6 );
  // t->SetFile("cross_section_data/data_5_35_0.5_6.0_6.0GeV.dat");
  /* --------- */
  t->PrepareOutput();
  t->BookHistos();
  t->SetLuminosity( pow( 10, 35 ) );
  t->ReadDataFile();
  t->PlotCrossSectionDependences();
  t->GeneratePseudoData( 100000, 1000, "lund_files");
  t->CompareGenEventsToOriginal();
  t->SaveOutput();
  delete t;

}


In order to produce a set of LUND files, one has to set

1) kinematic limits of the explored range. They are set through the method

 t->SetKinematicLimits( 5, 35, 0.5, 0.5, 10.9, 0.005, 75674, 11.0 );
 
where different arguments correspond to

theta_min, theta_max, theta_step, energy_min, energy_max, energy_step, _data_file_n_lines, beam_energy

More in details,

_data_file_n_lines

corresponds to the number of lines of the data files ( e.g. cross_section_data/data_5_35_0.5_10.9_11.0GeV.dat) that gives the number of kinematics points used to map the cross section dependence.

n. b. It can be obtained from the shell through

wc -l _data_file_name

e.g.:

[pisanos@jlabmn software_to_git]$ wc -l cross_section_data/data_5_35_0.5_6.0_6.0GeV.dat 

46872 cross_section_data/data_5_35_0.5_6.0_6.0GeV.dat

[pisanos@jlabmn software_to_git]$ 


2) File name containing the cross-section values generated in the previous step:

  t->SetFile("cross_section_data/data_5_35_0.5_10.9_11.0GeV.dat");


3) Default output rootfile - with some monitoring histograms - is named, in default, as the data file name, where the extension ".dat" is replaced by ".root".

t->PrepareOutput();

It will produce the file

data_5_35_0.5_10.9_11.0GeV.root

A further argument can be passed to add suffix before the .root extension:

 t->PrepareOutput("_test");

for example, will produce the file

data_5_35_0.5_10.9_11.0GeV_test.root


4) Luminosity for rate estimate can be set through

 t->SetLuminosity( pow( 10, 35 ) );

that, in this case, corresponds to 10^35 cm^-2 s^-1.


5) The method

  t->GeneratePseudoData( 100000, 1000, "lund_files");

is the actual generator. Arguments are

n1 ( = 100000 in the example): number of inclusive events to be generated

n2 ( = 1000   in the example): number of events to store per lund files

n3 ( = "lund_files"         ): directory for storing the lund files - if any.


6) Finally, the (optional) method

  t->CompareGenEventsToOriginal();

adds in the output rootfile some monitoring histograms as the comparison of the generated distributions to the original ones.



p.s. To use the pre-existing cross-section tables, the settings

  /* 11 GeV beam energy */
  t->SetKinematicLimits( 5, 35, 0.5, 0.5, 10.9, 0.005, 75674, 11.0 );
  t->SetFile("cross_section_data/data_5_35_0.5_10.9_11.0GeV.dat");
  /* --------- */
  /* 6 GeV beam energy */
  // t->SetKinematicLimits( 5, 35, 0.5, 0.5, 6.0, 0.005, 46872, 6 );
  // t->SetFile("cross_section_data/data_5_35_0.5_6.0_6.0GeV.dat");
  /* --------- */


can be applied, for 11 or 6 GeV respectively.
