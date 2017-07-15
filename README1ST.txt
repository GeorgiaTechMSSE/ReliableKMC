===================================================
Reiliable Kinetic Monte Carlo Simulaiton Package
===================================================
Reliable kinetic Monte Carlo (rKMC) is an extension of SPPARKS, which can simulate kinetics of chemical reaction and phase transition when the kinetic rates as model input are not precisely known. That is, model-form uncertainty is present. Tranditional approach to assess the input uncertainty effect is sensitivity analysis, where a slightly different input parameters are applied for each simulation run. Many runs of simulation are needed to evaluate the impact on the output. In rKMC, interval-valued rates are used in simulation so that the sensitivity can be assessed with only one single simulation run.

by Prof. Yan Wang
Georgia Institute of Technology, USA

The package is under GNU GENERAL PUBLIC LICENSE and AS IT IS. 
===================================================
To cite:

Wang, Y. (2013) Reliable kinetic Monte Carlo simulation based on random set sampling. Soft Computing, 17(8), pp.1439-1451.

===================================================




===================================================
Controlled Kinetic Monte Carlo Simulaiton Package
===================================================
cKMC is an extension of SPPARKS, which can simulate both top-down (e.g. AFM lithography, 
focused ion beam lithography, etc.) and bottom-up (e.g. chemical vapor deposition, 
physical vapor deposition, etc.) nanomanufacturing processes.

by Prof. Yan Wang
Georgia Institute of Technology, USA

The package is under GNU GENERAL PUBLIC LICENSE and AS IT IS. 
===================================================
To cite:

Wang, Y. (2011) "Controlled kinetic Monte Carlo simulation of nanomanufacturing processes," 
2011 ASME International Design Engineering Technical Conferences & The Computer and 
Information in Engineering Conference (IDETC/CIE2011), Aug. 28-31, 2011, 
Washington, DC, Paper No.DETC2011-48570 

Wang, Y. (2016) Controlled kinetic Monte Carlo simulation for computer-aided nanomanufacturing. Journal of Micro and Nano-manufacturing, 4(1), 011001.

===================================================
INSTALLATION
----------------------------------------------
A. CYGWIN - Windows environment

1. download and install cygwin from http://cygwin.com/install.html
   make sure the following packages are installed (by selecting the install options , 
   not "skip") during the cygwin setup. Note that you can always change options by running 
   "setup.exe" or "setup_x86.exe" under cygwin root directory.
   
   1) bash (under category "Shells")
   2) csh/tcsh (under category "Shells")
   3) make (under category "Devel")
   4) gcc-g++ (under category "Devel")
   5) diffutils (under category "Utils")
 
   
2. Compile
   2a. start a cygwin Terminal, go to directory "spparks-30Apr10/src/STUBS", 
       for instance, if the package was unzipped at "c:\user\spparks-30Apr10",
       use command "cd /cygdrive/c/usr/spparks-30Apr10/src/STUBS" to go to the directory

   2b. compile the stub using command: "make"
   
       Make sure that "libmpi.a" is generated under the directory "spparks-30Apr10/src/STUBS" before next step. 
       
   2c. go to directory "spparks-30Apr10/src" using command "cd ..", 
       compile the package using command: "make cygwin",

   If successfully compiled, an executable file "spk_cygwin.exe" will be generated 
   under the directory "spparks-30Apr10/src".


-----------------------------------------------
B. Linux

   compile a serial version
   ========================
   a. go to directory "spparks-30Apr10/src/STUBS", type command: "make"
   
   Make sure that "libmpi.a" is generated under the directory "spparks-30Apr10/src/STUBS" before next step. 
       
   b. go to directory "spparks-30Apr10/src", type command: "make serial",

   If successfully compiled, an executable file "spk_serial" will be generated 
   under the directory "spparks-30Apr10/src".   


   compile a MPI parallel version
   ========================
   a. go to directory "spparks-30Apr10/src/STUBS", type command: "make"
   
   Make sure that "libmpi.a" is generated under the directory "spparks-30Apr10/src/STUBS" before next step. 
       
   b. go to directory "spparks-30Apr10/src", type command: "make linux",

   If successfully compiled, an executable file "spk_linux" will be generated 
   under the directory "spparks-30Apr10/src".   
   
===================================================
TO RUN 
----------------------------------------------
A. CYGWIN - Windows environment

go to directory of input files, such as "spparks-30Apr10/nanomfg_examples/"
type command "../src/spk_cygwin.exe < in.AFM" to run the AFM lithography simulation.

Other examples can be run similarly. Just need to write the input files properly. Please refer to
the following paper for the details of special commands of cKMC, in addition to 
the standard SPPARKS commands. 
   
Wang, Y. (2011) "Controlled kinetic Monte Carlo simulation of nanomanufacturing processes," 
2011 ASME International Design Engineering Technical Conferences & The Computer and 
Information in Engineering Conference (IDETC/CIE2011), Aug. 28-31, 2011, 
Washington, DC, Paper No.DETC2011-48570  
   
   
   
=====================================================
If there is any question, contact: yan.wang@me.gatech.edu
=====================================================

