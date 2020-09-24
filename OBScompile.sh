echo off
#To compile
#gcc ./Code/BVP_ObsExt.c ./Modified-Gravity-Files/Linear-Scalar-Gauss-Bonnet-Files/Replace-in-Funcs/Defn_PHYSout.c ./Modified-Gravity-Files/Linear-Scalar-Gauss-Bonnet-Files/Replace-in-Funcs/Defn_ICBCout.c -lm -o ./Execs/BVP_ObsExt.exe
#gcc ./Code/BVP_ObsExt_Pert.c ./Modified-Gravity-Files/Perturbative-sGB-Files/Replace-in-Funcs/Defn_PHYSout.c ./Modified-Gravity-Files/Perturbative-sGB-Files/Replace-in-Funcs/Defn_ICBCout.c ./Funcs/Defn_KERRout.c -lm -o ./Execs/BVP_ObsExt_Pert.exe
#Modified gravity Files location
gcc ./Code/BVP_ObsExt.c ./Funcs/Defn_ISCOLRout.c ./Modified-Gravity-Files/Linear-Scalar-Gauss-Bonnet-Files/Replace-in-Funcs/Defn_ICBCout.c ./Funcs/Defn_KERRout.c -lm -o ./Execs/BVP_ObsExt.exe
gcc ./Code/BVP_ObsExt_Pert.c ./Funcs/Defn_ISCOLRout.c ./Modified-Gravity-Files/Perturbative-sGB-Files/Replace-in-Funcs/Defn_ICBCout.c ./Funcs/Defn_KERRout.c -lm -o ./Execs/BVP_ObsExt_Pert.exe
#
echo OBS compiled.
