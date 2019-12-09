#!/bin/bash
echo off
nfields=10
#
maxderiv=2
echo Number of fields  : $nfields
echo Maximum derivative: $maxderiv
gcc ./Gens/Export_MapleGen.c -o ./Execs/mgen.exe
echo Maple generator compiled.
#valgrind --leak-check=full --track-origins=yes
./Execs/mgen.exe $nfields $maxderiv
echo Maple generator executed.
#
echo Run Maple file now. Then
read -p "Press [Enter] key to continue compilation..."
#
gcc ./Gens/Export_NNZGen.c -o ./Execs/nnzgen.exe
echo NNZ generator compiled.
#valgrind --leak-check=full --track-origins=yes
./Execs/nnzgen.exe $nfields $maxderiv
echo NNZ generator executed.
#
gcc ./Gens/Export_CGen.c ./Funcs/Defn_JSout.c -o ./Execs/cgen.exe
echo C generator compiled.
#valgrind --leak-check=full --track-origins=yes
./Execs/cgen.exe $nfields $maxderiv
echo C generator executed.
#
#gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPMIXED.exe
#gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPSCHW.exe
#
#gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPKERR.exe
gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPEDGB.exe
echo BVP Compiled.
#
#
