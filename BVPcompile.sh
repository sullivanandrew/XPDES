#!/bin/bash
echo off
#
#gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPMIXED.exe
#gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPSCHW.exe
#
#gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPKERR.exe
gcc ./Code/BVP.c ./Code/BVP_main.c ./Code/BVP_GENsys.c ./Code/BVP_GENICBC.c ./Code/BVP_LSsolver.c ./Code/BVP_out.c ./Code/BVP_grid.c ./Funcs/Defn_FEout.c ./Funcs/Defn_ICBCout.c ./Funcs/Defn_JSout.c -lm -lgsl -lgslcblas -o ./Execs/BVPEDGB.exe
echo BVP Compiled.
#
#
