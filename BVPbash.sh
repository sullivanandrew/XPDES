echo off
echo Watch me delete this shit...
rm -r "./Data/BVPout/BVPout_iterations"
rm -r "./Data/BVPout/BVPout_sols"
rm -r "./Data/BVPout/BVPout_sys"
mkdir "./Data/BVPout/BVPout_iterations"
mkdir "./Data/BVPout/BVPout_sols"
mkdir "./Data/BVPout/BVPout_sys"
echo Successfully deleted that shit.
#valgrind --leak-check=full --track-origins=yes
#valgrind --tool=massif --time-unit=B
#ms_print massif.out.#####
#
#./Execs/BVPMIXED.exe 2 10 10 4 -3 1 0 0.0000 0.0000 1.0 0.0
#./Execs/BVPSCHW.exe 4 31 26 10 -3 1 0 0.0000 0.0100 0.25 0.0
#./Execs/BVPKERR.exe 8 31 26 12 -3 1 0 0.0000 0.0005 0.25 0.1
#valgrind --tool=massif --time-unit=B
#valgrind --leak-check=full --track-origins=yes
#./Execs/BVPEDGB.exe 10 11 11 4 -1 0 0 0.0000 0.0001 0.25 0.0
./Execs/BVPEDGB.exe 10 21 21 10 -1 0 0 0.0000 0.0001 0.25 0.0
#
#
#
echo BVP Executed.
