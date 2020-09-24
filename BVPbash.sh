echo off
echo Watch me delete Data directories
#rm -r "./Data/BVPoutLIN/BVPout_iterations"
#rm -r "./Data/BVPoutLIN/BVPout_sols"
#rm -r "./Data/BVPoutLIN/BVPout_sys"
#mkdir "./Data/BVPoutLIN/BVPout_iterations"
#mkdir "./Data/BVPoutLIN/BVPout_sols"
#mkdir "./Data/BVPoutLIN/BVPout_sys"
#rm -r "./Data/BVPoutEXP/BVPout_iterations"
#rm -r "./Data/BVPoutEXP/BVPout_sols"
#rm -r "./Data/BVPoutEXP/BVPout_sys"
#mkdir "./Data/BVPoutEXP/BVPout_iterations"
#mkdir "./Data/BVPoutEXP/BVPout_sols"
#mkdir "./Data/BVPoutEXP/BVPout_sys"
#
#rm -r "./Data/BVPoutKERR/BVPout_iterations"
#rm -r "./Data/BVPoutKERR/BVPout_sols"
#rm -r "./Data/BVPoutKERR/BVPout_sys"
#mkdir "./Data/BVPoutKERR/BVPout_iterations"
#mkdir "./Data/BVPoutKERR/BVPout_sols"
#mkdir "./Data/BVPoutKERR/BVPout_sys"
echo Successfully deleted Data directories.
#
#./Execs/BVPKERR.exe 10 31 31 16 -5 1.0 0.000 0.9 -0.001 0.1 1
#./Execs/BVPKERR.exe 10 31 31 16 -5 1.0 0.000 0.0 -0.001 0.001 0
#
./Execs/BVPEDGB.exe 10 31 21 12 -5 1.0 0.100 0.1 -0.001 0.0 0
#
echo BVP Executed.
