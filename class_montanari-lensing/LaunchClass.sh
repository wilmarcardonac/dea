#! /bin/bash
# This script read the explanatory.ini file of CLASS code
# and replace a line corresponding to a parameter (e.g., Omega_m) with another
# one with a new vlue of that parameter.
# Then it runs class.
# This is repeated iteratively.

#Change input explanatory file
fin="explanatory_lens.ini"
#Default name of the output file according to CLASS
fout="output/test_cl.dat"
#Additional label
label="Nbin5_drL"

#Choose the values of the parameter, and change its name
#NB: case sensitive, this replace ALL line containing the string to be replaced
for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645 
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in  67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done  
          done
        done
      done
    done
  done
done


for i in 0.02193 0.02209 0.02241 0.02257
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645 
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in  67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1168 0.1183 0.1213 0.1228
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in  67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9547 0.9596 0.9694 0.9743
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in  67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.1757e-9 2.19111e-9 2.22193e-9 2.23734e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in  67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in  65.95 66.61 67.273 68.59
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in 67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.050 0.055 0.065 0.070
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  1.00
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

for i in 0.02225
do
  sed -e 's/^omega_b.*/omega_b='$i'/'  $fin > temp1
  for j in 0.1198
  do
    sed -e 's/^omega_cdm.*/omega_cdm='$j'/' temp1 > temp2
    for k in 0.9645
    do
      sed -e 's/^n_s.*/n_s='$k'/' temp2 > temp3
      for l in 2.20652e-9
      do
        sed -e 's/^A_s.*/A_s='$l'/' temp3 > temp4
        for m in 67.27
        do
          sed -e 's/^H0.*/H0='$m'/' temp4 > temp5
          for n in  0.060
          do
            sed -e 's/^m_ncdm.*/m_ncdm='$n'/' temp5 > temp6
            for o in  0.950 0.975 1.025 1.050
            do
              sed -e 's/^lensing convergence rescale.*/lensing convergence rescale='$o'/' temp6 > temp
              echo Computing omega_b=$i , omega_cdm=$j , ns=$k  and As=$l ,H0=$m, m_ncdm=$n, kappa_r=$o
              mv temp $fin \
              && ./class $fin \
              && cp $fout output/cl_omB_$i.omCDM_$j.ns_$k.As_$l.H0_$m.mncdm_$n.kappar_$o.$label.dat
            done
          done
        done
      done
    done
  done
done

echo "End of computation"
