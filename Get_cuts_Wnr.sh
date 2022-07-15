#grep "^12 " Wannier_functions_band0.txt



N=13
file_tag=_slices_${N}


#Cut (0,N) to (N,0)
rm Cut_2${file_tag}.txt
n1_=${N}
n2_=-1
for i in {0..12}
do
n1_=$(echo "${n1_}-1" | bc -l)
n2_=$(echo "${n2_}+1" | bc -l)

n_ind=$(echo "${n1_} + ${N}*${n2_}" | bc -l)
#echo "${n_ind}"
line=$(grep "^${n_ind} " Wannier_functions_band0.txt)

echo "${line}" >> Cut_2${file_tag}.txt
done

#Cut (0,(N-1)/2) to (N-1, (N-1)/2)
rm Cut_1${file_tag}.txt
n1_=-1
n2_=$(echo "(${N}-1)/2" | bc -l)
n2_=$(printf "%1.0f" ${n2_})
for i in {0..12}
do
n1_=$(echo "${n1_}+1" | bc -l)
n_ind=$(echo "${n1_} + ${N}*${n2_}" | bc -l)

#echo "${n1_} ${n2_}  ${n_ind}"
line=$(grep "^${n_ind} " Wannier_functions_band0.txt)

echo "${line}" >> Cut_1${file_tag}.txt
done

#Cut ((N-1)/2,0) to ((N-1)/2, N-1)
rm Cut_3${file_tag}.txt
n1_=$(echo "(${N}-1)/2" | bc -l)
n1_=$(printf "%1.0f" ${n1_})
n2_=-1
for i in {0..12}
do
n2_=$(echo "${n2_}+1" | bc -l)
n_ind=$(echo "${n1_} + ${N}*${n2_}" | bc -l)

line=$(grep "^${n_ind} " Wannier_functions_band0.txt)

echo "${line}" >> Cut_3${file_tag}.txt
done



