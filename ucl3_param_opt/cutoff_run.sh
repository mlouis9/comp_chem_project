#!/bin/bash

cutoffs="250 300 350 400 450 500 550 600 650 700 750 800 850 900"
cp2k_bin=cp2k.popt
input_file=ucl3.inp
output_file=ucl3.out
no_proc_per_calc=2
no_proc_to_use=16

counter=1
max_parallel_calcs=$((no_proc_to_use / no_proc_per_calc))

for ii in $cutoffs; do
    work_dir=cutoff_${ii}Ry

    if [ ! -d "$work_dir" ]; then
        echo "Directory $work_dir not found!"
        continue
    fi

    (
        cd "$work_dir" || exit 1
        [ -f "$output_file" ] && rm "$output_file"
        mpirun -np $no_proc_per_calc --bind-to none "$cp2k_bin" -o "$output_file" "$input_file"
    ) &

    mod_test=$((counter % max_parallel_calcs))
    if [ "$mod_test" -eq 0 ]; then
        wait
    fi
    counter=$((counter + 1))
done
wait

