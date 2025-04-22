#!/bin/bash
 
cutoffs="250 300 350 400 450 500 550 600 650 700 750 800 850 900"
 
basis_file=BASIS_SET
potential_file=GTH_POTENTIALS
template_file=template.inp
input_file=ucl3.inp
 
rel_cutoff=60
 
for ii in $cutoffs ; do
    work_dir=cutoff_${ii}Ry
    if [ ! -d $work_dir ] ; then
        mkdir $work_dir
    else
        rm -r $work_dir/*
    fi
    sed -e "s/LT_rel_cutoff/${rel_cutoff}/g" \
        -e "s/LT_cutoff/${ii}/g" \
        $template_file > $work_dir/$input_file
    cp $basis_file $work_dir
    cp $potential_file $work_dir
done
