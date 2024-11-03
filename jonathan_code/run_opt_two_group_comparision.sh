#!/bin/bash
if [ "$1" == "-h" ] ; then
        echo -e "A script that does a simulation to compare two sample typs

Usage: `basename $0` -t <sample_type1> -T <sample_type2>  -i <input_tables_directory> -o <out_directory> -c <chromosome> 
		      -p <p_value_threshold> -m <mark> -g <groups_column> 
	<sample_type1>: sample_type 1 
	<sample_type2>: sample_type 2
	<input_tables_directory>: The path to the input final tables directory
 	<chromosome>: The chromosome (e.g chr10)
	<out_directory>: The path to the output directory where the table will be written
	<p_value_threshold>: the p-value to filter final table
	<mark>: The mark to be compared
	<groups_column> The name of the column in the metdata table that contains sample_type1 and sample_type2"	

exit 0
fi

input_tables_directory=/projects/epigenomics_assembly/jsteif/IHEC/postgresql_tables/ALL/final_tables

while [ $# -gt 0 ]
do
        case "$1" in
        -t) sample_type1="$2"; shift;;
	-T) sample_type2="$2"; shift;;
        -i) input_tables_directory="$2"; shift;;
	-o) out_directory="$2"; shift;;
	-c) chromosome="$2";shift;;
	-p) p_value_threshold="$2";shift;;
	-m) mark="$2";shift;;
	-g) groups_column="$2";shift;;
        esac
        shift
done


PATH=/gsc/software/linux-x86_64-centos7/R-3.6.3/bin:$PATH
#export PATH=/projects/epigenomics3/epigenomics3_results/users/esu/python3.7/bin:$PATH

Rscript /projects/blueprint_CLL_dart2/jsteif/scripts/opt_two_group_comparision.R -t $sample_type1 -T $sample_type2 -i $input_tables_directory -o $out_directory -c $chromosome -p $p_value_threshold -m $mark -g $groups_column

