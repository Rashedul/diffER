#!/bin/bash
if [ "$1" == "-h" ] ; then
        echo -e "A script that creates an sbatch file for a two group comparision

Usage: `basename $0` -t <sample_type1> -T <sample_type2>  -i <input_tables_directory> -o <out_directory> 
                      -p <p_value_threshold> -m <mark> -g <groups_column>
        <sample_type1>: sample_type 1
        <sample_type2>: sample_type 2
        <input_tables_directory>: The path to the input final tables directory
        <out_directory>: The path to the output directory where the table will be written
        <p_value_threshold>: the p-value to filter final table
        <mark>: The mark to be compared
	<groups_column> The name of the column in the metdata table that contains sample_type1 and sample_type2"

exit 0
fi

input_tables_directory=/projects/blueprint_CLL_dart2/jsteif/symlinked_final_tables

while [ $# -gt 0 ]
do
        case "$1" in
        -t) sample_type1="$2"; shift;;
        -T) sample_type2="$2"; shift;;
        -i) input_tables_directory="$2"; shift;;
        -o) out_directory="$2"; shift;;
        -p) p_value_threshold="$2";shift;;
        -m) mark="$2";shift;;
	-g) groups_column="$2";shift;;
        esac
        shift
done

PATH=/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin:$PATH
export PATH=/projects/epigenomics3/epigenomics3_results/users/esu/python3.7/bin/:$PATH

echo "#!/bin/bash
#SBATCH --partition=upgrade
#SBATCH -D $out_directory
#SBATCH --mem 370000" > $out_directory/sbatch_"$mark"_"$sample_type1"_"$sample_type2"_"$p_value_threshold".sh

declare -a chromosomes=( "chr22" "chr21" "chr20" "chr19" "chr18" "chr17" "chr16" "chr15" "chr14" "chr13" "chr12" "chr11" "chr10" "chr9" "chr8" "chr7" "chr6" "chr5" "chr4" "chr3" "chr2" "chr1" )

for c in "${chromosomes[@]}"; do
echo "bash /projects/blueprint_CLL_dart2/jsteif/scripts/run_opt_two_group_comparision.sh -t $sample_type1 -T $sample_type2 -p $p_value_threshold -c $c -m $mark -i $input_tables_directory -g $groups_column -o $out_directory" >> $out_directory/sbatch_"$mark"_"$sample_type1"_"$sample_type2"_"$p_value_threshold".sh ; done
