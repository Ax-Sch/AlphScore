#!/bin/bash
#Script to concatenate all files that contain alphafold derived features
director_combined_files="data/combine2_protein"
out_dir="data/merge_all"
first_file="$director_combined_files/"$(ls $director_combined_files | head -n1)

mkdir -p $out_dir
mkdir -p $out_dir"/tmp"
zcat $first_file | head -n1 > $out_dir"/tmp/header.csv"
zcat_sting=""

for file in $(ls $director_combined_files)
do
f_size=$(stat --printf="%s" $director_combined_files"/"$file)
if [ $f_size -gt 100 ]
then
zcat_string+=" ${director_combined_files}/${file}"
fi
done

echo $zcat_string
zcat $zcat_string | grep -v "chr,pos(1-based),ref,alt,aaref,aaalt," | sort -t, -T $out_dir"/tmp/" -k2,2 -k3,3n - | cat $out_dir"/tmp/header.csv" - | bgzip -c > $out_dir"/all_possible_values_concat.csv.gz"
