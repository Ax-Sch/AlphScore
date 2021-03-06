#!/bin/bash

#Script to concatenate all files that contain alphafold derived features
director_combined_files=$1
out_dir=$2
first_file="$director_combined_files/"$(ls $director_combined_files | head -n1)

echo "using header of: $first_file"

mkdir -p $out_dir
mkdir -p $out_dir"/tmp"
zcat $first_file | head -n1 > $out_dir"/header.csv"
zcat_sting=""

# loop over directory content and extract collect non empty files with TRUE in filename 
for file in $(ls $director_combined_files)
do
if [[ $file == *"TRUE"* ]]; then
f_size=$(stat --printf="%s" $director_combined_files"/"$file)
if [ $f_size -gt 100 ]
then
zcat_string+=" ${director_combined_files}/${file}"
fi
fi
done

#relevant_cols=$(zcat results/validation_set/validation_set_w_AlphScore.csv.gz | awk -v RS=',' '/^b_factor$/{print NR} /^SOLVENT_ACCESSIBILITY_core$/{print NR} /^Uniprot_acc_split$/{print NR} /^gnomad_train$/{print NR} /AlphScore/{print NR; exit}' | tr "\n" ",")"2,3,4,5,6,7,8,9,10"


echo $zcat_string

cat $out_dir"/header.csv" > $out_dir"/header_short.csv"

zcat $zcat_string  | grep -v "chr,pos(1-based)," | sort -V -t, -T $out_dir"/tmp/" -k1 -k2 - | \
cat $out_dir"/header_short.csv" - | tr "," "\t" | bgzip -c > $out_dir"/all_possible_values_concat.csv.gz"

tabix -s 1 -b 2 -e 2 $out_dir"/all_possible_values_concat.csv.gz"

rm -rf $out_dir"/tmp"

#zcat $zcat_string | grep -v "chr,pos(1-based),ref,alt,aaref,aaalt," | sort -t, -T $out_dir"/tmp/" -k2,2 -k3,3n - | cat $out_dir"/header.csv" - | bgzip -c > $out_dir"/all_possible_values_concat.csv.gz"


