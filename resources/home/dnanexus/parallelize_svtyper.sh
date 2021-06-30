#!/bin/bash

input=$1
directory=$2
output=$3
input_bam=$4
max_reads=$5

input_lines=$(grep -v \# $input | wc -l)
echo "Number of calls to process: $input_lines"
if [[ $input_lines -eq 0 ]]; then
    echo "No calls. Moving on."
    exit 0
fi

procunits=$(nproc)

# We need to add IDs to the variants. svtyper-sso required that
IFS=$'\t' #Split by tab

while read v; do
    if [[ $v == "#"* ]]; then
        echo "$v" >> "$directory/tmp.vcf"
    else
        tmp=($v)  #The array assignment tmp=($line) splits the value on whatever characters IFS contains
        tmp[2]="${tmp[0]}_${tmp[1]}_r$RANDOM" # IDs of the form chr1_392433_r242, should be unique enough
        line="${tmp[*]}"
        echo "$line" >>  "$directory/tmp.vcf"
    fi
done < $input

svtyper-sso --core $procunits --batch_size 500 --max_reads "$max_reads" -i "$directory/tmp.vcf" -B $input_bam > "$directory/out.vcf" || exit 1

echo "Combine GTs"
python /home/dnanexus/combine_genotypes.py $input "$directory/out.vcf" > "$directory/out_gt_combined.vcf" || exit 1

# Assemble everything. We add the original input to the ouput (some entries already might have a GT)
grep \# $input > $output
#grep -v \# $input >> $output # Don't add the original

# Remove IDs again
grep -v \# "$directory/out_gt_combined.vcf" | while read v; do
    tmp=($v)  #The array assignment tmp=($line) splits the value on whatever characters IFS contains
    tmp[2]="."
    line="${tmp[*]}"
    echo "$line" >>  $output
done

## THE FOLLOWING IS USING THE ORIGINAL SVTYPER IMPLEMENTATION

# procunits=$(nproc)
# threads=$(nproc)
# threads=$(expr $threads \* 4) 

# if [[ $input_lines -ge $threads ]]; then
#     lines=$(expr $input_lines / $threads)
#     split -d -a 5 -l $lines $input $directory
# else
#     split -d -a 5 -l $input_lines $input $directory
# fi

# echo "Processing units: $procunits"
# echo "Number of chunks to process: $threads"


# i=0
# for item in $directory*; do
#     i=$(expr $i + 1)
#     grep \# $input > $directory/$i
#     grep -v \# $item >> $directory/$i
#     #echo "svtyper -B $input_bam -i $directory/$i >> $directory/$i" >> $output.cmds

#     num_variants_in_chunk=$(grep -v \# $directory/$i | wc -l)
#     if [[ $num_variants_in_chunk -gt 0 ]]; then
#         echo "$num_variants_in_chunk calls in chunk $i."
#         echo "$i" >> $output.cmds
#         cat $directory/$i
#         echo "---------------------------"
#     else 
#         echo "No calls in chunk $i. Not sending this one to SVTyper."
#     fi
# done

# # We don't have the memfree option is the Ubuntu 14.04 version  of parallel
# #parallel --memfree 5G --retries 2 --verbose -a $output.cmds eval 2> /dev/null
# #parallel --retries 1 --verbose -a $output.cmds eval 2> /dev/null

# cat $output.cmds
# echo "---------------------------"

# command="echo \"Running SVtyper on chunk {}\"; svtyper -B $input_bam -i $directory/{} >> $directory/{} || exit 1"
# cat $output.cmds | xargs -P $procunits -i bash -c "$command" || exit 1


# grep \# $input > $output
# for item in $directory/*; do
#     grep -v \# $item >> $output
# done
