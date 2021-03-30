#!/usr/bin/env bash

illumina_bam=$1 #The name of the Illumina BAM file for which to call structural variants containing mapped reads.
illumina_bai=$2 #The name of the corresponding index for the Illumina BAM file.
ref_fasta=$3 #The name of the reference file that matches the reference used to map the Illumina inputs.
ref_index=$4 #The name of the corresponding index for the reference genome file.
prefix="" #(Optional) If provided, all output files will start with this. If absent, the base of the BAM file name will be used. CWL need to be adjusted
filter_short_contigs=$5 #If selected, SV calls will not be generated on contigs shorter than 1 MB.
run_breakdancer=$6 #If selected, the program Breakdancer will be one of the SV callers run.
run_breakseq=$7 #If selected, the program BreakSeq2 will be one of the SV callers run.
run_manta=$8 #If selected, the program Manta will be one of the SV callers run.
run_cnvnator=${9} #If selected, the program CNVnator will be one of the SV callers run.
run_lumpy=${10} #If selected, the program Lumpy will be one of the SV callers run.
run_delly_deletion=${11} #If selected, the deletion module of the program Delly2  will be one of the SV callers run.
run_delly_insertion=${12} #If selected, the insertion module of the program Delly2  will be one of the SV callers run.
run_delly_inversion=${13} #If selected, the inversion module of the program Delly2  will be one of the SV callers run.
run_delly_duplication=${14} #If selected, the duplication module of the program Delly2  will be one of the SV callers run.
run_genotype_candidates=${15} #If selected, candidate events determined from the individual callers will be genotyped and merged to create a consensus output.
run_svviz="False"
svviz_only_validated_candidates="False"
dnanexus="False"

echo "illumina_bam: $illumina_bam"  >> report.txt
echo "illumina_bai: $illumina_bai"  >> report.txt
echo "ref_fasta: $ref_fasta"  >> report.txt
echo "ref_index: $ref_index"  >> report.txt
echo "filter_short_contigs: $filter_short_contigs"  >> report.txt
echo "run_breakdancer: $run_breakdancer"  >> report.txt
echo "run_breakseq: $run_breakseq"  >> report.txt
echo "run_manta: $run_manta"  >> report.txt
echo "run_cnvnator: $run_cnvnator"  >> report.txt
echo "run_lumpy: $run_lumpy"  >> report.txt
echo "run_delly_deletion: $run_delly_deletion"  >> report.txt
echo "run_delly_insertion: $run_delly_insertion"  >> report.txt
echo "run_delly_inversion: $run_delly_inversion"  >> report.txt
echo "run_delly_duplication: $run_delly_duplication"  >> report.txt
echo "run_genotype_candidates: $run_genotype_candidates"  >> report.txt

# Create output files
echo "breakseq_stdout" > breakseq_stdout.txt
echo "breakseq_err" > breakseq_err.txt
echo "breakdancer_err" > breakdancer_err.txt
echo "manta_stdout" > manta_stdout.txt
echo "manta_err" > manta_err.txt
echo "cnvnator_stdout" > cnvnator_stdout.txt
echo "cnvnator_err" > cnvnator_err.txt
echo "sambamba_stdout" > sambamba_stdout.txt
echo "sambamba_err" > sambamba_err.txt
echo "delly_del_stdout" > delly_del_stdout.txt
echo "delly_del_err" > delly_del_err.txt
echo "delly_inv_stdout" > delly_inv_stdout.txt
echo "delly_inv_err" > delly_inv_err.txt
echo "delly_dup_stdout" > delly_dup_stdout.txt
echo "delly_dup_err" > delly_dup_err.txt
echo "delly_ins_stdout" > delly_ins_stdout.txt
echo "delly_ins_err" > delly_ins_err.txt
echo "lumpy_stdout" > lumpy_stdout.txt
echo "lumpy_err" > lumpy_err.txt
touch result_lumpy.vcf result_manta.vcf result_breakdancer.vcf result_manta_alignment_stats_summary.txt result_cnvnator.vcf result_cnvnator.output result_delly_deletion.vcf result_delly_inversion.vcf result_delly_duplication.vcf result_delly_insertion.vcf result_breakseq.vcf result_breakdancer.svtyped.vcf result_breakseq.svtyped.vcf result_cnvnator.svtyped.vcf result_delly.svtyped.vcf result_lumpy.svtyped.vcf result_manta.svtyped.vcf result_survivor_sorted.vcf result_survivor_combined_genotyped.vcf



check_threads(){
    breakdancer_threads=$(top -n 1 -b -d 10 | grep -c breakdancer)
    cnvnator_threads=$(top -n 1 -b -d 10 | grep -c cnvnator)
    sambamba_processes=$(top -n 1 -b -d 10 | grep -c sambamba)
    manta_processes=$(top -n 1 -b -d 10 | grep -c manta)
    breakseq_processes=$(top -n 1 -b -d 10 | grep -c breakseq)
    delly_processes=$(top -n 1 -b -d 10 | grep -c delly)
    lumpy_processes=$(top -n 1 -b -d 10 | grep -c lumpy)
    active_threads=$(python /getThreads.py "${breakdancer_threads}" "${cnvnator_threads}" "${sambamba_processes}" "${manta_processes}" "${breakseq_processes}" "${delly_processes}" "${lumpy_processes}")
    
    while [[ $active_threads -ge $(nproc) ]]; do
        echo "Waiting for 60 seconds"
        breakdancer_threads=$(top -n 1 -b -d 10 | grep -c breakdancer)
        cnvnator_threads=$(top -n 1 -b -d 10 | grep -c cnvnator)
        sambamba_processes=$(top -n 1 -b -d 10 | grep -c sambamba)
        manta_processes=$(top -n 1 -b -d 10 | grep -c manta)
        breakseq_processes=$(top -n 1 -b -d 10 | grep -c breakseq)
        delly_processes=$(top -n 1 -b -d 10 | grep -c delly)
        lumpy_processes=$(top -n 1 -b -d 10 | grep -c lumpy)
        active_threads=$(python /getThreads.py "${breakdancer_threads}" "${cnvnator_threads}" "${sambamba_processes}" "${manta_processes}" "${breakseq_processes}" "${delly_processes}" "${lumpy_processes}")
        sleep 60
    done
}

if [[ ! -f "${illumina_bam}" ]] || [[ ! -f "${ref_fasta}" ]]; then
    if [[ "${dnanexus}" == "True" ]]; then
        dx-jobutil-report-error "ERROR: An invalid (nonexistent) input file has been specified."
    else
        echo "ERROR: An invalid (nonexistent) input file has been specified."
        echo "${illumina_bam}"
        echo "${ref_fasta}"
        exit 1
    fi
fi


cp "${ref_fasta}" /tmp/ref.fa
ln -s /tmp/ref.fa

if [[ "${run_breakdancer}" != "True" ]] && [[ "${run_breakseq}" != "True" ]] && [[ "${run_manta}" != "True" ]] && [[ "${run_cnvnator}" != "True" ]] && [[ "${run_lumpy}" != "True" ]] && [[ "${run_delly_deletion}" != "True" ]] && [[ "${run_delly_insertion}" != "True" ]] && [[ "${run_delly_inversion}" != "True" ]] && [[ "${run_delly_duplication}" != "True" ]]; then
    echo "WARNING: Did not detect any SV modules requested by the user through command-line flags."
    echo "Running with default SV modules: Breakdancer, Breakseq, Manta, CNVnator, Lumpy, and Delly Deletion"
    run_breakdancer="True"
    run_breakseq="True"
    run_manta="True"
    run_cnvnator="True"
    run_lumpy="True"
    run_delly_deletion="True"
fi

if [[ "${ref_index}" == "None" ]]; then
    samtools faidx ref.fa &
else
    cp "${ref_index}" /tmp/ref.fa.fai
    ln -s /tmp/ref.fa.fai
fi

ref_genome=$(python /home/dnanexus/get_reference.py)
lumpy_exclude_string=""
if [[ "${ref_genome}" == "b37" ]]; then
    lumpy_exclude_string="-x /home/dnanexus/b37.bed"
elif [[ "$ref_genome" == "hg19" ]]; then
    lumpy_exclude_string="-x /home/dnanexus/hg19.bed"
else
    lumpy_exclude_string="-x /home/dnanexus/hg38.bed"
fi

export lumpy_scripts="/home/dnanexus/lumpy-sv/scripts"

# Get extension and threads
extn=${illumina_bam##*.}
threads="$(nproc)"
threads=$((threads - 3))

echo "Set up and index BAM/CRAM"

# Check if BAM file has already been processed -- if so, continue
if [[ -f "/home/dnanexus/in/done.txt" ]]; then
    echo "BAM file and index both exist in the mounted volume; continuing"
else
    # Allow for CRAM files
    if [[ "${extn}" == "cram" ]] || [[ "${extn}" == "CRAM" ]]; then
        echo "CRAM file input"
        mkfifo tmp_input.bam
        samtools view "${illumina_bam}" -bh -@ "${threads}" -T ref.fa -o - | tee tmp_input.bam > input.bam & 
        samtools index tmp_input.bam
        wait
        mv tmp_input.bam.bai input.bam.bai
        rm tmp_input.bam

        mv input.bam /home/dnanexus/in/input.bam
        mv input.bam.bai /home/dnanexus/in/input.bam.bai

        rm "${illumina_bam}" && touch "${illumina_bam}"
    elif [[ "${illumina_bai}" == "None" ]]; then
        echo "BAM file input, no index exists"
        mv "${illumina_bam}" input.bam
        samtools index input.bam

        mv input.bam /home/dnanexus/in/input.bam
        mv input.bam.bai /home/dnanexus/in/input.bam.bai

        touch "${illumina_bam}"
    else
        echo "BAM file input, index exists"
        #cp "${illumina_bam}" /home/dnanexus/in/input.bam
        #cp "${illumina_bai}" /home/dnanexus/in/input.bam.bai
        cp "${illumina_bam}" /tmp/input.bam
        cp "${illumina_bai}" /tmp/input.bam.bai

        #touch "${illumina_bam}"
    fi

    #touch /home/dnanexus/in/done.txt
fi

#ln -s /home/dnanexus/in/input.bam
#ln -s /home/dnanexus/in/input.bam.bai
ln -s /tmp/input.bam
ln -s /tmp/input.bam.bai

wait

echo "Generate contigs"

samtools view -H input.bam | python /getContigs.py "${filter_short_contigs}" > contigs

#mkdir -p /home/dnanexus/out/log_files/

if [[ "${run_breakseq}" == "True" || "${run_manta}" == "True" ]]; then
    echo "Launching jobs that cannot be parallelized by contig"
fi



# JOBS THAT CANNOT BE PARALLELIZED BY CONTIG breakseq_stdout
# BREAKSEQ2
if [[ "${run_breakseq}" == "True" ]]; then
    echo "BreakSeq"
    #mkdir -p /home/dnanexus/out/log_files/breakseq_logs/
    bplib="/breakseq2_bplib_20150129/breakseq2_bplib_20150129.gff"
    work="breakseq2"
    timeout 6h /home/dnanexus/breakseq2-2.2/scripts/run_breakseq2.py --reference ref.fa \
        --bams input.bam --work "${work}" \
        --bwa /usr/local/bin/bwa --samtools /usr/local/bin/samtools \
        --bplib_gff "${bplib}" \
        --nthreads "$(nproc)" --bplib_gff "${bplib}" \
        --sample "${prefix}" 1>> breakseq_stdout.txt 2>> breakseq_err.txt &
fi

# MANTA

echo "manta_err" > manta_err.txt
if [[ "${run_manta}" == "True" ]]; then
    echo "Manta"
    timeout 6h runManta 1>> manta_stdout.txt 2>> manta_err.txt &
fi

# PREPARE FOR BREAKDANCER
if [[ "${run_breakdancer}" == "True" ]]; then
    timeout 2h /breakdancer/cpp/bam2cfg -o breakdancer.cfg input.bam
fi

concat_breakdancer_cmd=
concat_cnvnator_cmd=

count=0
delly_deletion_concat=""
delly_inversion_concat=""
delly_duplication_concat=""
delly_insertion_concat=""
lumpy_merge_command=""

if [[ "${run_delly_deletion}" == "True" ]] || [[ "${run_delly_insertion}" == "True" ]] || [[ "${run_delly_inversion}" == "True" ]] || [[ "${run_delly_duplication}" == "True" ]]; then
   run_delly="True"
fi

# Process management for launching jobs
if [[ "${run_cnvnator}" == "True" ]] || [[ "${run_delly}" == "True" ]] || [[ "${run_breakdancer}" == "True" ]] || [[ "${run_lumpy}" == "True" ]]; then
    echo "Launching jobs parallelized by contig"

    
    

    while read -r contig; do
        if [[ $(samtools view input.bam "${contig}" | head -n 20 | wc -l) -ge 10 ]]; then
            echo "Running on contig ${contig}"
            count=$((count + 1))
            
            if [[ "${run_breakdancer}" == "True" ]]; then
                echo "Running Breakdancer for contig ${contig}"
                timeout 4h /breakdancer/cpp/breakdancer-max breakdancer.cfg input.bam -o "${contig}" > breakdancer-"${count}".ctx 2> breakdancer_err.txt &
                concat_breakdancer_cmd="${concat_breakdancer_cmd} breakdancer-${count}.ctx"
            fi

            check_threads

            if [[ "$run_cnvnator" == "True" ]]; then
                echo "Running CNVnator for contig ${contig}"
                runCNVnator "${contig}" "${count}" 1>> cnvnator_stdout.txt 2>> cnvnator_err.txt &
                concat_cnvnator_cmd="${concat_cnvnator_cmd} output.cnvnator_calls-${count}"
            fi

            check_threads

            if [[ "${run_delly}" == "True" ]] || [[ "${run_lumpy}" == "True" ]]; then
                echo "Running sambamba view"
                timeout 2h sambamba view -h -f bam -t "$(nproc)" input.bam "${contig}" > chr."${count}".bam 2>> sambamba_err.txt
                echo "Running sambamba index"
                sambamba index -t "$(nproc)" chr."${count}".bam 1>> sambamba_stdout.txt 2>> sambamba_err.txt
                
                check_threads

                if [[ "${run_delly_deletion}" == "True" ]]; then  
                    echo "Running Delly (deletions) for contig $contig"
                    timeout 6h delly -t DEL -o "${count}".delly.deletion.vcf -g ref.fa chr."${count}".bam 1>> delly_del_stdout.txt 2>> delly_del_err.txt & 
                    delly_deletion_concat="${delly_deletion_concat} ${count}.delly.deletion.vcf"
                fi

                check_threads

                if [[ "${run_delly_inversion}" == "True" ]]; then 
                    echo "Running Delly (inversions) for contig $contig"
                    timeout 6h delly -t INV -o $count.delly.inversion.vcf -g ref.fa chr."${count}".bam 1>> delly_inv_stdout.txt 2>> delly_inv_err.txt & 
                    delly_inversion_concat="${delly_inversion_concat} ${count}.delly.inversion.vcf"
                fi

                check_threads

                if [[ "${run_delly_duplication}" == "True" ]]; then 
                    echo "Running Delly (duplications) for contig ${contig}"
                    timeout 6h delly -t DUP -o "${count}".delly.duplication.vcf -g ref.fa chr."${count}".bam 1>> delly_dup_stdout.txt 2>> delly_dup_err.txt & 
                    delly_duplication_concat="${delly_duplication_concat} ${count}.delly.duplication.vcf"
                fi

                check_threads

                if [[ "${run_delly_insertion}" == "True" ]]; then 
                    echo "Running Delly (insertions) for contig ${contig}"
                    timeout 6h delly -t INS -o "${count}".delly.insertion.vcf -g ref.fa chr."${count}".bam 1>> delly_ins_stdout.txt 2>> delly_ins_err.txt & 
                    delly_insertion_concat="$delly_insertion_concat $count.delly.insertion.vcf"
                fi
                
                check_threads

                if [[ "${run_lumpy}" == "True" ]]; then
                    echo "Running Lumpy for contig ${contig}"
                    timeout 6h /home/dnanexus/lumpy-sv/bin/lumpyexpress -B chr."${count}".bam -o lumpy."${count}".vcf ${lumpy_exclude_string} -k 1>> lumpy_stdout.txt 2>> lumpy_err.txt & 
                    lumpy_merge_command="$lumpy_merge_command lumpy.$count.vcf"
                fi
            fi

            check_threads

        fi
    done < contigs
fi


wait


echo "Converting results to VCF format"

(if [[ "${run_lumpy}" == "True" ]]; then
    echo "Convert Lumpy results to VCF format"
    ls -sh *.vcf
    python /convertHeader.py "${prefix}" "${lumpy_merge_command}" | vcf-sort -c | uniq > lumpy.vcf

    if [[ -f lumpy.vcf ]]; then
        cp lumpy.vcf result_lumpy.vcf

        python /vcf2bedpe.py -i lumpy.vcf -o lumpy.gff
        python /Lumpy2merge.py lumpy.gff "${prefix}" 1.0
    else
        echo "No outputs of Lumpy found. Continuing."
    fi
fi) &

(if [[ "${run_manta}" == "True" ]]; then
    echo "Convert Manta results to VCF format"
    if [[ ! -f manta/results/variants/diploidSV.vcf.gz && ! -f manta/results/stats/alignmentStatsSummary.txt ]]; then
        echo "No outputs of Manta found. Continuing."
    else  
        #cp manta/results/variants/diploidSV.vcf.gz /home/dnanexus/out/sv_caller_results/"${prefix}".manta.diploidSV.vcf.gz
        mv manta/results/variants/diploidSV.vcf.gz .
        gunzip diploidSV.vcf.gz
        cp diploidSV.vcf result_manta.vcf
        python /Manta2merge.py 1.0 diploidSV.vcf "${prefix}"

        #cp manta/results/stats/alignmentStatsSummary.txt /home/dnanexus/out/sv_caller_results/"${prefix}".manta.alignmentStatsSummary.txt
        cp manta/results/stats/alignmentStatsSummary.txt result_manta_alignment_stats_summary.txt
    fi

fi) &

(if [[ "${run_breakdancer}" == "True" ]] && [[ -n "${concat_breakdancer_cmd}" ]]; then
    echo "Convert Breakdancer results to VCF format"
    # cat contents of each file into output file: lack of quotes intentional
    cat $concat_breakdancer_cmd > breakdancer.output

    if [[ -f breakdancer.output ]]; then
        #cp breakdancer.output /home/dnanexus/out/sv_caller_results/"${prefix}".breakdancer.ctx

        python /BreakDancer2Merge.py 1.0 breakdancer.output "${prefix}"

        python /convert_breakdancer_vcf.py < breakdancer.output > breakdancer.vcf
        cp breakdancer.vcf result_breakdancer.vcf
    else
        echo "No outputs of Breakdancer found. Continuing."
    fi
fi) &

(if [[ "${run_cnvnator}" == "True" ]] && [[ -n "${concat_cnvnator_cmd}" ]]; then
    echo "Convert CNVnator results to VCF format"
    # cat contents of each file into output file: lack of quotes intentional
    cat $concat_cnvnator_cmd > cnvnator.output

    if [[ -f cnvnator.output ]]; then
        perl /usr/utils/cnvnator2VCF.pl cnvnator.output > cnvnator.vcf

        #cp cnvnator.vcf /home/dnanexus/out/sv_caller_results/"${prefix}".cnvnator.vcf
        #cp cnvnator.output /home/dnanexus/out/sv_caller_results/"${prefix}".cnvnator.output
        cp cnvnator.vcf result_cnvnator.vcf
        cp cnvnator.output result_cnvnator.output

    else
        echo "No outputs of CNVnator found. Continuing."
    fi
fi) &

(if [[ "${run_breakseq}" == "True" ]]; then
    echo "Convert Breakseq results to VCF format"
    if [[ ! -f breakseq2/breakseq_genotyped.gff && ! -f breakseq2/breakseq.vcf.gz && ! -f breakseq2/final.bam ]]; then
        echo "No outputs of Breakseq found. Continuing."
    else
        mv breakseq2/breakseq.vcf.gz .
        gunzip breakseq.vcf.gz

        #cp breakseq2/breakseq_genotyped.gff /home/dnanexus/out/sv_caller_results/"${prefix}".breakseq.gff
        #cp breakseq.vcf /home/dnanexus/out/sv_caller_results/"${prefix}".breakseq.vcf
        #cp breakseq2/final.bam /home/dnanexus/out/sv_caller_results/"${prefix}".breakseq.bam

        cp breakseq.vcf result_breakseq.vcf
    fi

    # Do the log files after we copy the output so that the 
    # cd /home/dnanexus command doesn't spoil singularity
    if [[ -z $(find "${work}" -name "*.log") ]]; then
        echo "No Breakseq log files found."
    else
        cd "${work}" || return
        find ./*.log | tar -zcvf log.tar.gz -T -
        rm -rf ./*.log
        mv log.tar.gz /home/dnanexus/out/log_files/breakseq_logs/"$prefix".breakseq.log.tar.gz
        cd /home/dnanexus || return
    fi


fi) &

(if [[ "${run_delly_deletion}" == "True" ]]; then 
    echo "Convert Delly deletion results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_deletion_concat}" | vcf-sort -c | uniq > delly.deletion.vcf

    if [[ -f delly.deletion.vcf ]]; then
        #cp delly.deletion.vcf /home/dnanexus/out/sv_caller_results/"${prefix}".delly.deletion.vcf
        cp delly.insertion.vcf result_delly_deletion.vcf
    else
        echo "No Delly deletion results found. Continuing."
    fi
fi) &

(if [[ "${run_delly_inversion}" == "True" ]]; then
    echo "Convert Delly inversion results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_inversion_concat}" | vcf-sort -c | uniq > delly.inversion.vcf

    if [[ -f delly.inversion.vcf ]]; then
        #cp delly.inversion.vcf /home/dnanexus/out/sv_caller_results/"${prefix}".delly.inversion.vcf
        cp delly.insertion.vcf result_delly_inversion.vcf
    else
        echo "No Delly inversion results found. Continuing."
    fi
fi) &

(if [[ "${run_delly_duplication}" == "True" ]]; then
    echo "Convert Delly duplication results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_duplication_concat}" | vcf-sort -c | uniq > delly.duplication.vcf

    if [[ -f delly.duplication.vcf ]]; then
        #cp delly.duplication.vcf /home/dnanexus/out/sv_caller_results/"${prefix}".delly.duplication.vcf
        cp delly.insertion.vcf result_delly_duplication.vcf
    else
        echo "No Delly duplication results found. Continuing."
    fi
fi) &

(if [[ "${run_delly_insertion}" == "True" ]]; then
    echo "Convert Delly insertion results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_insertion_concat}" | vcf-sort -c | uniq > delly.insertion.vcf

    if [[ -f delly.insertion.vcf ]]; then
        #cp delly.insertion.vcf /home/dnanexus/out/sv_caller_results/"${prefix}".delly.insertion.vcf
        cp delly.insertion.vcf result_delly_insertion.vcf
    else
        echo "No Delly insertion results found. Continuing."
    fi
fi) &

wait

set -e
# Verify that there are VCF files available
if [[ -z $(find . -name "*.vcf") ]]; then
    if [[ "${dnanexus}" == "True" ]]; then
        dx-jobutil-report-error "ERROR: SVTyper requested, but candidate VCF files required to genotype. No VCF files found."
    else
        echo "ERROR: SVTyper requested, but candidate VCF files required to genotype. No VCF files found."
        exit 1
    fi
fi
set +e




# Run SVtyper and SVviz
if [[ "${run_genotype_candidates}" == "True" ]]; then

    # Only install SVTyper if necessary
    #pip install git+https://github.com/hall-lab/svtyper.git -q &
    source /miniconda/etc/profile.d/conda.sh
    conda activate svtyper_env

    echo "Running SVTyper"
    #mkdir -p /home/dnanexus/out/svtyped_vcfs/

    i=0
    # Breakdancer
    if [[ "${run_breakdancer}" == "True" ]]; then
        echo "Running SVTyper on Breakdancer outputs"
        mkdir svtype_breakdancer
        if [[ -f breakdancer.vcf ]]; then
            bash /home/dnanexus/parallelize_svtyper.sh breakdancer.vcf svtype_breakdancer result_breakdancer.svtyped.vcf input.bam

            sed -i 's/SAMPLE/breakdancer/g' result_breakdancer.svtyped.vcf
        else
            "No Breakdancer VCF file found. Continuing."
        fi
    fi

    # Breakseq
    if [[ "${run_breakseq}" == "True" ]]; then
        echo "Running SVTyper on BreakSeq outputs"
        mkdir svtype_breakseq
        if [[ -f breakseq.vcf ]]; then
            bash /home/dnanexus/parallelize_svtyper.sh breakseq.vcf svtype_breakseq result_breakseq.svtyped.vcf input.bam
        else
            echo "No BreakSeq VCF file found. Continuing."
        fi
    fi

    # CNVnator
    if [[ "${run_cnvnator}" == "True" ]]; then
        echo "Running SVTyper on CNVnator outputs"
        mkdir svtype_cnvnator
        if [[ -f cnvnator.vcf ]]; then
            python /get_uncalled_cnvnator.py | python /add_ciend.py 1000 > cnvnator.ci.vcf < cnvnator.vcf
            bash /home/dnanexus/parallelize_svtyper.sh cnvnator.vcf svtype_cnvnator result_cnvnator.svtyped.vcf input.bam
        else
            echo "No CNVnator VCF file found. Continuing."
        fi
    fi

    # Delly
    if [[ "${run_delly}" == "True" ]]; then
        echo "Running SVTyper on Delly outputs"
        if [[ -z $(find . -name "delly*vcf") ]]; then
            echo "No Delly VCF file found. Continuing."
        else
            for item in delly*vcf; do
                mkdir svtype_delly_"${i}"
                bash /home/dnanexus/parallelize_svtyper.sh "${item}" svtype_delly_"${i}" delly.svtyper."${i}".vcf input.bam
                i=$((i + 1))
            done

            grep \# delly.svtyper.0.vcf > result_delly.svtyped.vcf

            for item in delly.svtyper*.vcf; do
                grep -v \# "${item}" >> result_delly.svtyped.vcf
            done
        fi
    fi

    # Lumpy
    if [[ "${run_lumpy}" == "True" ]]; then
        echo "Running SVTyper on Lumpy outputs"
        mkdir svtype_lumpy
        if [[ -f lumpy.vcf ]]; then
            bash /home/dnanexus/parallelize_svtyper.sh lumpy.vcf svtype_lumpy result_lumpy.svtyped.vcf input.bam
        else
            echo "No Lumpy VCF file found. Continuing."
        fi
    fi

    # Manta
    if [[ "${run_manta}" == "True" ]]; then
        echo "Running SVTyper on Manta outputs"
        if [[ -f diploidSV.vcf ]]; then
            mv diploidSV.vcf result_manta.svtyped.vcf
        else
            echo "No Manta VCF file found. Continuing."
        fi
    fi

    wait

    # deactivate svtyper
    source deactivate

    # Prepare inputs for SURVIVOR
    echo "Preparing inputs for SURVIVOR"
    for item in *svtyped.vcf; do
        python /adjust_svtyper_genotypes.py "${item}" > adjusted.vcf
        mv adjusted.vcf "${item}"
        echo "Adding ${item} to SURVIVOR inputs"
        echo "${item}" >> survivor_inputs
    done

    # Prepare SVtyped VCFs for upload
    #for item in *svtyped.vcf; do
    #    cp "${item}" /home/dnanexus/out/svtyped_vcfs/"${item}"
    #done

    # Run SURVIVOR
    echo "Running SURVIVOR"
    survivor merge survivor_inputs 1000 1 1 0 0 10 survivor.output.vcf

    # Prepare SURVIVOR outputs for upload
    vcf-sort -c > survivor_sorted.vcf < survivor.output.vcf
    sed -i 's/SAMPLE/breakdancer/g' survivor_sorted.vcf
    python /combine_combined.py survivor_sorted.vcf "${prefix}" survivor_inputs /all.phred.txt | python /correct_max_position.py > result_survivor_combined_genotyped.vcf
    cp survivor_sorted.vcf result_survivor_sorted.vcf

fi

#cd /home/dnanexus/out/log_files/ && find . -type d -empty -delete && find . -maxdepth 1 -mindepth 1 -type d -exec tar czf {}.tar.gz {} --remove-files \;







