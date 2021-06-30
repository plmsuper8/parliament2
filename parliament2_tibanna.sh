#!/usr/bin/env bash

illumina_bam=$1 #The name of the Illumina BAM file for which to call structural variants containing mapped reads.
illumina_bai="${illumina_bam}.bai" #The name of the corresponding index for the Illumina BAM file.
ref_fasta=$2 #The name of the reference file that matches the reference used to map the Illumina inputs.
ref_index="${ref_fasta}.fai" #The name of the corresponding index for the reference genome file.
prefix="result" #(Optional) If provided, all output files will start with this. If absent, the base of the BAM file name will be used. CWL need to be adjusted
filter_short_contigs=$3 #If selected, SV calls will not be generated on contigs shorter than 1 MB.
run_breakdancer=$4 #If selected, the program Breakdancer will be one of the SV callers run.
run_breakseq=$5 #If selected, the program BreakSeq2 will be one of the SV callers run.
run_manta=$6 #If selected, the program Manta will be one of the SV callers run.
run_cnvnator=${7} #If selected, the program CNVnator will be one of the SV callers run.
run_lumpy=${8} #If selected, the program Lumpy will be one of the SV callers run.
run_delly_deletion=${9} #If selected, the deletion module of the program Delly2  will be one of the SV callers run.
run_delly_insertion=${10} #If selected, the insertion module of the program Delly2  will be one of the SV callers run.
run_delly_inversion=${11} #If selected, the inversion module of the program Delly2  will be one of the SV callers run.
run_delly_duplication=${12} #If selected, the duplication module of the program Delly2  will be one of the SV callers run.
run_genotype_candidates=${13} #If selected, candidate events determined from the individual callers will be genotyped and merged to create a consensus output.
sample_name=${14} # Name of the last 'result' column in the final vcf
#min_qual=${17} # Hard filter
run_svviz="False"
svviz_only_validated_candidates="False"
dnanexus="False"

mkdir -p /tmp/output

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
    echo "ERROR: An invalid (nonexistent) input file has been specified."
    exit 1
fi

if [[ ! -f "${illumina_bai}" ]]; then
    echo "ERROR: The bam index file is missing"
    exit 1
fi

if [[ ! -f "${ref_index}" ]]; then
    echo "ERROR: The reference fasta file is missing"
    exit 1
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
    echo "FAI file not specified. Exiting"
    exit 1
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
        echo "ERROR: CRAM files not supported"
        exit 1
    elif [[ "${illumina_bai}" == "None" ]]; then
        echo "ERROR: Illumina index file missing"
        exit 1
    else
        echo "BAM file input, index exists"
    fi

fi

ln -s "${illumina_bam}" input.bam
ln -s "${illumina_bai}" input.bam.bai

wait

echo "Generate contigs"

samtools view -H input.bam | python /getContigs.py "${filter_short_contigs}" > contigs || exit 1

mkdir -p /tmp/output/log_files/

if [[ "${run_breakseq}" == "True" || "${run_manta}" == "True" ]]; then
    echo "Launching jobs that cannot be parallelized by contig"
fi

# JOBS THAT CANNOT BE PARALLELIZED BY CONTIG
# BREAKSEQ2
if [[ "${run_breakseq}" == "True" ]]; then
    echo "BreakSeq"
    #mkdir -p /home/dnanexus/out/log_files/breakseq_logs/
    mkdir -p /tmp/output/log_files/breakseq_logs/
    bplib="/breakseq2_bplib_20150129/breakseq2_bplib_20150129.gff"
    work="breakseq2"
    timeout 6h /home/dnanexus/breakseq2-2.2/scripts/run_breakseq2.py --reference ref.fa \
        --bams input.bam --work "${work}" \
        --bwa /usr/local/bin/bwa --samtools /usr/local/bin/samtools \
        --bplib_gff "${bplib}" \
        --nthreads "$(nproc)" --bplib_gff "${bplib}" \
        --sample "${sample_name}" 1> /tmp/output/log_files/breakseq_logs/"${prefix}".breakseq.stdout.log 2> /tmp/output/log_files/breakseq_logs/"${prefix}".breakseq.stderr.log || exit 1 &
fi

# MANTA
if [[ "${run_manta}" == "True" ]]; then
    echo "Manta"
    #mkdir -p /home/dnanexus/out/log_files/manta_logs/
    mkdir -p /tmp/output/log_files/manta_logs/

    timeout 6h runManta 1> /tmp/output/log_files/manta_logs/"${prefix}".manta.stdout.log 2> /tmp/output/log_files/manta_logs/"${prefix}".manta.stderr.log || exit 1 &
fi

# PREPARE FOR BREAKDANCER
if [[ "${run_breakdancer}" == "True" ]]; then
    echo "Prepare Breakdancer"
    timeout 2h /breakdancer/cpp/bam2cfg -o breakdancer.cfg input.bam || exit 1
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
    mkdir -p /tmp/output/log_files/breakdancer_logs/
    mkdir -p /tmp/output/log_files/cnvnator_logs/
    mkdir -p /tmp/output/log_files/delly_deletion_logs/
    mkdir -p /tmp/output/log_files/delly_duplication_logs/
    mkdir -p /tmp/output/log_files/delly_insertion_logs/
    mkdir -p /tmp/output/log_files/delly_inversion_logs/
    mkdir -p /tmp/output/log_files/lumpy_logs/
    mkdir -p /tmp/output/log_files/sambamba_logs/

    while read -r contig; do
        if [[ $(samtools view input.bam "${contig}" | head -n 20 | wc -l) -ge 10 ]]; then
            echo "Running on contig ${contig}"
            count=$((count + 1))
            
            if [[ "${run_breakdancer}" == "True" ]]; then
                echo "Running Breakdancer for contig ${contig}"
                timeout 4h /breakdancer/cpp/breakdancer-max breakdancer.cfg input.bam -o "${contig}" > breakdancer-"${count}".ctx 2> /tmp/output/log_files/breakdancer_logs/"${prefix}".breakdancer."${contig}".stderr.log || exit 1 &
                concat_breakdancer_cmd="${concat_breakdancer_cmd} breakdancer-${count}.ctx"
            fi

            check_threads

            if [[ "$run_cnvnator" == "True" ]]; then
                echo "Running CNVnator for contig ${contig}"
                runCNVnator "${contig}" "${count}" 1> /tmp/output/log_files/cnvnator_logs/"${prefix}".cnvnator."${contig}".stdout.log 2> /tmp/output/log_files/cnvnator_logs/"${prefix}".cnvnator."${contig}".stderr.log || exit 1 &
                concat_cnvnator_cmd="${concat_cnvnator_cmd} output.cnvnator_calls-${count}"
            fi

            check_threads

            if [[ "${run_delly}" == "True" ]] || [[ "${run_lumpy}" == "True" ]]; then
                echo "Running sambamba view"
                timeout 2h sambamba view -h -f bam -t "$(nproc)" input.bam "${contig}" > chr."${count}".bam 2> /tmp/output/log_files/sambamba_logs/"${prefix}".sambamba."${contig}".stderr.log || exit 1
                echo "Running sambamba index"
                sambamba index -t "$(nproc)" chr."${count}".bam 1> /tmp/output/log_files/sambamba_logs/"${prefix}".sambamba."${contig}".stdout.log 2> /tmp/output/log_files/sambamba_logs/"${prefix}".sambamba."${contig}".stderr.log || exit 1
                
                check_threads

                if [[ "${run_delly_deletion}" == "True" ]]; then  
                    echo "Running Delly (deletions) for contig $contig"
                    timeout 6h delly -t DEL -o "${count}".delly.deletion.vcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_deletion_logs/"${prefix}".delly_deletion."${contig}".stdout.log 2> /tmp/output/log_files/delly_deletion_logs/"${prefix}".delly_deletion."${contig}".stderr.log || exit 1 & 
                    #delly_v0_8_7 call -t DEL -o "${count}".delly.deletion.bcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_deletion_logs/"${prefix}".delly_deletion."${contig}".stdout.log 2> /tmp/output/log_files/delly_deletion_logs/"${prefix}".delly_deletion."${contig}".stderr.log || exit 1 & 
                    #bcftools view "${count}".delly.deletion.bcf > "${count}".delly.deletion.vcf
                    delly_deletion_concat="${delly_deletion_concat} ${count}.delly.deletion.vcf"
                fi

                check_threads

                if [[ "${run_delly_inversion}" == "True" ]]; then 
                    echo "Running Delly (inversions) for contig $contig"
                    timeout 6h delly -t INV -o $count.delly.inversion.vcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_inversion_logs/"${prefix}".delly_inversion."${contig}".stdout.log 2> /tmp/output/log_files/delly_inversion_logs/"${prefix}".delly_inversion."${contig}".stderr.log || exit 1 & 
                    #delly_v0_8_7 call -t INV -o "${count}".delly.inversion.bcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_inversion_logs/"${prefix}".delly_inversion."${contig}".stdout.log 2> /tmp/output/log_files/delly_inversion_logs/"${prefix}".delly_inversion."${contig}".stderr.log || exit 1 & 
                    #bcftools view "${count}".delly.inversion.bcf > "${count}".delly.inversion.vcf
                    delly_inversion_concat="${delly_inversion_concat} ${count}.delly.inversion.vcf"
                fi

                check_threads

                if [[ "${run_delly_duplication}" == "True" ]]; then 
                    echo "Running Delly (duplications) for contig ${contig}"
                    timeout 6h delly -t DUP -o "${count}".delly.duplication.vcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_duplication_logs/"${prefix}".delly_duplication.stdout.log 2> /tmp/output/log_files/delly_duplication_logs/"${prefix}".delly_duplication.stderr.log || exit 1 & 
                    #delly_v0_8_7 call -t DUP -o "${count}".delly.duplication.bcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_duplication_logs/"${prefix}".delly_duplication."${contig}".stdout.log 2> /tmp/output/log_files/delly_duplication_logs/"${prefix}".delly_duplication."${contig}".stderr.log || exit 1 & 
                    #bcftools view "${count}".delly.duplication.bcf > "${count}".delly.duplication.vcf
                    delly_duplication_concat="${delly_duplication_concat} ${count}.delly.duplication.vcf"
                fi

                check_threads

                if [[ "${run_delly_insertion}" == "True" ]]; then 
                    echo "Running Delly (insertions) for contig ${contig}"
                    timeout 6h delly -t INS -o "${count}".delly.insertion.vcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_insertion_logs/"${prefix}".delly_insertion."${count}".stdout.log 2> /tmp/output/log_files/delly_insertion_logs/"${prefix}".delly_insertion."${count}".stderr.log || exit 1 & 
                    #delly_v0_8_7 call -t INS -o "${count}".delly.insertion.bcf -g ref.fa chr."${count}".bam 1> /tmp/output/log_files/delly_insertion_logs/"${prefix}".delly_insertion."${contig}".stdout.log 2> /tmp/output/log_files/delly_insertion_logs/"${prefix}".delly_insertion."${contig}".stderr.log || exit 1 & 
                    #bcftools view "${count}".delly.insertion.bcf > "${count}".delly.insertion.vcf
                    delly_insertion_concat="$delly_insertion_concat $count.delly.insertion.vcf"
                fi
                
                check_threads

                if [[ "${run_lumpy}" == "True" ]]; then
                    echo "Running Lumpy for contig ${contig}"
                    timeout 6h /home/dnanexus/lumpy-sv/bin/lumpyexpress -B chr."${count}".bam -o lumpy."${count}".vcf ${lumpy_exclude_string} -k 1> /tmp/output/log_files/lumpy_logs/"${prefix}".lumpy."${count}".stdout.log 2> /tmp/output/log_files/lumpy_logs/"${prefix}".lumpy."${count}".stderr.log || exit 1 & 
                    lumpy_merge_command="$lumpy_merge_command lumpy.$count.vcf"
                fi
            fi

            check_threads

        fi
    done < contigs
fi

wait

echo "Converting results to VCF format"
mkdir -p /tmp/output/sv_caller_results/

(if [[ "${run_lumpy}" == "True" ]]; then
    echo "Convert Lumpy results to VCF format"
    ls -sh *.vcf
    python /convertHeader.py "${prefix}" "${lumpy_merge_command}" | vcf-sort -c | uniq > lumpy.vcf || exit 1

    if [[ -f lumpy.vcf ]]; then
        cp lumpy.vcf /tmp/output/sv_caller_results/"${prefix}".lumpy.vcf

        python /vcf2bedpe.py -i lumpy.vcf -o lumpy.gff || exit 1
        python /Lumpy2merge.py lumpy.gff "${prefix}" 1.0 || exit 1
    else
        echo "No outputs of Lumpy found. Continuing."
    fi
fi) &

(if [[ "${run_manta}" == "True" ]]; then
    echo "Convert Manta results to VCF format"
    if [[ ! -f manta/results/variants/diploidSV.vcf.gz && ! -f manta/results/stats/alignmentStatsSummary.txt ]]; then
        echo "No outputs of Manta found. Continuing."
    else  
        cp manta/results/variants/diploidSV.vcf.gz /tmp/output/sv_caller_results/"${prefix}".manta.diploidSV.vcf.gz
        mv manta/results/variants/diploidSV.vcf.gz .
        gunzip diploidSV.vcf.gz
        python /Manta2merge.py 1.0 diploidSV.vcf "${prefix}" || exit 1

        cp manta/results/stats/alignmentStatsSummary.txt /tmp/output/sv_caller_results/"${prefix}".manta.alignmentStatsSummary.txt
    fi

fi) &

(if [[ "${run_breakdancer}" == "True" ]] && [[ -n "${concat_breakdancer_cmd}" ]]; then
    echo "Convert Breakdancer results to VCF format"
    # cat contents of each file into output file: lack of quotes intentional
    cat $concat_breakdancer_cmd > breakdancer.output

    if [[ -f breakdancer.output ]]; then
        cp breakdancer.output /tmp/output/sv_caller_results/"${prefix}".breakdancer.ctx

        python /BreakDancer2Merge.py 1.0 breakdancer.output "${prefix}" || exit 1

        python /convert_breakdancer_vcf.py < breakdancer.output > breakdancer.vcf || exit 1
        cp breakdancer.vcf /tmp/output/sv_caller_results/"${prefix}".breakdancer.vcf
    else
        echo "No outputs of Breakdancer found. Continuing."
    fi
fi) &

(if [[ "${run_cnvnator}" == "True" ]] && [[ -n "${concat_cnvnator_cmd}" ]]; then
    echo "Convert CNVnator results to VCF format"
    # cat contents of each file into output file: lack of quotes intentional
    cat $concat_cnvnator_cmd > cnvnator.output

    if [[ -f cnvnator.output ]]; then
        perl /usr/utils/cnvnator2VCF.pl cnvnator.output > cnvnator.vcf || exit 1

        cp cnvnator.vcf /tmp/output/sv_caller_results/"${prefix}".cnvnator.vcf
        cp cnvnator.output /tmp/output/sv_caller_results/"${prefix}".cnvnator.output
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

        cp breakseq2/breakseq_genotyped.gff /tmp/output/sv_caller_results/"${prefix}".breakseq.gff
        cp breakseq.vcf /tmp/output/sv_caller_results/"${prefix}".breakseq.vcf
        cp breakseq2/final.bam /tmp/output/sv_caller_results/"${prefix}".breakseq.bam
    fi

    # Do the log files after we copy the output so that the 
    # cd /home/dnanexus command doesn't spoil singularity
    if [[ -z $(find "${work}" -name "*.log") ]]; then
        echo "No Breakseq log files found."
    else
        cd "${work}" || return
        find ./*.log | tar -zcvf log.tar.gz -T -
        rm -rf ./*.log
        mv log.tar.gz /tmp/output/log_files/breakseq_logs/"$prefix".breakseq.log.tar.gz
        cd /home/dnanexus || return
    fi


fi) &

(if [[ "${run_delly_deletion}" == "True" ]]; then 
    echo "Convert Delly deletion results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_deletion_concat}" | vcf-sort -c | uniq > delly.deletion.vcf || exit 1

    if [[ -f delly.deletion.vcf ]]; then
        cp delly.deletion.vcf /tmp/output/sv_caller_results/"${prefix}".delly.deletion.vcf
    else
        echo "No Delly deletion results found. Continuing."
    fi
fi) &

(if [[ "${run_delly_inversion}" == "True" ]]; then
    echo "Convert Delly inversion results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_inversion_concat}" | vcf-sort -c | uniq > delly.inversion.vcf || exit 1

    if [[ -f delly.inversion.vcf ]]; then
        cp delly.inversion.vcf /tmp/output/sv_caller_results/"${prefix}".delly.inversion.vcf
    else
        echo "No Delly inversion results found. Continuing."
    fi
fi) &

(if [[ "${run_delly_duplication}" == "True" ]]; then
    echo "Convert Delly duplication results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_duplication_concat}" | vcf-sort -c | uniq > delly.duplication.vcf || exit 1

    if [[ -f delly.duplication.vcf ]]; then
        cp delly.duplication.vcf /tmp/output/sv_caller_results/"${prefix}".delly.duplication.vcf
    else
        echo "No Delly duplication results found. Continuing."
    fi
fi) &

(if [[ "${run_delly_insertion}" == "True" ]]; then
    echo "Convert Delly insertion results to VCF format"
    python /convertHeader.py "${prefix}" "${delly_insertion_concat}" | vcf-sort -c | uniq > delly.insertion.vcf || exit 1

    if [[ -f delly.insertion.vcf ]]; then
        cp delly.insertion.vcf /tmp/output/sv_caller_results/"${prefix}".delly.insertion.vcf
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
    mkdir -p /tmp/output/svtyped_vcfs/

    svtyper_start=`date +%s`

    i=0
    # Breakdancer
    if [[ "${run_breakdancer}" == "True" ]]; then
        echo "Running SVTyper on Breakdancer outputs"
        mkdir svtype_breakdancer
        if [[ -f breakdancer.vcf ]]; then
            bash /home/dnanexus/parallelize_svtyper.sh breakdancer.vcf svtype_breakdancer "${prefix}".breakdancer.svtyped.vcf input.bam 1000 || exit 1

            sed -i 's/SAMPLE/breakdancer/g' "${prefix}".breakdancer.svtyped.vcf
        else
            "No Breakdancer VCF file found. Continuing."
        fi
    fi

    # Breakseq
    if [[ "${run_breakseq}" == "True" ]]; then
        echo "Running SVTyper on BreakSeq outputs"
        mkdir svtype_breakseq
        if [[ -f breakseq.vcf ]]; then
            #mv breakseq.vcf "${prefix}".breakseq.svtyped.vcf
            bash /home/dnanexus/parallelize_svtyper.sh breakseq.vcf svtype_breakseq "${prefix}".breakseq.svtyped.vcf input.bam 1000 || exit 1
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
            bash /home/dnanexus/parallelize_svtyper.sh cnvnator.vcf svtype_cnvnator "${prefix}".cnvnator.svtyped.vcf input.bam 1000 || exit 1
            #mv cnvnator.vcf "${prefix}".cnvnator.svtyped.vcf
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
                bash /home/dnanexus/parallelize_svtyper.sh "${item}" svtype_delly_"${i}" delly.svtyper."${i}".vcf input.bam 1000 || exit 1
                i=$((i + 1))
            done

            grep \# delly.svtyper.0.vcf > "${prefix}".delly.svtyped.vcf

            for item in delly.svtyper*.vcf; do
                grep -v \# "${item}" >> "${prefix}".delly.svtyped.vcf
            done

            ## Merge files and the run SVtype - slower
            # # create header of merged file
            # dellyfiles=(delly*vcf)
            # grep \# "${dellyfiles[0]}" > d.merged.vcf # don't call it delly, so that it's not in the loop below

            # touch d.merged.tmp.vcf
            # for item in delly*vcf; do
            #     grep -v \# "${item}" >> d.merged.tmp.vcf
            # done

            # # Shuffle the variants, since some are harder to genotype than others
            # #shuf d.merged.tmp.vcf > d.merged.tmp.shuffled.vcf
            # cat d.merged.tmp.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > d.merged.tmp.sorted.vcf
            # #grep -v \# d.merged.tmp.shuffled.vcf >> d.merged.vcf
            # grep -v \# d.merged.tmp.sorted.vcf >> d.merged.vcf

            # mkdir svtype_delly
            # bash /home/dnanexus/parallelize_svtyper.sh d.merged.vcf svtype_delly delly.svtyper.vcf input.bam

            # mv delly.svtyper.vcf "${prefix}".delly.svtyped.vcf

        fi
    fi

    # Lumpy
    if [[ "${run_lumpy}" == "True" ]]; then
        echo "Running SVTyper on Lumpy outputs"
        mkdir svtype_lumpy
        if [[ -f lumpy.vcf ]]; then
            bash /home/dnanexus/parallelize_svtyper.sh lumpy.vcf svtype_lumpy "${prefix}".lumpy.svtyped.vcf input.bam 1000 || exit 1
        else
            echo "No Lumpy VCF file found. Continuing."
        fi
    fi

    # Manta
    if [[ "${run_manta}" == "True" ]]; then
        echo "Running SVTyper on Manta outputs"
        if [[ -f diploidSV.vcf ]]; then
            mv diploidSV.vcf "${prefix}".manta.svtyped.vcf
        else
            echo "No Manta VCF file found. Continuing."
        fi
    fi

    wait

    svtyper_end=`date +%s`
    svtyper_runtime=$((svtyper_end-svtyper_start))
    echo "SVTyper runtime: $svtyper_runtime seconds"

    # deactivate svtyper
    source deactivate

    # Prepare inputs for SURVIVOR
    echo "Preparing inputs for SURVIVOR"
    for item in *svtyped.vcf; do
        echo "Adding ${item} to SURVIVOR inputs"
        python /adjust_svtyper_genotypes.py "${item}" > adjusted.vcf || exit 1
        mv adjusted.vcf "${item}"
        echo "${item}" >> survivor_inputs
        echo "DONE"
    done

    # Prepare SVtyped VCFs for upload
    for item in *svtyped.vcf; do
        cp "${item}" /tmp/output/svtyped_vcfs/"${item}"
    done

    # Run SURVIVOR
    echo "Running SURVIVOR"
    survivor merge survivor_inputs 1000 1 1 0 0 10 survivor.output.vcf || exit 1

    # Prepare SURVIVOR outputs for upload
    echo "Sorting SURVIVOR output"
    vcf-sort -c > survivor_sorted.vcf < survivor.output.vcf || exit 1
    sed -i 's/SAMPLE/breakdancer/g' survivor_sorted.vcf
    echo "Run combine script"
    python /combine_combined.py survivor_sorted.vcf "${sample_name}" survivor_inputs /all.phred.txt | python /correct_max_position.py > /tmp/output/"${prefix}".combined.genotyped.vcf || exit 1
    cp survivor_sorted.vcf /tmp/output/"${prefix}".survivor_sorted.vcf
    echo "Sorting again"
    cat /tmp/output/"${prefix}".combined.genotyped.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > "${prefix}".combined.genotyped.sorted.vcf || exit 1
    bgzip -c "${prefix}".combined.genotyped.sorted.vcf > "${prefix}".combined.genotyped.sorted.vcf.gz || exit 1
    tabix -p vcf "${prefix}".combined.genotyped.sorted.vcf.gz || exit 1

    #cp "${prefix}".combined.genotyped.sorted.vcf /tmp/output/"${prefix}".combined.genotyped.unfiltered.vcf
    cp "${prefix}".combined.genotyped.sorted.vcf.gz /tmp/output/"${prefix}".combined.genotyped.sorted.vcf.gz
    cp "${prefix}".combined.genotyped.sorted.vcf.gz.tbi /tmp/output/"${prefix}".combined.genotyped.sorted.vcf.gz.tbi

    # Hard filtering - NOT DOING THIS FOR NOW
    # bcftools view -i 'QUAL>'${min_qual} "${prefix}".combined.genotyped.sorted.vcf.gz > "${prefix}".combined.genotyped.filtered.vcf
    # cp "${prefix}".combined.genotyped.filtered.vcf /tmp/output/"${prefix}".combined.genotyped.filtered.vcf
    # bgzip -c "${prefix}".combined.genotyped.filtered.vcf > "${prefix}".combined.genotyped.filtered.vcf.gz
    # tabix -p vcf "${prefix}".combined.genotyped.filtered.vcf.gz

    # cp "${prefix}".combined.genotyped.filtered.vcf.gz /tmp/output/"${prefix}".combined.genotyped.filtered.vcf.gz
    # cp "${prefix}".combined.genotyped.filtered.vcf.gz.tbi /tmp/output/"${prefix}".combined.genotyped.filtered.vcf.gz.tbi

    
fi

WRITEABLEDIR=`pwd`

cd /tmp/output/log_files/ && find . -type d -empty -delete && find . -maxdepth 1 -mindepth 1 -type d -exec tar czf {}.tar.gz {} --remove-files \;

cd /tmp/output
zip -r result_tmp.zip .
cd $WRITEABLEDIR
cp /tmp/output/result_tmp.zip result.zip
cp /tmp/output/"${prefix}".combined.genotyped.sorted.vcf.gz variants.vcf.gz
cp /tmp/output/"${prefix}".combined.genotyped.sorted.vcf.gz.tbi variants.vcf.gz.tbi


