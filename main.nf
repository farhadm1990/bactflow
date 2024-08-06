nextflow.enable.dsl=2
NXF_CONDA_ENABLED=true
nextflow.preview.output=true
params.help = false

def helpMessage = """

Usage: nextflow run bactflow [options]

Options:
   
    --setup_only            If true, only runs envSetUp(), default false
    --fastq_dir             Absolute path to the fastq_pass directory (required). 
    --concat_reads          Default true, it concatenates all your ONT basecaller 4000-chunk reads into one fastq file. Set it to false if it is already concatenated.
    --extension             String; extention of basecalled fastq files; default '.fastq.gz'
    --cpus                  Number of available cpus; default 1.
    --coverage_filter       If you want to normalize all your genomes to a certain coverage (default false).
    --coverage              Only if '--coverage_filter true'; default is 50.
    --genome_size           Genome size for coverage normalizaiton. Only if '--coverage_filter true'; default is 6.
    --out_dir               Output directory of your final results. Default "genebrosh_output"
    --tensor_batch          Medaka tensorflow batch size. Lower it in low coverage genomes. Default 200.
    --nanofilter            Filtering reads for length and quality; default true.
    --min_length            If '--nanofilter' true, filter reads below a certain read length (default 1000). 
    --min_quality           If '--nanofilter' true, filter reads below a certain read quality (default 16 for R10.4.1 flowcells). 
    --medaka_polish         If true, it will polish assembled genomes by medaka (dfault false).
    --basecaller_model      Basecaller model for medaka polishing step. 'r1041_e82_400bps_hac_v4.2.0'
    --checkm_lineag_check   If true, the genomes will be checked for their lineage completeness in one bin (default false).
    --genome_extension      Required if '--checkm_lineag_check true'; default fasta.
    --run_flye              If true, it runs Flye assembler; default true.
    --run_unicycler         If true, it runs Unicycler hybrid assemlber, default false.
    --run_megahit           If true, it runs Megahit assembler, default false.
    --run_spades            If true, it runs Spades assembler, default false.
    --tax_class             If true, it runs GTBtk taxonomic classification, default true.
    --prok_annot            If true, it runs gene annotaiton by Prokka, default false. 
    --run_checkm            If true, it runs checmk lineage and phylogenetic tree workflow.
    --checkm_db             An absolute path to the Checkm database.  
    --gtdbtk_data_path      Absolute path to the GDBtk database. 
    --run_quast             Post-assembly stats by Quaset, default true.
    --genome_dir            Path to already assembled genomes, only to run post-assembly tasks, e.g. taxonomy classification, gene annotations and quast or checkm 

"""


workflow {

    if (params.help) {
    log.info helpMessage
    exit 0
   }


    main:

        
        if(params.setup_only){
            def env_check =  envSetUP()
            testify(env_check)
        } else {
            def env_check =  envSetUP()
            testify(env_check)

        def pooled_out
        def dedup_fastq
        def filt_fastqs
        def cov_fastqs
        def fastas_fold
        def quast_out
        def circ_fasta
         // Get the concatenated fastq files
        if (params.run_flye) {
            if (params.concat_reads){
                pooled_out = fastqConcater(
                params.fastq_dir, 
                params.extension
            )
            } else {
                pooled_out = Channel.fromPath(params.fastq_dir + "/*.fastq")
                                .collect() // creates a value channel as such that it puts all the fastq files in one list channel.
            }
           

            //deduplication
            dedup_fastq = deduper(
            env_check,
            pooled_out,
            params.cpus
            )
            
            if(params.nanofilter) {
              filt_fastqs = nano_read_filt(
                    env_check,
                    dedup_fastq,
                    params.min_quality,
                    params.min_length)
                   

            } else {
                filt_fastqs = dedup_fastq
            }
            
            
            if(params.coverage_filter) {
                cov_fastqs = coverage_filt(
                    env_check,
                    filt_fastqs,
                    params.coverage,
                    params.genome_size
                )

                fastas_fold = assembly_flye1(
                env_check,
                cov_fastqs,
                params.cpus,
                params.coverage,
                params.genome_size,
                params.min_length,
                params.min_quality,
                params.basecaller_model,
                params.tensor_batch,
                params.medaka_polish
            )
                
            } else {
                fastas_fold = assembly_flye2(
                env_check,
                filt_fastqs,
                params.cpus,
                params.coverage,
                params.genome_size,
                params.min_length,
                params.min_quality,
                params.basecaller_model,
                params.tensor_batch,
                params.medaka_polish
            )
            }

            
            

     
            
        } else if (params.run_unicycler) {
            pooled_out = fastqConcater(
            env_check,
            params.fastq_dir, 
            params.extension
            )
           fastas_fold = assembly_unicycler(
                env_check,
                pooled_out,
                params.cpus,
                params.output_dir
            )
        } else if (params.run_spades) {
          fastas_fold =   assembly_spades(
                env_check,
                params.cpus,
                params.output_dir
            )
        } else if (params.run_megahit) {
           fastas_fold =  assembly_megahit(
                env_check,
                params.cpus,
                params.output_dir
            )
        } else {
            fastas_fold = Channel.fromPath(params.genome_dir)
                                    .collect() // 
        }
         
        // circulator
        
        if (params.circle_genome){
            circ_fasta = circulator(
                env_check,
                fastas_fold
            )
        } else {
            circ_fasta = fastas_fold
        }
        // Annotate the genomes
        if (params.prok_annot) {
            prokAnnot(
            env_check,
            circ_fasta,
            params.cpus
        )
        }
        if (params.tax_class) {
            taxonomyGTDBTK(
                env_check,
                circ_fasta,    
                params.cpus,
                params.genome_extension,
                params.gtdbtk_data_path
        )
        }
        if (params.run_checkm) {
            checkm_lineage(
            env_check,
            circ_fasta,
            params.checkm_db,
            params.cpus,
            params.genome_extension

        )
        }

        // quast stats
        if (params.run_quast) {
            quast_out = quast_check(
            env_check,
            circ_fasta,
            params.cpus
            )
        }
        }
        
        
        
        
    
}

output {
    directory params.out_dir
}


process envSetUP {
    output:
    path('environment_created'), emit: env_check//, optional: true //to prevent stoping if it failed

    script:

    """
    
    source \$(conda info --base)/etc/profile.d/conda.sh
    
    conda_con=\$(conda env list | grep "bactflow" | awk '{print \$1}' | grep -o '[a-zA-Z]' | wc -l)
    if [ \$conda_con -eq 0 ]
    then
        if command -v mamba
        then
            mamba env create -f ${baseDir}/config.yml
        else
            conda env create -f ${baseDir}/config.yml
        fi

        conda activate bactflow
        
        
        
      


    else
        
        conda activate bactflow
        
    fi
    
    bash ${baseDir}/r_installer_pkg.sh
    
    touch environment_created
    """
}

process testify {
    
    input:
    path env_check
    
    // conda 'cofig.yml'

    output:
   
    path("checked.txt"), emit: env_testify

    publish:
    env_testify >> 'env_testify'

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow
    conda env list | grep bactflow > checked.txt
    
    """
}

process fastqConcater {
    cpus params.cpus

    input:
    
    path fastq_dir
    val extension

    when:
    params.concat_reads

    output:
   // file('concatenated_fq_are_ready') 
    path("${fastq_dir}/pooled/*.fastq"), emit: pooled_out


    script:
    """
    
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

    num_file=\$(ls "${fastq_dir}"| wc -l)

    if [ -d "${fastq_dir}"/pooled ] 
    then 
        num_pooled=\$(ls "${fastq_dir}"/pooled | wc -l)
        num_dif=\$(expr "\${num_file}" - "\${num_pooled}" )
        if [ \$num_dif -eq 1 ]
        then
            touch concatenated_fq_are_ready
        else 
            if [ "${extension}" = ".fastq.gz" ]
            then

                for i in "${fastq_dir}"/*
                do 
                    if [ -d \$i ] && [ "\$(basename \$i)" != "pooled" ]
                    then

                        name=\$(basename \$i)
                        zcat \$i/*fastq.gz >> "${fastq_dir}"/"\${name}"_pooled.fastq
                        mv "${fastq_dir}"/"\${name}"_pooled.fastq "${fastq_dir}"/pooled
                    fi
                    
                    
                done
                touch concatenated_fq_are_ready

            elif [ "${extension}" = ".fastq" ] 
            then	
                for i in "${fastq_dir}"/*    
                do 
                    if [ -d \$i ] && [ "\$(basename \$i)" != "pooled" ]
                    then

                        name=\$(basename \$i)
                        cat \$i/*fastq.gz >> "${fastq_dir}"/"\${name}"_pooled.fastq
                        mv "${fastq_dir}"/"\${name}"_pooled.fastq "${fastq_dir}"/pooled
                    fi
                done
                touch concatenated_fq_are_ready
            else
                echo "Your extention is not recognized!"

            fi
        fi



    else 
        
        mkdir -p "${fastq_dir}"/pooled

        if [ "${extension}" = ".fastq.gz" ]
        then

            for i in "${fastq_dir}"/*
            do 
                if [ -d \$i ] && [ "\$(basename \$i)" != "pooled" ]
                then

                    name=\$(basename \$i)
                    zcat \$i/*fastq.gz >> "${fastq_dir}"/"\${name}"_pooled.fastq
                    mv "${fastq_dir}"/"\${name}"_pooled.fastq "${fastq_dir}"/pooled
                fi
                
                
            done
            touch concatenated_fq_are_ready

        elif [ "${extension}" = ".fastq" ] 
        then	
            for i in "${fastq_dir}"/*    
            do 
                if [ -d \$i ] && [ "\$(basename \$i)" != "pooled" ]
                then

                    name=\$(basename \$i)
                    cat \$i/*fastq.gz >> "${fastq_dir}"/"\${name}"_pooled.fastq
                    mv "${fastq_dir}"/"\${name}"_pooled.fastq "${fastq_dir}"/pooled
                fi
            done
            touch concatenated_fq_are_ready
        else
            echo "Your extention is not recognized!"

        fi

    fi
    """
}

// read deduplication
process deduper {
    cpus params.cpus

    input:
    path env_check
    path pooled_out
    val cpus

    output:
    path('dedup/*_dedup.fastq'), emit: dedup_fastq

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

    if [ ! -d dedup ]
    then 
        mkdir -p dedup 
    fi

    for i in ${pooled_out}
    do 
        name=\$(basename \$i | cut -f 1 -d'.')
        seqkit rmdup  \$i -s -j ${cpus} > dedup/"\${name}"_dedup.fastq
    done
    """
}

// Filtering nano pore reads
process nano_read_filt {
    cpus 1

    input:
    path env_check
    path dedup_fastq
    val min_quality
    val min_length

    when:
    params.nanofilter

    output:
    path('asm_out_dir/*_filt.fastq'), emit: filt_fastqs


    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

    if [ ! -d asm_out_dir ]
    then 
        mkdir -p asm_out_dir 
    fi

    for i in ${dedup_fastq}
    do 
        name=\$(basename \$i | cut -f 1 -d'.')
        NanoFilt -q ${min_quality} -l ${min_length} \$i > asm_out_dir/"\${name}"_filt.fastq

        if [ -f asm_out_dir/"\${name}"_filt.fastq ]
        then 
            echo "File \$i is now filtered and ready to be used :) "
        fi
    done

    """
}

// coverage filtering
process coverage_filt {
    cpus 3

    input:
    path env_check
    path filt_fastqs
    val coverage
    val genome_size

    when:
    params.coverage_filter

    output:
    path('asm_out_dir/cov_filt/*.fastq'), emit: cov_fastqs

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

    for i in ${filt_fastqs}
    do 
        if [ ! -d asm_out_dir/cov_filt ]
        then 
            mkdir -p asm_out_dir/cov_filt
        fi 

        name=\$(basename \$i | cut -f1 -d".")
        rasusa reads -c ${coverage} -g ${genome_size}mb  \$i > asm_out_dir/cov_filt/\$name.fastq
    done
    """
}

// assembling: with coverage filter
process assembly_flye1 {
    cpus params.cpus
    debug true
  //  errorStrategy 'ignore'
    label 'Assemlby'
    tag "Assembling ${cov_fastqs}"
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false

    input:
    path env_check
    path cov_fastqs
    val cpus
    val coverage
    val genome_size
    val min_length
    val min_quality
    val basecaller_model
    val tensor_batch
    val medaka_polish

    when:
    params.coverage_filter

    output:
  
    path('asm_out_dir/polished_fasta'), emit: fastas_fold

    script:
    
    """
    
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow


    medaka tools download_models --quiet
    if [ ! -d asm_out_dir ]
    then
        mkdir -p asm_out_dir
    fi 

    for i in ${cov_fastqs}
    do
    
        out_name=\$(basename \$i | cut -f 1 -d'_')

        echo "running flye..."
        
        flye --nano-raw \$i -t ${cpus} -i 2 --out-dir asm_out_dir/"\${out_name}"_flye  #--asm-coverage ${coverage} -g ${genome_size}m

        
        if [ '${medaka_polish}' == "true" ]
        then 
            echo "polishing fasta reads..."
            mini_align -i \$i -r asm_out_dir/"\${out_name}"_flye/assembly.fasta -P -m -p asm_out_dir/"\${out_name}"_flye/read_to_draft_\$out_name -t ${cpus} 

            medaka consensus asm_out_dir/"\${out_name}"_flye/read_to_draft_\$out_name.bam asm_out_dir/"\${out_name}"_flye/\$out_name.hdf --batch ${tensor_batch} --threads ${cpus} --model '${basecaller_model}'

            medaka stitch asm_out_dir/"\${out_name}"_flye/\$out_name.hdf  asm_out_dir/"\${out_name}"_flye/assembly.fasta asm_out_dir/"\${out_name}"_flye/"\${out_name}"_polished.fasta
            rm -rf asm_out_dir/"\${out_name}"_flye/*bam* asm_out_dir/"\${out_name}"_flye/*.hdf asm_out_dir/"\${out_name}"_flye/*.fai asm_out_dir/"\${out_name}"_flye/*.mmi asm_out_dir/"\${out_name}"_flye/*.bed
        fi

        
    done

    if [ ! -d asm_out_dir/polished_fasta  ]
    then
        mkdir -p asm_out_dir/polished_fasta 
    fi

    cp asm_out_dir/*_flye/*_polished.fasta  asm_out_dir/polished_fasta 

    # Final message 

    echo "your polished fasta files are ready in asm_out_dir/polished_fasta."

   


     
    """
    // important: don't pass numeric values between quotes. 
}


// assembling: with no coverage filter
process assembly_flye2 {
    cpus params.cpus
    debug true
   // errorStrategy 'ignore'
    label 'Assemlby'
    tag "Assembling ${filt_fastqs}"
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false

    input:
    path env_check
    path filt_fastqs
    val cpus
    val coverage
    val genome_size
    val min_length
    val min_quality
    val basecaller_model
    val tensor_batch
    val medaka_polish

    when:
    ! params.coverage_filter

    output:
   
    path('asm_out_dir/polished_fasta'), emit: fastas_fold

    script:
    
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

   

    medaka tools download_models --quiet
    if [ ! -d asm_out_dir ]
    then
        mkdir -p asm_out_dir
    fi 

    for i in ${filt_fastqs}
    do
    
        out_name=\$(basename \$i | cut -f 1 -d'_')

        echo "running flye..."
        
        flye --nano-raw \$i -t ${cpus} -i 2 --out-dir asm_out_dir/"\${out_name}"_flye  #--asm-coverage ${coverage} -g ${genome_size}m

        if [ '${medaka_polish}' == "true" ]
        then 
            echo "polishing fasta reads..."
            mini_align -i \$i -r asm_out_dir/"\${out_name}"_flye/assembly.fasta -P -m -p asm_out_dir/"\${out_name}"_flye/read_to_draft_\$out_name -t ${cpus} 

            medaka consensus asm_out_dir/"\${out_name}"_flye/read_to_draft_\$out_name.bam asm_out_dir/"\${out_name}"_flye/\$out_name.hdf --batch ${tensor_batch} --threads ${cpus} --model '${basecaller_model}'

            medaka stitch asm_out_dir/"\${out_name}"_flye/\$out_name.hdf  asm_out_dir/"\${out_name}"_flye/assembly.fasta asm_out_dir/"\${out_name}"_flye/"\${out_name}"_polished.fasta
            rm -rf asm_out_dir/"\${out_name}"_flye/*bam* asm_out_dir/"\${out_name}"_flye/*.hdf asm_out_dir/"\${out_name}"_flye/*.fai asm_out_dir/"\${out_name}"_flye/*.mmi asm_out_dir/"\${out_name}"_flye/*.bed
        fi

        
    done

    if [ ! -d asm_out_dir/polished_fasta  ]
    then
        mkdir -p asm_out_dir/polished_fasta 
    fi

    cp asm_out_dir/*_flye/*_polished.fasta  asm_out_dir/polished_fasta 

    # Final message 

    echo "your polished fasta files are ready in asm_out_dir/polished_fasta."

    

     
    """
    // important: don't pass numeric values between quotes. 
}

// circulating the genomes
process circulator {
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false

    input:
    path env_check
    path fastas_fold
    
    when:
    params.circle_genome

    output:
    path("circulated_fasta"), emit: circ_fasta 

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

    echo "Running circlator"

    for i in "${fastas_fold}"/*.fasta
    do  
        prefix=\$(basename \$i | cut -f1 -d'.')

        if [ ! -d circulatd_"\${prefix}" ]
        then 
            mkdir -p circulatd_"\${prefix}" 
        fi 

        if [ ! -d circulated_fasta ]
        then 
            mkdir -p circulated_fasta
        fi 

        circlator fixstart \$i circulatd_"\${prefix}" 

        cp circulatd_"\${prefix}".fasta circulated_fasta && rm -rf circulatd_*
    done

    """
}

// gene annotations

process prokAnnot {
    cpus params.cpus
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    
    input:
    path env_check
    path circ_fasta
    val cpus 
    
    when:
    params.prok_annot

    output:
    path('prokk_out'), optional: true //so that it deons't stop upon failing

    script:
    
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

  

    bash ${projectDir}/prokka_annot.sh -g "${circ_fasta}" -c ${cpus}
    
    """
}

// taxonomy classification by gtdbtk
process taxonomyGTDBTK {
    cpus params.cpus
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false

    input:
    path env_check
    path circ_fasta
    val cpus
    val genome_extension
    val gtdbtk_data_path
    
    output:
    path('gtdbtk_out'),  optional: true //so that it deons't stop upon failing

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow
    # Upgrade for gtdbtk
    python -m pip install gtdbtk --upgrade

    bash ${projectDir}/gtdbtk.sh -g '${circ_fasta}' -c ${cpus} -e '${genome_extension}' -d '${gtdbtk_data_path}'
    """
}

process checkm_lineage {
    if (params.cpus > 1) {
        cpus params.cpus -1 
    }

    publishDir "${params.out_dir}", mode: 'copy', overwrite: false

    input:
    path env_check
    path circ_fasta
    val checkm_db
    val cpus
    val genome_extension

    when:
    params.run_checkm

    output:
    tuple path('checkm_lineage.txt'), path('taxon_tree.newick'), path('genome_tree.newick'), path('genome_tree.tree'), emit: checkm_out,  optional: true //so that it deons't stop upon failing

    script:
    """
    #!/usr/bin/bash
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow
    
    pip install --upgrade checkm-genome
    checkm data setRoot '${checkm_db}'
    checkm lineage_wf -t ${cpus} --pplacer_threads ${cpus} -x '${genome_extension}' '${circ_fasta}' checkm_lineage && \
    checkm qa  -t ${cpus} checkm_lineage/lineage.ms checkm_lineage/  > checkm_lineage.txt 

    checkm tree -r --nt -t ${cpus}  -x '${genome_extension}' --pplacer_threads ${cpus}  '${circ_fasta}' checkm_tree && checkm tree_qa -o 4 --tab_table -f taxon_tree.newick checkm_tree && checkm tree_qa -o 3 --tab_table -f genome_tree.newick checkm_tree

    


    # Building the genome-based tree
    Rscript -e "
    library(Biostrings)
    library(msa)
    library(ape)
    library(tidyverse)
    library(readr)
    library(seqinr)

    seqs <- Biostrings::readDNAStringSet('checkm_tree/storage/tree/concatenated.fasta', format = 'fasta')
    als <- msa(seqs)

    als_seqinr <- msaConvert(als, type = 'seqinr::alignment')
    
    dis <- dist.alignment(als_seqinr, 'identity')
    tr <- nj(dis)

    write.tree(phy = tr, file = 'genome_tree.tree')
    "

    




    """


}


// quast assembly stats
process quast_check {
    cpus params.cpus
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    errorStrategy 'ignore'

    input:
    path env_check
    path circ_fasta
    val cpus

    when:
    params.run_quast

    output:
    path('quast_report'), emit: quast_out, optional: true

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate bactflow

    #Update numpy 
    pip install --upgrade numpy

    
    if [ ! -d quast_stat ]
    then 
        mkdir -p quast_stat

    fi 

    quast '${circ_fasta}'/*.fasta -o quast_stat -t ${cpus}
    """
}









