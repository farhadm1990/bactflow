nextflow.enable.dsl=2
NXF_CONDA_ENABLED=true


workflow {
    main:

        

        envSetUP()
       
        def pooled_out
        def dedup_fastq
        def filt_fastqs
        def cov_fastqs
        def fastas_fold
         // Get the concatenated fastq files
        if (params.run_flye) {
            pooled_out = fastqConcater(
            params.fastq_dir, 
            params.extension
            )

            //deduplication
            dedup_fastq = deduper(
            pooled_out,
            params.cpus
            )
            
            if(params.nanofilter) {
              filt_fastqs = nano_read_filt(
                    dedup_fastq,
                    params.min_quality,
                    params.min_length)
                   

            } else {
                filt_fastqs = dedup_fastq
            }
            
            
            if(params.coverage_filter) {
                cov_fastqs = coverage_filt(
                    filt_fastqs,
                    params.coverage,
                    params.genome_size
                )

                fastas_fold = assembly_flye1(
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
            pooled_fastq_files = fastqConcater(
            params.fastq_dir, 
            params.extension
            )
            assembly_unicycler(
                pooled_fastq_files,
                params.cpus,
                params.output_dir
            )
        } else if (params.run_spades) {
            assembly_spades(
                params.cpus,
                params.output_dir
            )
        } else if (params.run_megahit) {
            assembly_megahit(
                params.cpus,
                params.output_dir
            )
        }
         
        // Annotate the genomes
        if (params.prok_annot) {
            prokAnnot(
            fastas_fold,
            params.cpus
        )
        }
        if (params.tax_class) {
            taxonomyGTDBTK(
                fastas_fold,    
                params.cpus,
                params.genome_extension,
                params.gtdbtk_data_path
        )
        }
        if (params.checkm_lineag_check) {
            checkm_lineage(
            fastas_fold,
            params.checkm_db,
            params.cpus,
            params.genome_extension

        )
        }
        
    
}




process envSetUP {
    output:
    file('environment_created')//, optional: true //to prevent stoping if it failed

    script:
    """
    source \$(conda info --base)/etc/profile.d/conda.sh
    # source '/home/hackion/miniforge-pypy3/etc/profile.d/conda.sh'
    conda_con=\$(conda env list | grep "ont_helper" | awk '{print \$1}' | grep -o '[a-zA-Z]' | wc -l)
    if [ \$conda_con -eq 0 ]
    then
        mamba env create -f ${projectDir}/config.yml
        
        conda activate ont_helper
        bash ${projectDir}/r_installer_pkg.sh
        touch environment_created
        
      


    else
        touch environment_created
        conda activate ont_helper
    fi
    """
}

process fastqConcater {
    cpus params.cpus
    input:
    path fastq_dir
    val extension

    output:
   // file('concatenated_fq_are_ready') 
    path("${fastq_dir}/pooled/*.fastq"), emit: pooled_out


    script:
    """
    
    
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
    path pooled_out
    val cpus

    output:
    path('dedup/*_dedup.fastq'), emit: dedup_fastq

    script:
    """
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
    path dedup_fastq
    val min_quality
    val min_length

    output:
    path('asm_out_dir/*_filt.fastq'), emit: filt_fastqs


    script:
    """
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
    path filt_fastqs
    val coverage
    val genome_size

    output:
    path('asm_out_dir/cov_filt/*_filt.fastq'), emit: cov_fastqs

    script:
    """
    for i in ${filt_fastqs}
    do 
        if [ ! -d ./asm_out_dir/cov_filt ]
        then 
            mkdir -p ./asm_out_dir/cov_filt
        fi 

        name=\$(basename \$i | cut -f1 -d".")
        rasusa reads -c ${coverage} -g ${genome_size}mb  \$i > ./asm_out_dir/cov_filt/\$name.fastq
    done
    """
}

// assembling: with coverage filter
process assembly_flye1 {
    cpus params.cpus
    debug true
    errorStrategy 'ignore'
    label 'Assemlby'
    tag "Assembling ${cov_fastqs}"
    publishDir "${out_dir}/circulated_fasta", mode: 'copy', overwrite: false

    input:
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
    // tupl path("asm_out_dir/circulated_fasta/*.fasta"), emit: fastas
    path("circulated_fasta"), emit: fastas_fold

    script:
    
    """
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

    ######## Running circlator #########

    for i in asm_out_dir/polished_fasta/*.fasta
    do  
        prefix=\$(basename \$i | cut -f1 -d'.')

        if [ ! -d asm_out_dir/polished_fasta/circulatd_"\${prefix}" ]
        then 
            mkdir -p asm_out_dir/polished_fasta/circulatd_"\${prefix}" 
        fi 

        if [ ! -d asm_out_dir/circulated_fasta ]
        then 
            mkdir -p asm_out_dir/circulated_fasta
        fi 

        circlator fixstart \$i asm_out_dir/polished_fasta/circulatd_"\${prefix}" 

        cp asm_out_dir/polished_fasta/circulatd_"\${prefix}".fasta asm_out_dir/circulated_fasta && rm -rf asm_out_dir/polished_fasta/circulatd_*
    done


     
    """
    // important: don't pass numeric values between quotes. 
}


// assembling: with no coverage filter
process assembly_flye2 {
    cpus params.cpus
    debug true
    errorStrategy 'ignore'
    label 'Assemlby'
    tag "Assembling ${filt_fastqs}"
    publishDir "${out_dir}/circulated_fasta", mode: 'copy', overwrite: false

    input:
    path filt_fastqs
    val cpus
    val coverage
    val genome_size
    val min_length
    val min_quality
    val basecaller_model
    val tensor_batch

    when:
    ! params.coverage_filter

    output:
    // tupl path("asm_out_dir/circulated_fasta/*.fasta"), emit: fastas
    path("circulated_fasta"), emit: fastas_fold

    script:
    
    """
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

    ######## Running circlator #########

    for i in asm_out_dir/polished_fasta/*.fasta
    do  
        prefix=\$(basename \$i | cut -f1 -d'.')

        if [ ! -d asm_out_dir/polished_fasta/circulatd_"\${prefix}" ]
        then 
            mkdir -p asm_out_dir/polished_fasta/circulatd_"\${prefix}" 
        fi 

        if [ ! -d asm_out_dir/circulated_fasta ]
        then 
            mkdir -p asm_out_dir/circulated_fasta
        fi 

        circlator fixstart \$i asm_out_dir/polished_fasta/circulatd_"\${prefix}" 

        cp asm_out_dir/polished_fasta/circulatd_"\${prefix}".fasta asm_out_dir/circulated_fasta && rm -rf asm_out_dir/polished_fasta/circulatd_*
    done

     
    """
    // important: don't pass numeric values between quotes. 
}

// gene annotations

process prokAnnot {
    cpus params.cpus
    publishDir "${out_dir}/prokk_out", mode: 'copy', overwrite: false
    
    input:
    path fastas_fold
    val cpus 
    

    output:
    path("prokk_out"), optional: true //so that it deons't stop upon failing

    script:
    
    """
    echo "This is the ${cpus} number." > cpus_prok.txt

    source \$(conda info --base)/etc/profile.d/conda.sh
    conda activate ont_helper
    bash ${projectDir}/prokka_annot.sh -g "${fastas_fold}"/*.fasta -c ${cpus}
    
    """
}

// taxonomy classification by gtdbtk
process taxonomyGTDBTK {
    cpus params.cpus
    publishDir "final_output/gtdbtk_out", mode: 'copy', overwrite: false
    input:
    path fastas_fold
    val cpus
    val genome_extension
    val gtdbtk_data_path
    
    output:
    path("gtdbtk_out"),  optional: true //so that it deons't stop upon failing

    script:
    """
    echo "This is the ${cpus} number." 
    bash ${projectDir}/gtdbtk.sh -g '${fastas_fold}' -c ${cpus} -e '${genome_extension}' -d '${gtdbtk_data_path}'
    """
}

process checkm_lineage {
    cpus params.cpus

    publishDir "${out_dir}", mode: copy, overwrite: false

    input:
    path fastas_fold
    val checkm_db
    val cpus
    val genome_extension

    output:
    tupl path('checkm_lineage.txt'), path('taxon_tree.newick'), path('genome_tree.newick'),  optional: true //so that it deons't stop upon failing

    script:
    """
    echo "This is the ${cpus} number." 
    checkm data setRoot '${dcheckm_db}'
    checkm lineage_wf -t ${cpus} --pplacer_threads ${cpus} -x '${genome_extension}' '${fastas_fold}' checkm_lineage && \
    checkm qa  -t ${cpus} checkm_lineage/lineage.ms checkm_lineage/  > checkm_lineage.txt 

    checkm tree -r --nt -t ${cpus}  -x '${genome_extension}' --pplacer_threads ${cpus}  '${fastas_fold}' checkm_tree && checkm tree_qa -o 4 --tab_table -f taxon_tree.newick checkm_tree && checkm tree_qa -o 3 --tab_table -f genome_tree.newick checkm_tree
    """

}







// nextflow.bak run genebrosh --fastq_dir /home/projects/cu_10168/people/farpan/data/nadjia2/fastq_pass --extension '.fastq.gz' --cpus 40 -profile conda --tensor_batch 200 --genome_size 5 --coverage 80 --tax_class  true --checkm_lineage_check true --prok_annot true  --coverage_filter false --nanofilter false -w nadjia2 -process.echo -resume 

