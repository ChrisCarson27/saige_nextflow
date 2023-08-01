/*
Chris Carson modifying the work of 
Florian Wuennemann for LPC at Penn

module load lapack/3.8.0
*/

process saige_step1_bin_binary {
    tag "step1_${phenoFile.baseName}" //tagging the processes with the corresponding pheno file
    publishDir "${params.outdir}/${phenoFile.baseName}/SAIGE_out_step1", mode: 'copy' //publish the stage 1 output files to outdir/phenofile/saige_out_step1

    input:
    tuple val(grm_plink_input), path(bed), path(bim), path(fam) //this is a tuple of all the bim/fam/bed files necessary for the plinkfile input
    each path(phenoFile) //Execute the process for each element in the input phenofiles
    val phenoCol
    val covarColList
    val qcovarColList
    val sampleIDcol

    output:
    val(phenoFile.baseName), emit: phenotype
    tuple val(phenoFile.baseName), path("step1_${phenoFile.baseName}_out.rda"), path("step1_${phenoFile.baseName}.varianceRatio.txt"), emit: step1_out
    
    script:
    """
    step1_fitNULLGLMM.R \
      --plinkFile=${grm_plink_input} \
      --phenoFile="${phenoFile}" \
      --phenoCol=${phenoCol} \
      --covarColList=${covarColList} \
      --qCovarColList=${qcovarColList} \
      --sampleIDColinphenoFile=${sampleIDcol} \
      --traitType=binary \
      --outputPrefix="step1_${phenoFile.baseName}_out" \
      --nThreads=${task.cpus} \
      --IsOverwriteVarianceRatioFile=${params.saige_step1_extra_flags}
    """
  }

  process saige_step1_bin_quantitative {
    tag "step1_${phenoFile.baseName}" //tagging the processes with the corresponding pheno file
    publishDir "${params.outdir}/${phenoFile.baseName}/SAIGE_out_step1", mode: 'copy' //publish the stage 1 output files to outdir/phenofile/saige_out_step1

    input:
    tuple val(grm_plink_input), path(bed), path(bim), path(fam) //this is a tuple of all the bim/fam/bed files necessary for the plinkfile input
    each path(phenoFile) //Execute the process for each element in the input phenofiles
    val phenoCol
    val covarColList
    val qcovarColList
    val sampleIDcol
    val invNormalize

    output:
    val(phenoFile.baseName), emit: phenotype
    tuple val(phenoFile.baseName), path("step1_${phenoFile.baseName}_out.rda"), path("step1_${phenoFile.baseName}.varianceRatio.txt"), emit: step1_out
    
    script:
    """
    step1_fitNULLGLMM.R \
      --plinkFile=${grm_plink_input} \
      --phenoFile="${phenoFile}" \
      --phenoCol=${phenoCol} \
      --covarColList=${covarColList} \
      --qCovarColList=${qcovarColList} \
      --sampleIDColinphenoFile=${sampleIDcol} \
      --invNormalize=${invNormalize}
      --traitType=quantitative \
      --outputPrefix="step1_${phenoFile.baseName}_out" \
      --nThreads=${task.cpus} \
      --IsOverwriteVarianceRatioFile=${params.saige_step1_extra_flags}
    """
  }


/*
 --sparseGRMFile=${TODO} \
        --sparseGRMSampleIDFile=${TODO} \ 
*/
process saige_step2_spa_bgen_binary {
  tag "step2_chr${chrom}.${phenotype}"
  publishDir "${params.outdir}/${phenotype}/SAIGE_out_step2", mode: 'copy'

  input:
  tuple val(phenotype), val(rda), val(varRatio)
  each chrom
  val(bgen_prefix)
  val(bgen_suffix)
  val(bgen_path)
  path(sampleFile)
  val minMAC
  val minMAF

  output:
  tuple val(phenotype), path("${phenotype}.chr${chrom}.SAIGE.gwas.txt"), emit: assoc_res

  script:
  """
  step2_SPAtests.R \
        --bgenFile=${bgen_path}/${bgen_prefix}${chrom}${bgen_suffix} \
        --bgenFileIndex=${bgen_path}/${bgen_prefix}${chrom}${bgen_suffix}.bgi \
        --sampleFile=${sampleFile} \
        --AlleleOrder=ref-first \
        --chrom=${chrom} \
        --minMAC=${minMAC} \
        --minMAF=${minMAF} \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${varRatio} \
        --pCutoffforFirth=${pcutoffFirth} \
        --is_Firth_beta=${firthbetastatus} \
        --LOCO=TRUE
  """
}

process saige_step2_spa_PLINK_binary {
  tag "step2_chr${chrom}.${phenotype}"
  publishDir "${params.outdir}/${phenotype}/SAIGE_out_step2", mode: 'copy'

  input:
  path(plinkfiles) 
  tuple val(phenotype), val(rda), val(varRatio)
  each chrom
  path(sampleFile)
  val minMAC
  val minMAF
  val loco

  output:
  tuple val(phenotype), path("${phenotype}.chr${chrom}.SAIGE.gwas.txt"), emit: assoc_res

  script:
  """
  step2_SPAtests.R \
        --bedFile=${plinkfiles}.bed      \
        --bimFile=${plinkfiles}.bim       \
        --famFile=${plinkfiles}.fam       \
        --sampleFile=${sampleFile} \
        --AlleleOrder=ref-first \
        --chrom=${chrom} \
        --minMAC=${minMAC} \
        --minMAF=${minMAF} \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${varRatio} \
        --pCutoffforFirth=${pcutoffFirth} \
        --is_Firth_beta=${firthbetastatus} \
        --LOCO=${loco}
  """
}

process saige_step2_spa_VCF_binary {
  tag "step2_chr${chrom}.${phenotype}"
  publishDir "${params.outdir}/${phenotype}/SAIGE_out_step2", mode: 'copy'

  input:
  path(vcfFile) 
  path(vcfFileIndex)
  path(vcfField)
  tuple val(phenotype), val(rda), val(varRatio)
  each chrom
  path(sampleFile)
  val minMAC
  val minMAF
  val loco

  output:
  tuple val(phenotype), path("${phenotype}.chr${chrom}.SAIGE.gwas.txt"), emit: assoc_res

  script:
  """
  step2_SPAtests.R \
        --vcfFile=${vcfFile}      \
        --vcfFileIndex=${vcfFileIndex}       \
        --vcfField=${vcfField}       \
        --sampleFile=${sampleFile} \
        --AlleleOrder=ref-first \
        --chrom=${chrom} \
        --minMAC=${minMAC} \
        --minMAF=${minMAF} \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${varRatio} \
        --pCutoffforFirth=${pcutoffFirth} \
        --is_Firth_beta=${firthbetastatus} \
        --LOCO=${loco}
  """
}



/*
Need to implement eventually. R files not included in script
process prepare_files {
  tag "${phenotype}"
  publishDir "${params.outdir}/${phenotype}", mode: 'copy'

  input:
  tuple val(phenotype), path(merged_assoc)

  output:
  tuple val(phenotype), path("saige_results_${phenotype}_top_n.csv"), emit: top_hits
  tuple val(phenotype), path("saige_results_${phenotype}.csv"), emit: merged_out
  

  script:
  """
  # creates 2 .csv files, saige_results_<params.output_tag>.csv, saige_results_top_n.csv
  concat_chroms.R \
    --saige_output_name='saige_results' \
    --filename_pattern='${phenotype}.*' \
    --output_tag='${phenotype}' \
    --top_n_sites=${params.top_n_sites} \
    --max_top_n_sites=${params.max_top_n_sites}
    mv saige_results_top_n.csv saige_results_${phenotype}_top_n.csv
  """
}
*/


