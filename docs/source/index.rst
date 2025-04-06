

vivaxGEN Microhaplotypes Pipeline Documentation
===============================================
This repository contains 2 pipelines for handling microhaplotype amplicon
sequencing data:

* DADA2-based MicroHaplotype Analysis Pipeline; Original software can be found at `Dada2 github <https://benjjneb.github.io/dada2/>`_ but is modified for this pipeline.

* SNP-based variant calling which output ordinary VCF file, based on
  `vivaxGEN NGS-Pipeline <https://github.com/vivaxgen/ngs-pipeline>`_.


Installation on local devices and HPCs
---------------------------------------

.. code-block:: console

	"${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/MicroHaps/main/install.sh)

Copy and paste command into the command line in the folder you want the install to be saved to (we recommend you create a specific tools or software folder). 
When prompted for ``"Pipeline base directory? [./vvg-MicroHaps]"`` press enter again for the install to proceed.

The installation requires ~ 20-45 minutes as most of R packages need to be recompiled
during installation.

Once the installation finished, it will show the command to activate the
pipeline, as such:

.. code-block:: console

	/path/to/vvg-MicroHaps/bin/activate

.. warning::

    **For Conda-based users!**

    Be sure you are not in a conda environment or in the (base) conda environment prior to installing. 
    To deactivate your conda environment or (base) environment, enter:

    .. code-block:: console

        conda deactivate

This activation command has to be executed before all commands of the pipeline
can be run. When activated, the terminal will show the ``(µhaps)`` prompt.

The installation process also performs indexing of the reference files.
However, in case that the indexing fails, please perform manual indexeing
using the command:

.. code-block:: console

	ngs-pl initialize --target wgs

To test your install, and read about programme specifications / options:

.. code-block:: console

 	ngs-pl run-discovery-variant-caller --help


Updating the pipeline
----------------------

To update the pipeline line, assuming that the environment has been activated,
run the following command:

.. code-block:: console

	$VVGBIN/update-pipeline.sh

-----


Introduction to DADA2-based MicroHaplotype Analysis
----------------------------------------------------

This pipeline contains a wrapper for the MIT-Broad team `DADA2 software <https://github.com/broadinstitute/malaria-amplicon-pipeline>`_.


Quick Start 
^^^^^^^^^^^

#. Running the pipeline is as easy executing the following command:

.. code-block:: console

    ngs-pl run-microhaplotype-caller -u 4 --panel pvvvg-mhap -o output path_to_fastq/*.fastq.gz

.. note::
    The FASTQ files should be in fastq.gz format (gzip-compressed), and the
    filenames should reflect the sample name, eg: ``sample_1_date_batch_pool_R1.fastq.gz``.

The ``-u`` option is used to specify the number of underscores to remove (counted in reverse order) to obtain the actual sample name.
For example, if the sample name is ``sample_1`` and the fastq file is named ``sample_1_date_batch_pool_R1.fastq.gz``, then the ``-u`` argument should should be 4.

The ``--panel`` option is used to specify the panel that will be used for the analysis. Currently, only ``pvvvg-mhap`` is available.

The ``-o`` option is used to specify the output directory. 

.. warning::

    **For laptop users!**    
    It is essential that you specify ``-j 1`` as this limits the number of jobs running at one time. Without this argument, 
    the pipeline will utilise too much system memory and crash.


#. When the command finishes, examine the content of ``output`` directory

.. code-block:: console

    output/
        alignments/
            marker_1.fasta
            marker_1.msa
            marker_2.fasta
            marker_2.msa
            ...
        samples/
            Sample_1/
                reads/
                    raw-0_R1.fastq.gz
                    raw-0_R2.fastq.gz
                maps/
                    final.bam
                    final.bam.bai
                    Sample_1-0.bam
                mhaps-reads/
                    primer-trimmed_R1.fastq.gz
                    primer-trimmed_R2.fastq.gz
                    target_R1.fastq.gz
                    target_R2.fastq.gz
                logs/
                    ...
            Sample_2/
                ...
            ...
        malamp/
            dada2/
                ...
            ASVSeqs.fasta
            ASVTable.txt
            asv_to_cigar
            depths.tsv
            outputCIGAR.tsv
            marker_missingness.png
            marker_missingness.tsv
            sample_missingness.png
            sample_missingness.tsv
            meta
        final.coverages.tsv
        final.depths.tsv
        mylog.txt
        stats.tsv


The primary output file of interest is the ``outputCIGAR.tsv`` which contains the haplotype and their frequencies across the samples.


Introduction to SNP-based variant calling
------------------------------------------

This pipeline contains a wrapper for `vivaxGEN NGS-Pipeline <https://github.com/vivaxgen/ngs-pipeline/>`_, a SNP-based variant
calling pipeline, which has been setup to process *P. vivax* sequence data.
This pipeline will produce ordinary VCF file that can be used for further
downstream analysis.

Currently, this pipeline was setup to use multi-step mode of vivaxGEN
NGS-Pipeline in a single command line with GATK-based workflow to generate VCF
file
There is also Freebayes-based workflow, but it requires modification of the
setting, which will not be covered in this documentation.

Quick Start
^^^^^^^^^^^

#. Running the pipeline is as easy executing the following command:

.. code-block:: console

    ngs-pl run-discovery-variant-caller -o output_dir path_to_fastq/*.fastq.gz

#. When the command finishes, examine the content of ``output_dir`` directory:

The layout of the output directory is:

.. code-block:: console

    output_dir/
        metafile/
                manifest.tsv
        analysis/
                SAMPLE-1/
                    gvcf/
                        SAMPLE-1-PvP01_all_v2.g.vcf.gz
                        SAMPLE-1-PvP01_all_v2.g.vcf.gz.tbi
                    maps/
                        mapped-final.bam
                        mapped-final.bam.bai
                    reads/
                        raw-0_R1.fastq.gz
                        raw-0_R2.fastq.gz
                    logs/
                        ...
                SAMPLE-2/
                ...
        joint/
            concatenated.vcf.gz
            vcfs/
        reports/
            ...
        completed_samples/
            SAMPLE-1/
                ...
            SAMPLE-2/
                ...
        failed_samples/
            ...
        
The primary output file of interest is the ``concatenated.vcf.gz`` which contains the SNPs and the genotype calls for each of the samples.