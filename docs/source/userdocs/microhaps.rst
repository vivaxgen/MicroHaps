

Introduction to DADA2-based MicroHaplotype Analysis
=====================================================

This pipeline contains a wrapper for the MIT-Broad team `DADA2 software <https://github.com/broadinstitute/malaria-amplicon-pipeline>`_.


Quick Start 
------------

#. Running the pipeline is as easy executing the following command:

.. code-block:: console

    $ ngs-pl run-microhaplotype-caller -u 4 --panel pvvvg-mhap -o output path_to_fastq/*.fastq.gz

.. note::
    The FASTQ files should be in fastq.gz format (gzip-compressed), and the
    filenames should reflect the sample name, eg: ``sample_1_date_batch_pool_R1.fastq.gz``.

The ``-u`` option is used to specify the number of underscores to remove (counted in reverse order) to obtain the actual sample name.
For example, if the sample name is ``sample_1`` and the fastq file is named ``sample_1_date_batch_pool_R1.fastq.gz``, then the ``-u`` argument should be 4.

The ``--panel`` option is used to specify the panel that will be used for the analysis.
The available panels are:
    * ``pvvvg-mhap``: *P. Vivax* panel
    * ``pfspotmal-drug``: SpotMalaria (*P. falciparum*) drugs-resistance panel
    * ``pfspotmal-mhap``: SpotMalaria (*P. falciparum*) microhaplotype panel
Additional panels can be added, for more details please refer to the :doc:`developer's documentation <../develdocs/implementation>`.

The ``-o`` option is used to specify the output directory. 

.. warning::

    **For laptop users!**    
    It is essential that you specify ``-j <n>``, where ``n`` is a small number (dependent on the available system memory & cpu), as this limits the number of jobs running at one time. Without this argument, 
    the pipeline will utilise too much system memory and crash.


#. When the command finishes, examine the content of ``output`` directory

.. code-block:: console
    
    $ tree output
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
