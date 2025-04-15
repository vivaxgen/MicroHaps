

vivaxGEN Microhaplotypes Pipeline Documentation
===============================================
This repository contains 2 pipelines for handling microhaplotype amplicon
sequencing data:

* DADA2-based MicroHaplotype Analysis Pipeline; Original software can be found at `Dada2 github <https://benjjneb.github.io/dada2/>`_ but is modified for this pipeline.

* SNP-based variant calling which output ordinary VCF file, based on
  `vivaxGEN NGS-Pipeline <https://github.com/vivaxgen/ngs-pipeline>`_.


.. toctree::
    :maxdepth: 2
    :caption: User Documentation

    userdocs/installation.rst
    userdocs/microhaps.rst
    userdocs/discovery_calling.rst


.. toctree::
    :maxdepth: 2
    :caption: Developer Documentation

    develdocs/implementation.rst
    