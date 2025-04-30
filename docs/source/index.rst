

vivaxGEN Microhaplotypes Pipeline Documentation
===============================================

VivaxGEN-Microhaps is an open-source pipeline designed for processing targeted
amplicon sequencing data.
It offers two distinct sub-pipelines:

* **Microhaplotype Calling:** This sub-pipeline, adapted from the
  `malaria-amplicon-pipeline <https://github.com/broadinstitute/malaria-amplicon-pipeline>`_
  and converted to use Snakemake, generates microhaplotype allele tables.

* **SNP-based Variant Calling:** Leveraging the
  `vivaxGEN NGS-Pipeline <https://github.com/vivaxgen/ngs-pipeline>`_,
  this sub-pipeline generates standard VCF files for SNP analysis.


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
    