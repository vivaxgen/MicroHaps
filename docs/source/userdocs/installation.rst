Installation on local devices and HPCs
=========================================

.. code-block:: console

	$ "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/MicroHaps/main/install.sh)

Copy and paste command into the command line in the folder you want the install to be saved to (we recommend you create a specific tools or software folder). 
When prompted for ``"Pipeline base directory? [./vvg-MicroHaps]"`` press enter again for the install to proceed.

The installation requires ~ 20-45 minutes as most of R packages need to be recompiled
during installation.

Once the installation finished, it will show the command to activate the
pipeline, as such:

.. code-block:: console

	$ /path/to/vvg-MicroHaps/bin/activate

.. warning::

    **For Conda-based users!**

    Be sure you are not in a conda environment or in the (base) conda environment prior to installing. 
    To deactivate your conda environment or (base) environment, enter:

    .. code-block:: console

        $ conda deactivate

This activation command has to be executed before all commands of the pipeline
can be run. When activated, the terminal will show the ``(µhaps)`` prompt.

The installation process also performs indexing of the reference files.
However, in case that the indexing fails, please perform manual indexing
using the command:

.. code-block:: console

	$ ngs-pl initialize --target wgs

To test your install, and read about programme specifications / options:

.. code-block:: console

	$ ngs-pl run-discovery-variant-caller --help


Updating the pipeline
----------------------

To update the pipeline including its dependencies, assuming that the environment has been activated,
run the following command:

.. code-block:: console

	$ $VVGBIN/update-box

To only update the pipeline code (so not updating its dependencies), run the following instead:

.. code-block:: console

        $ $VVGBIN/update-box --pull-repo-only


