import setuptools
import os


setuptools.setup(

	name="Menzies_MicroHaps",

	version="0.0.1",
	packages=["microhap_seq"],
	license="MIT",
	long_description="MicroHap sequencing command line tool",
	scripts= ["scripts/%s" % x for x in os.listdir("scripts")],
)
