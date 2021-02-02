#!/usr/bin/env python3
from setuptools import setup,find_packages
import textwrap

#The function for Python packing
setup(
        name="IsoSplitter",
        version="1.2",
        description=("identification of alternative splicing sites using long-read transcriptome without a reference genome"),
        author="Hengfu-Yin",
        author_email="hfyin@caf.ac.cn",
        license='University of Oxford Academic Use Licence',
        keywords="Splitter",
        url="https://github.com/Hengfu-Yin/IsoSplitter",
        packages=['Example','scripts'],
        package_dir={'Example':'Example','scripts':'scripts'},
        #scripts=['scripts'],
        package_data = {
                'Example':['*.txt'],
        },
        entry_points={
                'console_scripts': [
                        'ShortReadsAligner=scripts.ShortReadsAligner:Main',
                        'IsoSplittingAnchor=scripts.IsoSplittingAnchor:Main'
                ]
        },
        classifiers=textwrap.dedent("""
        	Development Status :: 3 - Alpha
        	Intended Audience :: Science/Research
        	License :: Other/Proprietary License
        	Operating System :: OS Independent
        	Programming Language :: Python :: 3.5
        	Programming Language :: Python :: 3.6
        	Programming Language :: Python :: 3.7
        	Topic :: Software Development :: Build Tools
        	Topic :: System :: Archiving :: Packaging
        	""").strip().splitlines(),
        #python_requires='>=3.5',
)
