import rematch

from setuptools import setup

VERSION = rematch.__version__

with open('README.md') as fh:
    README = fh.read()

setup(
    name='ReMatCh',
    version='{}'.format(VERSION),
    packages=['ReMatCh',
              'ReMatCh.modules'],
    package_dir={'ReMatCh': 'ReMatCh'},
    package_data={'ReMatCh': ['utils/*',
                              'modules/mlst_schemas/*'
                              'src/*']},
    data_files=[('', ['LICENSE'])],
    install_requires=[
        'biopython',
    ],
    description='Reads mapping against target sequences, checking mapping and consensus sequences production',
    long_description=README,
    long_description_content_type='text/markdown',
    keywords=['reference mapping', 'allele variant calling', 'consensus sequence', 'sequence presence/absence',
                'MLST multilocus sequence typing'],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Environment :: Console',
        'Operating System :: Unix',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    url='https://github.com/B-UMMI/ReMatCh',
    author='Miguel P. Machado',
    author_email='mpmachado@medicina.ulisboa.pt',
    license='GPL3',
    entry_points={
        'console_scripts': [
            'rematch.py = ReMatCh.rematch:main',
            'rematch = ReMatCh.rematch:main'
        ]
    },
    python_requires='>=3.4'
)
