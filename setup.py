from setuptools import setup,find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='stecfinder',
      version='1.1.2',
      description='In silico clustering and serotyping of Shigatoxin producing E. coli',
      long_description=readme(),
      classifiers=[
          'License :: OSI Approved :: GPLv3',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics STEC serotyping',
      url='https://github.com/LanLab/STECFinder',
      author='Michael Payne',
      author_email='michael.payne@unsw.edu.au',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['stecfinder=stecfinder.stecfinder:main'],
      },
      zip_safe=False)
