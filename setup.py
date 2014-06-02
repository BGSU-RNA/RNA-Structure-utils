from distutils.core import setup

setup(
    name='RNAStructure',
    version='0.1.1',
    author='Blake A. Sweeney',
    author_email='bsweene@bgsu.edu',
    packages=['rnastructure', 'rnastructure.primary', 'rnastructure.secondary',
              'rnastructure.tertiary', 'rnastructure.util'],
    url='',
    license='LICENSE.txt',
    description='Some tools to parse RNA Structure',
    long_description="""
This is a series of tools and wrappers for things related to RNA structure.
This primarly focuses on wrapping other programs in a nice python interface,
such as RNAalifold and UNAfold. However, this also provides tools for parsing
secondary structure. There is also a small wrapper around the CIF file reader
provided by NDB to create a more pythonic interface.
    """
)
