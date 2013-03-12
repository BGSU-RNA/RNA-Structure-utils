from distutils.core import setup

setup(
    name='RNAStructure',
    version='0.0.12',
    author='Blake A. Sweeney',
    author_email='bsweene@bgsu.edu',
    packages=['rnastructure', 'rnastructure.primary', 'rnastructure.secondary',
              'rnastructure.util'],
    url='',
    license='LICENSE.txt',
    description='Some tools to parse RNA Structure',
    long_description=open('README.mkd').read(),
)
