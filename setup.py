from setuptools import setup, find_packages

exec(open('apples/version.py').read())
setup(
        name='apples',    # This is the name of your PyPI-package.
        version=__version__,    # Update the version number for new releases
        # The name of your script, and also the command you'll be using for calling it
        # Also other executables needed
        scripts=['run_apples.py','build_applesdtb.py',],
        description='APPLES: a distance-based phylogenetic placement tool',
        long_description='APPLES stands for Accurate Phylogenetic Placement\
            with LEast Squares and addresses the problem of phylogenetic placement\
            of DNA and protein sequences into an already existing reference tree.',
        long_description_content_type='text/plain',
        url='https://github.com/balabanmetin/apples',
        author='Metin Balaban',
        author_email='balaban@ucsd.edu',
        packages=find_packages(),
        zip_safe = False,
        install_requires=['numpy','treeswift','treecluster'],
        include_package_data=True
)
