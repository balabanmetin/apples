from setuptools import setup,find_packages

setup(
        name='apples',    # This is the name of your PyPI-package.
        version='1.3.0',    # Update the version number for new releases
        scripts=['run_apples.py',], # The name of your scipt, and also the command you'll be using for calling it
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
        install_requires=['numpy','treeswift'],
        include_package_data=True
)
