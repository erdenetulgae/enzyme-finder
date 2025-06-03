from setuptools import setup, find_packages

setup(
    name='enzyme_finder',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'requests',
        'brendapy',
    ],
    entry_points={
        'console_scripts': [
            'get-taxid=laminarinase.get_taxid:main',
            'fetch-refs=laminarinase.fetch_refs_uniprot:main',
            'extract-laminarinase=laminarinase.extract_ref_laminarinase:main',
        ]
    }
)
