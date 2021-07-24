from setuptools import find_packages, setup

setup(
    name='prppi-cli',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'pyyaml',
        'biopython',
        'osprey'
    ],
    entry_points={
        'console_scripts': [
            'prppi=prppi_cli.cli: cli']
    },
)
