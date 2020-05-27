from setuptools import setup

import versioneer

requirements = [
    # package requirements go here
    'cutadapt',
    'pysam',
    'scipy>=1.3',
    'numpy',
    'pandas',
    'jinja2',
    'matplotlib',
    'click==7.1.2',
    'scanpy==1.5.0',
    'leidenalg',
    'louvain'
]

setup(
    name='SCOPE-tools',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Single Cell Omics Preparation Entity tools",
    license="Apache",
    author="Singleron Biotechnologies",
    author_email='luyang@singleronbio.com',
    url='https://github.com/SingleronBio/SCOPE-tools',
    packages=['scopetools'],
    entry_points={
        'console_scripts': [
            'scope=scopetools.cli:cli'
        ]
    },
    install_requires=requirements,
    include_package_data=True,
    keywords='SCOPE-tools',
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ]
)
