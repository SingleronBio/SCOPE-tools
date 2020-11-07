from setuptools import setup

import versioneer

requirements = [
    # package requirements go here
    'cutadapt==2.10',
    'pysam==0.16.0.1',
    'scipy==1.5.3',
    'numpy==1.19.4',
    'pandas==1.1.4',
    'jinja2==2.11.2',
    'matplotlib==3.3.2',
    'click==7.1.2',
    'scanpy==1.6.0',
    'leidenalg==0.8.2',
    'louvain==0.6.1'
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
