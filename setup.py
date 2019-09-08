from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='scopetools',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Single Cell Omics Preparation Entity Tools",
    author="luyang",
    author_email='luyang@singleronbio.com',
    url='https://bitbucket.org/luyangatsingleronbio/scopetools',
    packages=['scopetools'],
    entry_points={
        'console_scripts': [
            'scopetools=scopetools.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='scopetools',
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)
