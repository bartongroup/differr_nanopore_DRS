from setuptools import setup


setup(
    name='differr',
    version='0.1',
    description=(
        'identify sites with differential mismatches caused by m6A in ONT DRS data'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'differr = differr.main:differr',
            'phaserr = differr.phasing:identify_three_prime_phasing',
            'polya_test = differr.cluster:find_polya_clusters',
        ]
    },
    packages=[
        'differr',
    ],
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'click',
        'joblib',
        'pysam',
        'statsmodels'
    ],
)
