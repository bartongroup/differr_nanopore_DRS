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
