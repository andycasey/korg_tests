from distutils.core import setup
from setuptools import find_packages

setup(
    name='grok',
    version='0.1',
    description='Compare spectral synthesis codes',
    author='Andy Casey',
    author_email='andrew.casey@monash.edu',
    url='https://github.com/andycasey/grok',
    packages=find_packages(where="src"),
    package_dir={
        "grok": "src/grok",
    },
    package_data={
        "grok": ["data/photospheres/sun.mod", "data/transitions/turbospec.20180901.*"],
    },
    data_files=[
        ("grok/synthesis/korg", 
            [
                "src/grok/synthesis/korg/template.jl"
            ]
        ),
        ("grok/synthesis/moog", ["src/grok/synthesis/moog/moog_synth.template"]),
        ("grok/synthesis/turbospectrum", 
            [
            "src/grok/synthesis/turbospectrum/bsyn.template",
            "src/grok/synthesis/turbospectrum/babsma.template",
            ]
        )
    ],
    scripts=['bin/grok']
)