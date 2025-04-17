from setuptools import setup, find_packages

setup(
    name="MultiRIN",
    version="0.1",
    packages=find_packages(),
    scripts=[
        'generate_single.py',
        'generate_multi.py'
    ]
)