import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
setup(
    name='cardio-model', version='0.1',
    author='Mateus Felismino, Hugo Santos',
    license="MIT",
    description='Testing',
    packages=find_packages())