from setuptools import find_packages, setup
import versioneer

setup(
    name='wgs_qc_utils',
    description='Whole genome qc',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
