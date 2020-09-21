from setuptools import find_packages, setup
import versioneer

package_data={'ideogram_data': ['wgs_qc_utils/wgs_qc_utils/reader/ideogram/ideogram.txt']},


setup(
    name='wgs_qc_utils',
    description='Whole genome qc',
    packages=find_packages(),
    package_data={'': ['ideogram.txt']},
    include_package_data=True,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
