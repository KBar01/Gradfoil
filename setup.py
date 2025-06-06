import os
import subprocess
import shutil
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

class BuildWithCMake(build_py):
    def run(self):
        
        # 1. Build C++ executables using CMake
        build_dir = os.path.abspath("build")
        os.makedirs(build_dir, exist_ok=True)
        
        # builds the executables using cmake
        subprocess.check_call(["cmake", ".."], cwd=build_dir)
        subprocess.check_call(["cmake", "--build", "."], cwd=build_dir)

        # 2. Copy executables into the Python package
        output_bin_dir = os.path.join("gradfoil", "bin")
        os.makedirs(output_bin_dir, exist_ok=True)
        shutil.copy(os.path.join(build_dir, "CFoil_fwd"), output_bin_dir)
        shutil.copy(os.path.join(build_dir, "CFoil_AD"), output_bin_dir)

        # 3. Continue with the normal Python build
        super().run()


setup(
    name="gradfoil",
    version="0.1.0",
    packages=find_packages(),
    cmdclass={
        'build_py': BuildWithCMake,
    },
    include_package_data=True,
    package_data={
        'gradfoil': ['bin/*'],
    },
    install_requires=[],
    python_requires=">=3.7",
    author="Your Name",
    description="Python interface to CFoil tools",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
)