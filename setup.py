from pathlib import Path
import glob
from setuptools import setup, find_packages

# The directory containing this file
LOCAL_PATH = Path(__file__).parent

# The text of the README file
README_FILE = (LOCAL_PATH / "README.rst").read_text()

# Build (recursive) list of resource files
resource_files = []
for directory in glob.glob("resources/*/"):
    files = glob.glob(directory + "*")
    resource_files.append((directory, files))

# Load requirements, so they are listed in a single place
with open("requirements.txt") as fp:
    install_requires = [dep.strip() for dep in fp.readlines()]

# This call to setup() does all the work
setup(
    author_email="tresoldi@shh.mpg.de",
    author="Tiago Tresoldi",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Software Development :: Libraries",
    ],
    data_files=[("docs", glob.glob("docs/*"))] + resource_files,
    description="Library for multiple asymmetric alignments on different alphabets",
    entry_points={"console_scripts": ["malign=malign.__main__:main"]},
    include_package_data=True,
    install_requires=install_requires,
    keywords=["alignment", "sequence alignment", "multiple alphabet"],
    license="MIT",
    long_description_content_type="text/x-rst",
    long_description=README_FILE,
    name="malign",
    packages=find_packages(where="src"),  # ["malign", "resources", "docs"],
    package_dir={"": "src"},  # , "resources":"..", "docs":".."},
    project_urls={"Documentation": "https://malign.readthedocs.io"},
    test_suite="tests",
    tests_require=[],
    url="https://github.com/tresoldi/malign",
    version="0.2",  # remember to sync with __init__.py
    zip_safe=False,
)
