from pathlib import Path
from setuptools import setup


# Loads _version.py module without importing the whole package.
def get_version_and_cmdclass(pkg_path):
    import os
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location(
        "version",
        os.path.join(pkg_path, "_version.py"),
    )
    module = module_from_spec(spec)  # type: ignore
    spec.loader.exec_module(module)  # type: ignore
    return module.__version__, module.get_cmdclass(pkg_path)


version, cmdclass = get_version_and_cmdclass("xenopict")

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="xenopict",
    version=version,
    cmdclass=cmdclass,
    description="Library for publication quality depictions of small molecules.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="S. Joshua Swamidass",
    author_email="swamidass@wustl.edu",
    url="https://github.com/swamidasslab/xenopict/",
    packages=["xenopict"],
    install_requires=[
        "matplotlib>=3.5",
        "colorcet",
        "numpy",
        "simplejson",
        "decorator",
        "rdkit",
        "shapely"
    ],
    extra_requires={
        "pdf": ["cairosvg"],
        "pandas": ["pandas"],
        "ipython": ["ipython"],
        "all": ["cairosvg", "pandas", "ipython", "colorcet"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Framework :: IPython",
        "Framework :: Jupyter",
        "Framework :: Jupyter :: JupyterLab :: Extensions",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Multimedia :: Graphics",
        "Topic :: Multimedia :: Graphics :: Editors :: Vector-Based",
        "Topic :: Multimedia :: Graphics :: Presentation",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization",
        "Typing :: Typed",
    ],
)
