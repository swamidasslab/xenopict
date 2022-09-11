
from setuptools import setup


from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='PlotDot',
    version='0.0.5',
    description='Utility for plotting perceputally uniform dots.',
    long_description=long_description,
    long_description_content_type='text/markdown',     
	author='S. Joshua Swamidass',
    author_email='swamidass@wustl.edu',
    url='https://github.com/swamidasslab/plotdot/',
    packages=['plotdot'],
    install_requires=["six>=1.13.0", "colorio", "numpy"],
    extra_requires={"rdk": [
        "rdk", "shapely"
    ]}
)
