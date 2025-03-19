from setuptools import setup, find_packages

setup(
    name="libraries x Compton",  # Name of your package
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",  # Add any dependencies your project requires
        "matplotlib"
    ],
    include_package_data=True,
    description="Some lib for Compton data analysis in university laboratory",
    author="Andrea Morandi",
    author_email="morandi.andrea2002@gmail.com",
)