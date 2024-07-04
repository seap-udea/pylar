##################################################################
# License: http://github.com/seap-udea/pylar                     #
##################################################################
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='ipylar',
    author="Jorge I. Zuluaga, Ruben D. Molina, Juan F. Salazar and Jesus D. Gomez-Velez",
    author_email="jorge.zuluaga@udea.edu.co",
    description="Python utilities for the LAR model (Land Atmospheric Reservoir)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pypi.org/project/ipylar",
    keywords='hydrology, climate change',
    license='MIT',

    # ######################################################################
    # CLASSIFIER
    # ######################################################################
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
    version='1.0.6',

    # ######################################################################
    # FILES
    # ######################################################################
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    
    # ######################################################################
    # ENTRY POINTS
    # ######################################################################
    entry_points={
        'console_scripts': ['install=ipylar.install:main'],
    },

    # ######################################################################
    # TESTS
    # ######################################################################
    test_suite='nose.collector',
    tests_require=['nose'],

    # ######################################################################
    # DEPENDENCIES
    # ######################################################################
    install_requires=['matplotlib','pandas'],
    # 'mpmath','sympy','scikit-learn','openpyxl','statsmodels','ipywidgets'

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={"": ["data/*"]},
)
