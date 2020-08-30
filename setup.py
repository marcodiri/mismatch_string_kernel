import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mismatch_string_kernel",
    version="0.0.2",
    author="Marco Di Rienzo",
    author_email="diridevelops@gmail.com",
    description="Python3 implementation of mismatch string kernel",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/marcodiri/mismatch_string_kernel",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
