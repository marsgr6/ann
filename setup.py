import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="attractor_network",
    version="1.3",
    author="Mario Gonzalez Rodriguez",
    author_email="marsgr6@gmail.com",
    description="Attractor Network Implementation",
    long_description=long_description,
    url="https://gitlab.com/kbjimenes/periodic_neural_activity",
    install_requires=['random', 'scipy', 'numpy'], 
    python_requires='>=3',
    packages=setuptools.find_packages(),
)