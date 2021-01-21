from setuptools import setup
setup(
    name='prism',
    version='0.0.2',
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    packages=["src", "src.prism"],
    include_package_data=True,
    install_requires=["click"],
    entry_points="""
    [console_scripts]
    prism=src.prism.main:cli
    """,
)

