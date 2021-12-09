from distutils.core import setup

setup(
    name="ANT",
    version="0.1",
    description="Ambiguous Nucleotide Tool (ANT), software for generating and evaluating degenerate codons for natural and expanded genetic codes.",
    author="Martin Engqvist",
    author_email="martin_engqvist@hotmail.com",
    py_modules=["dna", "protein", "ANT", "base_class", "pyperclip", "colcol"],
    scripts=["ANT.py", "ANT_GUI.py"],
)
