#ANT
###Ambiguous Nucleotide Tool (ANT) is a software for generating and evaluating degenerate codons for natural and expanded genetic codes.

=====

ANT comes with a GUI and a defined API. To use the API a working installation of Python is needed. To use the GUI a working installation of wxPython is additionally needed.
To use the API simply import ANT.py into your python shell.
To use the gui simply run ANT_GUI.py.

GUI:
![ANT GUI](/Screenshot.png?raw=true "ANT GUI")


API Commands:

Instantiation can be done using a list of amino acids:
```
>>> import ANT
>>> codon_object = ANT.DegenerateCodon(input=['S', 'T', 'A', 'G'], table=1)
```


Alternatively, instantiation can be done using a degerate codon:
```
>>> import ANT
>>> codon_object = ANT.DegenerateCodon(input='RSC', table=1)
```

To retreive the degenerate triplet:
```
>>>codon_object.getTriplet()
‘RSC’
```

To get all encoded amino acids (targets and off-targets):
```
>>> codon_object.getEncoded()
['S', 'T', 'A', 'G']
```


To get target amino acids:
```
>>> codon_object.getTarget()
['S', 'T', 'A', 'G']
```

To get off-target amino acids:
```
>>> codon_object.getOffTarget()
[]
```

To get which amino acids may still be chosen without further off-targets:
```
>>> codon_object.getPossible()
['C', 'R']
```

To get which genetic code was used:
```
>>> codon_object.getTable()
1
```

To get which codons are specified by the degenerate codon:
```
>>> codon_object.getCodons()
['ACC', 'GCC', 'AGC', 'GGC'] 
```


To get how many times each amino acid is encoded:
```
>>> codon_object.getCodonsPerAA()
{'*': 0, 'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': 0, 'I': 0, 'H': 0, 'K': 0,
 'M': 0, 'L': 0, 'N': 0, 'Q': 0, 'P': 0, 'S': 1, 'R': 0, 'U': 0, 'T': 1, 'W': 0,
 'V': 0, 'Y': 0}
```


To get altenative codons with the same number of off-target amino acids:
```
>>> codon_object.getAlternatives()
[['RSC', 'A', 'S', 'T', 'G'], ['RST', 'A', 'S', 'T', 'G']]
```


To get a more extensive list of altenative codons (some with more off-target amino acids):
```
>>> codon_object.getExtendedAlternatives()
[['RSC', 'A', 'S', 'T', 'G'], ['RST', 'A', 'S', 'T', 'G'], ['DSC', 'A', 'S', 'T'
, 'G', 'C'], ['DST', 'A', 'S', 'T', 'G', 'C'], ['DSA', 'A', 'S', 'T', 'G', '*',
'R'], ['DSG', 'A', 'S', 'T', 'G', 'R', 'W']]
```


To get a full report:
```
>>> codon_object.getReport()
Degenerate codon: RSC
Genetic code: 1
Real codons: ['ACC', 'GCC', 'AGC', 'GGC']
Target amino acids: ['A', 'S', 'T', 'G']
Off-target amino acids: []
Amino acids that can be added w/o further off-targets: ['C', 'R']
Codons for each amino acid: {'*': 0, 'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': 0, 'I': 0, 'H': 0, 'K': 0, 'M': 0, 'L': 0, 'N': 0, 'Q': 0, 'P': 0, 'S': 1, 'R': 0, 'U': 0, 'T': 1, 'W': 0, 'V': 0, 'Y': 0}
Library size (number of codons): 4
Clones to screen for 95% confidence: 11
Alternate codons with same number of off-target amino acids: [['RSC', 'A', 'S', 'T', 'G'], ['RST', 'A', 'S', 'T', 'G']]
```

