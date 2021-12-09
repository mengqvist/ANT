# The ambiguous nucleotide tool (ANT) is a free and open source tool aimed at
# generating and analysing degenerate codons to support research in protein engineering, directed evolution and synthetic biology.

# Copyright (C) 2015  Martin Engqvist |
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LICENSE:
# This file is part of ANT.
#
# ANT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# ANT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Get source code at: https://github.com/mengqvist/ANT
#


def one_to_three(one_letter):
    """
    Convert a one letter code amino acid to a three letter code.
    """
    assert one_letter.upper() in "FLSYCWPHERIMTNKVADQG*U", "Error, %s is not a valid amino acid" % one_letter

    AA = {
        "I": "Ile",
        "V": "Val",
        "L": "Leu",
        "F": "Phe",
        "C": "Cys",
        "M": "Met",
        "A": "Ala",
        "G": "Gly",
        "T": "Thr",
        "W": "Trp",
        "S": "Ser",
        "Y": "Tyr",
        "P": "Pro",
        "H": "His",
        "E": "Glu",
        "Q": "Gln",
        "D": "Asp",
        "N": "Asn",
        "K": "Lys",
        "R": "Arg",
        "*": "***",
        "U": "Uaa",
    }  # Uaa = unnatural amino acid
    return AA[one_letter.upper()]


def three_to_one(three_letter):
    """
    Convert a three letter code amino acid to a one letter code.
    """
    assert three_letter.upper() in [
        "ILE",
        "VAL",
        "LEU",
        "PHE",
        "CYS",
        "MET",
        "ALA",
        "GLY",
        "THR",
        "TRP",
        "SER",
        "TYR",
        "PRO",
        "HIS",
        "GLU",
        "GLN",
        "ASP",
        "ASN",
        "LYS",
        "ARG",
        "***",
        "UAA",
    ], (
        "Error, %s is not a valid amino acid" % three_letter
    )

    AA = {
        "ILE": "I",
        "VAL": "V",
        "LEU": "L",
        "PHE": "F",
        "CYS": "C",
        "MET": "M",
        "ALA": "A",
        "GLY": "G",
        "THR": "T",
        "TRP": "W",
        "SER": "S",
        "TYR": "Y",
        "PRO": "P",
        "HIS": "H",
        "GLU": "E",
        "GLN": "Q",
        "ASP": "D",
        "ASN": "N",
        "LYS": "K",
        "ARG": "R",
        "***": "*",
        "UAA": "U",
    }
    return AA[three_letter.upper()]


def one_to_full(one_letter):
    """
    Convert one-letter amino acid code to full amino acid name.
    """
    assert one_letter.upper() in "FLSYCWPHERIMTNKVADQG*U", "Error, %s is not a valid amino acid" % one_letter
    AA = {
        "F": "Phenylalanine",
        "L": "Leucine",
        "S": "Serine",
        "Y": "Tyrosine",
        "*": "Stop",
        "C": "Cysteine",
        "W": "Tryptophan",
        "P": "Proline",
        "H": "Histidine",
        "Q": "Glutamine",
        "R": "Arginine",
        "I": "Isoleucine",
        "M": "Methionine",
        "T": "Threonine",
        "N": "Asparagine",
        "K": "Lysine",
        "V": "Valine",
        "A": "Alanine",
        "D": "Aspartic acid",
        "E": "Glutamic acid",
        "G": "Glycine",
        "U": "Unnatural AA",
    }
    return AA[one_letter.upper()]


def full_to_one(full):
    """
    Convert full amino acid name to one-letter amino acid code.
    """
    assert full.lower() in [
        "phenylalanine",
        "leucine",
        "serine",
        "tyrosine",
        "stop",
        "cysteine",
        "tryptophan",
        "proline",
        "histidine",
        "glutamine",
        "arginine",
        "isoleucine",
        "methionine",
        "threonine",
        "asparagine",
        "lysine",
        "valine",
        "alanine",
        "aspartic acid",
        "glutamic acid",
        "glycine",
        "unnatural aa",
    ], (
        "Error, %s is not a valid amino acid" % full
    )

    AA = {
        "phenylalanine": "F",
        "leucine": "L",
        "serine": "S",
        "tyrosine": "Y",
        "stop": "*",
        "cysteine": "C",
        "tryptophan": "W",
        "proline": "P",
        "histidine": "H",
        "glutamine": "Q",
        "arginine": "R",
        "isoleucine": "I",
        "methionine": "M",
        "threonine": "T",
        "asparagine": "N",
        "lysine": "K",
        "valine": "V",
        "alanine": "A",
        "aspartic acid": "D",
        "glutamic acid": "E",
        "glycine": "G",
        "unnatural aa": "U",
    }
    return AA[full.lower()]


def three_to_full(three_letter):
    """
    Convert amino acid three letter code to full amino acid names.
    """
    return one_to_full(three_to_one(three_letter))


def full_to_three(full):
    """
    Convert full amino acid names to three letter code.
    """
    return one_to_three(full_to_one(full))


def count_aa(seq):
    """
    Count occurrences of all amino acids in sequence. Return as dictionary.
    """
    seq = seq.upper()
    assert (
        all([s in "IVLFCMAGTWSYPHEQDNKR*U" for s in seq]) is True
    ), "Error, unknown amino acids %s in sequence: %s" % (
        str([s for s in seq if s not in "XIVLFCMAGTWSYPHEQDNKR*"]),
        seq,
    )

    AA = {
        "I": seq.count("I"),
        "V": seq.count("V"),
        "L": seq.count("L"),
        "F": seq.count("F"),
        "C": seq.count("C"),
        "M": seq.count("M"),
        "A": seq.count("A"),
        "G": seq.count("G"),
        "T": seq.count("T"),
        "W": seq.count("W"),
        "S": seq.count("S"),
        "Y": seq.count("Y"),
        "P": seq.count("P"),
        "H": seq.count("H"),
        "E": seq.count("E"),
        "Q": seq.count("Q"),
        "D": seq.count("D"),
        "N": seq.count("N"),
        "K": seq.count("K"),
        "R": seq.count("R"),
        "*": seq.count("*"),
        "U": seq.count("U"),
    }
    return AA
