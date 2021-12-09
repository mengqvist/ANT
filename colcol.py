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


# This file mainly deals with color conversions, color transformations, and generating color scales.

import re


def is_rgb(input):
    """
    Check whether input is a valid RGB color.
    Return True if it is, otherwise False.
    """
    if len(input) == 3 and type(input) == tuple:
        if 0 <= input[0] <= 255 and 0 <= input[1] <= 255 and 0 <= input[2] <= 255:
            return True
        else:
            return False
    else:
        return False


def rgb_to_hex(rgb):
    """
    Convert RGB colors to hex.
    Input should be a tuple of integers (R, G, B).
    Output is a string.
    """
    # make sure input is ok
    assert is_rgb(rgb) is True, "Error, %s is not a valid RGB color." % str(rgb)

    # make conversion
    return "#%02x%02x%02x" % rgb


def is_hex(input):
    """
    Check whether an input string is a valid hex value.
    Return True if it is, otherwise False.
    """
    if type(input) is not str:
        return False

    regular_expression = re.compile(
        r"""^							#match beginning of string
										[#]{1} 						#exactly one hash
										[0-9a-fA-F]{6}				#exactly six of the hex symbols  0 to 9, a to f (big or small)
										$							#match end of string
										""",
        re.VERBOSE,
    )

    if regular_expression.match(input) == None:
        return False
    else:
        return True


def hex_to_rgb(hex):
    """
    Convert a hex color to RGB.
    Input should be a string.
    Output is a tuple of integers (R, G, B).
    """
    # make sure input is ok
    assert is_hex(hex) is True, "Error, %s is not a valid hex color." % hex

    # make the conversion
    hex = hex.lstrip("#")
    return tuple([int(hex[s : s + 2], 16) for s in range(0, len(hex), 2)])


def scale(col1, col2, white_mid=False):
    """
    Function makes a color scale from 0 to 100 using the supplied colors.
    The variable white_mid is boolean and determines whether the colors
    should transition over white in the middle or just smoothly run into each other.
    The function returns a dictionary of integer keys and hex color values corresponding to the scores 0 to 100.
    """
    assert is_hex(col1) or is_rgb(col1), "Error, the first color is not valid."
    assert is_hex(col2) or is_rgb(col2), "Error, the second color is not valid."
    assert type(white_mid) is bool, "Error, white_mid must be a boolean."

    if is_hex(col1):
        col1 = hex_to_rgb(col1)
    if is_hex(col2):
        col2 = hex_to_rgb(col2)

    color_dict = {}
    if white_mid is False:
        # determine how many points R, G and B should change with each score
        R_diff = (col2[0] - col1[0]) / 100.0
        G_diff = (col2[1] - col1[1]) / 100.0
        B_diff = (col2[2] - col1[2]) / 100.0

        # set starting values
        R, G, B = col1
        for i in range(101):
            color_dict[i] = rgb_to_hex((R + int(R_diff * i), G + int(G_diff * i), B + int(B_diff * i)))

    elif white_mid is True:
        first_half = scale((col1), (255, 255, 255))
        for key in range(0, 100, 2):
            color_dict[key / 2] = first_half[key]

        second_half = scale((255, 255, 255), col2)
        for key in range(0, 102, 2):
            color_dict[50 + key / 2] = second_half[key]

    return color_dict


def mix_colors(col1, col2):
    """
    Mix two colors and return the result.
    The input colors can be either RGB or hex.
    It is also possible to use one RGB value and one hex value.
    """
    color_dict = scale(col1, col2, white_mid=False)
    return color_dict[50]


def NextRGB(color=(0, 0, 0)):
    """
    Function for generating unique RGB colors.
    The input is a tuple of RGB colors (for example (124,1,34) and the method returns the "next" color.
    When R reaches 255 one is added to G and R is reset.
    When R and G both reach 255 one is added to B and R and G are reset.
    This should generate over 1.6 million colors (255*255*255)
    """
    assert is_rgb(color), "Error, the input must be a tuple of three integers between 0 and 255"
    R, G, B = color

    if R == 255 and G == 255 and B == 255:
        raise ValueError("R, G and B all have the value 255, no further colors are available.")
    elif R == 255 and G == 255:
        R = 0
        G = 0
        B += 1
    elif R == 255:
        R = 0
        G += 1
    else:
        R += 1
    return (R, G, B)
