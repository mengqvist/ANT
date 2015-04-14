#!/usr/bin/env python


#The ambiguous nucleotide tool (ANT) is a free and open source tool aimed at
#generating and analysing degenerate codons to support research in protein engineering, directed evolution and synthetic biology.

#Copyright (C) 2015  Martin Engqvist | 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LICENSE:
#This file is part of ANT.
#
#ANT is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
# 
#ANT is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Library General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software Foundation,
#Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#Get source code at: https://github.com/mengqvist/ANT
#


import random
import re


def Translate(DNA, table=1):
	"""
	Returns protein sequence from DNA string input.
	The table variable specifies which codon table should be used.
	table defaults to the standard codon table 1
	"""
	assert type(DNA) == str or type(DNA) == unicode, 'Error, input sequence must be a string or unicode'
	codons = CodonTable(table).getCodons()

	protein = []
	DNA = DNA.upper()
	DNA = DNA.replace('\n', '') #get rid of row breaks
	DNA = DNA.replace('U', 'T')

	for i in range(len(DNA)):
		if i%3==0:
			if i+3>len(DNA):
				pass
			elif any(DNA[i:(i+3)] in s for s in codons['F']):
				protein.append('F')
			elif any(DNA[i:(i+3)] in s for s in codons['L']):
				protein.append('L')
			elif any(DNA[i:(i+3)] in s for s in codons['S']):
				protein.append('S')
			elif any(DNA[i:(i+3)] in s for s in codons['Y']):
				protein.append('Y')
			elif any(DNA[i:(i+3)] in s for s in codons['*']):
				protein.append('*')
			elif any(DNA[i:(i+3)] in s for s in codons['C']):
				protein.append('C')
			elif any(DNA[i:(i+3)] in s for s in codons['W']):
				protein.append('W')
			elif any(DNA[i:(i+3)] in s for s in codons['P']):
				protein.append('P')
			elif any(DNA[i:(i+3)] in s for s in codons['H']):
				protein.append('H')
			elif any(DNA[i:(i+3)] in s for s in codons['E']):
				protein.append('E')
			elif any(DNA[i:(i+3)] in s for s in codons['R']):
				protein.append('R')
			elif any(DNA[i:(i+3)] in s for s in codons['I']):
				protein.append('I')
			elif any(DNA[i:(i+3)] in s for s in codons['M']):
				protein.append('M')
			elif any(DNA[i:(i+3)] in s for s in codons['T']):
				protein.append('T')
			elif any(DNA[i:(i+3)] in s for s in codons['N']):
				protein.append('N')
			elif any(DNA[i:(i+3)] in s for s in codons['K']):
				protein.append('K')
			elif any(DNA[i:(i+3)] in s for s in codons['V']):
				protein.append('V')
			elif any(DNA[i:(i+3)] in s for s in codons['A']):
				protein.append('A')
			elif any(DNA[i:(i+3)] in s for s in codons['D']):
				protein.append('D')
			elif any(DNA[i:(i+3)] in s for s in codons['Q']):
				protein.append('Q')
			elif any(DNA[i:(i+3)] in s for s in codons['G']):
				protein.append('G')
			elif any(DNA[i:(i+3)] in s for s in codons['U']): #special case allowing for unnatural AA
				protein.append('U')
			else:
				raise ValueError, '"%s" is not a valid codon' % DNA[i:(i+3)]
	return ''.join(protein)	


	
def GetCodons(AA, table=1, separate=False, exclude=False):
	'''
	Get the codons for a specified AA. Returns a list of strings.
	The variable table specifies which codon table should be used.
	table defaults to the standard codon table 1
	The separate variable determines whether codons for one amino acid with dissimilar first two nucleotides should be seperated out.
	For example if separate=False the codons for L are 	['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'].
	If separate=True they are split up as L = ['TTA', 'TTG'] and L2 = ['CTT', 'CTC', 'CTA', 'CTG'].
	exclude deterimes whether user-defined codons should be excluded or not. Valid values are True and False.
	'''
	AA = AA.upper()
	assert len(AA) == 1, 'Error, function takes a single amino acid as input'
	assert AA in 'FLSYCWPHERIMTNKVADQG*U', 'Error, %s is not a valid amino acid' % str(AA)

	codons = CodonTable(table, exclude).getCodons(separate)

	return codons[AA]	
	




def combine(input_list, max_total=50000):
	'''
	Takes a list of lists (nucleotides or amino acids) and makes every combination of these while retaining the internal order.
	For example [['A'], ['T','A','G'], ['C', 'A']] will result in six sequences: ['ATC', 'AAC', 'AGC', 'ATA', 'AAA', 'AGA']
	max_total puts a limit on the maximum number of sequences that can be computed. The number tends to increase explosively....
	The output is a list of the resulting sequences.
	'''

	#make sure the input will not result in too many sequences
	list_of_nums = [len(x) for x in input_list]
	total = reduce(lambda x, y: x*y, list_of_nums)
	assert total < max_total, 'The sequence "%s" would result in %s total sequences, this is above the set maximum of %s sequences.' % (string, total, max_total)

	#now combine them
	#first copy the output list the number of times that there are nucleotides in the current position
	output = []
	for pos in input_list:
		output_len = len(output)
		if output_len == 0:
			output.extend(pos)
			output_len = len(output)
		else:
			output.extend(output*(len(pos)-1)) #duplicate output the number of times that there are new nucleotides to add
			for i in range(0, len(pos)): #for every nucleotide to be added
				for j in range(0, output_len): #add that nucleotide the number of times that the prevous output was long (before the last duplication)
					output[j+i*output_len] += pos[i]
	return output #return the list

	
	
def listupper(t):
	'''
	Capitalizes all strings in nested lists. 
	'''
	if isinstance(t, list):
		return [listupper(s) for s in t]
	elif isinstance(t, str):
		return t.upper()
	else:
		return t
		
	
	
def UnAmb(input_string):
	'''
	Converts an ambiguous nucleotide sequence to a list of sequences containing only A, T, C and G (as appropriate).
	'''
	assert type(input_string) is str, 'Error, the input has to be a string.'
	input_string = input_string.upper()	

	pos_list = []
	for letter in input_string:
		assert letter in 'NMRWSYKVHDBGATC', 'Error, "%s" is not a valid ambigous nucleotide.' 
		if 'T' == letter:
			pos_list.append(['T'])

		elif 'C' == letter:
			pos_list.append(['C'])
			
		elif 'A' == letter:
			pos_list.append(['A'])

		elif 'G' == letter:
			pos_list.append(['G'])

		elif 'M' == letter:
			pos_list.append(['C','A'])

		elif 'Y' == letter:
			pos_list.append(['T','C'])

		elif 'K' == letter:
			pos_list.append(['T','G'])

		elif 'S' == letter:
			pos_list.append(['C','G'])

		elif 'W' == letter:
			pos_list.append(['T','A'])

		elif 'R' == letter:
			pos_list.append(['A','G'])

		elif 'H' == letter:
			pos_list.append(['T','C','A'])

		elif 'V' == letter:
			pos_list.append(['C','A','G'])

		elif 'D' == letter:
			pos_list.append(['T','A','G'])

		elif 'B' == letter:
			pos_list.append(['T','C','G'])

		elif 'N' == letter: 
			pos_list.append(['T','C','A','G'])

	return combine(pos_list) #call combine function and return the result as a list of strings

	
	
	
def commonNuc(nuc_list, greedy=False): 
	"""
	This function takes a list of lists and finds all degenerate symbols that represent 
	at least one nucleotide from each of the lists.
	The variable "greedy" determines whether the algorithm is greedy or not.
	
	With greedy=False	
	An example input is: [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']].
	T and C are both present in all lists, therefore, both 'T' and 'C' are acceptable returned as ['T', 'C'].

	With greedy=False
	Another example input is: [['G'], ['T'], ['T']].
	In this case either G or T is present in all lists, therefore the only acceptable output is ['K'] (ambiguous nucleotide for G and T). 

	
	With greedy=True 
	For the input: [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']]
	The greedy output includes all degenerate nucleotides that contain the desired regular nucleotides: 
	['C', 'T', 'Y', 'K', 'M', 'S', 'W', 'H', 'V', 'D', 'B', 'N']
	
	With greedy=True
	For the input: [['G'], ['T'], ['T']]
	The greedy output is: ['K', 'D', 'B', 'N']
	"""
	nuc_list = listupper(nuc_list)
	output = []
	
	if all(['A' in s for s in nuc_list]):
		output.append('A')
		
	if all(['G' in s for s in nuc_list]):
		output.append('G')

	if all(['C' in s for s in nuc_list]):
		output.append('C')

	if all(['T' in s for s in nuc_list]):
		output.append('T')

	if greedy is False and len(output)>0:
		return output
	
	
	if all(['C' in s or 'T' in s for s in nuc_list]):
		output.append('Y')
		
	if all(['G' in s or 'T' in s for s in nuc_list]):
		output.append('K')

	if all(['A' in s or 'C' in s for s in nuc_list]):
		output.append('M')

	if all(['C' in s or 'G' in s for s in nuc_list]):
		output.append('S')
		
	if all(['A' in s or 'T' in s for s in nuc_list]):
		output.append('W')
		
	if all(['A' in s or 'G' in s for s in nuc_list]):
		output.append('R')
		
	if greedy is False and len(output)>0:
		return output
		
		
	if all(['C' in s or 'T' in s or 'A' in s for s in nuc_list]):
		output.append('H')
		
	if all(['C' in s or 'A' in s or 'G' in s for s in nuc_list]):
		output.append('V')

	if all(['T' in s or 'A' in s or 'G' in s for s in nuc_list]):
		output.append('D')

	if all(['C' in s or 'T' in s or 'G' in s for s in nuc_list]):
		output.append('B')
		
	if greedy is False and len(output)>0:
		return output
		

	if all(['C' in s or 'T' in s or 'A' in s or 'G' in s for s in nuc_list]):
		output.append('N')
		
	return output

		
				
class CodonTable:
	'''
	A DNA codon object.
	Used to retrieve codon tables and codons for specified codon tables.	
	Pass a valid integer value when instantiating to choose which codon table to use.
	If exclude=True then certain codons will be excluded from the lists.
	'''
	def __init__(self, number, exclude=False):
		self.code = False
		self.code_num = False
		self.table = False
		self.codons = False

		#variable to hold user-defined data (from settings file)
		self.settings = dict()

		self.readSettings() #read settings file to get user-defined codon table and list of codons to exclude
		self.setTable(number) #get the specified codon table (returned as list of strings)
		self.setCodons(exclude) #convert the codon table information to codons
		

	def readSettings(self):
		'''
		Method which reads the settings file to get the user-defined codon table and which (if any) codons should get excluded.
		These are stored and used when computing degenerate codons (self.remove) or if codon table 1001 is chosen (the other variables). 
		'''

		#get the settings
		execfile('./settings.txt', self.settings)

		#make sure the settings are ok
		assert type(self.settings['code']) is str, 'Error, the Review the settings.txt file.'

		assert type(self.settings['codons_to_exclude']) is list, 'Error, the codons for exclusion must be in a list. Review the settings.txt file.' #make sure it is a list

		if len(self.settings['codons_to_exclude']) != 0: #if there are items in it
			self.settings['codons_to_exclude'] = [s.upper() for s in self.settings['codons_to_exclude']] #make uppercase
			for item in self.settings['codons_to_exclude']:
				assert re.match('^[ATCG]{3}$', item) != None, 'Error, %s is not a valid DNA codon to exclude. Please review the settings.txt file.' % item

		for aa in 'FLSYCWPHERIMTNKVADQG*U':
			assert aa in self.settings['AAs'], 'Error, the amino acid %s has not been specified. Review the settings.txt file.'
 
		assert self.settings['Base1'] == 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG', 'Error, the Base1 field is not correct. Review the settings.txt file.'
		assert self.settings['Base2'] == 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG', 'Error, the Base2 field is not correct. Review the settings.txt file.'
		assert self.settings['Base3'] == 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG', 'Error, the Base3 field is not correct. Review the settings.txt file.'


	def setTable(self, number):
		'''
		Find information for specified genetic code and use for downstream methods.
		Method is not intended for direct use.
		'''
		number = int(number)
		
		if number == 1:
			#Genetic Code [1]
			code = "Standard Code (transl_table=1)"
			AAs  = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "---M---------------M---------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 2:
			#Genetic Code [2]
			code = "Vertebrate Mitochondrial Code (transl_table=2)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
			Starts = "--------------------------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 3:
			#Genetic Code [3]
			code = "Yeast Mitochondrial Code (transl_table=3)"
			AAs  = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "----------------------------------MM----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 4:
			#Genetic Code [4]
			code = "Mold, Protozoan, Coelenterate Mitochondrial Code & Mycoplasma/Spiroplasma  Code (transl_table=4)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "--MM---------------M------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 5:
			#Genetic Code [5] 
			code = "Invertebrate Mitochondrial Code (transl_table=5)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
			Starts = "---M----------------------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 6:
			#Genetic Code [6]
			code = "Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)"
			AAs  = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 9:
			#Genetic Code [9]  
			code = "Echinoderm and Flatworm Mitochondrial Code (transl_table=9)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 10:
			#Genetic Code [10]   
			code = "Euplotid Nuclear Code (transl_table=10)"
			AAs  = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 11:
			#Genetic Code [11]
			code = "Bacterial, Archaeal and Plant Plastid Code (transl_table=11)"
			AAs  = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "---M---------------M------------MMMM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 12:
			#Genetic Code [12]
			code = "Alternative Yeast Nuclear Code (transl_table=12)"
			AAs  = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-------------------M---------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
																		   
		elif number == 13:
			#Genetic Code [13]
			code = "Ascidian Mitochondrial Code (transl_table=13)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
			Starts = "---M------------------------------MM---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 14:
			#Genetic Code [14]
			code = "Alternative Flatworm Mitochondrial Code (transl_table=14)"
			AAs  = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 15:
			#Genetic Code [15]
			code = "Blepharisma Nuclear Code (transl_table=15)"
			AAs  = "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 16:
			#Genetic Code [16]
			code = "Chlorophycean Mitochondrial Code (transl_table=16)"
			AAs  = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 21:
			#Genetic Code [21]
			code = "Trematode Mitochondrial Code (transl_table=21)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 22:
			#Genetic Code [22]
			code = "Scenedesmus obliquus mitochondrial (transl_table=22)"
			AAs  = "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "-----------------------------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 23:
			#Genetic Code [23]
			code = "Thraustochytrium Mitochondrial Code (transl_table=23)"
			AAs  = "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
			Starts = "--------------------------------M--M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 24:			
			#Genetic Code [24]
			code = "Pterobranchia mitochondrial code (transl_table=24)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
			Starts = "---M---------------M---------------M---------------M------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

		elif number == 25:			
			#Genetic Code [25]
			code = "Candidate Division SR1 and Gracilibacteria Code (transl_table=25)"
			AAs  = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
			Starts = "---M---------------M---------------M----------------------------"
			Base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
			Base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
			Base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
			
		elif number == 1001:
			#User-Defined Genetic Code [1001] is loaded from settings file.
			code = self.settings['code']
			AAs  =   self.settings['AAs']
			Starts = self.settings['Starts']
			Base1  = self.settings['Base1']
			Base2  = self.settings['Base2']
			Base3  = self.settings['Base3']
		
		else:
			raise ValueError, '%s is not a valid genetic code number' % number
		self.code_num = number
		self.code = code
		self.table = [code, AAs, Starts, Base1, Base2, Base3]
	
	
	def setCodons(self, exclude):
		'''
		Use a predefined codon table to generate a dictionary of amino acids with their codons.
		Method is not intended for direct use.
		If exclude=True then certain codons will be excluded from the lists. 
		The codons targeted for removal is read from the settings file.
		'''

		#check whether certain codons should be excluded
		if exclude is True:
			remove = self.settings['codons_to_exclude']
		else:
			remove = []

		#now get the amino acids
		code, AAs, Starts, Base1, Base2, Base3 = self.getTable()
		codons = {'start':[], 'F':[], 'L':[], 'S':[], 'Y':[], 'C':[], 'W':[], 'P':[], 'H':[], 'E':[], 'R':[], 'I':[], 'M':[], 'T':[], 'N':[], 'K':[], 'V':[], 'A':[], 'D':[], 'Q':[], 'G':[], '*':[], 'U':[]}
		for aa, s, b1, b2, b3 in zip(AAs, Starts, Base1, Base2, Base3):
			codon = b1+b2+b3

			if codon in remove:
				continue
			elif aa in 'FLSYCWPHERIMTNKVADQG*U':
				codons[aa].append(codon)
			else:
				raise ValueError, '"%s" is not a valid amino acid' % aa
				
			if s != '-': #if the codon is start
				codons['start'].append(codon)
			
		for key in codons.keys():
			if key != 'U':
				assert codons[key] != [], 'Error, there is no codon assigned to amino acid %s. Revise the user-edited codon table and the codon "exclusion list" in settings.txt' % key
			elif key == 'U' and self.code_num == 1001:
				assert codons[key] != [], 'Error, there is no codon assigned to amino acid %s. Revise the user-edited codon table and the codon "exclusion list" in settings.txt' % key


		self.codons = codons

	######## API intended for use #########
	def getExcluded(self):
		'''
		Return which codons were excluded from the computation.
		'''
		return self.settings['codons_to_exclude']

	def getCode(self):
		'''
		Return which genetic code is represented.
		The output is a string which specifies the code'
		Method is not intended for direct use.
		'''
		return self.code

		
	def getTable(self):
		'''
		Return the codon table data for the specified genetic code.
		The output is a list of strings.
		'''
		return self.table

		
	def getCodons(self, separate=False):
		'''
		Returns a dictionary of amino acids with their codons for the specified codon table.
		The separate variable determines whether codons for one amino acid with dissimilar first two nucleotides should be separated out.
		For example if separate=False the codons for L are 	['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'].
		If separate=True they are split up as L = [['TTA', 'TTG'], ['CTT', 'CTC', 'CTA', 'CTG']]
		'''
		assert type(separate) is bool, 'Error, "separate" must be True or False'
		if separate is False: #returns a dictionary containing the codon table
			return self.codons 
		elif separate is True:
			newdict = {}
			for aa in 'FLSYCWPHERIMTNKVADQG*U':
				f = lambda x: [codon[0:2] for codon in x] #function to get all first two nucleotides for an aa
				firsttwolist = list(set(f(self.codons[aa]))) #list of all unique first two nucleotides for an aa. For example ['TT', 'CT'] for leucine
#				print('aa', aa)
#				print('ftl', firsttwolist)
				if len(firsttwolist) > 1: #if there is more than one set of the first two nucleotides for this amino acid
					newdict[aa] = []
					for i in range(len(firsttwolist)):
						newdict[aa].append([x for x in self.codons[aa] if x[0:2] in firsttwolist[i]]) #add all the codons that match the first two
				else:
					newdict[aa] = self.codons[aa] #
			return newdict				

					
				
	def printTable(self):
		'''
		Print specified codon table.
		'''
		code, AAs, Starts, Base1, Base2, Base3 = self.table
		print('\nCode   = %s' % code)
		print('AAs    = %s' % AAs)
		print('Starts = %s' % Starts)
		print('Base1  = %s' % Base1)
		print('Base2  = %s' % Base2)
		print('Base3  = %s' % Base3)

	def printCodons(self):	
		'''
		Print codons specified by codon table.
		'''
		codons = self.getCodons()
		print('start = %s' %codons['start'])
		print('stop  = %s' %codons['stop'])
		for aa in 'FLSYCWPHERIMTNKVADQG*U':
			print('%s     = %s' % (aa, codons[aa]))
		
		
		
