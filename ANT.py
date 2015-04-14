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

from copy import deepcopy
import dna
import protein
import re
import math

	



class DegenerateCodon:
	'''
	Class that holds methods and values for computing the degenerate codon for a list of amino acids. 
	Alternatively, the class can evaluate an degenerate codon which is provided by the user. 
	Required input is an integer that determines the genetic code to use and either a list of desired amino acids in single letter code 
	OR a three-letter codon using the IUPAC Nucleotide ambiguity code (G, A, T, C, R, Y, W, S, M, K, H, B, V, D, N).
	
	
	The algorithm works as follows:
	
	* If a list of amino acids (as single letter code) is passed to the algorithm;
	All possible regular codons for those amino acids are looked up and returned as a list of lists.
	All nucleotides for first, second and third position are extracted into separate lists while retaining list structure.
	Separately, for pos 1, 2 and 3, all degenerate nucleotides that matches at least one nucleotide in each list are found.
	Every combination of these is then made to generate a set of possible degenerate codons.
	Each of these are scored by converting back to "real" codons, translating to amino acids and checking against the user-defined amino acid selection.
	The scoring function retains the degenerate codons which have the fewest off-target amino acids (encoded amino acids that were not chosen by the user).
	This can be one or several, but they all represent the minimum number of off-target amino acids.
	The scoring function then goes through the degenerate codons and selects that which encodes the fewest number of "real" codons, to decrease redundancy. 
		
	* If an degenerate codon is passed to the algorithm;
	The first, second and third degenerate nucleotide is converted to their actual nucleotide counterparts.
	Every combination of these is then made while maintaining the position of each base to generate the "real" codons encoded by the degenerate codons. 
	The resulting codons are then translated to amino acids using the specified genetic code.
	
	
	
	Regardless of which of the two ways is used to generate the codon object the methods to retrieve information are the same:

	To get the degenerate codon (as a three-character string):
	codon_object.getTriplet()
	
	To get all amino acids encoded by the current triplet (as a list of upper case amino acids in single letter code):
	codon_object.getEncoded()
	
	To get the target amino acids (as a list of upper case amino acids in single letter code):
	codon_object.getTarget()
	
	To get the off-target amino acids (as a list of upper case amino acids in single letter code):
	codon_object.getOffTarget()

	To get which amino acids can still be chosen without further off-target amino acids (as a list of upper case amino acids in single letter code):	
	codon_object.getPossible()
	
	To get which genetic code was used for the computation (as an integer):
	codon_object.getTable()
	
	To get which "real" codons are encoded by the degenerate triplet (as a list of upper case triplets using GATC code):	
	codon_object.getCodons()
	
	To get how often each amino acid is encoded by the degenerate triplet (as a dictionary with amino acid upper case single letter keys and integer values): 
	codon_object.getCodonsPerAA()
	
	To get which alternative degenerate triples (with the same number of off-target amino acids) might be used (as a list of three-letter strings of upper-case characters):
	codon_object.getAlternatives()
	
	To get an extended list of alternative triplets, and the amino acids they encode (this list is no longer limited to codons with the same number of off-target amino acids):
	codon_object.getExtendedAlternatives()
	
	To get a report containing all the above information plus the library size and screening burden (as a string):
	codon_object.getReport()
	'''
	
	def __init__(self, input, table=1):
		self.setTable(table)
		
		#input can be either a three-nucleotide string or a list of amino acids
		if len(input) == 3 and type(input) == str: #if string i.e. an degenerate codon
			input = input.upper()
			self.evaluateTriplet(input)
		elif type(input) == list: #if list, i.e. a list of amino acids to evaluate
			input = [s.upper() for s in input]
			self.setTarget(input)
		else:
			raise ValueError, 'The input is not valid'
		return

		
	###############################################################
	#######  Public methods intended for user interaction  ########
	###############################################################
	def getEncoded(self):
		'''
		Retrieves a list of the target+offtarget amino acids.
		I.e. all the amino acids encoded by the current triplet.
		Output is a list of amino acids in upper case single letter code. 
		'''
		return self.getTarget() + self.getOffTarget()
		
	def getTarget(self):
		'''
		Retrieves a list of the target amino acids.
		Output is a list of upper case amino acids in single letter code.
		'''
		return self.target
		
	def getOffTarget(self):
		'''
		Retrieves a list of the off-target amino acids.
		Output is a list of upper case amino acids in single letter code.
		'''
		return self.offtarget
		
	def getPossible(self):
		'''
		Retrieves a list of amino acids still possible without further off-targets.
		Output is a list of upper case amino acids in single letter code.
		'''
		return self.possible
		
	def getTriplet(self):
		'''
		Retrieves the degenerate codon.
		Output is a three-letter string of upper case characters.
		'''
		return self.triplet	
		
	def getTable(self):
		'''
		Retrieves which genetic code was used.
		Output is an integer.
		'''
		return self.table
		
	def getCodons(self):
		'''
		Retrieves a list of all the "real" codons encoded by the degenerate codon.
		Output is a list of three-letter strings of ACTG upper case characters.
		'''
		return dna.UnAmb(self.getTriplet())
	
	def getCodonsPerAA(self):
		'''
		Retrieves a dictionary specifying how many times each amino acid is coded for by the ambiguous codon.
		Output is a dictionary with amino acid upper case single letter keys and integer values.
		'''
		return protein.count_aa(dna.Translate(''.join(dna.UnAmb(self.getTriplet())), self.getTable()))
		
	def getExtendedAlternatives(self):
		'''
		To get an extended list of alternative triplets, and the amino acids they encode.
		This list is no longer limited to codons with the same number of off-target amino acids.	
		'''
		return self.extendedalternatives
		
	def getAlternatives(self):
		'''
		Retrieves alternative degenerate codons which have same number of off-target amino acids as best degenerate codon,
		but may differ in the number of codons per amino acid.
		Output is a list of degenerate codons as three-letter strings of upper-case characters.
		'''
		return self.alternatives
	
	def getReport(self):
		'''
		Retrieve a report containing all available data.
		Output is a string.
		'''
		#load data from settings file		
		temp = dict()
		execfile('./settings.txt', temp)
		assert type(temp['library_coverage']) is int, 'Error, the library coverage must be an integer between 1 and 99. Please review the settings.txt file.'
		assert 1 <= temp['library_coverage'] <= 99, 'Error, the library coverage must be an integer between 1 and 99. Please review the settings.txt file.'


		triplet = self.getTriplet()
		codons = self.getCodons()
		num_codons = len(codons)
		output = 'Degenerate codon: %s\n' % triplet
		output += 'genetic code: %s\n' % self.getTable()
		output += 'Codons which were excluded from the computation: %s\n' % temp['codons_to_exclude']
		output += 'Real codons encoded by the degenerate codon: %s\n' % codons
		output += 'Target amino acids: %s\n' % self.getTarget()
		output += 'Encoded amino acids: %s\n' % self.getEncoded()
		output += 'Off-target amino acids: %s\n' % self.getOffTarget()
		output += 'Amino acids that can be added w/o further off-targets: %s\n' % self.getPossible()
		output += 'Codons for each amino acid: %s\n' % self.getCodonsPerAA()
		output += 'Library size (number of codons): %s\n' % num_codons
		output += 'Clones to screen for %s%% library coverage: %s\n' % (temp['library_coverage'], int(-math.log(1-temp['library_coverage']/100.0)/(1/float(num_codons))))    #T=-ln(1-Pi)/Fi
		output += 'Alternate codons with same number of off-target amino acids: %s\n' % self.getAlternatives()
		output += 'Alternate codons with same number or more off-target amino acids: %s\n' % self.getExtendedAlternatives()
		return output
		
	################################################################

	
		
		
		
	################################################################
	######## Methods NOT intended for direct user interaction ######
	################################################################
	
	def sumupcodons(self, DesiredCodons):
		'''
		Takes a list of regular codon lists and does two things; first it splits them up based on the first, second and third position of the triplet,
		then it checks which nucleotide (degenerate allowed) that matches at least one of the nucleotides from each amino acid.	
		
		Example input, where alanine, cysteine and tyrosine are desired is: [['GCT', 'GCC', 'GCA', 'GCG'], ['TGT', 'TGC'], ['TAT', 'TAC']]
		The objective is to keep the list structure but to make three separate list where each 
		holds all the unique nucleotides for either the first, second or third base of the codons.
		The correct output for the first position would be: [['G'], ['T'], ['T']],
		for the second [['C'], ['G'], ['A']],
		and for the third [['T', 'C', 'A', 'G'], ['T', 'C'], ['T', 'C']].
		These lists can then be passed on to another function for finding a nucleotide that matches at least one of the nucleotides from each amino acid.	
		'''
		allcodon1 = []
		allcodon2 = []
		allcodon3 = []
		for entry in DesiredCodons:
			codon1 = []
			codon2 = []
			codon3 = []
			for codon in entry: #splits up the codons of each AA into triplets
				if codon[0] not in codon1:
					codon1.extend(codon[0]) #Takes first base in the triplet
				if codon[1] not in codon2:
					codon2.extend(codon[1]) #Takes the second base in triplet
				if codon[2] not in codon3:
					codon3.extend(codon[2]) #Takes third base in triplet
			allcodon1.append(codon1) #Adds up all the first bases
			allcodon2.append(codon2) #Adds up all the second bases
			allcodon3.append(codon3) #Adds up all the third bases
		return (allcodon1, allcodon2, allcodon3)


	def extra_list_elements(self, list_A, list_B): 
		'''
		Method for comparing two lists to find which elements are not present in both.
		The elements which are not common to both lists are returned as a list.
		Duplicates are ignored such that if list A has the integer 15 twice and list B only once, this will NOT result in 15 being returned.
		
		This method is used to check which off-target amino acids there are.
		'''
		not_in_B = [s for s in list_A if s not in list_B]
		not_in_A = [s for s in list_B if s not in list_A]
		not_in_both = list(set(not_in_A + not_in_B))
		return not_in_both
	
	
	def flatten_codon_list(self, codon_list):
		'''
		Some of the codons are list of lists (happens when the amino acid has codons at different parts of the codon circle).
		This method flattens the structure one level such that a list of codons remains.
		'''
		output = []
		#print(codon_list)
		for pos in codon_list:
			if pos == []:
				pass
			else:
				output_len = len(output)
				if isinstance(pos[0], str):
					if output_len == 0:
						output.append([pos])
					else:
						for o in range(output_len):
							output[o].append(pos)
						
				elif isinstance(pos[0], list):
					if output_len == 0:
						for p in pos:
							output.append([p])
					else:
						output.extend(deepcopy(output * (len(pos)-1)))
						for i in range(len(pos)):
							for j in range(output_len):
								output[j+i*output_len].append(pos[i])
		return output

	
	def find_degenerate(self, AA_list):
		'''
		Method for finding an degenerate codon encoding a list of desired amino acids.
		The method finds the codon(s) with fewest off-target amino acids.
		To reduce redundancy, the method then goes through all the best codons 
		(they all have the same number of off-target amino acids) and finds the one with the lowest number of codons. 
		If there are still more than one which are equivalent, the method then picks one WITHOUT a stop codon.
		
		The input is a list of upper case amino acids in single-letter code.
		The valid values are: FLSYCWPHERIMTNKVADQG*		
		
		The output is a tuple of the best degenerate codon and the off-target amino acids.
		The degenerate codon is a string of three of the following characters: GATCRYWSMKHBVDN
		The off-target amino acids is a list of upper case amino acids in single letter code.
		'''
		#make sure input is OK
		assert all([s in 'FLSYCWPHERIMTNKVADQG*U' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list
		
		#get all codons for chosen amino acids
		regular_triplets = [dna.GetCodons(aa, table=self.getTable(), separate=True, exclude=True) for aa in AA_list]

		#some of the codons are list of lists (happens when the amino acid has codons at different parts of the codon circle)
		#I need to flatten this into separate lists with which go on further
		regular_triplets = self.flatten_codon_list(regular_triplets)

		best_score = None
		all_alternatives = [] #to save the result of all possible triplets
		for codon_list in regular_triplets:
			#get all nucleotides for first, second and third position while retaining list structure		
			first, second, third = self.sumupcodons(codon_list) 
			
			#check which degenerate nucleotide can be used to find at least one match in each of the lists
			possible_triplets = dna.combine([dna.commonNuc(first), dna.commonNuc(second), dna.commonNuc(third)])
						
			#now go through them and see which is best
			for triplet in possible_triplets:
				#convert the triplet back to a list of real codons 
				Realcodons = dna.combine([dna.UnAmb(triplet[0]), dna.UnAmb(triplet[1]), dna.UnAmb(triplet[2])]) #condense the different codons for position 1, 2, and 3 to a list of triplets
			
				#Check which AA these codons code for
				ResultingAA = [dna.Translate(codon, table=self.getTable()) for codon in Realcodons]

				#compare which amino acids were desired with the ones resulting from the degenerate codon
				offtarget = sorted(self.extra_list_elements(AA_list, ResultingAA))

				#add to all options
				if any([True for s in all_alternatives if s[0] == triplet]) is False:
					all_alternatives.append([triplet]+AA_list+offtarget)
				
				#if there are fewer off-target amino acids with the new codon, keep it
				if len(offtarget) < best_score or best_score == None:
					best_score = len(offtarget)
					good_triplets = []
					good_triplets.append(triplet)
				elif len(offtarget) == best_score:
					good_triplets.append(triplet)
		
		#the saved triplets all have the same number of off-target amino acids, now keep the one with the lowest number of codons (to reduce ambiguity)
		best_triplet = None #for storing best degenerate triplet
		best_offtarget = None #for storing the off-target AA of the best triplet
		best_score = None #for storing the length of the off-target list
		alternatives = [] #for saving alternative triplets and their encoded amino acids
		for triplet in good_triplets:
			#convert the triplet back to a list of real codons 
			Realcodons = dna.combine([dna.UnAmb(triplet[0]), dna.UnAmb(triplet[1]), dna.UnAmb(triplet[2])]) #condense the different codons for position 1, 2, and 3 to a list of triplets

			#Check which AA these codons code for
			ResultingAA = [dna.Translate(codon, table=self.getTable()) for codon in Realcodons]

			#compare which amino acids were desired with the ones resulting from the degenerate codon
			offtarget = sorted(self.extra_list_elements(AA_list, ResultingAA))
				
			#save alternatives stats
			if any([True for s in alternatives if s[0] == triplet]) is False:
					alternatives.append([triplet]+AA_list+offtarget)
			
			#save the stats in case there are fewer codons
			if len(Realcodons) < best_score or best_score == None:
				
				#save stats
				best_score = len(Realcodons)
				best_triplet = triplet
				best_offtarget = offtarget
		
			#if another codon has same stats as the previous best one, replace the previous codon if it has an off-target stop
			elif len(Realcodons) == best_score and '*' in best_offtarget: 

				#save stats
				best_score = len(Realcodons)
				best_triplet = triplet
				best_offtarget = offtarget
			
			
		return best_triplet, best_offtarget, alternatives, all_alternatives

	
		
	
	def next_steps(self):
		'''
		Method for finding which other amino acids can be selected without introducing 
		further (more in number) off-target ones.
		The output is a list of upper case amino acids in single letter code.
		'''

		possibleAA = []
		targetAA = self.getTarget()
		
		unusedAA = self.extra_list_elements(targetAA, list('FLSYCWPHERIMTNKVADQG*U'))
		
		if len(unusedAA)>0:
			for AA in unusedAA:
				temptargetAA = targetAA[:]
				temptargetAA.append(AA)
				triplet, offtarget, alternatives, all_options = self.find_degenerate(temptargetAA)
				if len(offtarget) <= len(self.getOffTarget()):
					possibleAA.append(AA)

		return sorted(possibleAA)


	def evaluateTriplet(self, amb_codon):
		'''
		Evaluate the degenerate codon by computing which amino acids it codes for.
		The input is a string, three letters long and comprising only IUPAC Nucleotide ambiguity code.
		The valid values is any combination of three of the following: GATCRYWSMKHBVDN		
		'''
		#make sure input is OK
		assert type(amb_codon) is str and len(amb_codon) == 3, 'Error, the degenerate codon must be a string three characters long.'
		m = re.match('^[GATCRYWSMKHBVDN]{3}$', amb_codon)
		assert m != None, 'Error, the codon %s is not valid. It may only use the chracters GATCRYWSMKHBVDN.' % amb_codon
		
		#compute target amino acids and set variables
		self.target = list(set([dna.Translate(s, self.getTable()) for s in dna.UnAmb(amb_codon)]))
		self.setTriplet(amb_codon)
		self.setOffTarget([])
		
		#Now let's get the alternative codons
		triplet, offtarget, alternatives, all_options = self.find_degenerate(self.getTarget())
		
		self.setAlternatives(alternatives)
		self.setExtendedAlternatives(sorted(all_options, key=len))
		
		#see which other amino acids are possible without further off-target
		self.setPossible(self.next_steps())


	
	def setTarget(self, AA_list):
		'''
		Set which target amino acids are desired.
		The input is a list of upper case amino acids in the single-letter code.	
		The valid values are: FLSYCWPHERIMTNKVADQG*	
		'''
		#make sure input is OK
		assert all([s in 'FLSYCWPHERIMTNKVADQG*U' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list
		self.target = AA_list

		#compute the single triplet and the off-target AA
		triplet, offtarget, alternatives, all_options = self.find_degenerate(self.getTarget())
		
		self.setTriplet(triplet)
		self.setOffTarget(offtarget)
	
		self.setAlternatives(alternatives)
		self.setExtendedAlternatives(sorted(all_options, key=len))
		
		#see which other amino acids are possible without further off-target
		self.setPossible(self.next_steps())
		

	
	def setOffTarget(self, AA_list):
		'''
		Set which amino acids can still be chosen in the next step without further off-target amino acids.
		The input is a list of upper case amino acids in the single-letter code.	
		The valid values are: FLSYCWPHERIMTNKVADQG*
		'''
		assert all([s in 'FLSYCWPHERIMTNKVADQG*U' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list
		self.offtarget = AA_list

		
	def setPossible(self, AA_list):
		'''
		Set which amino acids can still be chosen in the next step without further off-target amino acids.
		The input is a list of upper case amino acids in the single-letter code.
		The valid values are: FLSYCWPHERIMTNKVADQG*
		'''
		assert all([s in 'FLSYCWPHERIMTNKVADQG*U' for s in AA_list]), 'Error, one or more of the amino acids %s are not valid.' % AA_list		
		self.possible = AA_list
	
	
	def setTriplet(self, amb_codon):
		'''
		Set the degenerate codon.
		The input is a string, three letters long and comprising only IUPAC Nucleotide ambiguity code.
		The valid values is any combination of three of the following: GATCRYWSMKHBVDN
		'''
		assert type(amb_codon) is str and len(amb_codon) == 3, 'Error, the degenerate codon must be a string three characters long.'
		m = re.match('^[GATCRYWSMKHBVDN]{3}$', amb_codon)
		assert m != None, 'Error, the codon %s is not valid. It may only use the chracters GATCRYWSMKHBVDN.' % amb_codon
		self.triplet = amb_codon
	
	def setExtendedAlternatives(self, options_list):
		'''
		Sets a list of lists containing alternative codons and all the amino acids encoded by them.
		'''
		self.extendedalternatives = options_list
		
	def setAlternatives(self, alternatives_list):
		'''
		Set alternative degenerate codons which have same number of off-target amino acids as best degenerate codon,
		but may differ in the number of codons per amino acid.
		Valid input is a list of upper-case three-letter degenerate codons composed of three of the following: GATCRYWSMKHBVDN
		'''
#		assert all(re.match('^[GATCRYWSMKHBVDN]{3}$', s) != None for s in alternatives_list), 'Error, the degenerate codon must be a string three characters long.'
		self.alternatives = alternatives_list
		
		
	def setTable(self, table):
		'''
		Set which genetic code to use.
		The input is an integer.
		The valid values are: 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25 and 1001
		'''
		table = int(table)
		assert table in [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 1001], 'Error, %s is an invalid genetic code.' % table
		self.table = table
	
	################################################################		

	
		
if __name__ == '__main__':

	#specify how to parse the arguments
	import argparse
	parser = argparse.ArgumentParser() 
	parser.add_argument('--codon')
	parser.add_argument('--aa', nargs='*')
	parser.add_argument('--table')
	args = parser.parse_args()

	assert (args.codon != None and args.aa != None) is not True, 'Error, you cannot specify codon and a set of amino acids at the same time.'
	assert (args.codon == None and args.aa == None) is False, 'Error, you must specify a codon or a set of amino acids.'

	#If a genetic code was specified, use it. Otherwise use 1.
	if args.table == None:
		table = 1
	else:
		table = args.table
	
	#Now use the codon, or amino acids depending on what was given.
	if args.codon == None: #If a set of amino acids were specified. 
		AA = args.aa
		codon_object = DegenerateCodon(AA, table)

	elif args.aa == None: #If a codon was specified. 
		codon = args.codon
		codon_object = DegenerateCodon(codon, table)
	else:
		raise ValueError


	print(codon_object.getReport())

