#!/usr/bin/env python
import fileinput
import re
import sys
import copy
import pymongo
import multiprocessing
from subprocess import call
import glob
from functools import partial

class hapmap_load:

	def __init__(self, logging = False):
		self.logging = logging

	def __load_header(self, header_line, structure):
		# remove the < >
		header_line = re.sub(r'<|>', '', header_line)
		# format the header line so that we have an array of list elements
		header_line = re.split(r',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', header_line)
		# reference element in the info structure
		elm_ref = ''
		if len(header_line) > 1:
			# each header line
			for line in header_line:
				val = line.split(r'=')
				# add the id as the key
				if val[0] == 'ID':
					# set the 
					elm_ref = structure[val[1]] = {}
				else:
					# for this id, add the values
					elm_ref[val[0].lower()] = re.sub('^\"|\"?$', '', val[1])

	def load_vcf(self, vcf_file, database):
		try:
			# Initialize the client for OrientDB
			client = pyorient.OrientDB(database.ip, database.port) 
			# open our database for writing
			client.db_open(database.name, database.user, database.password )
		except Exception as e :
			print "Database %s does not exist" % database.name
			sys.exit()
		# contains the information in the VCF's info field 
		info_structure = {}
		# contains the information in the VCF's alt field
		alt_structure = {}
		# list of the table columns (samples)
		tablecolumns = []
		# open and walk through our VCF file
		for num, line in enumerate(fileinput.input(vcf_header), 1):
			# remove them end lines!
			line = line.rstrip('\r\n') 
			# file header
			if line[:2] == '##':
				if self.logging:
					print "Loading header" 
				# remove the starting ## and split on the first =
				header_line = line[2:].split('=', 1)
				# write to the appropriate hash
				if header_line[0] == "INFO":
					self.__load_header(header_line[1], info_structure)
				elif header_line[0] == "ALT":
					self.__load_header(header_line[1], alt_structure)
			# sample header
			elif line[:1] == '#':
				if self.logging:
					print "Creating Samples"
				# split the line on the tab and drop the first 9 elements (they are not samples)
				tablecolumns = (re.split(r'\t+', line))[9:]
				# iterate through all the table columns 
				for sample in tablecolumns:
					try:
						client.record_create({ @Sample: { 'id' : sample } })
					except Exception as e:
						print "%s could not be added, skipping" % sample
			# calls
			else
				# split the line on the tabs
				variant_calls = re.split(r'\t+', line)
				# create the variant object
				variant = { 
							'id': variant_calls[2],
							'names': variant_calls[2].split(','),
							'chromosome': int(variant_calls[0]),
							'position': variant_calls[1],
							'filter': variant_calls[6],
							'reference_base': variant_calls[3],
							'alternate_bases': variant_calls[4].strip(r'<|>').split(','),
							'alternate_structure': alt_structure,
							'info': {}
						  }
				# walk over each info value
				for info_line in variant_calls[7].split(';'):
					# cut apart at the ='s
					info_line = info_line.split('=')
					# copy the variant info 
					variant['info'][info_line[0]] = copy.deepcopy(info_structure[info_line[0]])
					# find out how to store each value for the infos
					values = []
					if len(info_line) > 1:
						# get each value
						for elm in info_line[1].split(','):
							raw_type = variant['info'][info_line[0]]['type'].lower()
							# determin the type of the value
							if raw_type == 'integer':
								values.append( int(elm) )
							elif raw_type == 'double' or raw_type == 'float':
								values.append( float(elm) )
							else:
								values.append( elm )
					# take the value, comma seperated, and break it into an integer array (faster to search on and compare)
					variant['info'][info_line[0]]['value'] = values
				# try to insert the variant
				try:
					client.record_create({ @Variant: variant })
				except Exception as e:
					print "%s already exists in variants, skipping" % variant_calls[2]
				# add each call to the appropriate sample
				for ind, call in enumerate(variant_calls[9:]):
					# unless the call is ./. or .|.
					if re.search(r'(\d(\/))+\d', call) != None:
						# add the variant, phase (if it's | then phased, if / unphased), and genotype as an integer array
						# add a reference to the variant by getting it's object ID
						variant_call = { 
							'sample': tablecolumns[ind]['id'],
							'variant': variant_calls[2], 
							'phased': re.search(r'\|', call) != None, 
							'genotype': map(lambda x: int(x), re.split(r'\||\/', call))
						}
						# add the call
						try:
							client.record_create({ @Call: variant_call })
						except Exception as e:
							print "%s already exists in variants, skipping" % variant_calls[2]



if len(sys.argv) > 2:
	# create our obj
	loader = hapmap_load()
	# start loading variants
	loader.load_vcf(sys.argv[1], sys.argv[2])
else:
	print "please specify a VCF and OrientDB database"