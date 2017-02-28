#!/usr/bin/env ruby
require 'pp'
require 'json'

class HapmapVcfLoader
	attr_reader :samples, :variants
	def initialize()
		# contains the information in the VCF's info field 
		@info_structure = {}
		# contains the information in the VCF's alt field
		@alt_structure = {}
		# contain the list of samples columns ( > 9 ) in the VCF
		@samples = []
		# contain the list of variants in the VCF (col < 9)
		@variants = []
	end
	
	def load_vcf(vcf_file)
		# start reading the VCF file
		File.open(vcf_file, "r") do |f|
			f.each_line do |line|
				line = line.chomp
				# file header
				if line[0..1] == '##'
					# remove the starting ##
					line = line.slice(2, line.size).split(/=/, 2)
					# write to the appropriate hash
					case line[0]
						when "INFO"
							load_header(line[1], @info_structure)
						when "ALT"
							load_header(line[1], @alt_structure)
					end
				# table header
				elsif line[0] == '#'
					# remove the starting #
					line = line.slice(1, line.size)
					# split the line into an array and remove the first 9 (the variant info)
					table_cols = line.split(/\t/).drop(9)
					# for each column (sample), we'll create an entry
					table_cols.each do |col|
						@samples << { "_key" => col, "calls" => [] }
					end
				# table row
				else
					# break the table row into an array
					line = line.split(/\t/)
					# add the row's variant 
					@variants << { 
							"_key" => line[2],
							"names" => line[2].split(/,/),
							"chromosome"=> line[0].to_i,
							"position"=> line[1],
							"filter"=> line[6],
							"reference_base"=> line[3],
							"alternate_bases"=> line[4].split(/,/),
							"info" => @info_structure.clone
						  }
					# get each of the variant's info values and add them to our newly created variant
					line[7].split(/;/).each do |e|
						# cut apart at the ='s
						info = e.split(/=/)
						# take the value, comma seperated, and break it into an integer array (faster to search on and compare)
						@variants.last["info"][info[0]]['value'] = info[1].to_s.split(/\,/).map(&:to_i)
					end
					# get each of the samples's call and add it to the appropriate sample
					line[9..line.size].each_with_index do |call, ind|
						# unless it's ./. (meaning not recorded)
						unless call =~ /\.(\||\/)\./
							# add the variant, phase (if it's | then phased, if / unphased), and genotype as an integer array
							@samples[ind]["calls"] << { 
								"variant" => line[2], 
								"phased" => "#{!!(call =~ /\|/)}", 
								"genotype" => call.gsub(/(\||\/)/, '').scan(/\d/).map(&:to_i)
							}
						end
					end

				end
			end
		end
	end

	private
	
	def load_header(header_line, structure)
		# format the header line so that we have an array of list elements
		header_line = header_line.gsub(/<|>/, '').split(/,(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/)
		# reference element in the info structure
		elm_ref = ''
		header_line.each do |e|
			val = e.split(/=/, 2)
			# add the id as the key
			if val[0] == "ID"
				elm_ref = val[1]
				structure[elm_ref] = {}
			else
				# for this id, add the values
				structure[elm_ref][val[0].downcase] = val[1].to_s.gsub(/^\"|\"?$/, '')
			end
		end
		
	end
	
end



if ARGV.empty?
	puts "Please supply a VCF file" 
	exit
end

loader = HapmapVcfLoader.new
loader.load_vcf(ARGV.shift)
loader.variants.each { |var| puts "db._query('INSERT @document INTO variants'," + var.to_json + "})" }
loader.samples.each { |samp| puts "db._query('INSERT @document INTO variants'," + samp.to_json + "})" }
# loader.samples.each { |samp| pp "db._query('INSERT @document INTO samples', #{samp.to_json})"  }

