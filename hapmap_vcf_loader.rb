#!/usr/bin/env ruby
require 'json'

class HapmapVcfLoader
	attr_reader :samples, :variants
	def initialize()
		@table_cols = []
	end
	
	def load_variants(vcf_file, variant_output_file_name = "variants.json")
		# contains the information in the VCF's info field 
		info_structure = {}
		# contains the information in the VCF's alt field
		alt_structure = {}
		# file for outputting variants
		variant_output = File.open(variant_output_file_name, "w")
		# open and walk through our VCF file
		File.open(vcf_file, "r") do |f|
			f.each_line do |line|
				line = line.chomp
				# file header
				if line[0..1] == '##'
					# remove the starting ## and split on the first =
					line = line[2..line.size].split(/=/, 2)
					# write to the appropriate hash
					case line[0]
						when "INFO"
							load_header(line[1], info_structure)
						when "ALT"
							load_header(line[1], alt_structure)
					end
				elsif line[0] == '#'
					# remove the '#'
					line[0] = ''
					# split the line on tab and remove the first 9 elements (the variant info)
					@table_cols = line.split(/\t/).drop(9)
				else
					# break the table row into an array
					line = line.split(/\t/)
					# add the row's variant 
					variant = { 
							"_key" => line[2],
							"names" => line[2].split(/,/),
							"chromosome"=> line[0].to_i,
							"position"=> line[1],
							"filter"=> line[6],
							"reference_base"=> line[3],
							"alternate_bases"=> line[4].split(/,/),
							"alternate_structure"=> alt_structure,
							"info" => Marshal.load(Marshal.dump(info_structure))
						  }
					# get each of the variant's info values and add them to our newly created variant
					line[7].split(/;/).each do |e|
						# cut apart at the ='s
						info = e.split(/=/)
						# take the value, comma seperated, and break it into an integer array (faster to search on and compare)
						variant.last["info"][info[0]]['value'] = info[1].to_s.split(/\,/).map(&:to_i)
					end

					variant_output.puts variant.to_json
				end
			end
		end

		ensure
			variant_output.close unless sample_output.nil?
		end
	end

	def load_sample(vcf_file, sample_column, sample_output_file_name = "sample.json")
		load_table_header(vcf_file) if @table_cols.empty?
		# column info
		column = {}
		column[:name] = String(sample_column) rescue nill
		column[:number] = Integer(sample_column) rescue nill
		# add the name if sample_column is a number
		column[:name] = @table_cols[sample_column] if column[:name].nil?
		# add the number if the sample_column is a name
		column[:number] = @table_cols.find_index(sample_column) if column[:number].nil?

		# if we couldn't find it, get outta here
		raise "Could not find sample #{sample_column} in #{vcf_file}" if column[:name].nil? or column[:number].nil?

		# file for outputting variants
		sample_output = File.open(sample_output_file_name, "w")
		# construct the sample object
		sample_output.print "{ \"_key\" : #{column[:name]}, \"calls\" : ["

		# start reading the VCF file
		File.open(vcf_file, "r") do |f|
			f.each_line do |line|
				line = line.chomp
				# break the table row into an array
				line = line.split(/\t/)
				# get each of the samples's call and add it to the appropriate sample
				line[column[:number]].each do |call|
					# unless it's ./. (meaning not recorded)
					unless call =~ /\.(\||\/)\./
						# add the variant, phase (if it's | then phased, if / unphased), and genotype as an integer array
						sample = { 
							"variant" => line[2], 
							"phased" => "#{!!(call =~ /\|/)}", 
							"genotype" => call.split(/\||\//).map(&:to_i)
						}
						sample_output.print sample.to_json
					end
				end
			end
		end

		sample_output.print "]}"

		ensure
			sample_output.close unless sample_output.nil?
		end
	end

	def load_all_samples(vcf_file, sample_output_file_name = "sample.json")
		# populate the table_col if need be
		load_table_header(vcf_file) if @table_cols.empty?
		# from 9 - oblivion, load each sample
		@table_cols[9..@table_cols.size].each do | sample_column |
			load_sample(vcf_file, sample_column, "#{sample_column}_#{sample_output_file_name}")
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

	def load_table_header(vcf_file)
		# find the headding line in our VCF file
		File.open(vcf_file, "r").each_line do | line |
			# if we're on the #CROM line
			if line =~ /^#[^#]/
				# remove the '#'
				line[0] = ''
				# split the line on tab and remove the first 9 elements (the variant info)
				@table_cols = line.split(/\t/).drop(9)
				# get the heck out of this loop! it'd take forever otherwise
				next
			end
		end
	end
	
end



if ARGV.empty?
	puts "Please supply a VCF file" 
	exit
end

loader = HapmapVcfLoader.new
loader.load_variants(ARGV.shift)
loader.load_sample(ARGV.shift)




