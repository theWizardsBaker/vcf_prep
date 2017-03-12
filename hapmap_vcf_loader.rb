#!/usr/bin/env ruby
require 'json'
# require 'thread/pool'
require 'parallel'
require 'concurrent'

class HapmapVcfLoader
	attr_accessor :logging
	def initialize()
		@logging = false
		@table_cols = []
	end
	
	def load_variants(vcf_file, variant_output_file_name = "variants.json")
		begin 
			puts "Loading Variants" if @logging
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
						# split the line on tab 
						@table_cols = line.split(/\t/)
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
								"alternate_bases"=> line[4].gsub(/<|>/, '').split(/,/),
								"alternate_structure"=> alt_structure,
								"info" => {}
							  }
						# get each of the variant's info values and add them to our newly created variant
						line[7].split(/;/).each do |e|
							# cut apart at the ='s
							info_column = e.split(/=/)
							# copy the variant info 
							variant["info"][info_column[0]] = info_structure[info_column[0]].clone
							# find out how to store each value for the infos
							val = info_column[1].to_s.split(/\,/).map do |elm|
								case variant["info"][info_column[0]]["type"].downcase
									when "integer"
										elm.to_i
									when "double", "float"
										elm.to_f 
									else
										elm.to_s
								end
							end
							# take the value, comma seperated, and break it into an integer array (faster to search on and compare)
							variant["info"][info_column[0]]['value'] = val
						end

						variant_output.puts variant.to_json
					end
				end
			end
		ensure
			variant_output.close unless variant_output.nil?
		end
	end

	def load_sample(vcf_file, sample_column, sample_output_file_name = "sample.json", append_to_file)
		begin
			puts "Loading #{sample_column}" if @logging
			if @table_cols.empty?
				load_table_header(vcf_file) 
				puts "IS EMPTY"
			end
			# column info
			column = {}
			if (sample_column.is_a? Integer)
				# add the name if sample_column is a number, add 8 to our column so that we avoid the first 8 
				# variant rows
				column[:number] = sample_column + 8
				column[:name] = @table_cols[sample_column]
			else
				# add the number if the sample_column is a name
				column[:number] = @table_cols.find_index(sample_column)
				column[:name] = sample_column
			end

			# if we couldn't find it, get outta here
			raise "Could not find sample #{column} in #{vcf_file}" if column[:name].nil? or column[:number].nil?
			# file for outputting variants
			sample_output = File.open(sample_output_file_name, (append_to_file ? "a" : "w"))
			# construct the sample object
			sample_output.print "{ \"_key\" : \"#{column[:name]}\", \"calls\" : ["

			# start reading the VCF file
			File.open(vcf_file, "r") do |f|
				first = true
				f.each_line do |line, ind|
					if line =~ /^[^#]/
						line = line.chomp
						# break the table row into an array
						line = line.split(/\t/)
						# get each of the samples's call and add it to the appropriate sample
						call = line[column[:number]].to_s
						# unless it's ./. (meaning not recorded)
						unless call =~ /\.(\||\/)\./
							# add the variant, phase (if it's | then phased, if / unphased), and genotype as an integer array
							# this should save us a lot of space 
							genotype = call.split(/\||\//).map(&:to_i)
							if genotype.any? { |x| x > 0 }
								sample = { 
									"variant" => line[2], 
									"phased" => "#{!!(call =~ /\|/)}", 
									"genotype" => genotype
								}
							else 
								sample = { "variant" => line[2] }
							end
							# add the 'ole comma to seperate our json objects
							sample_output.print "," unless first
							sample_output.print "#{sample.to_json}"
							first = false
						end
					end
				end
			end
			sample_output.puts "]}"
			puts "#{sample_column} Done!" if @logging
		ensure
			sample_output.close unless sample_output.nil?
		end
	end

	def load_all_samples(vcf_file, sample_output_file_name = "samples.json", thread_count: Concurrent.processor_count)
		puts "Loading All Samples" if @logging
		# populate the table_col if need be
		load_table_header(vcf_file) if @table_cols.empty?
		# from 9 - oblivion, load each sample
		Parallel.each(@table_cols[9..@table_cols.size], in_processes: thread_count) do | sample_column |
			# t_pool.process do 
				load_sample(vcf_file, sample_column, "#{sample_column}_#{sample_output_file_name}", true)
			# end
		end
		# combine all the smaller files
		puts "Combinging samples" if @logging
		combine = %x{echo *_#{sample_output_file_name} | xargs cat > #{sample_output_file_name}}
		puts "Removing temp_samples" if @logging
		remove = %x{ rm *_#{sample_output_file_name}} if combine
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
				line = line.chomp
				# remove the '#'
				line[0] = ''
				# split the line on tab 
				@table_cols = line.split(/\t/)
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
vcf_file = ARGV.shift
loader.load_variants(vcf_file)
loader.load_all_samples(vcf_file)
# loader.load_sample(vcf_file, "ben7884")



