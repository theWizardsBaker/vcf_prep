#!/usr/bin/env ruby
require 'mongo'
require 'parallel'
require 'concurrent'
require 'pp'

class Hapmap_Load

	def initalize(logging = false)
		@logging = logging
	end

	def load_vcf(vcf_file, database_name, pool_count = Concurrent.processor_count)
		# select the database
		@client = Mongo::Client.new([ 'localhost:27017' ], :database => database_name)
		# split the inital file
		# system "csplit -f vcf_loader_head_tmp_ #{vcf_file} '/^[^#].*$/'"
		# system "csplit -f vcf_loader_head_tmp_ #{vcf_file} '/##HapMapVersion/'1 "
		# system "split -f vcf_loader_head_tmp_ #{vcf_file} -n $(`grep -n '#CHROM' example_csplit.txt | cut  -d : -f 1`)"
		header_end = %x(grep -n '^#CHROM' #{vcf_file} | cut  -d : -f 1)
		# load the header elements
		parse_header(vcf_file, header_end)
		# split the second file based on how many cores we have
		system "split -l$((`wc -l < vcf_loader_head_tmp_01`/#{pool_count})) #{vcf_file} vcf_loader_tmp_ "
		# remove the first tmps
		system "rm vcf_loader_head_tmp_*"
		# get all tmp files
		vcf_shards = %x(ls vcf_loader_tmp_*).split(/\n/)
		# run parallel processes 
		Parallel.each(vcf_shards, in_processes: pool_count) { |vcf_shard_file| load_rows(vcf_shard_file) }
		# add the indexs as a final step
		loader.add_index()
	end

	def add_index()
		# lets add the index if one doesn't exist
		@client[:variants].indexes.create_one({ '_key' => Mongo::TEXT }, :name => 'variant_search_index', :unique => true)
		# lets add the index if one doesn't exist
		@client[:samples].indexes.create_one({ '_key' => Mongo::TEXT }, :name => 'sample_search_index', :unique => true)
	end

	private

	def parse_header(vcf_header, header_end)
		puts "Begin load" if @logging
		# contains the information in the VCF's info field 
		@info_structure = {}
		# contains the information in the VCF's alt field
		@alt_structure = {}
		# open and walk through our VCF file
		File.open(vcf_header, "r").limit(header_end) do |f|
			f.each_line do |line|
				# remove them end lines!
				line = line.chomp
				# file header
				if line[0..1] == '##'
					puts "Loading header" if @logging
					# remove the starting ## and split on the first =
					header_line = line[2..line.size].split(/=/, 2)
					# write to the appropriate hash
					case header_line[0]
						when "INFO"
							load_header(header_line[1], @info_structure)
						when "ALT"
							load_header(header_line[1], @alt_structure)
					end
				# sample header
				elsif line[0] == '#'
					puts "Creating Samples" if @logging
					# split the line on the tab and drop the first 9 elements (they are not samples)
					@table_columns = line.split(/\t/)
					# walk over each element and add it to the mongo client
					@table_columns.each { |sample| @client[:samples].insert_one({ _key: sample, calls: [] }) }
				end
			end
		end
	end

	def load_header(header_line, structure)
		# format the header line so that we have an array of list elements
		header_line = header_line.gsub(/<|>/, '').split(/,(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/)
		# reference element in the info structure
		elm_ref = ''
		# each header line
		header_line.each do |line|
			val = line.split('=', 2)
			# add the id as the key
			if val[0] == 'ID'
				elm_ref = val[1]
				structure[elm_ref] = {}	
			else
				# for this id, add the values
				structure[elm_ref][val[0].downcase] = val[1].to_s.gsub(/^\"|\"?$/, '')
			end
		end
	end

	def load_rows(vcf_file)
		puts "Loading Variant" if @logging
		# walk over each line
		File.open(vcf_file, "r") do |f|
			f.each_line do |line|
				# don't grab the headers
				next if line =~ /^#/
				# remove the newline
				line = line.chomp
				# split the line on the tabs
				variant_calls = line.split(/\t/)
				# create the variant object
				variant = { 
						'_key' =>  variant_calls[2],
						'names' =>  variant_calls[2].split(/,/),
						'chromosome' => variant_calls[0].to_i,
						'position' => variant_calls[1],
						'filter' => variant_calls[6],
						'reference_base' => variant_calls[3],
						'alternate_bases' => variant_calls[4].gsub(/<|>/, '').split(/,/),
						'alternate_structure' => @alt_structure,
						'info' =>  {}
					  }
				# walk over each info value
				variant_calls[7].split(';').each do |info_line|
					# cut apart at the ='s
					info_column = info_line.split('=')
					# copy the variant info 
					variant['info'][info_column[0]] = @info_structure[info_column[0]].clone
					# find out how to store each value for the infos
					# break apart each info line by the ,  and map all values by type
					variant['info'][info_column[0]]['value'] = info_column[1].to_s.split(/\,/).map do |elm|
						case variant['info'][info_column[0]]['type'].downcase
							when 'integer'
								elm.to_i
							when 'double', 'float'
								elm.to_f 
							else
								elm.to_s
						end
					end
				end
				# try to insert the variant
				@client[:variants].insert_one(variant)
				# sample time
				variant_calls[9..-1].each_with_index do |call, ind|
					# unless the call is ./. or .|.
					unless call =~ /\.(\||\/)\./
						# add the variant, phase (if it's | then phased, if / unphased), and genotype as an integer array
						# add a reference to the variant by getting it's object ID
						var_ref = @client[:variants].find({ '_key' => variant_calls[2] })
						variant_call = { 
							'variant' => var_ref.first[:_id], 
							'phased' => !!(call =~ /\|/), 
							'genotype' => call.split(/\||\//).map(&:to_i)
						}
						puts "Inserting sample #{@table_columns[ind]}" if @logging
						# add the call to our sample
						@client[:samples].update_one({ '_key' => @table_columns[ind] }, { '$push' => { 'calls' => variant_call } })
					end
				end
			end
		end
	end
end



if ARGV.size >= 2
	# create our obj
	loader = Hapmap_Load.new
	# start loading variants
	loader.load_vcf(ARGV.shift, ARGV.shift)
else
	puts "please specify a VCF and mongo database"
end