#!/usr/bin/env ruby
require 'mongo'
require 'slop'
# require 'zlib'
# require 'parallel'
# require 'concurrent'
# require 'pp'

class Hapmap_Load

	def initalize(logging = false)
		@logging = logging
	end

	def load_vcf(vcf_file, database_options = {})
		puts "Begin load" if @logging
		begin
			Mongo::Logger.logger.level = ::Logger::FATAL
			# select the database
			client = Mongo::Client.new([ "#{database_options[:ip]}:27017" ], database_options)
			# contains the information in the VCF's info field 
			info_structure = {}
			# contains the information in the VCF's alt field
			alt_structure = {}
			# list of the table columns (samples)
			tablecolumns = []
			File.open(vcf_file, "r") do |f|
				f.each_line do |line|
					line = line.chomp
					# file header
					if line[0..1] == '##'
						puts "Loading header" if @logging
						# remove the starting ## and split on the first =
						header_line = line[2..line.size].split(/=/, 2)
						# write to the appropriate hash
						case header_line[0]
							when "INFO"
								load_header(header_line[1], info_structure)
							when "ALT"
								load_header(header_line[1], alt_structure)
						end
					# sample header
					elsif line[0] == '#'
						puts "Creating Samples" if @logging
						# split the line on the tab and drop the first 9 elements (they are not samples)
						# walk over each element and add it to the mongo client
						tablecolumns = line.split(/\t/)[9..-1].map do |sample| 
							sample_id = BSON::ObjectId.new
							{ _key: sample, _id: sample_id }
						end
						client[:samples].insert_many(tablecolumns)
						# Mongo::BulkWrite.get(client[:samples], tablecolumns, ordered:false, write_concern: Mongo::BulkWrite::UnorderedBulkWrite)
						# blk = Mongo::BulkWrite.new(client[:samples], tablecolumns, ordered:false, write_concern: 0)
					else
						puts "Creating Varians and Calls" if @logging
						# split the line on the tabs
						variant_calls = line.split(/\t/)
						# create the variant object id
						variant_id = BSON::ObjectId.new
						# create the variant object
						variant = { 
								_id: variant_id,
								_key: variant_calls[2],
								names: variant_calls[2].split(/,/),
								chromosome:variant_calls[0].to_i,
								position:variant_calls[1],
								filter:variant_calls[6],
								reference_base:variant_calls[3],
								alternate_bases:variant_calls[4].gsub(/<|>/, '').split(/,/),
								alternate_structure:alt_structure,
								info: {}
							  }
						# walk over each info value
						variant_calls[7].split(';').each do |info_line|
							# cut apart at the ='s
							info_column = info_line.split('=')
							# copy the variant info 
							variant[:info][info_column[0]] = info_structure[info_column[0]].clone
							# find out how to store each value for the infos
							# break apart each info line by the ,  and map all values by type
							val = info_column[1].to_s.split(/\,/).map do |elm|
								case variant[:info][info_column[0]][:type].downcase
									when 'integer'
										elm.to_i
									when 'double', 'float'
										elm.to_f 
									else
										elm.to_s
								end
							end
							variant[:info][info_column[0]][:value] = val if not val.empty?
						end
						# try to insert the variant
						client[:variants].insert_one(variant)

						calls = []
						# sample time
						variant_calls[9..-1].each_with_index do |call, ind|
							# unless the call is ./. or .|.
							unless call =~ /\.(\||\/)\./
								# add the variant, phase (if it's | then phased, if / unphased), and genotype as an integer array
								# store phase and what the genotype is if it differs from the reference
								calls.push({ 
									sample: tablecolumns[ind][:id],
									variant: variant_id, 
									phased: !!(call =~ /\|/), 
									genotype: call.split(/\||\//).map(&:to_i)
								})
							end
						end

						# add the call to our sample
						client[:calls].insert_many(calls)
						# Mongo::BulkWrite::UnorderedBulkWrite.new(:calls, calls, write_concern: 0)
					end
				end
			end
			# lets add the index if one doesn't exist
			client[:variants].indexes.create_one({ _key: 1 }, { name: 'variant_search_index', unique: true })
			# lets add the index if one doesn't exist
			client[:samples].indexes.create_one({ _key: 1 }, { name: 'samples_search_index', unique: true })
		rescue Exception => e
			puts e.message
			puts (e.backtrace or []).join("\n")
		ensure 
			# close our mongo connection
			client.close
		end
	end


	private

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
				structure[elm_ref][val[0].downcase.to_sym] = val[1].to_s.gsub(/^\"|\"?$/, '')
			end
		end
	end

end

opts = Slop.parse do |o|
  o.string '-d', '--database', 'database name'
  o.string '-i', '--ip', 'mongo address', default: 'localhost'
  o.string '-u', '--user', 'user name for database'
  o.string '-p', '--password', 'password for database'
  o.integer '-s', '--max-pool-size', 'maxium number of pools'
end

if ARGV.size >= 2
	# create our obj
	loader = Hapmap_Load.new
	# start loading variants
	loader.load_vcf(opts.arguments[0], opts.to_hash)
else
	puts "please specify a VCF and mongo database"
end