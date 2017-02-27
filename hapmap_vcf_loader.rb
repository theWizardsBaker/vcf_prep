#!/usr/bin/env ruby
require 'pp'

input_file = ARGV[0]

info_structure = {}

alt_structure = {}

table_columns = []

def load_header(header_line, structure)

	# format the header line so that we have an array of list elements
	header_line = header_line.gsub(/<|>/, '').split(/,(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)/)
	# reference element in the info structure
	elm_ref = ''
	header_line.each do |e|
		val = e.split(/=/, 2)
		if val[0] == "ID"
			elm_ref = val[1]
			structure[elm_ref] = {}
		else
			structure[elm_ref][val[0].downcase] = val[1].to_s.gsub(/^\"|\"?$/, '')
		end
	end
	
end

def load_table()

end

File.open(input_file, "r") do |f|
  f.each_line do |line|
  	line = line.chomp
  	# file header
    if line[0..1] == '##'
    	header_line = line.gsub(/^##/, '')
    	header_line = header_line.split(/=/, 2)
    	case header_line[0]
	    	when "INFO"
	    		load_header(header_line[1], info_structure)
	        when "ALT"
	    		load_header(header_line[1], alt_structure)
    	end

    # table header
    elsif line[0] == '#'
    	line = line.gsub(/^#/, '')
    	table_columns = line.split(/\t/)
    # table row
    else
    	line = line.split(/\t/)
    	row = { 
				"_key" => line[2],
				"names" => [line[2]],
				"chromosome"=> line[0],
			    "position"=> line[1],
			    "filter"=> line[6],
			    "reference_base"=> line[3],
			    "alternate_bases"=> line[4].split(/,/),
			    "info" => info_structure.clone
    		  }

    	line[7].split(/;/).each do |e|
    		info = e.split(/=/)
    		row['info'][info[0]]['value'] = info[1]
    	end

    	pp row
    end
  end
end
