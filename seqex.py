### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
	python seqex.py
	--in <FULL_PATH_TO_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	--contig <STRING, name of contig>
	--start <INT, start of region to extract>
	--end <INT, end of region to extract>
					"""

import re

# --- end of imports --- #


def load_multiple_fasta_file( multiple_fasta_file ):
	"""! @brief loads multiple fasta file into dict """
	
	content = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				content.update( { header: seq } )
				header = line.strip()[1:].split(' ')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		content.update( { header: seq } )
	
	return content


def get_seq_part_from_assembly( content, contig_name, start, end ):
	"""! @brief get seq from start to end; returns modified start and end position, if contig is shorter than requested sequence """
	
	if end > len( content[ contig_name ] ):
		if start < 0:
			return content[ contig_name ], 0, len( content[ contig_name ] )
		else:
			return content[ contig_name ][ start: ], start, len( content[ contig_name ] )
	else:
		if start < 0:
			return content[ contig_name ][ :end ], 0, end
		else:
			return content[ contig_name ][ start:end ], start, end


def collect_seq_parts( data, outputfile, content, flank_len ):
	"""! @brief write extracted seqs to file """
	
	with open( outputfile, "w" ) as out:
		for each in data:
			try:
				l_seq, l_start, l_end = get_seq_part_from_assembly( content, each['contig_name'], each['start']-flank_len, each['start'] )	#left is always upstream!!
				m_seq, m_start, m_end = get_seq_part_from_assembly( content, each['contig_name'], each['start'], each['end']+1 )
				r_seq, r_start, r_end = get_seq_part_from_assembly( content, each['contig_name'], each['end']+1, each['end']+flank_len+1 )		#right is always downstream!!
				out.write( '>' + each['contig_name']  +  '_' + str( l_start ) + '_' + str( r_end ) + '\n'  )
				seq = l_seq.lower() + m_seq.upper() + r_seq.lower()		#sequence of interest is highlighted by upper case
				while True:
					if len( seq ) > 80:
						out.write( seq[:81] + '\n' )
						seq = seq[81:]
					else:
						out.write( seq + '\n' )
						break
			except:
				print "failed to collect contig information: " + str( each )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	assembly_file = arguments[ arguments.index('--in')+1 ]
	outputfile = arguments[ arguments.index('--out')+1 ]
	contig = arguments[ arguments.index('--contig')+1 ]
	start = int( arguments[ arguments.index('--start')+1 ] )
	end = int( arguments[ arguments.index('--end')+1 ] )
	
	# --- execution of analysis --- #
	content = load_multiple_fasta_file( assembly_file )
    
	flank_len = 500
	
	collect_seq_parts( [ { 'contig_name': contig, 'start': start, 'end': end } ], outputfile, content, flank_len )


if __name__ == "__main__":
	if '--in' in sys.argv and '--out' in sys.argv and '--contig' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
