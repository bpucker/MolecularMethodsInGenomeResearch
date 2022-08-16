### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.43 ###

# position adjusted to 1-based system in v0.4

__usage__ = """
	python3 seqex3.py
	--in <FULL_PATH_TO_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_FILE>
	--contig <STRING, name of contig>
	--start <INT, start of region to extract>
	--end <INT, end of region to extract>
	
	optional:
	--revcomp <EXTRAC_REVERSE_COMPLEMENT>
					"""

import re, sys

# --- end of imports --- #


def load_multiple_fasta_file( multiple_fasta_file ):
	"""! @brief loads multiple fasta file into dict """
	
	content = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline().strip()[1:]
		if " " in header:
			header = header.split(' ')[0]
		if "\t" in header:
			header = header.split('\t')[0]
		if ":" in header:
			header = header.split(':')[0]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				content.update( { header: "".join( seq ) } )
				header = line.strip()[1:]
				if " " in header:
					header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
				if ":" in header:
					header = header.split(':')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		content.update( { header: "".join( seq ) } )
	
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


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c', 'n': 'n', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N' }
	for nt in seq:
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] )


def collect_seq_parts( data, outputfile, content, flank_len, revcomp_state ):
	"""! @brief write extracted seqs to file """
	
	with open( outputfile, "w" ) as out:
		for each in data:
			try:
				l_seq, l_start, l_end = get_seq_part_from_assembly( content, each['contig_name'], each['start']-flank_len, each['start']-1 )	#left is always upstream!!
				m_seq, m_start, m_end = get_seq_part_from_assembly( content, each['contig_name'], each['start']-1, each['end'] )
				r_seq, r_start, r_end = get_seq_part_from_assembly( content, each['contig_name'], each['end'], each['end']+flank_len+1 )		#right is always downstream!!
				out.write( '>' + each['contig_name']  +  '_' + str( l_start ) + '_' + str( r_end ) + '\n'  )
				seq = l_seq.lower() + m_seq.upper() + r_seq.lower()		#sequence of interest is highlighted by upper case
				if revcomp_state:
					seq = revcomp( seq )
				n=80
				chunks = [ seq[i:i+n] for i in range(0, len(seq), n)]
				for chunk in chunks:
					out.write( chunk + '\n' )
			except:
				print( "failed to collect contig information: " + str( each ) )


def main( arguments ):
	"""! @brief run all parts of this script """
	
	assembly_file = arguments[ arguments.index('--in')+1 ]
	outputfile = arguments[ arguments.index('--out')+1 ]
	contig = arguments[ arguments.index('--contig')+1 ]
	start = int( arguments[ arguments.index('--start')+1 ] )
	end = int( arguments[ arguments.index('--end')+1 ] )
	
	if "--revcomp" in arguments:
		revcomp_state = True
	else:
		revcomp_state = False
	
	# --- execution of analysis --- #
	content = load_multiple_fasta_file( assembly_file )
    
	flank_len = 500
	
	collect_seq_parts( [ { 'contig_name': contig, 'start': start, 'end': end } ], outputfile, content, flank_len, revcomp_state )


if '--in' in sys.argv and '--out' in sys.argv and '--contig' in sys.argv and '--start' in sys.argv and '--end' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
