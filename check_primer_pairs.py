### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
	python check_primer_pairs.py
	--in <FULL_PATH_TO_INPUT_FILE>
	--out <FULL_PATH_TO_OUTPUT_FOLDER>
	--ref <FULL_PATH_TO_FOLDER_WITH_REFERENCE_SEQUENCES>
					"""

import os, glob, sys, re
from operator import itemgetter

# --- end of imports --- #


def uniqueness_check( fw_info, rv_info ):
	"""! @brief run the uniqueness check for the sequences hit by BLASTn """
	
	# --- check for perfect match and uniqueness --- #
	fw_info = sorted( fw_info, key = itemgetter( 'score' ) )
	rv_info = sorted( rv_info, key = itemgetter( 'score' ) )
	if fw_info and rv_info:		#hits for both primers detected
		if fw_info[0]['chr'] == rv_info[0]['chr']:	#check
		
			unique = True
			# --- check unique --- #
			if len( fw_info ) > 1:
				#compare both best hits for fw primer
				if fw_info[0]['score']*0.9 < fw_info[1]['score']:
					unique = False
			if len( rv_info ) > 1:
				#compare both best hits for rv primer
				if rv_info[0]['score']*0.9 < rv_info[1]['score']:
					unique = False
						
			orientation = True
			# ---- check orientation towards each other and distance --- #
			if fw_info[0]['start'] < fw_info[0]['end']:
				if rv_info[0]['start'] < rv_info[0]['end']:
					orientation = False
			else:
				if rv_info[0]['start'] > rv_info[0]['end']:
					orientation = False
			# --- return information about expected PCR product --- #
			if unique and orientation:
				return { fw_info[0]['primer']:	{	'chr': fw_info[0]['chr'],
													'start': min( [ fw_info[0]['start'], fw_info[0]['end'], rv_info[0]['start'], rv_info[0]['end'] ] ),
													'end': max( [ fw_info[0]['start'], fw_info[0]['end'], rv_info[0]['start'], rv_info[0]['end'] ] ),
													'fw_primer': fw_info[0]['primer'],
													'rv_primer': rv_info[0]['primer']
												}
						}
			else:
				return False


def load_and_check_BLASTn_results( BLASTn_result_file, primer_lengths ):
	"""! @brief check uniqueness of both primers and perfect binding in each genome
	
		@return (dict) contains primers, expected PCR product size, position in the given accession
	 """
	 
	results = {}
	
	# --- load information from file --- #
	with open( BLASTn_result_file, "r" ) as f:
		fw_info = []
		rv_info = []
		valid_primer_pairs = {}
		line = f.readline()
		prev_locus = line.split('_%_')[0]
		prev_primer_pair = line.split('_%_')[1]
		while line:
			parts = line.strip().split('\t')
			current_locus = parts[0].split('_%_')[0]
			current_primer_pair = parts[0].split('_%_')[1]
			
			# --- saving results and reset temp variables, when next locus is reached --- #
			if current_locus != prev_locus:
				unique_check_results = uniqueness_check( fw_info, rv_info )
				if unique_check_results:
					valid_primer_pairs.update( unique_check_results )
				
				if valid_primer_pairs:
					results.update( { prev_locus: valid_primer_pairs } )
				
				fw_info = []
				rv_info = []
				valid_primer_pairs = {}
				prev_locus = current_locus
				prev_primer_pair = current_primer_pair
			
			# --- save results and reset temp variables, when next primer pair is reached --- #
			elif current_primer_pair != prev_primer_pair:
				unique_check_results = uniqueness_check( fw_info, rv_info )
				if unique_check_results:
					valid_primer_pairs.update( unique_check_results )
				
				if valid_primer_pairs:
					results.update( { prev_locus: valid_primer_pairs } )
				
				fw_info = []
				rv_info = []
				prev_primer_pair = current_primer_pair
				
				
			
			# --- do normal analysis until end of file is reached --- #
			# --- check hit quality --- #
			
			if int( parts[7] ) == primer_lengths[ parts[0] ]:	#check for binding at three prime end!
				if int( parts[7] ) - int( parts[6] ) > 12:
					# --- entry for hit --- #
					data = { 	'chr': parts[1],
								'start': int( parts[8] ),
								'end': int( parts[9] ),
								'score': float( parts[-1] ),
								'primer': parts[0]
							}
					
					if '_1'  == parts[0][-2:]:
						fw_info.append( data )
					else:
						rv_info.append( data )
			
			line = f.readline()
		# --- final addition of results --- #
		unique_check_results = uniqueness_check( fw_info, rv_info )
		if unique_check_results:
			valid_primer_pairs.update( unique_check_results )
				
		if valid_primer_pairs:
			results.update( { prev_locus: valid_primer_pairs } )
	
	return results


def check_primer_pair( BLASTn_query_file, reference_files, prefix, output_file, mapping_table ):
	"""! @brief checks primer pair in all three accessions """
	
	primer_lengths = get_all_primer_lengths( BLASTn_query_file )
	all_loci = []
	collected_infos = {}
	IDs = []
	for filename in reference_files[:10]:
		ID = filename.split('/')[-1].split('.')[0]
		IDs.append( ID )
				
		# --- run BLASTn --- #
		BLASTn_result_file = prefix + ID + "_BLASTn_result_file.txt"
		if not os.path.isfile( BLASTn_result_file ):
			blast_db = prefix + ID
			os.popen( "makeblastdb -in " + filename + " -out " + blast_db + " -dbtype nucl" )
			os.popen( "blastn -query " + BLASTn_query_file + " -db " + blast_db + " -out " + BLASTn_result_file + " -outfmt 6 -evalue 0.001 -word_size 6 -num_threads 8" )
		
		# --- analyze BLASTn results: check for binding and uniqueness --- #
		information = load_and_check_BLASTn_results( BLASTn_result_file, primer_lengths )
		working = information.keys()
		print "number of valid pairs in " + ID + ": " + str( len( working ) )
		all_loci += working
		collected_infos.update( { ID: information } )

	# --- write resuls into output file --- #
	IDs = sorted( IDs )
	with open( output_file, "w" ) as out:
		# --- generate header --- #
		header = [ "fw_primer", "rv_primer" ]
		for ID in IDs:
			header.append( ID + "_chr" )
			header.append( ID + "_start" )
			header.append( ID + "_end" )
		out.write( "\t".join( header ) + '\n' )
		# --- add entries --- #
		for locus in all_loci:
			new_line = [ ]
			for ID in IDs:
				try:
					locus_details = collected_infos[ ID ][ locus ].values()[0]
					print locus_details
					if len( new_line ) == 0:
						try:
							new_line.append( mapping_table[ locus_details['fw_primer'] ] )
							new_line.append( mapping_table[ locus_details['rv_primer'] ] )
						except KeyError:
							print locus_details
							print locus_details.keys()
							print locus_details['fw_primer']
							print locus_details['rv_primer']
							print "ERROR in mapping table"
					new_line.append( locus_details['chr'] )
					new_line.append( locus_details['start'] )
					new_line.append( locus_details['end'] )
				except KeyError:
					new_line.append( "n/a" )
					new_line.append( "n/a" )
					new_line.append( "n/a"  )
			out.write( '\t'.join( map( str, new_line ) ) + '\n' )


def get_all_primer_lengths( primer_file ):
	"""! @brief get lengthes of all primer sequences """
	
	primer_lengths = {}
	
	with open( primer_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == '>':
				header = line.strip()[1:]
				seq = f.readline().strip()
				primer_lengths.update( { header: len( seq ) } )
			line = f.readline()
	
	return primer_lengths


def modify_primer_file( infile, outfile ):
	"""! @brief add subfix number to all primers """
	
	counter = 0
	mapping_table = {}
	with open( outfile, "w" ) as out:
		with open( infile, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] == ">":
					if counter % 2 == 0:
						name = line.strip()
						out.write( ">primer" + str( counter / 2 ) + '_%_' +  str( counter / 2 )+ '_%_1\n'  )
						mapping_table.update( { "primer" + str( counter / 2 ) + '_%_' + str( counter / 2 ) + '_%_1': name[1:] } )
					else:
						name = line.strip()
						out.write( ">primer" + str( counter / 2 ) + '_%_' +  str( counter / 2 )+ '_%_2\n'  )
						mapping_table.update( { "primer" + str( counter / 2 ) + '_%_' +  str( counter / 2 )+ '_%_2': name[1:] } )
					counter += 1
				else:
					out.write( line )
				line = f.readline()
	return mapping_table


def main( arguments ):
	"""! @brief runs all parts of the primer pair analysis """
	
	prefix = arguments[ arguments.index('--out')+1 ]
	ref_dir = arguments[ arguments.index('--ref')+1 ]
	primer_file = arguments[ arguments.index('--in')+1 ]
	
	if prefix[-1] != '/':
		prefix += '/'
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	if ref_dir[-1] != '/':
		ref_dir += '/'
	
	reference_files = sorted( glob.glob( ref_dir + "*.fasta" )+glob.glob( ref_dir + "*.fas" )+glob.glob( ref_dir + "*.fa" )+glob.glob( ref_dir + "*.fna" ) +glob.glob( ref_dir + "*.FASTA" )+glob.glob( ref_dir + "*.FAS" )+glob.glob( ref_dir + "*.FA" )+glob.glob( ref_dir + "*.FNA" ))
	final_primer_file = prefix + "primers_for_investigation.fasta"
	output_file = prefix + "expected_PCR_products.txt"
	
	# --- modify primer file --- #
	mapping_table = modify_primer_file( primer_file, final_primer_file )
	
	# --- running primer check against all three accessions --- #
	check_primer_pair( final_primer_file, reference_files, prefix, output_file, mapping_table )


if __name__ == '__main__':
	
	if '--out' in sys.argv and '--in' in sys.argv and '--ref' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
