import os.path
import multiprocessing
import utils
import functools
import sys
import pickle


def index_fasta_samtools(fasta, region_None, region_outfile_none, print_comand_True):
	command = ['samtools', 'faidx', fasta, '', '', '']
	shell_true = False
	if region_None is not None:
		command[3] = region_None
	if region_outfile_none is not None:
		command[4] = '>'
		command[5] = region_outfile_none
		shell_true = True
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, shell_true, None, print_comand_True)
	return run_successfully, stdout


# Indexing reference file using Bowtie2
def indexSequenceBowtie2(referenceFile, threads):
	if os.path.isfile(str(referenceFile + '.1.bt2')):
		run_successfully = True
	else:
		command = ['bowtie2-build', '--threads', str(threads), referenceFile, referenceFile]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


# Mapping with Bowtie2
def mappingBowtie2(fastq_files, referenceFile, threads, outdir, conserved_True, numMapLoc):
	sam_file = os.path.join(outdir, str('alignment.sam'))

	# Index reference file
	run_successfully = indexSequenceBowtie2(referenceFile, threads)

	if run_successfully:
		command = ['bowtie2', '-k', str(numMapLoc), '-q', '', '--threads', str(threads), '-x', referenceFile, '', '--no-unal', '-S', sam_file]

		if len(fastq_files) == 1:
			command[9] = '-U ' + fastq_files[0]
		elif len(fastq_files) == 2:
			command[9] = '-1 ' + fastq_files[0] + ' -2 ' + fastq_files[1]
		else:
			return False, None

		if conserved_True:
			command[4] = '--sensitive'
		else:
			command[4] = '--very-sensitive-local'

		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)

	if not run_successfully:
		sam_file = None

	return run_successfully, sam_file


def split_cigar(cigar):
	cigars = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']

	splited_cigars = []
	numbers = ''
	for char in cigar:
		if char not in cigars:
			numbers += char
		else:
			splited_cigars.append([int(numbers), char])
			numbers = ''

	return splited_cigars


def recode_cigar_based_on_base_quality(cigar, bases_quality, softClip_baseQuality, mapping_position, direct_strand_true, softClip_cigarFlagRecode):
	cigar = split_cigar(cigar)
	soft_left = []
	soft_right = []
	cigar_flags_for_reads_length = ('M', 'I', 'S', '=', 'X')
	read_length_without_right_s = sum([cigar_part[0] for cigar_part in cigar if cigar_part[1] in cigar_flags_for_reads_length]) - (cigar[len(cigar) - 1][0] if cigar[len(cigar) - 1][1] == 'S' else 0)
	for x, base in enumerate(bases_quality):
		if ord(base) - 33 >= softClip_baseQuality:
			if x <= cigar[0][0] - 1:
				if cigar[0][1] == 'S':
					soft_left.append(x)
			elif x > read_length_without_right_s - 1:
				if cigar[len(cigar) - 1][1] == 'S':
					soft_right.append(x)

	left_changed = (False, 0)
	if len(soft_left) > 0:
		soft_left = min(soft_left) + 1
		if soft_left == 1:
			cigar = [[cigar[0][0], softClip_cigarFlagRecode]] + cigar[1:]
			left_changed = (True, cigar[0][0])
		elif cigar[0][0] - soft_left > 0:
			cigar = [[soft_left, 'S']] + [[cigar[0][0] - soft_left, softClip_cigarFlagRecode]] + cigar[1:]
			left_changed = (True, cigar[0][0] - soft_left)

	right_changed = (False, 0)
	if len(soft_right) > 0:
		soft_right = max(soft_right) + 1
		cigar = cigar[:-1]
		if soft_right - read_length_without_right_s > 0:
			cigar.append([soft_right - read_length_without_right_s, softClip_cigarFlagRecode])
			right_changed = (True, soft_right - read_length_without_right_s)
		if len(bases_quality) - soft_right > 0:
			cigar.append([len(bases_quality) - soft_right, 'S'])

	if left_changed[0]:
		# if direct_strand_true:
		mapping_position = mapping_position - left_changed[1]
	# if right_changed[0]:
	# 	if not direct_strand_true:
	# 		mapping_position = mapping_position + right_changed[1]

	return ''.join([str(cigar_part[0]) + cigar_part[1] for cigar_part in cigar]), str(mapping_position)


def verify_is_forward_read(sam_flag_bit):
	# 64 = 1000000
	forward_read = False
	bit = format(sam_flag_bit, 'b').zfill(7)
	if bit[-7] == '1':
		forward_read = True
	return forward_read


def verify_mapped_direct_strand(sam_flag_bit):
	# 16 = 10000 -> mapped in reverse strand
	direct_strand = False
	bit = format(sam_flag_bit, 'b').zfill(5)
	if bit[-5] == '0':
		direct_strand = True
	return direct_strand


def verify_mapped_tip(reference_length, mapping_position, read_length, cigar):
	tip = False
	if 'S' in cigar:
		cigar = split_cigar(cigar)
		if cigar[0][1] == 'S':
			if mapping_position - cigar[0][0] < 0:
				tip = True
		if cigar[len(cigar) - 1][1] == 'S':
			if mapping_position + cigar[len(cigar) - 1][0] > reference_length:
				tip = True
	return tip


def change_sam_flag_bit_mapped_reverse_strand_2_direct_strand(sam_flag_bit):
	bit = list(format(sam_flag_bit, 'b').zfill(5))
	bit[-5] = '0'
	return int(''.join(bit), 2)


def change_sam_flag_bit_mate_reverse_strand_2_direct_strand(sam_flag_bit):
	bit = list(format(sam_flag_bit, 'b').zfill(6))
	bit[-6] = '0'
	return int(''.join(bit), 2)


def move_read_mapped_reverse_strand_2_direct_strand(seq, bases_quality, sam_flag_bit, cigar):
	seq = utils.reverse_complement(seq)
	bases_quality = ''.join(reversed(list(bases_quality)))
	sam_flag_bit = change_sam_flag_bit_mapped_reverse_strand_2_direct_strand(sam_flag_bit)
	cigar = ''.join([str(cigar_part[0]) + cigar_part[1] for cigar_part in reversed(split_cigar(cigar))])
	return seq, bases_quality, str(sam_flag_bit), cigar


@utils.trace_unhandled_exceptions
def parallelized_recode_soft_clipping(line_collection, pickleFile, softClip_baseQuality, sequences_length, softClip_cigarFlagRecode):
	lines_sam = []
	for line in line_collection:
		line = line.splitlines()[0]
		if len(line) > 0:
			if line.startswith('@'):
				lines_sam.append(line)
			else:
				line = line.split('\t')
				if not verify_mapped_tip(sequences_length[line[2]], int(line[3]), len(line[9]), line[5]):
					line[5], line[3] = recode_cigar_based_on_base_quality(line[5], line[10], softClip_baseQuality, int(line[3]), verify_mapped_direct_strand(int(line[1])), softClip_cigarFlagRecode)
				lines_sam.append('\t'.join(line))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(lines_sam, writer)


def recode_soft_clipping_from_sam(sam_file, outdir, threads, softClip_baseQuality, reference_dict, softClip_cigarFlagRecode):
	pickle_files = []
	sequences_length = {}
	for x, seq_info in reference_dict.items():
		sequences_length[seq_info['header']] = seq_info['length']

	with open(sam_file, 'rtU') as reader:
		pool = multiprocessing.Pool(processes=threads)
		line_collection = []
		x = 0
		for x, line in enumerate(reader):
			line_collection.append(line)
			if x % 10000 == 0:
				pickleFile = os.path.join(outdir, 'remove_soft_clipping.' + str(x) + '.pkl')
				pickle_files.append(pickleFile)
				pool.apply_async(parallelized_recode_soft_clipping, args=(line_collection, pickleFile, softClip_baseQuality, sequences_length, softClip_cigarFlagRecode,))
				line_collection = []
		if len(line_collection) > 0:
			pickleFile = os.path.join(outdir, 'remove_soft_clipping.' + str(x) + '.pkl')
			pickle_files.append(pickleFile)
			pool.apply_async(parallelized_recode_soft_clipping, args=(line_collection, pickleFile, softClip_baseQuality, sequences_length, softClip_cigarFlagRecode,))
			line_collection = []
		pool.close()
		pool.join()

	os.remove(sam_file)

	new_sam_file = os.path.join(outdir, 'alignment_with_soft_clipping_recoded.sam')
	with open(new_sam_file, 'wt') as writer:
		for pickleFile in pickle_files:
			if os.path.isfile(pickleFile):
				lines_sam = None
				with open(pickleFile, 'rb') as reader:
					lines_sam = pickle.load(reader)
				if lines_sam is not None:
					for line in lines_sam:
						writer.write(line + '\n')
				os.remove(pickleFile)

	return new_sam_file


# Sort alignment file
def sortAlignment(alignment_file, output_file, sortByName_True, threads):
	outFormat_string = os.path.splitext(output_file)[1][1:].lower()
	command = ['samtools', 'sort', '-o', output_file, '-O', outFormat_string, '', '-@', str(threads), alignment_file]
	if sortByName_True:
		command[6] = '-n'
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	if not run_successfully:
		output_file = None
	return run_successfully, output_file


# Index alignment file
def indexAlignment(alignment_file):
	command = ['samtools', 'index', alignment_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, True)
	return run_successfully


def mapping_reads(fastq_files, reference_file, threads, outdir, conserved_True, numMapLoc, rematch_run, softClip_baseQuality, softClip_recodeRun, reference_dict, softClip_cigarFlagRecode):
	# Create a symbolic link to the reference_file
	reference_link = os.path.join(outdir, os.path.basename(reference_file))
	os.symlink(reference_file, reference_link)

	bam_file = None
	# Mapping reads using Bowtie2
	run_successfully, sam_file = mappingBowtie2(fastq_files, reference_link, threads, outdir, conserved_True, numMapLoc)

	if run_successfully:
		# Remove soft clipping
		if rematch_run == softClip_recodeRun or softClip_recodeRun == 'both':
			print 'Recoding soft clipped regions'
			sam_file = recode_soft_clipping_from_sam(sam_file, outdir, threads, softClip_baseQuality, reference_dict, softClip_cigarFlagRecode)

		# Convert sam to bam and sort bam
		run_successfully, bam_file = sortAlignment(sam_file, str(os.path.splitext(sam_file)[0] + '.bam'), False, threads)

		if run_successfully:
			os.remove(sam_file)
			# Index bam
			run_successfully = indexAlignment(bam_file)

	return run_successfully, bam_file, reference_link


def create_vcf(bam_file, sequence_to_analyse, outdir, counter, reference_file):
	gene_vcf = os.path.join(outdir, 'samtools_mpileup.sequence_' + str(counter) + '.vcf')

	command = ['samtools', 'mpileup', '--count-orphans', '--no-BAQ', '--min-BQ', '0', '--min-MQ', str(7), '--fasta-ref', reference_file, '--region', sequence_to_analyse, '--output', gene_vcf, '--VCF', '--uncompressed', '--output-tags', 'INFO/AD,AD,DP', bam_file]

	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
	if not run_successfully:
		gene_vcf = None
	return run_successfully, gene_vcf


# Read vcf file
class Vcf():
	def __init__(self, vcfFile):
		self.vcf = open(vcfFile, 'rtU')
		self.line_read = self.vcf.readline()
		while self.line_read.startswith('#'):
			self.line_read = self.vcf.readline()
		self.line = self.line_read

	def readline(self):
		self.line_stored = self.line
		self.line = self.vcf.readline()
		return self.line_stored

	def close(self):
		self.vcf.close()


def get_variants(gene_vcf):
	variants = {}

	vfc_file = Vcf(gene_vcf)
	line = vfc_file.readline()
	while len(line) > 0:
		fields = line.splitlines()[0].split('\t')
		if len(fields) > 0:
			fields[1] = int(fields[1])

			info_field = {}
			for i in fields[7].split(';'):
				i = i.split('=')
				if len(i) > 1:
					info_field[i[0]] = i[1]
				else:
					info_field[i[0]] = None

			format_field = {}
			format_field_name = fields[8].split(':')
			format_data = fields[9].split(':')

			for i in range(0, len(format_data)):
				format_field[format_field_name[i]] = format_data[i].split(',')

			fields_to_store = {'REF': fields[3], 'ALT': fields[4].split(','), 'info': info_field, 'format': format_field}
			if fields[1] in variants:
				variants[fields[1]][len(variants[fields[1]])] = fields_to_store
			else:
				variants[fields[1]] = {0: fields_to_store}

		line = vfc_file.readline()
	vfc_file.close()

	return variants


def indel_entry(variant_position):
	entry_with_indel = []
	entry_with_snp = None
	for i in variant_position:
		keys = variant_position[i]['info'].keys()
		if 'INDEL' in keys:
			entry_with_indel.append(i)
		else:
			entry_with_snp = i

	return entry_with_indel, entry_with_snp


def get_alt_noMatter(variant_position, indel_true):
	dp = sum(map(int, variant_position['format']['AD']))
	index_alleles_sorted_position = sorted(zip(map(int, variant_position['format']['AD']), range(0, len(variant_position['format']['AD']))), reverse=True)
	index_dominant_allele = None
	if not indel_true:
		ad_idv = index_alleles_sorted_position[0][0]

		if len([x for x in index_alleles_sorted_position if x[0] == ad_idv]) > 1:
			index_alleles_sorted_position = sorted([x for x in index_alleles_sorted_position if x[0] == ad_idv])

		index_dominant_allele = index_alleles_sorted_position[0][1]
		if index_dominant_allele == 0:
			alt = '.'
		else:
			alt = variant_position['ALT'][index_dominant_allele - 1]

	else:
		ad_idv = int(variant_position['info']['IDV'])

		if float(ad_idv) / float(dp) >= 0.5:
			if len([x for x in index_alleles_sorted_position if x[0] == index_alleles_sorted_position[0][0]]) > 1:
				index_alleles_sorted_position = sorted([x for x in index_alleles_sorted_position if x[0] == index_alleles_sorted_position[0][0]])

			index_dominant_allele = index_alleles_sorted_position[0][1]
			if index_dominant_allele == 0:
				alt = '.'
			else:
				alt = variant_position['ALT'][index_dominant_allele - 1]
		else:
			ad_idv = int(variant_position['format']['AD'][0])
			alt = '.'

	return alt, dp, ad_idv, index_dominant_allele


def count_number_diferences(ref, alt):
	number_diferences = 0

	if len(ref) != len(alt):
		number_diferences += 1

	for i in range(0, min(len(ref), len(alt))):
		if alt[i] != 'N' and ref[i] != alt[i]:
			number_diferences += 1

	return number_diferences


def get_alt_correct(variant_position, alt_noMatter, dp, ad_idv, index_dominant_allele, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	alt = None
	low_coverage = False
	multiple_alleles = False

	if dp >= minimum_depth_presence:
		if dp < minimum_depth_call:
			alt = 'N' * len(variant_position['REF'])
			low_coverage = True
		else:
			if ad_idv < minimum_depth_call:
				alt = 'N' * len(variant_position['REF'])
				low_coverage = True
			else:
				if float(ad_idv) / float(dp) < minimum_depth_frequency_dominant_allele:
					if index_dominant_allele is not None:
						variants_coverage = [int(variant_position['format']['AD'][i]) for i in range(0, len(variant_position['ALT']) + 1) if i != index_dominant_allele]
						if sum(variants_coverage) == 0:
							alt = alt_noMatter
						else:
							if float(max(variants_coverage)) / float(sum(variants_coverage)) > 0.5:
								multiple_alleles = True
								alt = 'N' * len(variant_position['REF'])
							elif float(max(variants_coverage)) / float(sum(variants_coverage)) == 0.5 and len(variants_coverage) > 2:
								multiple_alleles = True
								alt = 'N' * len(variant_position['REF'])
							else:
								alt = alt_noMatter
					else:
						multiple_alleles = True
						alt = 'N' * len(variant_position['REF'])
				else:
					alt = alt_noMatter
	else:
		low_coverage = True

	return alt, low_coverage, multiple_alleles


def get_alt_alignment(ref, alt):
	if alt is None:
		alt = 'N' * len(ref)
	else:
		if len(ref) != len(alt):
			if len(alt) < len(ref):
				if alt == '.':
					alt = ref
				alt += 'N' * (len(ref) - len(alt))
			else:
				if alt[:len(ref)] == ref:
					alt = '.'
				else:
					alt = alt[:len(ref)]

	return alt


def get_indel_more_likely(variant_position, indels_entry):
	indel_coverage = {}
	for i in indels_entry:
		indel_coverage[i] = int(variant_position['info']['IDV'])
	return indel_coverage.index(str(max(indel_coverage.values())))


def determine_variant(variant_position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, indel_true):
	alt_noMatter, dp, ad_idv, index_dominant_allele = get_alt_noMatter(variant_position, indel_true)

	alt_correct, low_coverage, multiple_alleles = get_alt_correct(variant_position, alt_noMatter, dp, ad_idv, index_dominant_allele, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)

	alt_alignment = get_alt_alignment(variant_position['REF'], alt_correct)

	return variant_position['REF'], alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment


def confirm_nucleotides_indel(ref, alt, variants, position_start_indel, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, alignment_true):
	alt = list(alt)

	for i in range(0, len(alt) - 1):
		if len(alt) < len(ref):
			new_position = position_start_indel + len(alt) - i - 1
			alt_position = len(alt) - i - 1
		else:
			if i + 1 > len(ref):
				break
			new_position = position_start_indel + 1 + i
			alt_position = 1 + i

		if alt[alt_position] != 'N':
			if new_position not in variants:
				if alignment_true:
					alt[alt_position] = 'N'
				else:
					alt = alt[: alt_position]
					break

			entry_with_indel, entry_with_snp = indel_entry(variants[new_position])
			new_ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[new_position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
			if alt_noMatter != '.' and alt[alt_position] != alt_noMatter:
				alt[alt_position] = alt_noMatter

	return ''.join(alt)


def snp_indel(variants, position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele):
	entry_with_indel, entry_with_snp = indel_entry(variants[position])

	if len(entry_with_indel) == 0:
		ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
	else:
		ref_snp, alt_correct_snp, low_coverage_snp, multiple_alleles_snp, alt_noMatter_snp, alt_alignment_snp = determine_variant(variants[position][entry_with_snp], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)

		indel_more_likely = entry_with_indel[0]
		if len(entry_with_indel) > 1:
			indel_more_likely = get_indel_more_likely(variants[position], entry_with_indel)

		ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = determine_variant(variants[position][indel_more_likely], minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, True)

		if alt_noMatter == '.':
			ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = ref_snp, alt_correct_snp, low_coverage_snp, multiple_alleles_snp, alt_noMatter_snp, alt_alignment_snp
		else:
			if alt_correct is None and alt_correct_snp is not None:
				alt_correct = alt_correct_snp
			elif alt_correct is not None and alt_correct_snp is not None:
				if alt_correct_snp != '.' and alt_correct[0] != alt_correct_snp:
					alt_correct = alt_correct_snp + alt_correct[1:] if len(alt_correct) > 1 else alt_correct_snp
			if alt_noMatter_snp != '.' and alt_noMatter[0] != alt_noMatter_snp:
				alt_noMatter = alt_noMatter_snp + alt_noMatter[1:] if len(alt_noMatter) > 1 else alt_noMatter_snp
			if alt_alignment_snp != '.' and alt_alignment[0] != alt_alignment_snp:
				alt_alignment = alt_alignment_snp + alt_alignment[1:] if len(alt_alignment) > 1 else alt_alignment_snp

			# if alt_noMatter != '.':
			# 	alt_noMatter = confirm_nucleotides_indel(ref, alt_noMatter, variants, position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
			# if alt_correct is not None and alt_correct != '.':
			# 	alt_correct = confirm_nucleotides_indel(ref, alt_correct, variants, position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, False)
			# if alt_alignment != '.':
			# 	alt_alignment = confirm_nucleotides_indel(ref, alt_alignment, variants, position, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, True)

	return ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment


def get_true_variants(variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence):
	variants_correct = {}
	variants_noMatter = {}
	variants_alignment = {}

	correct_absent_positions = {}
	correct_last_absent_position = ''

	noMatter_absent_positions = {}
	noMatter_last_absent_position = ''

	multiple_alleles_found = []

	counter = 1
	while counter <= len(sequence):
		if counter in variants:
			noMatter_last_absent_position = ''

			ref, alt_correct, low_coverage, multiple_alleles, alt_noMatter, alt_alignment = snp_indel(variants, counter, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele)

			if alt_alignment != '.':
				variants_alignment[counter] = {'REF': ref, 'ALT': alt_alignment}

			if alt_noMatter != '.':
				variants_noMatter[counter] = {'REF': ref, 'ALT': alt_noMatter}

			if alt_correct is None:
				if counter - len(correct_last_absent_position) in correct_absent_positions:
					correct_absent_positions[counter - len(correct_last_absent_position)]['REF'] += ref
				else:
					correct_absent_positions[counter] = {'REF': ref, 'ALT': ''}
				correct_last_absent_position += ref
			else:
				if alt_correct != '.':
					if len(alt_correct) < len(ref):
						if len(alt_correct) == 1:
							correct_absent_positions[counter + 1] = {'REF': ref[1:], 'ALT': ''}
						else:
							correct_absent_positions[counter + 1] = {'REF': ref[1:], 'ALT': alt_correct[1:]}

						correct_last_absent_position = ref[1:]
					else:
						variants_correct[counter] = {'REF': ref, 'ALT': alt_correct}
						correct_last_absent_position = ''
				else:
					correct_last_absent_position = ''

			if multiple_alleles:
				multiple_alleles_found.append(counter)

			counter += len(ref)
		else:
			variants_alignment[counter] = {'REF': sequence[counter - 1], 'ALT': 'N'}

			if counter - len(correct_last_absent_position) in correct_absent_positions:
				correct_absent_positions[counter - len(correct_last_absent_position)]['REF'] += sequence[counter - 1]
			else:
				correct_absent_positions[counter] = {'REF': sequence[counter - 1], 'ALT': ''}
			correct_last_absent_position += sequence[counter - 1]

			if counter - len(noMatter_last_absent_position) in noMatter_absent_positions:
				noMatter_absent_positions[counter - len(noMatter_last_absent_position)]['REF'] += sequence[counter - 1]
			else:
				noMatter_absent_positions[counter] = {'REF': sequence[counter - 1], 'ALT': ''}
			noMatter_last_absent_position += sequence[counter - 1]

			counter += 1

	for position in correct_absent_positions:
		if position == 1:
			variants_correct[position] = {'REF': correct_absent_positions[position]['REF'], 'ALT': 'N'}
		else:
			if position - 1 not in variants_correct:
				variants_correct[position - 1] = {'REF': sequence[position - 2] + correct_absent_positions[position]['REF'], 'ALT': sequence[position - 2] + correct_absent_positions[position]['ALT']}
			else:
				variants_correct[position - 1] = {'REF': variants_correct[position - 1]['REF'] + correct_absent_positions[position]['REF'][len(variants_correct[position - 1]['REF']) - 1:], 'ALT': variants_correct[position - 1]['ALT'] + correct_absent_positions[position]['ALT'][len(variants_correct[position - 1]['ALT']) - 1 if len(variants_correct[position - 1]['ALT']) > 0 else 0:]}

	for position in noMatter_absent_positions:
		if position == 1:
			variants_noMatter[position] = {'REF': noMatter_absent_positions[position]['REF'], 'ALT': 'N'}
		else:
			if position - 1 not in variants_noMatter:
				variants_noMatter[position - 1] = {'REF': sequence[position - 2] + noMatter_absent_positions[position]['REF'], 'ALT': sequence[position - 2] + noMatter_absent_positions[position]['ALT']}
			else:
				variants_noMatter[position - 1] = {'REF': variants_noMatter[position - 1]['REF'] + noMatter_absent_positions[position]['REF'][len(variants_noMatter[position - 1]['REF']) - 1:], 'ALT': variants_noMatter[position - 1]['ALT'] + noMatter_absent_positions[position]['ALT'][len(variants_noMatter[position - 1]['ALT']) - 1 if len(variants_noMatter[position - 1]['ALT']) > 0 else 0:]}

	return variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found


def clean_variant_in_extra_seq_left(variant_dict, position, length_extra_seq, multiple_alleles_found, number_multi_alleles):
	number_diferences = 0

	if position + len(variant_dict[position]['REF']) - 1 > length_extra_seq:
		if multiple_alleles_found is not None and position in multiple_alleles_found:
			number_multi_alleles += 1

		temp_variant = variant_dict[position]
		del variant_dict[position]
		variant_dict[length_extra_seq] = {}
		variant_dict[length_extra_seq]['REF'] = temp_variant['REF'][length_extra_seq - position:]
		variant_dict[length_extra_seq]['ALT'] = temp_variant['ALT'][length_extra_seq - position:] if len(temp_variant['ALT']) > length_extra_seq - position else temp_variant['REF'][length_extra_seq - position]
		number_diferences = count_number_diferences(variant_dict[length_extra_seq]['REF'], variant_dict[length_extra_seq]['ALT'])
	else:
		del variant_dict[position]

	return variant_dict, number_multi_alleles, number_diferences


def clean_variant_in_extra_seq_rigth(variant_dict, position, sequence_length, length_extra_seq):
	if position + len(variant_dict[position]['REF']) - 1 > sequence_length - length_extra_seq:
		variant_dict[position]['REF'] = variant_dict[position]['REF'][: - (position - (sequence_length - length_extra_seq)) + 1]
		variant_dict[position]['ALT'] = variant_dict[position]['ALT'][: - (position - (sequence_length - length_extra_seq)) + 1] if len(variant_dict[position]['ALT']) >= - (position - (sequence_length - length_extra_seq)) + 1 else variant_dict[position]['ALT']

	number_diferences = count_number_diferences(variant_dict[position]['REF'], variant_dict[position]['ALT'])

	return variant_dict, number_diferences


def cleanning_variants_extra_seq(variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found, length_extra_seq, sequence_length):
	number_multi_alleles = 0
	number_diferences = 0

	counter = 1
	while counter <= sequence_length:
		if counter <= length_extra_seq:
			if counter in variants_correct:
				variants_correct, number_multi_alleles, number_diferences = clean_variant_in_extra_seq_left(variants_correct, counter, length_extra_seq, multiple_alleles_found, number_multi_alleles)
			if counter in variants_noMatter:
				variants_noMatter, ignore, ignore = clean_variant_in_extra_seq_left(variants_noMatter, counter, length_extra_seq, None, None)
			if counter in variants_alignment:
				variants_alignment, ignore, ignore = clean_variant_in_extra_seq_left(variants_alignment, counter, length_extra_seq, None, None)
		elif counter > length_extra_seq and counter <= sequence_length - length_extra_seq:
			if counter in variants_correct:
				if counter in multiple_alleles_found:
					number_multi_alleles += 1
				variants_correct, number_diferences_found = clean_variant_in_extra_seq_rigth(variants_correct, counter, sequence_length, length_extra_seq)
				number_diferences += number_diferences_found
			if counter in variants_noMatter:
				variants_noMatter, ignore = clean_variant_in_extra_seq_rigth(variants_noMatter, counter, sequence_length, length_extra_seq)
			if counter in variants_alignment:
				variants_alignment, ignore = clean_variant_in_extra_seq_rigth(variants_alignment, counter, sequence_length, length_extra_seq)
		else:
			if counter in variants_correct:
				del variants_correct[counter]
			if counter in variants_noMatter:
				del variants_noMatter[counter]
			if counter in variants_alignment:
				del variants_alignment[counter]

		counter += 1

	return variants_correct, variants_noMatter, variants_alignment, number_multi_alleles, number_diferences


def get_coverage(gene_coverage):
	coverage = {}

	with open(gene_coverage, 'rtU') as reader:
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				line = line.split('\t')
				coverage[int(line[1])] = int(line[2])

	return coverage


def get_coverage_report(coverage, sequence_length, minimum_depth_presence, minimum_depth_call, length_extra_seq):
	if len(coverage) == 0:
		return sequence_length - 2 * length_extra_seq, 100.0, 0.0

	count_absent = 0
	count_lowCoverage = 0
	sum_coverage = 0

	counter = 1
	while counter <= sequence_length:
		if counter > length_extra_seq and counter <= sequence_length - length_extra_seq:
			if coverage[counter] < minimum_depth_presence:
				count_absent += 1
			else:
				if coverage[counter] < minimum_depth_call:
					count_lowCoverage += 1
				sum_coverage += coverage[counter]
		counter += 1

	mean_coverage = 0
	percentage_lowCoverage = 100
	if sequence_length - 2 * length_extra_seq - count_absent > 0:
		mean_coverage = float(sum_coverage) / float(sequence_length - 2 * length_extra_seq - count_absent)
		percentage_lowCoverage = float(count_lowCoverage) / float(sequence_length - 2 * length_extra_seq - count_absent) * 100

	return count_absent, percentage_lowCoverage, mean_coverage


# Get genome coverage data
def compute_genome_coverage_data(alignment_file, sequence_to_analyse, outdir, counter):
	genome_coverage_data_file = os.path.join(outdir, 'samtools_depth.sequence_' + str(counter) + '.tab')
	command = ['samtools', 'depth', '-a', '-q', '0', '-r', sequence_to_analyse, alignment_file, '>', genome_coverage_data_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, True, None, False)
	return run_successfully, genome_coverage_data_file


def write_variants_vcf(variants, outdir, sequence_to_analyse, sufix):
	vcf_file = os.path.join(outdir, str(sequence_to_analyse + '.' + sufix + '.vcf'))
	with open(vcf_file, 'wt') as writer:
		writer.write('##fileformat=VCFv4.2' + '\n')
		writer.write('#' + '\t'.join(['SEQUENCE', 'POSITION', 'ID_unused', 'REFERENCE_sequence', 'ALTERNATIVE_sequence', 'QUALITY_unused', 'FILTER_unused', 'INFO_unused', 'FORMAT_unused']) + '\n')
		for i in sorted(variants.keys()):
			writer.write('\t'.join([sequence_to_analyse, str(i), '.', variants[i]['REF'], variants[i]['ALT'], '.', '.', '.', '.']) + '\n')

	compressed_vcf_file = vcf_file + '.gz'
	command = ['bcftools', 'convert', '-o', compressed_vcf_file, '-O', 'z', vcf_file]
	run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
	if run_successfully:
		command = ['bcftools', 'index', compressed_vcf_file]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)

	if not run_successfully:
		compressed_vcf_file = None

	return run_successfully, compressed_vcf_file


def parse_fasta_inMemory(fasta_memory):
	fasta_memory = fasta_memory.splitlines()
	sequence_dict = {}
	for line in fasta_memory:
		if len(line) > 0:
			if line.startswith('>'):
				sequence_dict = {'header': line[1:], 'sequence': ''}
			else:
				sequence_dict['sequence'] += line

	return sequence_dict


def compute_consensus_sequence(reference_file, sequence_to_analyse, compressed_vcf_file, outdir, sufix):
	sequence_dict = None

	gene_fasta = os.path.join(outdir, str(sequence_to_analyse + '.fasta'))

	run_successfully, stdout = index_fasta_samtools(reference_file, sequence_to_analyse, gene_fasta, False)
	if run_successfully:
		command = ['bcftools', 'consensus', '-f', gene_fasta, compressed_vcf_file]
		run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
		if run_successfully:
			sequence_dict = parse_fasta_inMemory(stdout)

	return run_successfully, sequence_dict


def create_sample_consensus_sequence(outdir, sequence_to_analyse, reference_file, variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence, length_extra_seq):
	variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found = get_true_variants(variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence)

	variants_correct, variants_noMatter, variants_alignment, number_multi_alleles, number_diferences = cleanning_variants_extra_seq(variants_correct, variants_noMatter, variants_alignment, multiple_alleles_found, length_extra_seq, len(sequence))

	run_successfully = False
	consensus = {'correct': {}, 'noMatter': {}, 'alignment': {}}
	for variant_type in ['variants_correct', 'variants_noMatter', 'variants_alignment']:
		run_successfully, compressed_vcf_file = write_variants_vcf(eval(variant_type), outdir, sequence_to_analyse, variant_type.split('_', 1)[1])
		if run_successfully:
			run_successfully, sequence_dict = compute_consensus_sequence(reference_file, sequence_to_analyse, compressed_vcf_file, outdir, variant_type.split('_', 1)[1])
			if run_successfully:
				consensus[variant_type.split('_', 1)[1]] = {'header': sequence_dict['header'], 'sequence': sequence_dict['sequence'][length_extra_seq:len(sequence_dict['sequence']) - length_extra_seq]}

	return run_successfully, number_multi_alleles, consensus, number_diferences


@utils.trace_unhandled_exceptions
def analyse_sequence_data(bam_file, sequence_information, outdir, counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, threads):
	count_absent = None
	percentage_lowCoverage = None
	meanCoverage = None
	number_diferences = 0

	# Create vcf file (for multiple alleles check)
	run_successfully, gene_vcf = create_vcf(bam_file, sequence_information['header'], outdir, counter, reference_file)

	if run_successfully:
		# Create coverage tab file
		run_successfully, gene_coverage = compute_genome_coverage_data(bam_file, sequence_information['header'], outdir, counter)

		if run_successfully:
			variants = get_variants(gene_vcf)

			coverage = get_coverage(gene_coverage)

			run_successfully, number_multi_alleles, consensus_sequence, number_diferences = create_sample_consensus_sequence(outdir, sequence_information['header'], reference_file, variants, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, sequence_information['sequence'], length_extra_seq)

			count_absent, percentage_lowCoverage, meanCoverage = get_coverage_report(coverage, sequence_information['length'], minimum_depth_presence, minimum_depth_call, length_extra_seq)

	utils.saveVariableToPickle([run_successfully, counter, number_multi_alleles, count_absent, percentage_lowCoverage, meanCoverage, consensus_sequence, number_diferences], outdir, str('coverage_info.' + str(counter)))


def get_sequence_information(fasta_file, length_extra_seq):
	sequence_dict = {}
	headers = []

	with open(fasta_file, 'rtU') as reader:
		blank_line_found = False
		sequence_counter = 0
		temp_sequence_dict = {}
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				if not blank_line_found:
					if line.startswith('>'):
						if len(temp_sequence_dict) > 0:
							if temp_sequence_dict.values()[0]['length'] - 2 * length_extra_seq > 0:
								sequence_dict[temp_sequence_dict.keys()[0]] = temp_sequence_dict.values()[0]
								headers.append(temp_sequence_dict.values()[0]['header'])
							else:
								print temp_sequence_dict.values()[0]['header'] + ' sequence ignored due to length <= 0'
							temp_sequence_dict = {}

						if line[1:] in headers:
							sys.exit('Found duplicated sequence headers')

						sequence_counter += 1
						temp_sequence_dict[sequence_counter] = {'header': line[1:], 'sequence': '', 'length': 0}
					else:
						temp_sequence_dict[sequence_counter]['sequence'] += line
						temp_sequence_dict[sequence_counter]['length'] += len(line)
				else:
					sys.exit('It was found a blank line between the fasta file above line ' + line)
			else:
				blank_line_found = True

		if len(temp_sequence_dict) > 0:
			if temp_sequence_dict.values()[0]['length'] - 2 * length_extra_seq > 0:
				sequence_dict[temp_sequence_dict.keys()[0]] = temp_sequence_dict.values()[0]
				headers.append(temp_sequence_dict.values()[0]['header'])
			else:
				print temp_sequence_dict.values()[0]['header'] + ' sequence ignored due to length <= 0'

	return sequence_dict, headers


def determine_threads_2_use(number_sequences, threads):
	if number_sequences >= threads:
		return 1
	else:
		return threads / number_sequences


def sequence_data(sample, reference_file, bam_file, outdir, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, debug_mode_true):
	sequence_data_outdir = os.path.join(outdir, 'sequence_data', '')
	utils.removeDirectory(sequence_data_outdir)
	os.mkdir(sequence_data_outdir)

	sequences, headers = get_sequence_information(reference_file, length_extra_seq)

	threads_2_use = determine_threads_2_use(len(sequences), threads)

	pool = multiprocessing.Pool(processes=threads)
	for sequence_counter in sequences:
		sequence_dir = os.path.join(sequence_data_outdir, str(sequence_counter), '')
		utils.removeDirectory(sequence_dir)
		os.makedirs(sequence_dir)
		pool.apply_async(analyse_sequence_data, args=(bam_file, sequences[sequence_counter], sequence_dir, sequence_counter, reference_file, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, threads_2_use,))
	pool.close()
	pool.join()

	run_successfully, sample_data, consensus_files, consensus_sequences = gather_data_together(sample, sequence_data_outdir, sequences, outdir.rsplit('/', 2)[0], debug_mode_true, length_extra_seq)

	return run_successfully, sample_data, consensus_files, consensus_sequences


def chunkstring(string, length):
	return (string[0 + i:length + i] for i in range(0, len(string), length))


def write_consensus(outdir, sample, consensus_sequence):
	consensus_files = {}
	for consensus_type in ['correct', 'noMatter', 'alignment']:
		consensus_files[consensus_type] = os.path.join(outdir, str(sample + '.' + consensus_type + '.fasta'))
		with open(consensus_files[consensus_type], 'at') as writer:
			writer.write('>' + consensus_sequence[consensus_type]['header'] + '\n')
			fasta_sequence_lines = chunkstring(consensus_sequence[consensus_type]['sequence'], 80)
			for line in fasta_sequence_lines:
				writer.write(line + '\n')
	return consensus_files


def gather_data_together(sample, data_directory, sequences_information, outdir, debug_mode_true, length_extra_seq):
	run_successfully = True
	counter = 0
	sample_data = {}

	consensus_files = None
	consensus_sequences_together = {'correct': {}, 'noMatter': {}, 'alignment': {}}

	write_consensus_first_time = True

	genes_directories = [d for d in os.listdir(data_directory) if not d.startswith('.') and os.path.isdir(os.path.join(data_directory, d, ''))]
	for gene_dir in genes_directories:
		gene_dir_path = os.path.join(data_directory, gene_dir, '')

		files = [f for f in os.listdir(gene_dir_path) if not f.startswith('.') and os.path.isfile(os.path.join(gene_dir_path, f))]
		for file_found in files:
			if file_found.startswith('coverage_info.') and file_found.endswith('.pkl'):
				file_path = os.path.join(gene_dir_path, file_found)

				if run_successfully:
					run_successfully, sequence_counter, multiple_alleles_found, count_absent, percentage_lowCoverage, meanCoverage, consensus_sequence, number_diferences = utils.extractVariableFromPickle(file_path)

					for consensus_type in consensus_sequence:
						consensus_sequences_together[consensus_type][sequence_counter] = {'header': consensus_sequence[consensus_type]['header'], 'sequence': consensus_sequence[consensus_type]['sequence']}

					if write_consensus_first_time:
						for consensus_type in ['correct', 'noMatter', 'alignment']:
							file_to_remove = os.path.join(outdir, str(sample + '.' + consensus_type + '.fasta'))
							if os.path.isfile(file_to_remove):
								os.remove(file_to_remove)
						write_consensus_first_time = False
					consensus_files = write_consensus(outdir, sample, consensus_sequence)

					gene_identity = 0
					if sequences_information[sequence_counter]['length'] - 2 * length_extra_seq - count_absent > 0:
						gene_identity = 100 - (float(number_diferences) / (sequences_information[sequence_counter]['length'] - 2 * length_extra_seq - count_absent)) * 100

					sample_data[sequence_counter] = {'header': sequences_information[sequence_counter]['header'], 'gene_coverage': 100 - (float(count_absent) / (sequences_information[sequence_counter]['length'] - 2 * length_extra_seq)) * 100, 'gene_low_coverage': percentage_lowCoverage, 'gene_number_positions_multiple_alleles': multiple_alleles_found, 'gene_mean_read_coverage': meanCoverage, 'gene_identity': gene_identity}
					counter += 1

		if not debug_mode_true:
			utils.removeDirectory(gene_dir_path)

	if counter != len(sequences_information):
		run_successfully = False

	return run_successfully, sample_data, consensus_files, consensus_sequences_together


rematch_timer = functools.partial(utils.timer, name='ReMatCh module')


@rematch_timer
def runRematchModule(sample, fastq_files, reference_file, threads, outdir, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, minimum_gene_coverage, conserved_True, debug_mode_true, numMapLoc, minimum_gene_identity, rematch_run, softClip_baseQuality, softClip_recodeRun, reference_dict, softClip_cigarFlagRecode):
	rematch_folder = os.path.join(outdir, 'rematch_module', '')
	utils.removeDirectory(rematch_folder)
	os.mkdir(rematch_folder)

	# Map reads
	run_successfully, bam_file, reference_file = mapping_reads(fastq_files, reference_file, threads, rematch_folder, conserved_True, numMapLoc, rematch_run, softClip_baseQuality, softClip_recodeRun, reference_dict, softClip_cigarFlagRecode)

	if run_successfully:
		# Index reference file
		run_successfully, stdout = index_fasta_samtools(reference_file, None, None, True)
		if run_successfully:
			print 'Analysing alignment data'
			run_successfully, sample_data, consensus_files, consensus_sequences = sequence_data(sample, reference_file, bam_file, rematch_folder, threads, length_extra_seq, minimum_depth_presence, minimum_depth_call, minimum_depth_frequency_dominant_allele, debug_mode_true)

			if run_successfully:
				print 'Writing report file'
				number_absent_genes = 0
				number_genes_multiple_alleles = 0
				mean_sample_coverage = 0
				with open(os.path.join(outdir, 'rematchModule_report.txt'), 'wt') as writer:
					writer.write('\t'.join(['#gene', 'percentage_gene_coverage', 'gene_mean_read_coverage', 'percentage_gene_low_coverage', 'number_positions_multiple_alleles', 'percentage_gene_identity']) + '\n')
					for i in range(1, len(sample_data) + 1):
						writer.write('\t'.join([sample_data[i]['header'], str(round(sample_data[i]['gene_coverage'], 2)), str(round(sample_data[i]['gene_mean_read_coverage'], 2)), str(round(sample_data[i]['gene_low_coverage'], 2)), str(sample_data[i]['gene_number_positions_multiple_alleles']), str(round(sample_data[i]['gene_identity'], 2))]) + '\n')

						if sample_data[i]['gene_coverage'] < minimum_gene_coverage or sample_data[i]['gene_identity'] < minimum_gene_identity:
							number_absent_genes += 1
						else:
							mean_sample_coverage += sample_data[i]['gene_mean_read_coverage']
							if sample_data[i]['gene_number_positions_multiple_alleles'] > 0:
								number_genes_multiple_alleles += 1

					if len(sample_data) - number_absent_genes > 0:
						mean_sample_coverage = float(mean_sample_coverage) / float(len(sample_data) - number_absent_genes)
					else:
						mean_sample_coverage = 0

					writer.write('\n'.join(['#general', '>number_absent_genes', str(number_absent_genes), '>number_genes_multiple_alleles', str(number_genes_multiple_alleles), '>mean_sample_coverage', str(round(mean_sample_coverage, 2))]) + '\n')

					print '\n'.join([str('number_absent_genes: ' + str(number_absent_genes)), str('number_genes_multiple_alleles: ' + str(number_genes_multiple_alleles)), str('mean_sample_coverage: ' + str(round(mean_sample_coverage, 2)))])

	if not debug_mode_true:
		utils.removeDirectory(rematch_folder)

	return run_successfully, sample_data if 'sample_data' in locals() else None, {'number_absent_genes': number_absent_genes if 'number_absent_genes' in locals() else None, 'number_genes_multiple_alleles': number_genes_multiple_alleles if 'number_genes_multiple_alleles' in locals() else None, 'mean_sample_coverage': round(mean_sample_coverage, 2)} if 'mean_sample_coverage' in locals() else None, consensus_files if 'consensus_files' in locals() else None, consensus_sequences if 'consensus_sequences' in locals() else None
