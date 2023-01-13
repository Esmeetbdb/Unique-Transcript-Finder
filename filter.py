def filter_counts(gtf_to_filter, max_count):
	filtered_file = open(gtf_to_filter+'filtered_counts','w')
	for line in open(gtf_to_filter):
		if "\ttranscript\t" in line:
			content = line.strip().split('\t')
			info = content[8].split(';')
			count = int(info[-2].replace(' transcript_count \"','').replace('\"'))
			if count <= max_count:
				filtered_file.writelines(line)
				#maybe add way to get exons in filtered file as well?
	filtered_file.close()


def filter_frequency(gtf_to_filter, max_frequency):
	filtered_file = open(gtf_to_filter+'filtered_counts','w')
	for line in open(gtf_to_filter):
		if "\ttranscript\t" in line:
			content = line.strip().split('\t')
			info = content[8].split(';')
			frequency = float(info[-2].replace(' transcript_frequency \"','').replace('\"'))
			if frequency <= max_frequency:
				filtered_file.writelines(line)
				#maybe add way to get exons in filtered file as well?
	filtered_file.close()


def filter_coverage(gtf_to_filter, min_cov):
	filtered_file = open(gtf_to_filter+'filtered_cov','w')
	for line in open(gtf_to_filter):
		if "\ttranscript\t" in line:
			content = line.strip().split('\t')
			info = content[8].split(';')
			coverage = float(info[-5].replace(' cov \"','').replace('\"'))
			if coverage >= min_cov:
				filtered_file.writelines(line)
				#maybe add way to get exons in filtered file as well?
	filtered_file.close()

def filter_annotated(gtf_to_filter):
	filtered_file = open(gtf_to_filter+'filtered_annotated','w')
	for line in open(gtf_to_filter):
		if "\ttranscript\t" in line:
			if '; ref_gene_name \"' in line:
				filtered_file.writelines(line)

	filtered_file.close()
