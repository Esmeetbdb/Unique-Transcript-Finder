def make_dict(out_path, gtf_to_filter):
	file_dict = {}
	
	for line in open(out_path+gtf_to_filter):
		if line[0] == '#':
			continue
		if '\ttranscript\t' in line:
			k = line
			file_dict[k] = []
		else:
			file_dict[k].append(line)
	return file_dict



def filter_dict(file_dict, max_count, max_frequency, min_cov, min_tpm, min_fpkm, annotated, out_path, gtf_to_filter):
	filtered_file = open(out_path + gtf_to_filter + '_filtered', 'w')
	for line in file_dict:
		content = line.strip().split('\t')
		info = content[8].split(';')
		print(info)
		count = int(info[-2].replace(' transcript_count \"','').replace('\"',''))
		frequency = float(info[-1].replace(' transcript_frequency \"','').replace('\"',''))
		cov = float(info[-6].replace(' cov \"','').replace('\"', ''))
		tpm = float(info[-4].replace(' TPM \"', '').replace('\"', ''))
		fpkm = float(info[-5].replace(' FPKM \"', '').replace('\"', ''))
		if '; ref_gene_name \"' in line:
			in_ref_gen = True
		else:
			in_ref_gen = False

		if count <= max_count and frequency <= max_frequency and cov >= min_cov and tpm >= min_tpm and fpkm >= min_fpkm:
			if annotated == True:
				if in_ref_gen == True:
					filtered_file.writelines(line)
					for exon in file_dict[line]:
						filtered_file.writelines(exon)
			else:
				filtered_file.writelines(line)
				for exon in file_dict[line]:
					filtered_file.writelines(exon)
		
	

