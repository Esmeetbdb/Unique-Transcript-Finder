def run_ggfcomp(container_loc, reference_transcriptome, input_file, out_prefix, out_path):
	import os
	cmd = 'singularity exec {} gffcompare -i {} -r {} -o {}'.format(container_loc, input_file,reference_transcriptome,out_path+out_prefix)
	os.system(cmd)

def count_transcripts(out_prefix, out_path):
	transcript_counts = {}
	transcript_info = {}
	file = open("non_included_transcripts.txt","w")
	for line in open(out_path+out_prefix+'.tracking'):
		count = 0
		content = line.strip().split('\t')
		total_count = len(content)-4
		#print(line)
#		if content[3] == "u":
#			file.writelines(line+'\n')
#			continue
		transcript_id = content[0]
		transcript_info[transcript_id] = []
		for item in content[4:]:
			if item != '-':
				count+=1
				transcript_info[transcript_id].append(item.split(',')[0])
		
		transcript_counts[transcript_id] = (count, count/total_count)
	file.close()
	return transcript_counts,transcript_info
		

def extract_info(transcript_info, threshold, above_thresh):
	tpm_dict = {}
	fpkm_dict = {}
	len_dict = {}
	cov_dict = {}
	exon_dict = {}
	for transcript_id in transcript_info:
		if above_thresh == True: #to get transcript info of transcripts present in more than n individuals, where n is the threshold
			if len(transcript_info[transcript_id]) >= threshold:
				tpm_dict[transcript_id] = []
				fpkm_dict[transcript_id] = []
				len_dict[transcript_id] = []
				cov_dict[transcript_id] = []
				exon_dict[transcript_id] = []

				for ind in transcript_info[transcript_id]:
					content = ind.split('|')
					tpm_dict[transcript_id].append(float(content[4]))
					fpkm_dict[transcript_id].append(float(content[3]))
					len_dict[transcript_id].append(int(content[6]))
					cov_dict[transcript_id].append(float(content[5]))
					exon_dict[transcript_id].append(int(content[2]))
		else: #to get transcript info of transcripts present in less than n individuals, where n is the threshold
			if len(transcript_info[transcript_id]) <= threshold:
				tpm_dict[transcript_id] = []
				fpkm_dict[transcript_id] = []
				len_dict[transcript_id] = []
				cov_dict[transcript_id] = []
				exon_dict[transcript_id] = []			

				for ind in transcript_info[transcript_id]:
					content = ind.split('|')
					#print(content)
					tpm_dict[transcript_id].append(float(content[4]))
					fpkm_dict[transcript_id].append(float(content[3]))
					#print(content[6])
					len_dict[transcript_id].append(int(content[6]))
					cov_dict[transcript_id].append(float(content[5]))
					exon_dict[transcript_id].append(int(content[2]))
					
	return tpm_dict, fpkm_dict, len_dict, cov_dict, exon_dict

def counts_into_db(outprefix, transcript_counts, database_name, out_path):
	import pickle
	new_db = open(out_path+database_name, 'w')
	database = out_path+outprefix+'.combined.gtf'
	database_dict = {}
	for line in open(database):
		if '\ttranscript\t' in line:
			content = line.strip().split('\t')
			info = content[8].split(';')
			transcript_id = info[0].replace('transcript_id \"','').replace('\"','')
			count = transcript_counts[transcript_id][0]
			freq = transcript_counts[transcript_id][1]
			new_line = line.strip() + 'transcript_count \"{}\"; transcript_frequency \"{}\";\n'.format(count,freq)
			new_db.writelines(new_line)
			database_dict[transcript_id] = (count,freq)
		else:
			new_db.writelines(line)
	with open('{}.pickle'.format(out_path+database_name),'wb') as file:
		pickle.dump(database_dict, file)

	new_db.close()



