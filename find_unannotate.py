def find_unannotated(database, threshold):
	file = open('unannotated_common_transcripts.txt', 'w')
	for line in open(database):
		if "\ttranscript\t" not in line:
			continue
		if "gene_name" in line:
			continue
		content = line.strip().split('\t')
		info = content[8].split(';')
		count = int(info[-3].replace(' transcript_count \"', '').replace('transcript_count ','').replace('\"', ''))
		if count >= threshold:
			file.writelines('name: {}, count: {}, chr: {}, start: {}, stop: {}\n'.format(info[0], str(count), content[0],content[3],content[4]))
	file.close()
find_unannotated('/proj/sens2017106/esmee_txdb/other_scripts/transcript_database.gtf', 5)
		
