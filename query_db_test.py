import pickle
def query_ind(database, gtf, gffcompare_container, output_prefix, out_path):
	import os
	cmd = 'singularity exec {} gffcompare -r {} -R -o {} {}'.format(gffcompare_container,out_path+database,out_path+output_prefix, gtf)
	os.system(cmd)

def tmap_to_dict(gtf, outprefix):
	patient_gtf = gtf.split('/')[-1]
	tmap = gtf.replace(patient_gtf,'{}.{}.tmap'.format(outprefix,patient_gtf))
	tmap_dict = {}
	for line in open(tmap):
		content = line.strip().split('\t')		
		transcript_id = content[4]
		database_id = content[1]
		tmap_dict[transcript_id] = database_id
	return tmap_dict

def unpickle_db(database_pickle):
	unpickle_file = open(database_pickle, 'rb')
	database_dict = pickle.load(unpickle_file)	
	return database_dict

def annotate_counts(database_dict, tmap_dict, gtf, outprefix, out_path):
	import os
	import subprocess
	patient_gtf = gtf.split('/')[-1]
	annotated_gtf = open(out_path+patient_gtf+'counts_test','w')
	for line in open(gtf):
		if "\ttranscript\t" in line:
			content = line.strip().split('\t')
			t_id = content[8].split(';')[1].replace(' transcript_id \"','').replace('\"','')
			#t_id = content[8].split(';')[1].replace(' transcript_id \"','').replace('\"','').replace('.','\.')
			class_code = content[2]

			if class_code != 'u':
				db_id = tmap_dict[t_id]
				count = database_dict[db_id][0]
				freq = database_dict[db_id][1]
				new_line = line.strip() + '; transcript_count \"' + str(count) + '\" ' +'; transcript_frequency \"' + str(freq) + '\"' + '\n'
				annotated_gtf.writelines(new_line)

			else:
				annotated_gtf.writelines(line.strip + ';transcript_count \"0\" ; transcript_frequency \"0.0\"')
		else:
			annotated_gtf.writelines(line)
	annotated_gtf.close()						


