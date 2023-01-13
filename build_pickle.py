def counts_from_db(database_name, out_path):
	import pickle
	database_dict = {}
	for line in open(out_path+database_name):
		if '\ttranscript\t' in line:
			content = line.strip().split('\t')
			info = content[8].split(';')
			# extract count
			# extract freq
			# extract transcript id
			database_dict[transcript id] = (count, freq)
	with open('{}.pickle'.format(out_path+database_name), 'wb') as file:
		pickle.dump(database_dict, file) 
