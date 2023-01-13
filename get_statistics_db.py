def get_avg_per_transcript(info_dict_per_transcript):
	average_list = [] #list with averages per transcript for plotting and t-test
	avg_dict = {} #dict with average per transcript (key), will contain standard deviation of info per transcript as well
	std_list = []
	for transcript in info_dict_per_transcript:
		inf = info_dict_per_transcript[transcript]
		mean = sum(inf)/len(inf)
		variance = sum([((x-mean)**2) for x in inf])/len(inf)
		std = variance ** 0.5
		average_list.append(mean)
		avg_dict[transcript] = (mean,std)
		std_list.append(std)
	return average_list, avg_dict, std_list

def print_stats_to_file(output_file, avg_dict):
	with open(output_file) as file:
		file.writelines('Transcript_ID, mean, standard_deviation\n')
		for transcript in avg_dict:
			file.writelines('{}, {}, {}\n'.format(transcript,avg_dict[transcript][0],avg_dict[transcript][1]))		

def t_test(average_list_1, average_list_2):
	from scipy import stats as st
	from scipy.stats import levene

	stats_lev, p_val_lev = levene(average_list_1, average_list_2)

	if p_val_lev >= 0.05:
		equal_var_assumed = True
	else:
		equal_var_assumed = False

	stat, p_val = st.ttest_ind(average_list_1, average_list_2, equal_val = equal_var_assumed)

	return stat, p_val
