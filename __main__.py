import argparse
import sys
import os
import pickle

def make_db(args):
#	import merge_db as mdb
#	mdb.merge_db(args.input, args.database_name, args.min_len, args.min_cov, args.min_fpkm, args.min_tpm, args.min_iso)
	
	import count_transcripts as ct
	ct.run_ggfcomp(args.gffcompare_container,args.reference_transcriptome,args.input,args.out_prefix, args.target)
	transcript_counts,transcript_info = ct.count_transcripts(args.out_prefix, args.target)
	ct.counts_into_db(args.out_prefix, transcript_counts, args.database_name, args.target)

def query_db(args):
	if args.create_database_pickle == True:
		import build_pickle.py as bp
		bp.couts_from_db(args.database, args.target_dir)

	import query_db_test as qdb
	qdb.query_ind(args.database, args.patient_gtf, args.gffcompare_container, args.out_prefix, args.target_dir)
	tmap_dict = qdb.tmap_to_dict(args.patient_gtf, args.out_prefix)
	db_pickle = args.target_dir + args.database + '.pickle'
	db_dict = qdb.unpickle_db(db_pickle)

	qdb.annotate_counts(db_dict, tmap_dict, args.patient_gtf, args.out_prefix, args.target_dir)
	
def stats_db(args):
	import count_transcripts as ct
	import get_statistics as gs


	transcript_counts,transcript_info = ct.count_transcripts(args.out_prefix, args.target_dir)
	tpm_dict_rare, fpkm_dict_rare, len_dict_rare, cov_dict_rare, exon_dict_rare = ct.extract_info(transcript_info, args.rare_threshold, False)

	rare_tpm_avg_list, rare_tpm_avg_dict, rare_tpm_std_list = gs.get_avg_per_transcript(tpm_dict_rare)
	gs.print_stats_to_file(args.target_dir+'rare_tpm_info.csv',rare_tpm_avg_dict)

	rare_fpkm_avg_list, rare_fpkm_avg_dict, rare_fpkm_std_list = gs.get_avg_per_transcript(fpkm_dict_rare)
	gs.print_stats_to_file(args.target_dir+'rare_fpkm_info.csv',rare_fpkm_avg_dict)


	rare_len_avg_list, rare_len_avg_dict, rare_len_std_list = gs.get_avg_per_transcript(len_dict_rare)
	gs.print_stats_to_file(args.target_dir+'rare_len_info.csv',rare_len_avg_dict)


	rare_cov_avg_list, rare_cov_avg_dict, rare_cov_std_list = gs.get_avg_per_transcript(cov_dict_rare)
	gs.print_stats_to_file(args.target_dir+'rare_cov_info.csv',rare_cov_avg_dict)


	rare_exon_avg_list, rare_exon_avg_dict, rare_exon_std_list = gs.get_avg_per_transcript(exon_dict_rare)
	gs.print_stats_to_file(args.target_dir+'rare_exon_info.csv',rare_exon_avg_dict)

	
	tpm_dict_common, fpkm_dict_common, len_dict_common, cov_dict_common, exon_dict_common = ct.extract_info(transcript_info, args.common_threshold, True)
	
	common_tpm_avg_list, common_tpm_avg_dict, common_tpm_std_list = gs.get_avg_per_transcript(tpm_dict_common)
	gs.print_stats_to_file(args.target_dir+'common_tpm_info.csv',common_tpm_avg_dict)

	common_fpkm_avg_list, common_fpkm_avg_dict, common_fpkm_std_list = gs.get_avg_per_transcript(fpkm_dict_common)
	gs.print_stats_to_file(args.target_dir+'common_fpkm_info.csv',common_fpkm_avg_dict)

	common_len_avg_list, common_len_avg_dict, common_len_std_list = gs.get_avg_per_transcript(len_dict_common)
	gs.print_stats_to_file(args.target_dir+'common_len_info.csv',common_len_avg_dict)
	
	common_cov_avg_list, common_cov_avg_dict, common_cov_std_list = gs.get_avg_per_transcript(cov_dict_common)
	gs.print_stats_to_file(args.target_dir+'common_cov_info.csv',common_cov_avg_dict)

	common_exon_avg_list, common_exon_avg_dict, common_exon_std_list = gs.get_avg_per_transcript(exon_dict_common)
	gs.print_stats_to_file(args.target_dir+'common_exon_info.csv',common_exon_avg_dict)


	#t-tests
	stat_tpm, p_tpm, mean_rare_tpm, mean_common_tpm = gs.t_test(rare_tpm_avg_list, common_tpm_avg_list)
	stats_fpkm, p_fpkm,  mean_rare_fpkm, mean_common_fpkm = gs.t_test(rare_fpkm_avg_list, common_fpkm_avg_list)
	stats_len, p_len,  mean_rare_len, mean_common_len = gs.t_test(rare_len_avg_list, common_len_avg_list)
	stats_cov, p_cov,  mean_rare_cov, mean_common_cov = gs.t_test(rare_cov_avg_list, common_cov_avg_list)
	stats_exon, p_exon, mean_rare_exon, mean_common_exon = gs.t_test(rare_exon_avg_list, common_exon_avg_list)

	with open(args.out_prefix+'_statistics.csv', 'w') as file:
		file.writelines('Compared, Statistic, P-value, MeanRare, MeanCommon\n')
		file.writelines('TPM, {}, {}, {}, {}\n'.format(stat_tpm, p_tpm, mean_rare_tpm, mean_common_tpm))
		file.writelines('FPKM, {}, {}, {}, {}\n'.format(stats_fpkm, p_fpkm, mean_rare_fpkm, mean_common_fpkm))
		file.writelines('Length, {}, {}, {}, {}\n'.format(stats_len, p_len, mean_rare_len, mean_common_len))
		file.writelines('Coverage, {}, {}, {}, {}\n'.format(stats_cov, p_cov, mean_rare_cov, mean_common_cov))
		file.writelines('Exon_number, {}, {}, {}, {}\n'.format(stats_exon, p_exon, mean_rare_exon, mean_common_exon))


def filter_query(args):
	import filter_new as f
	query_dict = f.make_dict(args.output_folder, args.query_file)
	f.filter_dict(query_dict, args.max_count, args.max_frequency, args.min_cov, args.min_tpm, args.min_fpkm, args.annotated, args.output_folder, args.query_file)

def main():
	
	parser = argparse.ArgumentParser(description = "build, analyse and query a RNAseq transcript database")
	
	subparsers = parser.add_subparsers()
	
	parser_build = subparsers.add_parser("build", help = "build help")
	parser_build.add_argument("input", type=str, help = "path to input file that contains paths to all gtf files to be included in the database")
	parser_build.add_argument("gffcompare_container", type=str, help = "path to singularity container with gff compare")
	parser_build.add_argument("reference_transcriptome", type=str, help = "path to reference transcriptome/genome")
	parser_build.add_argument("--database_name","-n", type =str, default="transcript_database.gtf", help="name of the database. Database will be created in the current folder unless a path is specified. The final database with counts will have a _counts suffix.")
	parser_build.add_argument("--min_len","-m", type=int, default=50, help="minimum length of transcripts to be included in the database.")
	parser_build.add_argument("--min_cov","-c", type=int, default=0, help="minimum coverage of transcripts to be included in the database.")
	parser_build.add_argument("--min_fpkm","-F", type=int, default=0, help="minimum fpkm of transcripts to be included in the database.")
	parser_build.add_argument("--min_tpm","-T", type = int, default=0, help="minimum TPM of transcripts to be included in the database.")
	parser_build.add_argument("--min_iso","-f", type=float, default=0.0, help="minimum isoform frequency of transcripts to be included in the database")
	parser_build.add_argument("--out_prefix","-O",type=str, default = "transcript_database", help="out prefix for files generated by gff compare when counting transcript occurence in database")
	parser_build.add_argument("--target",type=str, default = "", help="Path where output files should be stored. Must end with \. If left empty files will be created in the current directory.")
	parser_build.set_defaults(func=make_db)

	parser_query = subparsers.add_parser("query", help = "build query")
	parser_query.add_argument("patient_gtf", type=str, help = "path to the GTF file with patient transcript information")
	parser_query.add_argument("gffcompare_container", type=str, help = "path to singularity container with gff compare")
	parser_query.add_argument("--database", type=str, default="transcript_database.gtf", help = "path to the database created in the build step or other database with transcript information in the correct format")
	parser_query.add_argument("--database_pickle", type=str, default = "transcript_database.gtf.pickle", help = "path to the pickle file containing count information obtained from the database")
	parser_query.add_argument("--create_database_pickle", type=bool, default = False, help = "If the pickle file with count information is lost this option can be used to make a new one prior to querying.")
	parser_build.add_argument("--target_dir",type=str, default = "", help="Path where output files should be stored. Must end with \. If left empty files will be created in the current directory.")
	parser_query.add_argument("--out_prefix", type=str, default = "query", help = "path to the GTF file with patient transcript information")
	parser_query.set_defaults(func=query_db)
	
	parser_statistics = subparsers.add_parser("statistics", help = "help statistics")
	parser_statistics.add_argument("rare_threshold",type=int, help = "The max number of times a transcript can be found in the database for it to be considered a rare transcript")
	parser_statistics.add_argument("common_threshold",type=int, help = "The min number of times a transcript must be found in the database for it to be considered a common transcript")
	parser_statistics.add_argument("--out_prefix","-O",type=str, default = "transcript_database", help="Outprefix used in the database build step.")
	parser_statistics.add_argument("--target_dir",type=str, default = "", help="Path to the folder that contains all files created during database building. If left empty the current folder is taken as the correct one")
	parser_statistics.set_defaults(func=stats_db)

	parser_filter = subparsers.add_parser("filter", help = "help filter")
	parser_filter.add_argument("query_file", type=str, help = "the gtf file that should be filtered. Created in the query step.")
	parser_filter.add_argument("--max_count", type = int, default = 1, help = "The max number of times a transcript can be found in the database to pass the filter.")
	parser_filter.add_argument("--max_frequency", type=float, default = 1.0, help = "The max frequecy at which a transcript can be included in the database to pass the filter. If this is changed it is adviced to set max_count arbitrarily high so that the count does not influence filtering.")
	parser_filter.add_argument("--min_cov",type=float, default = 0.0, help = "The min coverage a transcript should have to pass the filter")
	parser_filter.add_argument("--min_tpm", type = float, default = 0.0, help = "The min tpm a transcript should have to pass the filter")
	parser_filter.add_argument("--min_fpkm", type = float, default = 0.0, help = "The min fpkm a transcript should have to pass the filter")
	parser_filter.add_argument("--annotated", type = bool, default = True, help = "If True, only transcript that are annotated in the reference genome pass the filter,")
	parser_filter.add_argument("--output_folder", type=str, default = "", help="Path where output files should be stored. Must end with \. If left empty files will be created in the current directory.")
	parser_filter.set_defaults(func=filter_query)

	args = parser.parse_args(sys.argv[1:])
	args.func(args)
	
if __name__ == '__main__':
	main()
