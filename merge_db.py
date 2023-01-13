def merge_db(input, output, min_len, min_cov, min_fpkm, min_tpm, min_iso):
	import os
	cmd = 'stringtie --merge -o {} -m {} -c {} -F {} -T {} -f {} -g 50 -i {}'.format(output,min_len, min_cov, min_fpkm, min_tpm, min_iso, input)
	os.system(cmd)



