def get_STRGid(gene_name,gtf):
    for line in open(gtf):
        if '\ttranscript\t' in line:
            if "ref_gene_name \"{}\"".format(gene_name) in line:
                info = line.strip().split('\t')[-1]
                STRG_id = info.split(';')[0].replace('gene_id \"', '').replace('\"', '')
    return STRG_id

def get_unique_STRGid(STRGid, gtf, max_count):
    unique = []
    for line in open(gtf):
        if '\ttranscript\t' in line and "gene_id \"{}\"".format(STRGid) in line:
            info = line.strip().split('\t')[-1]
            count = int(info.split(';')[-2].replace(' transcript_count \"', '').replace('\"', ''))
            if count <= max_count:
                unique.append(info.split(';')[1].replace(' transcript_id \"', '').replace('\"', ''))

    return unique

def get_gene_list(gene_list_file):
    gene_list = []
    for line in open(gene_list_file):
        gene_list.append(line.strip())
    return gene_list

def get_all_unique(gene_list, gtf, max_count, prefix):
    with open('{}_unique.txt'.format(prefix), 'w') as file:
        for gene in gene_list:
            STRGid = get_STRGid(gene, gtf)
            unique = get_unique_STRGid(STRGid, gtf, max_count)
            line = '{}:'.format(gene)
            for transcript in unique:
                line += '{},'.format(transcript)


