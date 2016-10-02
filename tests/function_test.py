from Bio import SeqIO

def get_genbank_record(input_genome_path):
    with open(input_genome_path) as input_genome_handle:
        records = list(SeqIO.parse(input_genome_handle, "genbank"))
    print(len(records))
    print(type(records))
    return(records)

if __name__ == "__main__":
    recs = get_genbank_record("/home/nicholas/GitHub/pyutilsnrw/tests/references/n_equitans.gbk")
    print(len(recs))
    print(type(recs))
