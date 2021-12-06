import Bio
# from Bio import Seq
from Bio import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
# from Bio import SearchIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import re
import collections
from collections import Counter
import os



class Hosts_Finder():
    def __init__(self, user_input):
        self.user_input = user_input
        #self.integrase = ['Integrase (Y-int)', 'Integrase', 'Tyrosine Integrase', 'serine integrase', 'Integrase (S-int)',
                  #'Serine Integrase', 'tyrosine integrase', 'Integrase, (S-int)']
        Entrez.email = "efgjakubo@gmail.com"
        E_VALUE_THRESH = 1e-12
        self.attp_sequence = ''
        self.tax_class = {}
        self.consensus_seq = ''

    def __repr__(self):
        return f"<Hosts_Finder {self.name}>"

    def get_user_input(self) -> str:
        id = "Z18946"
        return id

    def fetch_data_from_entrez(self, id):
        handle = Entrez.efetch(db='nucleotide', id=id, rettype="gbwithparts", retmode="text")
        return handle

    def write_file(self, filename, handle):
        if os.path.exists(filename):
            os.remove(filename)
        with open(filename, 'w') as out_handle:
            out_handle.write(handle.read())
        handle.close()
        result_handle = open(filename)
        return result_handle

    def map_location(self, result_handle):
        integrase_gene = 'integrase'
        for record in SeqIO.parse(result_handle, "genbank"):
            for f in record.features:
                if f.type == "CDS":
                    # print(f.qualifiers["product"][0])
                    product = f.qualifiers["product"][0]
                    if integrase_gene.upper() in product.upper():
                        my_start = f.location._start.position
                        my_end = f.location._end.position + 400
                        return [my_start, my_end]

    def blast(self, id, positions):
        result_handle = NCBIWWW.qblast("blastn", "nt", id, entrez_query='txid201174[ORGN]', query_from=positions[0], query_to=positions[1])
        return result_handle

    def save_file(self, filename, result_handle):
        if os.path.exists(filename):
            os.remove(filename)
        with open(filename+'.xml', 'w') as saved_xml_file:
            blast_results = result_handle.read()
            saved_xml_file.write(blast_results)
        return saved_xml_file

    def get_attp(self, saved_xml_file):
        test=0
        if os.path.exists("accession_out.txt"):
            os.remove("accession_out.txt")
        accession_out = open("accession_out.txt", 'w')
        if os.path.exists("query_out.fasta"):
            os.remove("query_out.fasta")
        query_out = open("query_out.fasta", 'w')

        for record in NCBIXML.parse(open(saved_xml_file)):
            for align in record.alignments:
                accession_out.write(align.accession + '\n')
                for hsp in align.hsps:
                    # if hsp.expect < E_VALUE_THRESH:
                    if hsp.identities < 75:
                        #                print("\n\n*Alignment*")
                        #                print("sequence:", align.title)
                        #                print("identities:", hsp.identities)
                        #                print("e value:", hsp.expect)
                        #                print(hsp.query[0:75] + "...")
                        #                print(hsp.match[0:75] + "...")
                        #                print(hsp.sbjct[0:75] + "...")

                        test += 1
                        query_out.write('>Test' + str(test) + '\n' + hsp.query[0:75] + '\n')

                attp_sequence = hsp.query[0:75]  # replace with code to capture true sequence

        query_out.close()
        sequences = [s for s in SeqIO.parse('query_out.fasta', 'fasta')]
        max_len = max([len(s.seq) for s in sequences])
        GAPS = '-'
        for seq in sequences:
            padding = GAPS * (max_len - len(seq.seq))
            seq.seq += padding

        SeqIO.write(sequences, 'query_out.fasta', 'fasta')

        accession_out.close()
        # print('count HSP = ', count)
        # print('count total = ', count_total, '\n')
        return attp_sequence

    def get_tax_id(self, filename):
        #output_file = open('tax_ids_from_blasthits.out', 'w')
        file_ids = open('accession_out.txt', 'r')
        ids_list = file_ids.readlines()

        ids_parsed = []
        taxonomy = []

        for acc_id in ids_list:
            ids_parsed.append(acc_id.replace("\n", ""))

        for acc_id in ids_parsed:
            # print('\n{}'.format(acc_id))
            seqio = Entrez.esummary(db="nucleotide", id=acc_id, retmode="xml")
            seqio_read = Entrez.read(seqio)
            seqio.close()

            tax_id = seqio_read[0]['TaxId']
            # print("Tax Id = ", tax_id)
            #output_file.write('\n\n{} \t TaxID{} \n'.format(acc_id, tax_id))
            taxinfo = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            taxinfo_read = Entrez.read(taxinfo)
            taxinfo.close()

            tax_list = taxinfo_read[0]['LineageEx']

            for tax_element in tax_list:
                if tax_element['Rank'] == 'family':
                    tax_class = tax_element['ScientificName']
                    taxonomy.append(tax_class)

            taxonomy_dict = Counter(taxonomy)
            sorted_taxonomy = dict(sorted(taxonomy_dict.items(), key=lambda x: x[1], reverse=True))

        return sorted_taxonomy

    def get_consensus(self, filename):
        alignment = AlignIO.read(open(filename), "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()
        return consensus

    def search(self):
        handle = self.fetch_data_from_entrez(self.user_input)
        positions = self.map_location(handle)
        result_handle = self.blast(self.user_input, positions)
        self.save_file('blast_results', result_handle)
        self.attp_sequence = self.get_attp("blast_results.xml")
        self.tax_class = self.get_tax_id("accession_out.txt")
        self.consensus_seq = self.get_consensus('query_out.fasta')

    def cleanup(self):
        if os.path.exists("blast_results.xml"):
            os.remove("blast_results.xml")
        if os.path.exists("accession_out.txt"):
            os.remove("accession_out.txt")
        if os.path.exists("query_out.fasta"):
            os.remove("query_out.fasta")




# user_input = "MK660712"
# print('starting search...')
# hf = Hosts_Finder(user_input)
# hf.search()
# #
# print(hf.attp_sequence)
# print(hf.tax_class)
# print(hf.consensus_seq)


# handle = hf.fetch_data_from_entrez(user_input)
# positions = hf.map_location(handle)
# # print(positions)
# # positions = [24253, 25769]
# result_handle = hf.blast(user_input, positions)
# blast_results_file = hf.save_file('blast_results', result_handle)
# hf.attp_sequence = hf.get_attp("blast_results.xml")
# # print(hf.attp_sequence)
# hf.tax_class = hf.get_tax_id("accession_out.txt")
# # print(hf.tax_class)
# hf.consensus_seq = hf.get_consensus('query_out_2.fasta')
# # print(hf.consensus_seq)

