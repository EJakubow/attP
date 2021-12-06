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



class Hosts_Finder():
    def __init__(self, user_input):
        self.user_input = user_input
        self.integrase = ['Integrase (Y-int)', 'Integrase', 'Tyrosine Integrase', 'serine integrase', 'Integrase (S-int)',
                  'Serine Integrase', 'tyrosine integrase']
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
        with open(filename, 'w') as out_handle:
            out_handle.write(handle.read())
        handle.close()
        result_handle = open(filename)
        return result_handle

    def map_location(self, result_handle):
        for record in SeqIO.parse(result_handle, "genbank"):
            for f in record.features:
                if f.type == "CDS" and "gene" in f.qualifiers:
                    gene = f.qualifiers['gene'][0]
                    for i in self.integrase:
                        if f.qualifiers["product"][0] == i:
                            my_start = f.location._start.position
                            my_end = f.location._end.position + 400
                            return [my_start, my_end]

    def blast(self, id, positions):
        result_handle = NCBIWWW.qblast("blastn", "nt", id, entrez_query='txid201174[ORGN]', query_from=positions[0], query_to=positions[1])
        return result_handle

    # def save_file(self, filename, result_handle):
    #     with open(filename+'.xml', 'w') as saved_xml_file:
    #         blast_results = result_handle.read()
    #         saved_xml_file.write(blast_results)
    #     return saved_xml_file

    def get_attp(self, saved_xml_file):  #change variable to filename?
        test=0
        accession_out = open("templates/accession_out.txt", 'w')
        query_out = open("templates/query_out.fasta", 'w')

        for record in NCBIXML.parse(open(saved_xml_file)):
            for align in record.alignments:
                # print("Accession number:", align.accession)
                accession_out.write(align.accession + '\n')
                # print("match: %s " % align.title[:100])
                for hsp in align.hsps:
                    # if hsp.expect < E_VALUE_THRESH:
                    if hsp.identities < 75:
                        # print("\n\n*Alignment*")
                        # print("sequence:", align.title)
                        # print("identities:", hsp.identities)
                        # print("e value:", hsp.expect)
                        # print(hsp.query[0:75] + "...")
                        # print(hsp.match[0:75] + "...")
                        # print(hsp.sbjct[0:75] + "...")
                        test += 1
                        query_out.write('>Test' + str(test) + '\n' + hsp.query[0:75] + '\n')

                    attp_sequence = hsp.query[0:75]  # replace with code to capture true sequence
        accession_out.close()
        query_out.close()
        # print('count HSP = ', count)
        # print('count total = ', count_total, '\n')
        return attp_sequence

    def get_tax_id(self, filename):
        #output_file = open('tax_ids_from_blasthits.out', 'w')
        file_ids = open('templates/accession_out.txt', 'r')
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

        return(Counter(taxonomy))

    def get_consensus(self, filename):
        alignment = AlignIO.read(open(filename), "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus()
        return consensus

    def search(self):
        handle = hf.fetch_data_from_entrez(self.user_input)
        positions = hf.map_location(handle)
        result_handle = hf.blast(user_input, positions)
        hf.save_file('blast_results', result_handle)
        hf.attp_sequence = hf.get_attp("blast_results.xml")
        hf.tax_class = hf.get_tax_id("accession_out.txt")
        hf.consensus_seq = hf.get_consensus('query_out_2.fasta')



user_input = "Z18946"
hf = Hosts_Finder(user_input)
hf.search()

print(hf.attp_sequence)
print(hf.tax_class)
print(hf.consensus_seq)


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

