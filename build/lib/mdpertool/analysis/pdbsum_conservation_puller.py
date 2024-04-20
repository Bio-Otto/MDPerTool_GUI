from urllib import request
import pandas as pd


# url_root = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetText.pl?pdb="
# pdb_id = '1pvs'
# pdb_chain = 'A'
#
# html = url_root + pdb_id + '&chain=' + pdb_chain
#
# file = request.urlopen(url_root + pdb_id + '&chain=' + pdb_chain)
#
# for cnt, line in enumerate(file):
#     if cnt > 17:
#         decoded_line = line.decode("utf-8")
#         conserv_score = float(decoded_line[53:57].replace(',', '.'))
#         conserv_res = decoded_line[12:18]
#         if conserv_score == 9.9:
#             print(conserv_res, conserv_score)
from Bio.PDB import PDBParser
from Bio.SeqIO.PdbIO import PdbSeqresIterator


def get_residue_list(pdb_file_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file_path)
    structure_residue_list = []

    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                structure_residue_list.append(str(residue.resname)+str(residue.id[1]))
    return structure_residue_list


def get_conservation_scores(pdb_id, chain_id, cutoff, bound_pdb):
    actual_res_IDs = []
    actual_scores = []
    url_root = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetText.pl?pdb="

    html_address = url_root + pdb_id.lower() + '&chain=' + chain_id
    file = request.urlopen(html_address)

    for count, line in enumerate(file):
        if count > 17:
            decoded_line = line.decode("utf-8")

            try:
                conserv_score = float(decoded_line[53:57].replace(',', '.'))
                conserv_res = decoded_line[12:18].strip()

                if conserv_score >= cutoff and conserv_res != '':
                    actual_res_IDs.append(conserv_res)
                    actual_scores.append(conserv_score)
            except:
                pass

    actual_res_IDs = list(filter(None, actual_res_IDs))

    if len(actual_scores) == len(actual_res_IDs):
        print("Conservation Scores pulled succesfully")
        try:
            filtered_conserv_scores = []
            structure_res_ls = get_residue_list(pdb_file_path=bound_pdb)
            filtered_resIDs = sorted(set(structure_res_ls).intersection(actual_res_IDs), key=lambda x: structure_res_ls.index(x))

            for cnt, res_i in enumerate(filtered_resIDs):
                idx = actual_res_IDs.index(res_i)
                filtered_conserv_scores.append(actual_scores[idx])

            return filtered_resIDs, filtered_conserv_scores
        except Exception as err:
            print("There are problems in getting conservation scores\n", err)
            return actual_res_IDs, actual_scores

    else:
        print("Pulling Data from conservation DB has failed")
