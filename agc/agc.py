#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
import numpy as np
import nwalign3 as nw

# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Maxime TISON"
__copyright__ = "Universite de Paris"
__credits__ = ["Maxime TISON"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Maxime TISON"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as  monfich:
        seq = ""
        for line in monfich:
            if line.startswith(">"):
                if len(seq) > minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line[:-1]
        yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    list_read = [read for read in read_fasta(amplicon_file, minseqlen)]
    set_read = list(set(list_read))
    count_read = []
    for i in range(0, len(set_read), 1):
        if list_read.count(set_read[i]) >= mincount:
            count_read.append([set_read[i], list_read.count(set_read[i])])
    count_read.sort(key= lambda x: x[1], reverse=True)
    for i in range(0, len(count_read), 1):
        yield(count_read[i])


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            if id_seq in kmer_dict[kmer]:
                continue
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    kmer_list = list(cut_kmer(sequence, kmer_size))
    id_seqs = [id for kmer in kmer_list if kmer in kmer_dict for id in kmer_dict[kmer]]
    best_mate = Counter(id_seqs).most_common(2)
    return [best_id[0] for best_id in best_mate]


def detect_chimera(perc_identity_matrix):
    seq_similar = []
    list_std =[np.std(elem) for elem in perc_identity_matrix]
    mean_std = statistics.mean(list_std)
    for i in range(0, len(perc_identity_matrix), 1):
        if perc_identity_matrix[i][0] > perc_identity_matrix[i][1]:
            seq_similar.append(0)
        else:
            seq_similar.append(1)
    if mean_std > 5:
        if seq_similar.count(0) >= 1 and seq_similar.count(1) >= 1:
            return(True)
        else:
            return(False)
    else:
	    return(False)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    match = os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH"))
    seq_len = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    OTU_list = [seq_len[0]]
    
    for elem in seq_len:
        for otu in OTU_list:
            alignment_list = nw.global_align(otu[0], elem[0], gap_open=-1,
                                             gap_extend=-1, matrix= match)
            
            if get_identity(alignment_list) < 97:
                OTU_list.append(elem)
    return(OTU_list)


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, 'wt') as my_out_file:
        for i, sequence in enumerate(OTU_list):
            my_out_file.write(f">OTU_{i+1} occurrence:{sequence[1]}\n")
            my_out_file.write(fill(sequence[0])+"\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici


if __name__ == '__main__':
    main()
