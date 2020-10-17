#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
#import statistics
#from random import randint
#import random
#from operator import itemgetter
import networkx as nx
import matplotlib.pyplot as plt
#random.seed(9001)



__author__ = "HOLLIER Laëtitia"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["HOLLIER Laëtitia"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "HOLLIER Laëtitia"
__email__ = "laetitia-hollier@outlook.fr"
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
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

##########################################################
########### 1. Création du graphe de de Bruijn ###########
##########################################################

def read_fastq(nom):
    """
    La fonction read_fastq prend en argument
    nom : un fichier fastq (str)
    renvoit les reads de ce fichier
    sous forme d'un argument de sequences
    """
    with open(nom) as filin:
        for line in enumerate(filin):
            yield next(filin)[:-1]
            next(filin)
            next(filin)

def cut_kmer(seq, k):
    """
    La fonction cut_mer prend 2 arguments:
    seq : une sequence sous forme d'une chaine de caractères (str)
    k : la taille du k-mer (integer)
    renvoit les k-mers uniques présents dans les reads
    """
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]


def build_kmer_dict(nom,k):
    """
    La fonction build_kmer_dict prend en argument
    nom : un fichier fastq (str)
    k : la taille du k-mer (integer)
    renvoit un dictionnaire comportant le k-mer (str) et
    la valeur du nombre d'occurence de ce k-mer (int)
    """
    dict_kmer = {}
    for i in read_fastq(nom):
        for j in cut_kmer(i,k):
            if j not in dict_kmer:
                dict_kmer[j] = 1
            else:
                dict_kmer[j] += 1
    return dict_kmer


def build_graph(dico):
    """
    La fonction build_graph prend en argument
    dico : dictionnaire possédant
    les k-mers et leur nombre d'occurence
    et renvoit un arbre orienté et pondéré
    """
    arbre_kmer = nx.DiGraph()
    for mot in dico:
        arbre_kmer.add_edge(mot[0:len(mot)-1],
            mot[1:],weight = dico[mot])
    print(nx.info(arbre_kmer))
    nx.draw(arbre_kmer)
    plt.show()
    plt.savefig("Figure1: 1er Graphe Orienté.png")
    return arbre_kmer



#================================================
# Main program
#================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    kmer = build_kmer_dict(sys.argv[2],int(sys.argv[4]))
    graphe = build_graph(kmer)