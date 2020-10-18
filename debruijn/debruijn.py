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
import statistics
import random
import networkx as nx
import matplotlib.pyplot as plt
random.seed(9001)



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
    La fonction read_fastq prend en entrée
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
    La fonction cut_mer prend 2 entrée:
    seq : une sequence sous forme d'une chaine de caractères (str)
    k : la taille du k-mer (integer)
    renvoit les k-mers uniques présents dans les reads
    """
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]


def build_kmer_dict(nom,k):
    """
    La fonction build_kmer_dict prend en entrée
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
    La fonction build_graph prend en entrée
    dico : dictionnaire de k-mers (str) et leur nombre d'occurence (int)
    renvoit un arbre orienté et pondéré
    """
    arbre_kmer = nx.DiGraph()
    for mot in dico:
        arbre_kmer.add_edge(mot[0:len(mot)-1],
            mot[1:],weight = dico[mot])
    print(nx.info(arbre_kmer))
    options = {'node_color': "red", "node_size":500}
    nx.draw(arbre_kmer,with_labels=True, font_weight='bold', **options)
    plt.show()
    plt.savefig("Figure1: 1er Graphe Orienté de Bruijn.png")
    return arbre_kmer


##########################################################
########### 2. Parcours du graphe de de Bruijn ###########
##########################################################


def get_starting_nodes(arbre):
    """
    La fonction get_starting_nodes prend en entrée
    arbre: object networkx DiGraph()
    renvoit une liste de noeuds d'entrée (str ou int)
    """
    nodes_in = []
    for node in arbre.nodes:
        if len(list(arbre.predecessors(node))) == 0:
            nodes_in.append(node)
    return nodes_in


def get_sink_nodes(arbre):
    """
    La fonction get_sink_nodes prend en entrée
    arbre: object networkx DiGraph()
    renvoit une liste de noeuds de sortie (str ou int)
    """
    nodes_out = []
    for node in arbre.nodes:
        if len(list(arbre.successors(node))) == 0:
            nodes_out.append(node)
    return nodes_out


def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def get_contigs(arbre, entree, sortie):
    """
    La fonction get_sink_nodes prend en entrée
    arbre: object networkx DiGraph()
    entree: une liste de noeurds d'entree (str)
    sortie: une liste de noeurds de sortie (str)
    renvoit une liste de tuple avec le contig (str) et sa taille (int)
    """
    contigs = []
    for debut in entree:
        for fin in sortie:
            if list(nx.all_simple_paths(arbre, debut, sortie)) != []:
                path = nx.shortest_path(arbre,debut,fin)
                seq =""
                for i in range(len(path)-1):
                    seq += path[i][-1]
                contigs.append((seq, len(seq)))
    return contigs



def save_contigs(contig,fichier):
    """
    La fonction save_contigs prend en entrée
    contig: une liste de tuple
    fichier: nom de fichier de sortie
    renvoit un fichier fasta
    """
    with open(fichier, "w") as filout:
        for i,seq in enumerate(contig):
            filout.write(">contig_{} len={}\n".format(i+1,seq[1]))
            filout.write(fill(seq[0] + "\n"))
            filout.write("\n")
        return filout


##########################################################
###### 3. Simplification du graphe de de Bruijn ##########
##########################################################

def std(liste):
    """
    calcule de lecart type de la liste
    """
    return statistics.stdev(liste)


def path_average_weight(arbre,chemin):
    """
    prend un graphe (arbre): object networkx DiGraph()
    un chemin: liste de chaine de caractères (str)
    renvoit le poids moyen du chemin (float)
    """
    liste_weights = []
    for i,j,poids in arbre.subgraph(chemin).edges(data=True):
        liste_weights.append(poids["weight"])
    return statistics.mean(liste_weights)



def remove_paths(arbre,liste_chemin, delete_entry_node, delete_sink_node):
    """
    prend un graphe (arbre): object networkx DiGraph()
    liste_chemin: liste de chaine de caractères (str)
    delete_entry_node: variable booléenne (True/False)
    delete_sink_node: variable booléenne (True/False)
    retourne l'arbre nettoyé des chemins indésirables
    """
    for i in range(len(liste_chemin)):
        if delete_entry_node == True:
            arbre.remove_node(liste_chemin[i][0])
        if delete_sink_node == True:
            arbre.remove_node(liste_chemin[i][-1])
    return arbre


def select_best_path(arbre, liste_chemin, liste_taille,liste_weights,
    delete_entry_node = False, delete_sink_node =False):
    """
    prend un graphe (arbre): object networkx DiGraph()
    liste_chemin: liste de chaine de caractères (str)
    liste_taille: liste de la taille de chaque chemin (int)
    liste_weights: liste de poids moyen de chaque chemin (float)
    delete_entry_node: variable booléenne (False)
    delete_sink_node: variable booléenne (False)
    retourne l'arbre nettoyé des chemins indésirables

    Meilleur chemin: celui qui a le meilleur poids
    si 2 chemins avec même poids => on prend celui le plus long
    si même taille => on prend random
    """
    # initialisation des paramètres
    path_max_weight = 0
    path_max_length = 0
    path_max_indice = -1
    for i in range(len(liste_chemin)):
        # on prend la path avec le meilleur poids
        if liste_weights[i] > path_max_weight[i]:
            path_max_weight = liste_weights[i]
            path_max_length = liste_taille[i]
            path_max_indice = i
        # si on a même poids, recherche de la plus longue séquence
        elif liste_weights[i] == path_max_weight[i]:
            if liste_taille[i] > path_max_length[i]:
                path_max_weight = liste_weights[i]
                path_max_length = liste_taille[i]
                path_max_indice = i
            # si même taille de sequence, random entre les 2 seqs
            elif liste_taille[i] == path_max_length[i]:
                choix_path = random.randrange(0,1)
                if choix_path:
                    path_max_weight = liste_weights[i]
                    path_max_length = liste_taille[i]
                    path_max_indice = i

    arbre = remove_paths(arbre, liste_chemin[:path_max_indice]+
        liste_chemin[path_max_indice+1:], delete_entry_node,
        delete_sink_node)
    return arbre



def solve_bubble(arbre, ancetre, descendant):
    """
    prend un graphe (arbre): object networkx DiGraph()
    ancetre: noeud ancêtre, chaine de caractères (str)
    descendant: noeud descendant, chaine de caractères (str)
    renvoit un graphe (arbre) nettoyé de la bulle
    """
    solve_path = []
    solve_length = []
    solve_weight = []
    paths =nx.all_simple_paths(arbre, ancetre, descendant)
    for path in paths:
        solve_path.append(path)
        solve_length.append(len(path))
        solve_weight.append(path_average_weight(arbre,path))
    arbre = select_best_path(arbre, solve_path, solve_length,solve_weight)
    return arbre



#================================================
#================ Main program ==================
#================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # lecture du fichier
    file = sys.argv[2]

    # construction du graphe grace au dictionnaire kmer
    kmer = build_kmer_dict(file,int(sys.argv[4]))
    graphe = build_graph(kmer)

    # ecriture du/des contigs
    entry = get_starting_nodes(graphe)
    ending = get_sink_nodes(graphe)
    sequences = get_contigs(graphe,entry,ending)
    save_contigs(sequences, "path_contig.txt")


if __name__ == '__main__':
    main()