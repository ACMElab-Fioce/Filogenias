from Bio import Phylo
import pandas as pd
from dendropy import Tree
from collections import Counter

# Carrega a árvore a partir de um arquivo em formato NEXUS
tree = Tree.get_from_path("reroot_input_nexus.nexus.color.tree", "nexus")

dic_clado_taxons = {}
# Adiciona um identificador único em cada nó interno
clade_id = 0
for node in tree.postorder_node_iter():
    if not node.is_leaf():
        temp_list = []
        node.annotations.add_new(name= "clado",value=f"clado_{clade_id}")

        for taxon in node.leaf_nodes():
            nome_folha = str(taxon.taxon)
            nome_folha = nome_folha.strip("'")
            temp_list.append(nome_folha)
        dic_clado_taxons[f"clado_{clade_id}"] = temp_list

        clade_id += 1

# Escreve a árvore com os identificadores únicos em um arquivo em formato NEXUS
tree.write(path="arvore_identificadores.nexus", schema="nexus")



