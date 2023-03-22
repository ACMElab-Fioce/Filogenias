import pandas as pd
import ast
from dendropy import Tree
from collections import Counter


def diferenca(A,B):
    A = set(A)
    B = set(B)
    # return(A.difference(B), B.difference(A))
    return(list(A.difference(B)))


def get_synapomorphy(taxons):

    d = [dic_amostras_mut[taxon] for taxon in taxons]

    # Calcula a frequência das mutações
    todas_as_mutacoes = []
    num_amostras = len(d)
    for lista in d:
        todas_as_mutacoes += lista

    dic_mut_freq = {}
    dic_mut = Counter(todas_as_mutacoes)

    for mut, contagem in dic_mut.items():
        dic_mut_freq[mut] = contagem/num_amostras
    # print(dic_mut_freq)
    mutacoes_clado = [mut for mut, freq in dic_mut_freq.items() if freq >= 0.7]
    # print(mutacoes_clado)
    return (mutacoes_clado, dic_mut_freq)
    # return(clado, taxons, list(sinapomorfia), dic_mut_freq)

# Carrega a árvore a partir de um arquivo em formato NEXUS
tree = Tree.get_from_path("arvore_identificadores.nexus", "nexus")

df = pd.read_csv("nextclade.tsv", sep='\t')

amostras = df["seqName"].tolist()
mutacoes = df["substitutions"].tolist()

dic_amostras_mut = {amostra:mutacao.split(',') for amostra, mutacao in zip(amostras, mutacoes)}

with open("input_clados.tsv",'r') as input:
    linhas = input.readlines()
with open("anotacoes_sinapomorfias.tsv","w") as output:
    output.write(f"clado_basal\tclado_filho\tdiferenca_mut_clados\tmut_clado_basal\tmut_clado_filho\tfreq_clado_basal\tfreq_clado_filho\n")
    for linha in linhas[1:]:
        clados = linha.rstrip('\n').split(',')

        dic_clado_taxons = {}
        # Adiciona um identificador único em cada nó interno

        for clado in clados:
            for node in tree.postorder_node_iter():
                if not node.is_leaf():
                    temp_list = []
                    if node.annotations.get_value("clado") == clado:

                        for taxon in node.leaf_nodes():
                            # print(dir(taxon))
                            nome_folha = str(taxon.taxon)
                            nome_folha = nome_folha.strip("'")
                            temp_list.append(nome_folha)
                        dic_clado_taxons[f"{clado}"] = temp_list
        taxons_clado_basal = diferenca(dic_clado_taxons[clados[0]], dic_clado_taxons[clados[1]])
        if len(clados) > 2:
            taxons_clado_filho = diferenca(dic_clado_taxons[clados[1]], dic_clado_taxons[clados[2]])
        else:
            taxons_clado_filho = dic_clado_taxons[clados[1]]



        mutacoes_clado_basal, freq_mutacoes_clado_basal = get_synapomorphy(taxons_clado_basal)

        mutacoes_clado_filho, freq_mutacoes_clado_filho = get_synapomorphy(taxons_clado_filho)

        diferenca_mut_clados = diferenca(mutacoes_clado_filho, mutacoes_clado_basal)
        # print(diferenca_mut_clados)
        # exit(0)
        output.write(f"{clados[0]}\t{clados[1]}\t{diferenca_mut_clados}\t{mutacoes_clado_basal}\t{mutacoes_clado_filho}\t{freq_mutacoes_clado_basal}\t{freq_mutacoes_clado_filho}\n")
