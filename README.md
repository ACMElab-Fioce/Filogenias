# Filogenias

Esse conjunto de scripts serve para identificar as sinapomorfias que foram usadas para agrupar sequências em uma dada árvore filogenética.

O script "nomeia_clado.py" recebe uma árvore filogenética em formato nexus e escreve um nome genérico para cada clado e escreve a árvore modificada.

O script "checar_sinapomorfias.py" recebe a árvore modificada do script anterior, o arquivo "input_clados.tsv" (onde terão os clados de interesse) e um arquivo com as mutações de cada amostra.

Caso alguém queira analisar não o clado inteiro, mas apenas as sequências da politomia do clado, deve-se passar no arquivo "input_clados.tsv" o clado imediatamente posterior ao clado filho.
