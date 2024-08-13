#!/usr/bin/env python3
"""
promer-deid-simplified-customized
Miguel F. Romero, 2023
github.com/miferg
Format genomic mapping information to create depth and identity plots.
Edited by Gerardo Mejia, 2024
github.com/GerardoMej
Usage:
promer-deid-simplified-customized.py INPUT.show-coords.tsv OUTPUT_PREFIX

Make sure to previously run: show-coords -c -d -k -l -r -T INPUT.delta > INPUT.show-coords.tsv
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt

def get_promer(filename):
    """Abre el archivo delta y recupera los alineamientos, luego los guarda en un DataFrame"""
    promer = pd.read_csv(filename, sep='\t', skiprows=4, header=None)
    promer.columns = ['ref-start', 'ref-end', 'query-start', 'query-end', 'ref-alilen', 'query-alilen', 'identity', 'similarity', 'stop-codons-%', 'ref-length', 'query-length',
                      'ref-coverage', 'query-coverage', 'ref-frame', 'query-frame', 'ref-name', 'query-name']
    # Agregamos una nueva columna 'index' para el índice numérico
    promer['index'] = range(1, len(promer) + 1)
    return promer

def get_depth_d(promer):
    depth_d = {}
    for refname in promer['ref-name'].unique():
        depth_d[refname] = promer[promer['ref-name'] == refname].shape[0]
    return depth_d

def get_depth(depth_d):
    depth = pd.DataFrame.from_dict(depth_d, orient='index', columns=['depth']).reset_index()
    depth.columns = ['seqname', 'depth']
    # Agrega la columna 'index' para que coincida con la gráfica de identity
    depth['index'] = range(1, len(depth) + 1)
    return depth

def plot_promer(depth, promer, outprefix):
    plt.rcParams["figure.figsize"] = [12, 6]  # Ajusta el tamaño de la figura
    f, ax = plt.subplots(2, 1)

    depth.plot.line(x='index', y='depth', color='silver', style='-', alpha=1, ax=ax[0])
    ax[0].set_ylabel('Depth')
    ax[0].set_xlabel('')
    ax[0].get_legend().remove()
    ax[0].xaxis.set_tick_params(labelbottom=False)

    # Calcula el número de secuencias por cada 3,000 secuencias y gráficalo
    x_ticks = range(0, len(promer) + 1, 3000)
    x_labels = [str(x) for x in x_ticks]
    ax[1].set_xticks(x_ticks)
    ax[1].set_xticklabels(x_labels, rotation=0)
    ax[1].set_xlim(0, len(promer) + 1)
    ax[1].set_ylim(0, 100)  # Ajusta el límite del eje y

    # Ajusta el tamaño de los puntos en la gráfica de identity
    promer.plot.scatter(x='index', y='identity', c=promer.apply(lambda row: 'orchid' if row['query-frame'] * row['ref-frame'] < 0 else 'darkturquoise', axis=1), s=10, alpha=0.6, ax=ax[1])
    ax[1].set_ylabel('Identity %')
    ax[1].set_xlabel('Index')  # Etiqueta del eje x

    plt.savefig(outprefix + '.deid.png', dpi=300)

def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    coordsfile = sys.argv[1]
    outprefix = sys.argv[2]

    promer = get_promer(coordsfile)
    depth_d = get_depth_d(promer)
    depth = get_depth(depth_d)
    plot_promer(depth, promer, outprefix)

    print('Plot saved as ' + outprefix + '.deid.png')

    depth.to_csv(outprefix + '.depth.tsv', sep='\t', index=False)
    print('Depth table saved as ' + outprefix + '.depth.tsv')

if __name__ == '__main__':
    main()

