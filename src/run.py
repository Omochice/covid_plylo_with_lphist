import json
import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
from Bio import SeqIO
from Bio.Phylo import TreeConstruction, draw
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

import plot
import visualize
from local_pattern import make_lphist

project_dir = Path(__file__).resolve().parents[1]
fastafile = os.path.join(project_dir, "data/sequences.fasta")
use_file = os.path.join(project_dir, "data/use_sequences.fasta")
# use_file = os.path.join(project_dir, "data/med.fasta")


def extraction(infh: str, outfh: str) -> None:
    print("extraction...")
    with open(outfh, "w") as f:
        for entry in SeqIO.parse(infh, "fasta"):
            desc = entry.description
            if desc.split("|")[-1].strip() == "complete genome":
                SeqIO.write(entry, f, "fasta")


def make_lower_triangular_matrix(names: list, targets: list) -> tuple:
    pass


def get_cmap(countries: pd.Series, ) -> dict:
    cmap = [
        "red", "blue", "yellow", "green", "pink", "aqua", "maroon", "black",
        "navy", "lime", "purple", "teal", "gold", "orange"
    ]
    countries = set(countries.values.tolist())
    print((countries))
    d = {}
    for i, c in zip(countries, cmap):
        d[i] = generate_random_color()
    print(d, cmap)
    return d


def main():
    df = pd.read_csv(os.path.join(project_dir, "data/sequences.csv"),
                     index_col=1)

    extraction(fastafile, use_file)
    weight = plot.calc_weights(use_file)
    results = {}
    lphist_table = []
    names = []
    # country_to_color = get_cmap(df["Geo_Location"])
    cmap = {}

    for entry in SeqIO.parse(use_file, "fasta"):
        accession = entry.id
        country = df.at[accession, "Geo_Location"]

        x_coordinates, y_coordinates = plot.calculate_coordinate(
            entry.seq, weight)
        lphist = make_lphist(x_coordinates, y_coordinates)
        results[accession] = {
            "accession": accession,
            "country": country,
            "x_coordinates": x_coordinates,
            "y_coordinates": y_coordinates,
            "norm_lphist": lphist.tolist()
        }
        show_name = f"{accession}-{country}"

        lphist_table.append(lphist)
        names.append(show_name)
        if country == "Japan":
            cmap[show_name] = "red"

    with open(os.path.join(project_dir, "data/result.json"), "w") as f:
        json.dump(results, f, indent=4)

    with open("triplet_weight.json", "w") as f:
        json.dump(weight, f)
    visualize.make_codon(weight,
                         "triplet_table",
                         os.path.join(project_dir, "data"),
                         save_xlxs=False)

    for entry in results.values():
        save_path = os.path.join(project_dir, "data/img",
                                 entry["accession"] + ".png")
        figure = plot.generate_image(entry, save_path)

    # 距離行列（下三角）
    l_name = len(names)
    dm = np.zeros((l_name, l_name))
    for i in range(l_name):
        for j in range(i, l_name):
            dm[j, i] = np.sum((lphist_table[i] - lphist_table[j])**2)**0.5
            # dm[j, i] = np.sum(abs(lphist_table[i] - lphist_table[j]))

    # plylo用のlist
    dm_list = []
    for i in range(l_name):
        dm_list.append(list(dm[i, :i + 1]))

    #============================
    #Neighbor Joining法
    #============================
    #Distance matrix型の用意
    DM = TreeConstruction._DistanceMatrix(names, dm_list)

    #nj木に変換
    tree = DistanceTreeConstructor().nj(DM)

    #クレード名か'Inner#'のものを''に改名する（見づらいので）
    for clade in tree.find_clades():
        if 'Inner' in clade.name:
            clade.name = ''

    fig = plt.figure(figsize=(30, 30), dpi=150)
    ax = fig.add_subplot(1, 1, 1)
    draw(tree, axes=ax, label_colors=cmap)
    # for i, (key, value) in enumerate(cmap.items()):
    #     fig.text(1, 1 - 0.05 * i, key, color=value, ha="right")
    # plt.savefig('aaa.png')
    plt.title("corona virus")
    plt.show()


if __name__ == "__main__":
    main()
