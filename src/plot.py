import glob
import json
import os
from collections import Counter, OrderedDict
from itertools import product

import matplotlib.pyplot as plt
from Bio import SeqIO

import utils


def calc_weights(infh: str, alignment: bool = True) -> dict:
    print("Calc weights ...")
    # classes = os.listdir(root_dir)
    result = Counter()

    for entry in SeqIO.parse(infh, "fasta"):
        result += Counter(utils.window_search(entry.seq, each=3))

    # for class_name in classes:
    #     print(f"{class_name}...", end="")
    #     filenames = glob.glob(os.path.join(root_dir, class_name, "*.gbk"))
    #     if len(filenames) >= lower_limit:
    #         for filename in filenames:
    #             for entry in SeqIO.parse(filename, "genbank"):
    #                 result += Counter(utils.window_search(entry.seq, each=3))
    #     print("finished.")
    self_info = utils.information_content(result)

    if alignment:
        output = OrderedDict()
        for base_comb in product(("a", "t", "g", "c"), repeat=3):
            triplet = "".join(base_comb)
            output[triplet] = self_info[triplet]
        return output
    else:
        return self_info


def calculate_coordinate(seq: str, weight: dict) -> tuple:
    VECTORS = {"a": (1, 1), "t": (-1, 1), "g": (-1, -1), "c": (1, -1)}
    x_coordinates = [0]
    y_coordinates = [0]

    for triplet in utils.window_search(seq, overhang="before"):
        x_coordinates.append(x_coordinates[-1] +
                             VECTORS[triplet[-1]][0] * weight.get(triplet, 1))
        y_coordinates.append(y_coordinates[-1] +
                             VECTORS[triplet[-1]][1] * weight.get(triplet, 1))

    return x_coordinates, y_coordinates


def generate_image(entry: dict, savepath: str) -> None:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x_coo = entry["x_coordinates"]
    y_coo = entry["y_coordinates"]
    ax.plot(x_coo, y_coo, color="black", linewidth=1)

    ax.set_xlim(min(x_coo), max(x_coo))
    ax.set_ylim(min(y_coo), max(y_coo))
    ax.set_aspect('equal')
    ax.set_title(f"{entry['accession']} in {entry['country']}")
    plt.savefig(savepath)
    plt.close()


if __name__ == "__main__":
    gbk_data = "/home/mochi/human.gbk"

    ### test
    with open(
            "/home/mochi/workspace/senior_thesis/generate_weight/weight_limit_5.json",
            "r") as f:
        weights = json.load(f)
    for entry in SeqIO.parse(gbk_data, "genbank"):
        generate_image(entry.seq, "/home/mochi/testfig.png", weights)
