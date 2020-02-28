import re
from functools import cmp_to_key
from itertools import product

import numpy as np
from Bio import SeqIO


def get_taxonID(record) -> str:
    for feature in record.features:
        if feature.type == "source":
            db_xref = feature.qualifiers["db_xref"][0]  # 2回以上sourceが現れるのはいないはず
            taxonomyID = db_xref.split(":")[1]
            return taxonomyID


def taxon_cmp(a: str, b: str) -> int:
    prefixes = ("super", "", "sub", "infra", "parv")
    suffixes = ("kingdom", "phylum", "class", "cohort", "order", "family",
                "tribe", "genus", "species")
    sortlist = {}
    for i, (suf, pre) in enumerate(product(suffixes, prefixes)):
        sortlist[pre + suf] = i

    if a == b:
        return 0
    elif sortlist.get(a, -1) < sortlist.get(b, -1):
        return -1
    else:
        return 1


def window_search(sequence: str, each=3, overhang=None) -> str:
    pretty_sequence = re.findall(r"[atgc]", str(sequence.lower()))

    if overhang == "before":
        for i in range(1, each):
            yield "".join(pretty_sequence[:i])

    for i in range(0, len(pretty_sequence) - each + 1):
        yield "".join(pretty_sequence[i:i + each])

    if overhang == "after":
        for i in range(-1 * each + 1, 0, 1):
            yield "".join(pretty_sequence[i:])


def information_content(counter: dict) -> dict:
    bases = {"a", "t", "g", "c"}

    weight_dict = {}
    for triplet, content in counter.items():
        first_2_base = triplet[:2]  # ここがATなら
        # denominator = 0
        # for base in bases:
        #     denominator += counter[first_2_base +
        #                            base]  # ここはATA+ATT+ATG+ATCになる
        denominator = sum([counter[first_2_base + base] for base in bases])
        weight_dict[triplet] = -1 * np.log2(content / denominator)

    return weight_dict


if __name__ == "__main__":
    temp = [
        'parvorder', 'genus', 'superkingdom', 'kingdom', 'phylum', 'subphylum',
        'superclass', 'class', 'subclass', 'infraclass', 'cohort', 'subcohort',
        'superorder', 'order', 'suborder', 'infraorder', 'superfamily',
        'family', 'subfamily', 'species'
    ]
    print(sorted(temp, key=cmp_to_key(taxon_cmp)))
