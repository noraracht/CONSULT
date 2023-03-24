import csv
import argparse
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Creates a taxonomic ancestor lookup table."
)
parser.add_argument(
    "--input-path",
    type=str,
    required=True,
)
parser.add_argument(
    "--output-dir", type=str, required=False, default="./taxonomy_lookup/"
)
parser.add_argument(
    "--from-taxonomy-database",
    action=argparse.BooleanOptionalAction,
)
args = parser.parse_args()


PATH_RANKS = Path(args.input_path)
DIR_LOOKUPS = Path(args.output_dir)
DIR_LOOKUPS.mkdir(parents=True, exist_ok=True)

TAXA_LEVELS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

if __name__ == "__main__":
    taxID_lookup = defaultdict(list)
    genome_lookup = defaultdict(int)
    level_lookup = defaultdict(str)

    if args.from_taxonomy_database:
        taxID_reldict = {}
        with open(PATH_RANKS, "r", newline="") as csvfile:
            reader = list(csv.reader(csvfile, delimiter="|"))
            for row in reader:
                tid1 = int(str(row[0]).strip())
                tid2 = int(str(row[1]).strip())
                taxID_reldict[tid1] = tid2
                level_lookup[tid1] = str(row[2]).strip()

        for tid1, tid2 in taxID_reldict.items():
            if level_lookup[tid1] in TAXA_LEVELS:
                taxID_lookup[tid1].append(tid1)
                parent = tid2
                while parent != 1:
                    if level_lookup[parent] in TAXA_LEVELS:
                        taxID_lookup[tid1].append(parent)
                    parent = taxID_reldict[parent]
    else:
        with open(PATH_RANKS, "r", newline="") as csvfile:
            reader = list(csv.DictReader(csvfile, delimiter="\t"))
        for row_ranks in reader:
            genome_lookup[row_ranks["genome"]] = row_ranks["species"]
            for ix, taxa1 in enumerate(TAXA_LEVELS[:]):
                tid1 = row_ranks[taxa1]
                if tid1 not in taxID_lookup.keys() and tid1 != 0:
                    for taxa2 in TAXA_LEVELS[ix:]:
                        taxID_lookup[tid1].append(row_ranks[taxa2])
                    level_lookup[row_ranks[taxa1]] = taxa1
                else:
                    break

        with open(DIR_LOOKUPS / "genome_lookup", "w") as f:
            lookup_str = ""
            for key, val in genome_lookup.items():
                lookup_str += f"{key} {val}"
                lookup_str += "\n"
            f.write(lookup_str)

    with open(DIR_LOOKUPS / "taxonomy_lookup", "w") as f:
        lookup_str = ""
        for taxa, ancestors in taxID_lookup.items():
            lookup_str += f"{taxa} {','.join([str(a_taxa) for a_taxa in ancestors])}"
            lookup_str += "\n"
        f.write(lookup_str)

    with open(DIR_LOOKUPS / "level_lookup", "w") as f:
        lookup_str = ""
        for key, val in level_lookup.items():
            lookup_str += f"{key} {val}"
            lookup_str += "\n"
        f.write(lookup_str)
