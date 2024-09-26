from glob import glob
from functions import bash


def fetch_research_data():
    bash("rm -rf Research-Data")
    bash("git clone https://github.com/noashachar/Research-Data")
    bash("rm -rf Research-Data/.git")
    # this item is too similar to CR/3, oversaturating similarites between other items
    bash("rm Research-Data/CR/2.*")


def get_paths_pairs():
    return [
        *zip(
            sorted(glob("Research-Data/*/*.desc.txt")),
            sorted(glob("Research-Data/*/*.genes.txt")),
        )
    ]


paths_pairs = get_paths_pairs()

if not paths_pairs:
    fetch_research_data()
    paths_pairs = get_paths_pairs()
