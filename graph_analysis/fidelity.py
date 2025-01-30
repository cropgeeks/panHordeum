"""Fidelity check between the input and graph versions of a genome.

It will compare homologous sequences from the both versions of the genome
and check if their sequence are exactly identical or not. The homology
is assumed to be given by the order of the sequences in both files.

For each genome, reports which of the sequences are identical as a 
dictionary with the name of the sequence as key and the boolean result as value.
If all the sequences are identical, returns a global True instead.

This script assumes that:

- The original genomes are in folder './query_genomes' in fasta format.
- The graph version are in folder './graph_genomes' in fasta format.
- There is a 1 to 1 correspondence between files in both folders.
- The files are named in a way that both versions of the same genome conincide
  in their order (1st <-> 1st, 2nd <-> 2nd, and so on). This is guaranteed, 
  for instance, is the files have the same name (although it is not the only way).
  
Dependencies:
- Biopython
"""

# =============================================================================
#               IMPORTS
# =============================================================================

from pathlib import Path
from Bio import SeqIO


# =============================================================================
#                FUNCTIONS
# =============================================================================


def comp_genomes(query, graph) -> bool | dict[str, bool]:
    comp_dict = {
        recq.id: recq.seq.upper() == recg.seq.upper()
        or recq.seq.upper() == recg.seq.reverse_complement().upper()
        for recq, recg in zip(query, graph)
    }
    return all(comp_dict.values()) or comp_dict


def main():

    query_dir = Path("query_genomes")
    graph_dir = Path("graph_genomes")

    results = {}
    for query_path, graph_path in zip(
        sorted(query_dir.glob("*.fasta")), sorted(graph_dir.glob("*.fasta"))
    ):
        query = SeqIO.parse(query_path, "fasta")
        graph = SeqIO.parse(graph_path, "fasta")
        results[query_path.stem] = comp_genomes(query, graph)

    for genome, result in results.items():
        if isinstance(result, dict):
            print(f"{genome} is NOT identical in all their sequences")
            for seqname, val in result.items():
                res = "identical" if val else "different"
                print(f"\t - {seqname}: {res}")
        else:
            print(f"{genome} IS identical in all their sequences")


if __name__ == "__main__":
    main()
