# initialize
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import pandas as pd

# email to use the Entrez service
Entrez.email = "max@alignbio.org"

# Function to fetch the sequence from NCBI
def fetch_sequence(accession_number):
    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

# Function to extract coordinates and sequences for the genes of interest
def extract_gene_info(genome_record, verbose=False):
    gene_info = []
    for feature in genome_record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", [""])[0]
            start = int(feature.location.start)  # Convert ExactPosition to int
            end = int(feature.location.end)      # Convert ExactPosition to int
            strand = "+" if feature.location.strand == 1 else "-"
            sequence = feature.extract(genome_record.seq)
            
            if verbose:
                print(f"Gene: {gene_name}, Start: {start}, End: {end}, Strand: {strand}")
            
            gene_info.append({
                "Gene": gene_name,
                "Start": start,
                "End": end,
                "Strand": strand,
                "Sequence": str(sequence)
            })
    
    # Return as a DataFrame for easier manipulation
    return pd.DataFrame(gene_info)

