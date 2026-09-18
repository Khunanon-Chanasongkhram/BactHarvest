#!/usr/bin/env python3
"""
GenoTy.py v1.1
==============
Scans NCBI Genomes (Bacteria, Archaea, Eukaryota, Viruses) for specific genetic markers.

Features:
- Multi-Taxon Support: Bacteria, Archaea, Eukaryota, Viruses.
- Assembly Level Control: Reference, Annotated, or Total.
- Smart Filtering: Ignores regulators/repressors.
- API Key Support: High-speed scanning.

Usage:
------
# 1. Check Bacteria for nifH (Reference Genomes only)
python GenoTy.py --gene "nifH" --email "your@email.com" --api_key "YOUR_KEY_HERE" --taxon Bacteria --level Reference

# 2. Check Archaea for specific enzymes (All Annotated Genomes)
python GenoTy.py --gene "mcrA" --email "your@email.com" --api_key "YOUR_KEY_HERE" --taxon Archaea --level Annotated

# 3. Check Viruses (e.g., Phages) for a gene
python GenoTy.py --gene "terminase" --email "your@email.com" --api_key "YOUR_KEY_HERE" --taxon Viruses --level Total

Output:
-------
Saves a CSV file containing:
Organism Name | Genome Accession | GeneID | Description | Genotype Status
==============
Programmer: Khunanon Chanasongkhram
"""

import argparse
import sys
import time
import subprocess
import importlib.util
import re
import os

# ==========================================
# 0. COLOR SETTINGS & DEPENDENCIES
# ==========================================
class Col:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    CYAN = '\033[96m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

# Enables color support in Windows legacy command prompts
os.system("")

def check_and_install(package_name, import_name=None):
    """Checks if a Python library is installed. If not, installs it automatically."""
    if import_name is None:
        import_name = package_name
    if importlib.util.find_spec(import_name) is None:
        print(f"{Col.YELLOW}[*] Library '{package_name}' not found. Installing...{Col.RESET}")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
            print(f"{Col.GREEN}[+] '{package_name}' installed successfully!{Col.RESET}")
        except subprocess.CalledProcessError:
            print(f"{Col.RED}[!] Failed to install '{package_name}'. Please install manually.{Col.RESET}")
            sys.exit(1)

print(f"{Col.CYAN}[-] Checking system dependencies...{Col.RESET}")
check_and_install("biopython", "Bio")
check_and_install("pandas")
check_and_install("tqdm")

import pandas as pd
from tqdm import tqdm
from Bio import Entrez

# Keywords that indicate a gene is a regulator/repressor, not a functional marker
BAD_KEYWORDS = ["regulator", "repressor", "inhibitor", "pseudo"]

# ==========================================
# 1. THE GENOTYPER CLASS
# ==========================================

class GenoTyper:
    def __init__(self, email, api_key=None, taxon="Bacteria", level="Reference"):
        Entrez.email = email

        if api_key:
            Entrez.api_key = api_key
            print(f"{Col.CYAN}[-] API Key detected. High-speed mode enabled (10 req/sec).{Col.RESET}")
        else:
            print(f"{Col.CYAN}[-] No API Key. Running in standard mode (3 req/sec).{Col.RESET}")

        self.taxon = taxon
        self.level = level
        self.ref_genomes = {}

    def _build_genome_query(self):
        """Constructs the NCBI search string for the Genome Assembly database."""
        taxon_q = "all[filter]" if self.taxon == "All" else f'"{self.taxon}"[Organism]'

        if self.level == "Reference":
            level_q = 'AND "reference genome"[Filter]'
        elif self.level == "Annotated":
            level_q = 'AND "refseq has annotation"[Properties]'
        else:
            level_q = ""

        return f'({taxon_q} {level_q} AND "latest"[Filter])'

    def _build_gene_query(self, gene_name):
        """Constructs the search string for the Gene database."""
        taxon_q = "" if self.taxon == "All" else f'AND "{self.taxon}"[Organism]'
        return f'("{gene_name}"[Gene/Protein Name] OR "{gene_name}"[Symbol]) {taxon_q} AND srcdb_refseq[PROP]'

    def _fetch_batches(self, db, webenv, query_key, count, batch_size, desc):
        """
        Generator: fetches Entrez summaries in batches with retry logic.
        Yields a list of document summaries per batch.
        """
        for start in tqdm(range(0, count, batch_size), desc=desc):
            summaries = None
            for _ in range(3):
                try:
                    handle = Entrez.esummary(db=db, retstart=start, retmax=batch_size,
                                             webenv=webenv, query_key=query_key)
                    summaries = Entrez.read(handle)
                    break
                except Exception:
                    time.sleep(2)

            if summaries is None:
                continue

            if 'DocumentSummarySet' in summaries:
                yield summaries['DocumentSummarySet']['DocumentSummary']
            elif isinstance(summaries, list):
                yield summaries

    def _extract_gene_id(self, gene_doc):
        """Extracts GeneID from a gene document, handling varying NCBI XML formats."""
        if hasattr(gene_doc, 'attributes') and 'uid' in gene_doc.attributes:
            return gene_doc.attributes['uid']
        return gene_doc.get('Id', '') or gene_doc.get('uid', '')

    def fetch_reference_list(self):
        """
        STEP 1: Fetches the list of valid genomes (Assemblies) to build an allow-list
        for fast TaxID lookup during gene scanning.
        """
        query = self._build_genome_query()
        print(f"\n{Col.CYAN}[-] GenoTy: Mapping Genomes...{Col.RESET}")
        print(f"    Target: {Col.BOLD}{self.taxon}{Col.RESET} | Level: {Col.BOLD}{self.level}{Col.RESET}")

        try:
            handle = Entrez.esearch(db="assembly", term=query, retmax=500000, usehistory="y")
            record = Entrez.read(handle)
            count = int(record["Count"])
            webenv, query_key = record["WebEnv"], record["QueryKey"]

            print(f"{Col.YELLOW}[*] Found {count} matching genomes. Indexing details...{Col.RESET}")

            for doc_list in self._fetch_batches("assembly", webenv, query_key, count, 300,
                                                f"{Col.CYAN}Indexing{Col.RESET}"):
                for doc in doc_list:
                    self.ref_genomes[str(doc['Taxid'])] = {
                        "accession": doc['AssemblyAccession'],
                        "organism": doc['Organism']
                    }

            print(f"{Col.GREEN}[+] Successfully indexed {len(self.ref_genomes)} genomes.{Col.RESET}\n")

        except Exception as e:
            print(f"{Col.RED}[!] Critical Error fetching genome list: {e}{Col.RESET}")
            sys.exit(1)

    def search_gene(self, gene_name):
        """
        STEP 2: Searches for the gene, filters out false positives,
        and matches results to the genome allow-list.
        """
        print(f"{Col.CYAN}[-] GenoTy: Genotyping for marker '{gene_name}'...{Col.RESET}")

        try:
            handle = Entrez.esearch(db="gene", term=self._build_gene_query(gene_name),
                                    retmax=100000, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])
            webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]

            print(f"{Col.YELLOW}[*] Found {count} gene candidates. Filtering and Matching...{Col.RESET}")

            hits = []
            for doc_list in self._fetch_batches("gene", webenv, query_key, count, 300,
                                                f"{Col.CYAN}Scanning{Col.RESET}"):
                for gene_doc in doc_list:
                    try:
                        organism = gene_doc.get('Organism', {})
                        if 'TaxID' not in organism:
                            continue

                        gene_taxid = str(organism['TaxID'])
                        if gene_taxid not in self.ref_genomes:
                            continue

                        gene_symbol = gene_doc.get('Name', gene_name)
                        if gene_symbol.lower() != gene_name.lower():
                            continue

                        desc = gene_doc.get('Description', '')
                        if any(kw in desc.lower() for kw in BAD_KEYWORDS):
                            continue

                        ref_data = self.ref_genomes[gene_taxid]
                        clean_organism = re.sub(r'\s*\(.*?\)', '', ref_data['organism']).strip()
                        hits.append({
                            "Organism Name": clean_organism,
                            "Genome Accession": ref_data['accession'],
                            "Genotype Status": "Positive",
                            "Gene Symbol": gene_symbol,
                            "Gene Description": desc,
                            "GeneID": self._extract_gene_id(gene_doc),
                            "TaxID": gene_taxid
                        })
                    except KeyError:
                        continue

            return hits

        except Exception as e:
            print(f"{Col.RED}[!] Error during gene search: {e}{Col.RESET}")
            return []

    def run(self, gene_name, output_file):
        """Orchestrates the entire process: Indexing -> Scanning -> Saving."""
        self.fetch_reference_list()
        results = self.search_gene(gene_name)

        if results:
            df = pd.DataFrame(results)
            df = df.drop_duplicates(subset=['Genome Accession'])
            df = df.sort_values(by="Organism Name")
            df.to_csv(output_file, index=False)
            print(f"\n{Col.GREEN}[+] SUCCESS: Genotyped {len(df)} strains positive for '{gene_name}'.{Col.RESET}")
            print(f"{Col.GREEN}[+] Data saved to: {output_file}{Col.RESET}")
        else:
            print(f"\n{Col.RED}[-] No genomes found containing '{gene_name}' in the selected group.{Col.RESET}")

# ==========================================
# 2. MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Scans NCBI Genomes (Bacteria, Archaea, Eukaryota, Viruses) for specific genetic markers.")

    parser.add_argument("--gene", required=True, help="Target gene marker (e.g. 'nifH', 'lacZ')")
    parser.add_argument("--email", required=True, help="Your email for NCBI")
    parser.add_argument("--taxon", default="Bacteria",
                        choices=["Bacteria", "Archaea", "Eukaryota", "Viruses", "All"],
                        help="Taxonomic group to scan (Default: Bacteria)")
    parser.add_argument("--level", default="Reference",
                        choices=["Reference", "Annotated", "Total"],
                        help="Genome Assembly Level (Default: Reference)")
    parser.add_argument("--api_key", default=None, help="NCBI API Key (Speeds up processing)")
    parser.add_argument("--out", default=None, help="Output CSV filename")

    args = parser.parse_args()

    if not args.out:
        clean_gene = args.gene.replace(" ", "_").replace("/", "-")
        args.out = f"GenoTy_{args.taxon}_{clean_gene}.csv"

    typer = GenoTyper(args.email, args.api_key, args.taxon, args.level)
    try:
        typer.run(args.gene, args.out)
    except KeyboardInterrupt:
        print(f"\n{Col.RED}[!] Interrupted by user. Exiting...{Col.RESET}")
