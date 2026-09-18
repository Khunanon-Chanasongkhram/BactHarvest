import pandas as pd
import glob

OUTPUT_FILENAME = "GenoTy_Merged_Output.csv"

COLUMN_ORDER = ['Organism Name', 'Genome Accession', 'Gene Symbol', 'Gene Count',
                'Confidence', 'GeneID', 'Gene Description', 'Genotype Status', 'TaxID']

def merge_genoty_files():
    # Exclude the output file itself to prevent re-merging on repeated runs
    input_files = [f for f in glob.glob("*.csv") if f != OUTPUT_FILENAME]
    print(f"[-] Found {len(input_files)} CSV files to merge.")

    dfs = []
    for filename in input_files:
        try:
            df = pd.read_csv(filename)
            if 'Genome Accession' in df.columns:
                dfs.append(df)
            else:
                print(f"[!] Skipping {filename}: 'Genome Accession' column not found.")
        except Exception as e:
            print(f"[!] Error reading {filename}: {e}")

    if not dfs:
        print("[!] No valid GenoTy CSV files found to merge.")
        return

    combined_df = pd.concat(dfs, ignore_index=True)

    agg_rules = {
        'Gene Symbol':      lambda x: ', '.join(sorted(set(x.dropna().astype(str)))),
        'Gene Description': lambda x: ' | '.join(x.dropna().astype(str)),
        'GeneID':           lambda x: ', '.join(x.dropna().astype(str)),
        'Organism Name':    'first',
        'Genotype Status':  'first',
        'TaxID':            'first',
    }
    # Only aggregate columns that are present in the data
    final_agg_rules = {k: v for k, v in agg_rules.items() if k in combined_df.columns}

    print("[-] Merging data...")
    merged_df = combined_df.groupby('Genome Accession', as_index=False).agg(final_agg_rules)

    merged_df['Gene Count'] = merged_df['Gene Symbol'].apply(
        lambda x: len([g for g in x.split(',') if g.strip()])
    )
    merged_df['Confidence'] = merged_df['Gene Count'].apply(
        lambda n: 'High' if n == 8 else 'Low'
    )

    cols_to_use = [c for c in COLUMN_ORDER if c in merged_df.columns]
    merged_df = merged_df[cols_to_use]

    merged_df.to_csv(OUTPUT_FILENAME, index=False)
    print(f"[+] Success! Merged {len(combined_df)} rows into {len(merged_df)} unique genomes.")
    print(f"[+] Output saved to: {OUTPUT_FILENAME}")

if __name__ == "__main__":
    merge_genoty_files()
