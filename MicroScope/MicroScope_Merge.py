import pandas as pd
import glob
import os
import sys

folder_path = os.path.dirname(os.path.abspath(__file__))
all_files = sorted(glob.glob(os.path.join(folder_path, "results_*.csv")))

if not all_files:
    print("No result CSV files found (expected files matching results_*.csv).")
    sys.exit(1)

print(f"Found {len(all_files)} CSV file(s):")

df_list = []
for filename in all_files:
    print(f"  - Reading: {os.path.basename(filename)}")
    try:
        df = pd.read_csv(filename, index_col=None, header=0)
        df_list.append(df)
    except Exception as e:
        print(f"    WARNING: Skipped {os.path.basename(filename)} — {e}")

if not df_list:
    print("No data could be read from any file.")
    sys.exit(1)

df_combined = pd.concat(df_list, axis=0, ignore_index=True)
print(f"\nTotal rows before deduplication: {len(df_combined)}")

# Determine correct column names (MicroScope outputs capitalised headers)
pmid_col = next((c for c in df_combined.columns if c.strip().upper() == "PMID"), None)
org_col = next((c for c in df_combined.columns if c.strip().lower() == "organism"), None)

if pmid_col and org_col:
    df_cleaned = df_combined.drop_duplicates(subset=[pmid_col, org_col], keep='first')
else:
    print("WARNING: Could not find PMID/Organism columns — deduplicating on all columns.")
    df_cleaned = df_combined.drop_duplicates(keep='first')

print(f"Unique rows after deduplication: {len(df_cleaned)}")

# Sort by performance_value descending (matches MicroScope output order)
val_col = next((c for c in df_cleaned.columns if c.strip().lower() == "performance_value"), None)
if val_col:
    df_cleaned = df_cleaned.sort_values(val_col, ascending=False, na_position='last')

output_path = os.path.join(folder_path, 'Siderophore_Combined_Cleaned.csv')
df_cleaned.to_csv(output_path, index=False, encoding='utf-8')
print(f"\nMerged file saved to: {output_path}")

print("\nSample data (first 5 rows):")
print(df_cleaned.head().to_string(index=False))
