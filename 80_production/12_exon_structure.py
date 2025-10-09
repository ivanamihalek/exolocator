#!/usr/bin/env python3
from sys import stdout

from MySQLdb.constants.FIELD_TYPE import DATETIME

from el_utils.processes import *
from el_utils.special_gene_sets import *
from el_utils.known_exon_utils import exons_for_gene

from config import Config

import numpy as np
import pandas as pd

def add_cds_protein_coords(
    df,
    gene_col="gene_id",
    start_col="start_in_gene",
    end_col="end_in_gene",
    strand_col="strand",
    is_coding_col="is_coding",
    is_canonical_col="is_canonical",
    canon_transl_start_col="canon_transl_start",
    canon_transl_end_col="canon_transl_end",
):
    """
    Add CDS and protein coordinates per canonical coding exon.
    - cds_start/cds_end are 1-based positions along the CDS (transcript order).
    - protein_start/protein_end are amino-acid indices = ceil(cds_pos/3).
    - coding_genomic_start/end are the genomic interval of the coding portion inside the exon,
      trimmed by the transcript-level canonical translation interval if present.

    Only exons with is_canonical==1 AND is_coding==1 are assigned; others remain NA.
    """
    df = df.copy()

    # new columns
    for c in (
        "coding_genomic_start",
        "coding_genomic_end",
        "cds_start",
        "cds_end",
        "protein_start",
        "protein_end",
    ):
        df[c] = pd.NA

    # coerce numeric fields (robust to float/int mix in input)
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")
    df[canon_transl_start_col] = pd.to_numeric(df[canon_transl_start_col], errors="coerce")
    df[canon_transl_end_col] = pd.to_numeric(df[canon_transl_end_col], errors="coerce")

    # process gene by gene
    for gene_id, gdf in df.groupby(gene_col):
        # select canonical coding exons (we only assign CDS/protein for these)
        canon_mask = (gdf[is_canonical_col] == 1) & (gdf[is_coding_col] == 1)
        canon = gdf[canon_mask].copy()
        if canon.empty:
            continue

        # determine transcript-level canonical translation interval (if any)
        transl_values = pd.concat(
            [gdf[canon_transl_start_col].dropna(), gdf[canon_transl_end_col].dropna()]
        )
        if len(transl_values) > 0:
            transl_min = int(transl_values.min())
            transl_max = int(transl_values.max())
            use_transl_bounds = True
        else:
            transl_min = transl_max = None
            use_transl_bounds = False

        # normalize strand token (allow '-' or '-1' etc)
        strand_token = str(canon[strand_col].iloc[0]).strip()
        neg_strand = strand_token in ("-", "-1")

        # sort exons in transcript order
        # For + strand: ascending genomic coordinates; for - strand: descending
        canon_sorted = canon.sort_values(by=start_col, ascending=not neg_strand)

        # accumulate CDS positions along transcript (1-based)
        cds_cursor = 1
        for idx, row in canon_sorted.iterrows():
            s = int(row[start_col])
            e = int(row[end_col])
            exon_lo = min(s, e)
            exon_hi = max(s, e)

            # intersect exon with transcript-level translation bounds if available
            if use_transl_bounds:
                coding_lo = max(exon_lo, transl_min)
                coding_hi = min(exon_hi, transl_max)
            else:
                coding_lo = exon_lo
                coding_hi = exon_hi

            if coding_hi < coding_lo:
                # this canonical exon has no coding bases after clipping -> skip
                continue

            coding_len = coding_hi - coding_lo + 1
            cds_start = cds_cursor
            cds_end = cds_cursor + coding_len - 1

            # store
            df.at[idx, "coding_genomic_start"] = int(coding_lo)
            df.at[idx, "coding_genomic_end"] = int(coding_hi)
            df.at[idx, "cds_start"] = int(cds_start)
            df.at[idx, "cds_end"] = int(cds_end)
            # protein indices: ceil(cds_pos/3)  -> (cds_pos + 2) // 3
            df.at[idx, "protein_start"] = int((cds_start + 2) // 3)
            df.at[idx, "protein_end"] = int((cds_end + 2) // 3)

            cds_cursor = cds_end + 1

    return df

import json
import pandas as pd

def add_protein_features(
    df,
    json_path,
    protein_start_col="protein_start",
    protein_end_col="protein_end",
    min_overlap_fraction=0.25,
):
    """
    Annotate exons with overlapping UniProt features.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing protein_start and protein_end columns.
    json_path : str
        Path to UniProt JSON file with 'features' list.
    protein_start_col, protein_end_col : str
        Column names for exon protein coordinate range.
    min_overlap_fraction : float
        Minimum fraction (0-1) of exon protein range overlapping a feature
        for it to count as "largely overlapping".

    Returns
    -------
    pd.DataFrame
        Input dataframe with added columns for feature overlaps.
    """

    df = df.copy()

    # --- Load UniProt feature JSON ---
    with open(json_path, "r") as f:
        data = json.load(f)

    features = data.get("features", [])
    parsed_features = []
    for feat in features:
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        desc = feat.get("description", "")
        ftype = feat.get("type", "")
        if start is None or end is None:
            continue
        parsed_features.append(
            {"start": int(start), "end": int(end), "desc": desc, "type": ftype}
        )

    # --- Annotate exons ---
    for i, row in df.iterrows():
        exon_start = row[protein_start_col]
        exon_end = row[protein_end_col]

        if pd.isna(exon_start) or pd.isna(exon_end):
            continue

        overlaps = []
        exon_len = exon_end - exon_start + 1

        for feat in parsed_features:
            fstart, fend = feat["start"], feat["end"]

            # compute overlap
            overlap_start = max(exon_start, fstart)
            overlap_end = min(exon_end, fend)
            overlap = max(0, overlap_end - overlap_start + 1)
            frac = overlap / exon_len if exon_len > 0 else 0

            if frac >= min_overlap_fraction:
                overlaps.append(feat)

        # Add triplets for each overlapping feature
        for j, feat in enumerate(overlaps, start=1):
            df.loc[i, f"feature{j}_start"] = feat["start"]
            df.loc[i, f"feature{j}_end"] = feat["end"]
            df.loc[i, f"feature{j}_desc"] = feat["desc"] or feat["type"]

    return df



#########################################
def main():
    gene_name = "USH2A"
    outdir = "/home/ivana/scratch"
    species = 'homo_sapiens'
    db = connect_to_mysql(Config.mysql_conf_file)
    cursor = db.cursor()
    qry = f"select ensembl_gene_id  from identifier_maps.hgnc where approved_symbol='{gene_name}'"
    ensembl_stable_gene_id = hard_landing_search(cursor, qry)[0][0]
    print(ensembl_stable_gene_id)

    db_name ="homo_sapiens_core_113_38"
    gene_id = stable2gene(cursor, ensembl_stable_gene_id, db_name=db_name)
    print(gene_id)
    outfnm = f"{gene_name}.exons.tsv"
    with open(outfnm, "w") as outfile:
        table_header = ['exon no',"gene_id", "exon_id", "start_in_gene", "end_in_gene",
                      "canon_transl_start", "canon_transl_end", "exon_seq_id", "strand", "phase",
                      "provenance", "is_coding", "is_canonical", "is_constitutive", "covering_exon",
                      "covering_provenance", "analysis_id"]
        print("\t".join(table_header), file=outfile)
        exons_for_gene(cursor, gene_id, db_name, 0, outfile, stdout)

    df = pd.read_csv(outfnm, sep="\t")
    df['strand'] = df['strand'].replace(-1, '-')
    df['strand'] = df['strand'].replace(1, '+')
    # normalize missing translation markers (your data used -1)
    df['canon_transl_start'] = df['canon_transl_start'].replace(-1, np.nan)
    df['canon_transl_end']   = df['canon_transl_end'].replace(-1, np.nan)

    print(df)
    df2 =  add_cds_protein_coords(df)
    print(df2)
    df3 = add_protein_features(df2, "USH2A.domains.json", min_overlap_fraction=0.25)

    outfnm2 = f"{gene_name}.exons.annotated.tsv"
    df3.to_csv(outfnm2, sep="\t", index=False)



#########################################
if __name__ == '__main__':
    main()

'''
In v 101 a bunch of canonical transcript coordinates were missing for mus caroli
not sure if I should worry about that - there ar 75 cases like that
exmaple MGP_CAROLIEiJ_G0027698 Ppp1cc protein phosphatase 1 catalytic subunit gamma
The smae for mus pahari and mus spretus
'''