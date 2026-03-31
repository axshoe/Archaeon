"""
src/data/ncbi.py
----------------
NCBI Entrez API wrapper for Archaeon.

Handles two tasks:
  1. Searching the NCBI protein database for sequences from extreme-environment
     biomes using Entrez esearch + efetch (via Biopython).
  2. Running BLASTp similarity searches against known thermostable enzyme
     families to identify which reference enzyme a candidate most resembles
     and how similar it is (expressed as percent identity).

Why NCBI: The NCBI protein database (nr) is the most comprehensive public
protein database on the planet. It contains sequences from every sequenced
organism including thousands of extremophiles from hot springs, hydrothermal
vents, hypersaline lakes, and acidic environments. Free, documented, REST API.

Why BLAST: Sequence homology search is the standard first-pass method for
predicting function. If a new sequence shares 60%+ identity with a known
lipase from Sulfolobus acidocaldarius, it almost certainly IS a lipase. We
use BLAST identity as one signal in the IAS scoring pipeline.

Design note: This module deliberately avoids writing to disk during fetching.
All records are returned as in-memory dictionaries. The main CLI decides what
to persist.
"""

import os
import time
import logging
from typing import Optional
from io import StringIO

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------
# Known thermostable enzyme families to BLAST against.
# These are well-characterized extremophile enzyme reference sequences.
# Format: {family_name: NCBI accession of canonical reference sequence}
# Extending this dict adds new enzyme families to the pipeline with zero
# code changes.
# -------------------------------------------------------------------------
ENZYME_FAMILIES = {
    "lipase":     "AAC69568.1",   # Thermus thermophilus lipase
    "protease":   "AAA27706.1",   # Thermolysin (Bacillus thermoproteolyticus)
    "amylase":    "AAA22222.1",   # Alpha-amylase (Bacillus stearothermophilus)
    "cellulase":  "AAB42154.1",   # Endoglucanase (Thermotoga maritima)
    "xylanase":   "CAA93248.1",   # Xylanase (Thermomyces lanuginosus)
    "laccase":    "AAW33585.1",   # Laccase (Trametes versicolor thermostable)
    "peroxidase": "CAA32235.1",   # Manganese peroxidase (lignin-degrading)
    "isomerase":  "AAC43587.1",   # Glucose isomerase (Thermoanaerobacterium)
}

# Biome search terms that reliably return extremophile sequences from NCBI.
# These use NCBI organism/environment metadata fields.
EXTREME_BIOME_QUERIES = [
    "hot spring metagenome[organism]",
    "hydrothermal vent metagenome[organism]",
    "Sulfolobus[organism] AND enzyme[title]",
    "Thermus thermophilus[organism]",
    "Pyrococcus[organism] AND enzyme[title]",
]


def configure_entrez(email: str, api_key: Optional[str] = None) -> None:
    """
    Set Entrez credentials. Must be called before any NCBI queries.

    NCBI rate limits:
      - Without API key: 3 requests/second
      - With API key:    10 requests/second

    Always provide an email; NCBI uses it to contact you if your script
    is overloading their servers. Non-optional for production use.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    logger.info(f"Entrez configured for {email}")


def search_ncbi_proteins(
    query: str,
    database: str = "protein",
    max_results: int = 50,
    min_length: int = 100,
    max_length: int = 2000,
) -> list:
    """
    Search NCBI protein database and fetch sequence records.

    Args:
        query:       Entrez search string (see EXTREME_BIOME_QUERIES above)
        database:    NCBI database name; 'protein' is the right choice here
        max_results: Upper bound on returned records
        min_length:  Discard sequences shorter than this (aa). Very short
                     fragments are usually assembly artifacts, not real enzymes.
        max_length:  Discard sequences longer than this. Enzymes >2000 aa are
                     rare and often multi-domain complexes that behave poorly
                     as discrete industrial enzymes.

    Returns:
        List of dicts with keys:
          id, accession, description, sequence, length, organism, source_biome
    """
    if not Entrez.email:
        raise RuntimeError("Call configure_entrez() before searching NCBI.")

    logger.info(f"Searching NCBI protein: {query}")

    # esearch returns a list of GI numbers / accessions matching the query.
    # usehistory='y' stores results server-side for more efficient efetch.
    try:
        search_handle = Entrez.esearch(
            db=database,
            term=query,
            retmax=max_results,
            usehistory="y",
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
    except Exception as e:
        logger.error(f"NCBI esearch failed for '{query}': {e}")
        return []

    id_list = search_results.get("IdList", [])
    if not id_list:
        logger.warning(f"No results for query: {query}")
        return []

    logger.info(f"Found {len(id_list)} hits for '{query}'")

    # efetch retrieves actual sequence records in GenBank format.
    # GenBank contains the sequence plus rich metadata (organism, description,
    # taxonomy, isolation source). We parse with BioPython SeqIO.
    try:
        time.sleep(0.4)  # Respect NCBI rate limits (3 req/s without API key)
        fetch_handle = Entrez.efetch(
            db=database,
            id=",".join(id_list),
            rettype="gb",
            retmode="text",
        )
        records_raw = fetch_handle.read()
        fetch_handle.close()
    except Exception as e:
        logger.error(f"NCBI efetch failed: {e}")
        return []

    # Parse GenBank text into SeqRecord objects, then into plain dicts.
    records = []
    for record in SeqIO.parse(StringIO(records_raw), "genbank"):
        seq_str = str(record.seq)

        # Filter by length
        if len(seq_str) < min_length or len(seq_str) > max_length:
            continue

        # Filter out low-quality sequences with too many ambiguous residues
        if len(seq_str) > 0 and seq_str.count("X") / len(seq_str) > 0.1:
            continue

        # Extract organism and isolation source from feature table
        organism = "unknown"
        source_biome = "unknown"
        for feature in record.features:
            if feature.type == "source":
                organism = feature.qualifiers.get("organism", ["unknown"])[0]
                env_sample = feature.qualifiers.get("environmental_sample", [])
                isolation_source = feature.qualifiers.get("isolation_source", [""])[0]
                if env_sample or "metagenome" in organism.lower():
                    source_biome = isolation_source or "environmental metagenome"
                else:
                    source_biome = organism

        records.append({
            "id":           record.id,
            "accession":    record.name,
            "description":  record.description,
            "sequence":     seq_str,
            "length":       len(seq_str),
            "organism":     organism,
            "source_biome": source_biome,
        })

    logger.info(f"Parsed {len(records)} valid records")
    return records


def search_all_biomes(max_per_biome: int = 20) -> list:
    """
    Convenience wrapper: run all EXTREME_BIOME_QUERIES and merge results.
    Deduplicates by accession so overlapping queries don't double-count.
    The 0.5s sleep between queries satisfies NCBI TOS rate limits.
    """
    all_records = {}
    for query in EXTREME_BIOME_QUERIES:
        batch = search_ncbi_proteins(query, max_results=max_per_biome)
        for rec in batch:
            acc = rec["accession"]
            if acc not in all_records:
                all_records[acc] = rec
        time.sleep(0.5)

    result = list(all_records.values())
    logger.info(f"Total unique NCBI records after merging all biomes: {len(result)}")
    return result


def blast_against_enzyme_families(
    sequence: str,
    top_n_families: int = 3,
) -> list:
    """
    Run BLASTp on a query sequence against the NCBI nr database.

    Uses NCBI's remote BLAST service (NCBIWWW), which is free but slow
    (typically 20-60 seconds per query). Only run this on top IAS candidates.

    Returns hits sorted by percent identity (descending):
      [{family, hit_title, identity_pct, e_value, alignment_length, score}, ...]

    Why percent identity and not e-value?
    E-value measures statistical significance; identity measures functional
    similarity. An enzyme with 70% identity to a known laccase is very likely
    a laccase. 30% identity might be a distant homolog with different specificity.
    We report both; the IAS scorer uses identity as its BLAST signal.
    """
    logger.info(f"Running BLAST for sequence of length {len(sequence)}")
    hits = []

    try:
        # BLOSUM62 is the standard substitution matrix for protein homology.
        # hitlist_size=5 keeps results manageable.
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="nr",
            sequence=sequence,
            hitlist_size=5,
            matrix_name="BLOSUM62",
        )
        blast_records = list(NCBIXML.parse(result_handle))
    except Exception as e:
        logger.warning(f"BLAST failed: {e}. Returning empty hits.")
        return []

    if not blast_records or not blast_records[0].alignments:
        return []

    seen_families = set()
    for alignment in blast_records[0].alignments:
        hsp = alignment.hsps[0]  # Best HSP (high-scoring segment pair)

        # Percent identity: fraction of aligned positions that are identical
        identity_pct = (hsp.identities / hsp.align_length) * 100

        # Match to a known enzyme family by keyword search in the hit title
        title_lower = alignment.title.lower()
        matched_family = "unknown"
        for family_name in ENZYME_FAMILIES:
            if family_name in title_lower:
                matched_family = family_name
                break

        if matched_family in seen_families:
            continue
        seen_families.add(matched_family)

        hits.append({
            "family":           matched_family,
            "hit_title":        alignment.title[:80],
            "identity_pct":     round(identity_pct, 1),
            "e_value":          hsp.expect,
            "alignment_length": hsp.align_length,
            "score":            hsp.score,
        })

        if len(hits) >= top_n_families:
            break

    hits.sort(key=lambda x: x["identity_pct"], reverse=True)
    return hits


def fetch_reference_sequence(accession: str) -> Optional[str]:
    """
    Fetch a single protein sequence by accession number.
    Used to retrieve reference enzyme sequences for comparison panels.
    """
    try:
        time.sleep(0.4)
        handle = Entrez.efetch(
            db="protein",
            id=accession,
            rettype="fasta",
            retmode="text",
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except Exception as e:
        logger.warning(f"Could not fetch {accession}: {e}")
        return None
