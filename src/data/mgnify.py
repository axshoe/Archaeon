"""
src/data/mgnify.py
------------------
MGnify API wrapper for Archaeon.

MGnify (https://www.ebi.ac.uk/metagenomics) is the European Bioinformatics
Institute's metagenomic analysis platform. It processes raw environmental
sequencing data submitted by research groups worldwide and makes the annotated
results freely available via a JSON REST API.

Why MGnify vs NCBI:
  NCBI protein database is curated but skews toward cultured organisms.
  MGnify is exclusively metagenomic: it captures protein-coding sequences
  from uncultured organisms that have never been grown in a lab. Since most
  extremophile diversity is uncultured, MGnify gives access to novel enzyme
  candidates that simply don't exist in NCBI nr.

The MGnify API follows JSON:API specification. Responses have a 'data' key
containing a list of resource objects, each with an 'id', 'type', and
'attributes' dict. Pagination is handled via 'links.next'.

This module handles:
  1. Searching environmental samples by biome (e.g., "root:Environmental:Aquatic:
     Hydrothermal vents")
  2. Downloading protein sequences associated with those samples (from MGnify's
     protein database, which aggregates across all samples)
  3. Filtering by functional annotation when available (EC number, InterPro family)

API docs: https://www.ebi.ac.uk/metagenomics/api/v1/
"""

import time
import logging
from typing import Optional

import requests

logger = logging.getLogger(__name__)

MGNIFY_BASE_URL = "https://www.ebi.ac.uk/metagenomics/api/v1"

# MGnify biome lineage strings for extreme environments.
# These correspond to the controlled vocabulary used by MGnify to classify
# environmental samples. 'root:Environmental:Aquatic:Hydrothermal vents' is
# the literal path in their biome ontology.
EXTREME_BIOMES = {
    "hydrothermal_vent":     "root:Environmental:Aquatic:Hydrothermal vents",
    "hot_spring":            "root:Environmental:Terrestrial:Hot springs",
    "hypersaline":           "root:Environmental:Aquatic:Saline",
    "acid_mine":             "root:Environmental:Terrestrial:Mine",
    "deep_subsurface":       "root:Environmental:Terrestrial:Subsurface",
    "polar":                 "root:Environmental:Aquatic:Polar",
    "volcanic":              "root:Environmental:Terrestrial:Volcanic",
}

# Functional annotation filters: InterPro family IDs for major enzyme classes.
# If a protein sequence has an InterPro annotation matching one of these, it
# is a member of that structural/functional family.
ENZYME_INTERPRO = {
    "lipase":    "IPR013818",  # AB hydrolase superfamily
    "protease":  "IPR001995",  # Peptidase M4 thermolysin-like
    "amylase":   "IPR006048",  # Alpha-amylase, catalytic domain
    "cellulase": "IPR001547",  # Glycoside hydrolase superfamily
    "xylanase":  "IPR000254",  # Xylanase/chitin deacetylase
}


def _get(url: str, params: Optional[dict] = None, retries: int = 3) -> Optional[dict]:
    """
    Thin wrapper around requests.get with retry logic.

    MGnify's API occasionally times out or returns 503 during high load.
    3 retries with exponential backoff handles transient failures gracefully.
    """
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            elif resp.status_code == 429:
                # Rate limited: back off for 10 seconds
                logger.warning("MGnify rate limit hit; sleeping 10s")
                time.sleep(10)
            else:
                logger.warning(f"MGnify HTTP {resp.status_code} for {url}")
                return None
        except requests.RequestException as e:
            wait = 2 ** attempt
            logger.warning(f"MGnify request failed (attempt {attempt+1}): {e}. Retrying in {wait}s")
            time.sleep(wait)
    return None


def search_samples_by_biome(biome_key: str, max_samples: int = 20) -> list:
    """
    Retrieve environmental sample metadata for a given extreme biome.

    Args:
        biome_key:   Key from EXTREME_BIOMES dict
        max_samples: Max samples to retrieve (pagination stops here)

    Returns:
        List of dicts: {sample_id, biome_lineage, description, geo_loc}
    """
    biome_lineage = EXTREME_BIOMES.get(biome_key)
    if not biome_lineage:
        logger.error(f"Unknown biome key: {biome_key}. Options: {list(EXTREME_BIOMES.keys())}")
        return []

    logger.info(f"Searching MGnify samples for biome: {biome_lineage}")

    url = f"{MGNIFY_BASE_URL}/samples"
    params = {
        "lineage": biome_lineage,
        "page_size": min(max_samples, 100),
        "ordering": "-last-update",
    }

    data = _get(url, params)
    if not data:
        return []

    samples = []
    for item in data.get("data", []):
        attrs = item.get("attributes", {})
        samples.append({
            "sample_id":      item["id"],
            "biome_lineage":  attrs.get("biome-lineage", biome_lineage),
            "description":    attrs.get("sample-name", ""),
            "geo_loc":        attrs.get("geo-loc-name", "unknown"),
            "collection_date": attrs.get("collection-date", ""),
        })

    logger.info(f"Retrieved {len(samples)} samples for biome '{biome_key}'")
    return samples


def search_proteins_by_biome(
    biome_key: str,
    max_proteins: int = 50,
    min_length: int = 100,
    max_length: int = 2000,
) -> list:
    """
    Search MGnify's protein database for sequences from a given extreme biome.

    MGnify maintains a searchable protein database (MGnify Proteins) separate
    from its sample database. This endpoint returns protein sequences with
    their associated metadata including GO terms, InterPro annotations, and
    biome of origin.

    Args:
        biome_key:    Key from EXTREME_BIOMES dict
        max_proteins: Maximum proteins to retrieve
        min_length:   Minimum sequence length filter
        max_length:   Maximum sequence length filter

    Returns:
        List of dicts compatible with the NCBI record format so the scorer
        can treat both data sources uniformly.
    """
    biome_lineage = EXTREME_BIOMES.get(biome_key)
    if not biome_lineage:
        logger.error(f"Unknown biome: {biome_key}")
        return []

    logger.info(f"Searching MGnify proteins for biome: {biome_key}")

    url = f"{MGNIFY_BASE_URL}/proteins"
    params = {
        "biome_name": biome_lineage,
        "page_size": min(max_proteins, 100),
    }

    data = _get(url, params)
    if not data or "data" not in data:
        logger.warning(f"No MGnify protein data for biome: {biome_key}")
        return []

    records = []
    for item in data.get("data", []):
        attrs = item.get("attributes", {})
        sequence = attrs.get("sequence", "")

        if not sequence:
            continue
        if len(sequence) < min_length or len(sequence) > max_length:
            continue
        if sequence.count("X") / max(len(sequence), 1) > 0.1:
            continue

        # Extract functional annotation if available
        go_terms = attrs.get("go-slim", [])
        interpro_ids = [
            ann.get("accession", "")
            for ann in attrs.get("interpro-entries", [])
        ]

        # Map InterPro IDs to human-readable enzyme family names
        enzyme_family = "unknown"
        for family, ipr_id in ENZYME_INTERPRO.items():
            if ipr_id in interpro_ids:
                enzyme_family = family
                break

        records.append({
            "id":               item["id"],
            "accession":        item["id"],
            "description":      attrs.get("function-calls", enzyme_family),
            "sequence":         sequence,
            "length":           len(sequence),
            "organism":         "environmental metagenome",
            "source_biome":     biome_key,
            "source":           "MGnify",
            "go_terms":         go_terms,
            "interpro_ids":     interpro_ids,
            "enzyme_family":    enzyme_family,
        })

    logger.info(f"Parsed {len(records)} MGnify protein records for biome '{biome_key}'")
    return records


def search_all_extreme_biomes(max_per_biome: int = 25) -> list:
    """
    Convenience wrapper: search all extreme biomes and merge results.
    Deduplicates by MGnify protein ID.
    """
    all_records = {}
    for biome_key in EXTREME_BIOMES:
        batch = search_proteins_by_biome(biome_key, max_proteins=max_per_biome)
        for rec in batch:
            pid = rec["id"]
            if pid not in all_records:
                all_records[pid] = rec
        time.sleep(0.3)  # Be polite to EBI's servers

    result = list(all_records.values())
    logger.info(f"Total unique MGnify records across all biomes: {len(result)}")
    return result


def get_sample_runs(sample_id: str) -> list:
    """
    Retrieve analysis runs associated with a given MGnify sample ID.

    A 'run' in MGnify corresponds to a single sequencing run associated with
    the sample. Each run has its own assembly, annotation, and statistics.
    This is useful if you want to trace a candidate sequence back to its
    original raw sequencing run for provenance.

    Returns list of run accessions (strings).
    """
    url = f"{MGNIFY_BASE_URL}/samples/{sample_id}/runs"
    data = _get(url)
    if not data:
        return []
    return [item["id"] for item in data.get("data", [])]
