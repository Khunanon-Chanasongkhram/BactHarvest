#!/usr/bin/env python3
"""
MicroScope.py v1.1.0
=====================
A precision text-mining tool for extracting microbial phenotypes from PubMed.

Dependencies (auto-installed on first run):
-------------------------------------------
- biopython (Bio.Entrez)
- pandas
- tqdm
- certifi

Usage Examples:
---------------
# 1. Standard Run
python MicroScope.py --phenotype "nitrogen fixation" --email "your@email.com" --min_score 5

# 2. High Strictness
python MicroScope.py --phenotype "plastic degradation" --email "your@email.com" --min_score 8 --max_articles 500

# 3. Fast Run with API key and JSON output
python MicroScope.py --phenotype "phosphate solubilization" --email "your@email.com" \\
    --api_key "YOUR_KEY" --skip_genome --output_format both

# 4. Resume an interrupted run
python MicroScope.py --phenotype "nitrogen fixation" --email "your@email.com" \\
    --resume --cache_file nitrogen_cache.json

# 5. Find organisms wet-lab-proven NOT to have the phenotype (negative controls)
python MicroScope.py --phenotype "nitrogen fixation" --email "your@email.com" \\
    --negative_evidence

# 6. Restrict to a publication-year window
python MicroScope.py --phenotype "siderophore" --email "your@email.com" \\
    --year_from 2015 --year_to 2024

# 7. Look beyond the built-in genus whitelist (every new genus is checked
#    against NCBI Taxonomy, so --verify_taxonomy is required)
python MicroScope.py --phenotype "siderophore" --email "your@email.com" \\
    --discover_genera --verify_taxonomy

Notes:
------
- PubMed results are fetched newest-first (sort=pub_date), so --max_articles
  prioritizes recent literature before older papers. Requests larger than one
  NCBI page (10,000 PMIDs) are paged automatically via the Entrez history
  server, so --max_articles is not capped at 10,000.
- Evidence keywords are matched by word stem, so nominalised and British
  spellings score the same as the base verb: "IAA production", "phosphate
  solubilisation" and "biosynthesis of siderophores" all count as strong
  evidence, where earlier versions only credited "produce"/"solubilize".
- Organism recognition is normally limited to the curated BACTERIAL_GENERA
  whitelist. --discover_genera lifts that limit and accepts any plausible
  binomial instead, but each newly seen name must then resolve inside the
  Bacteria/Archaea subtree of NCBI Taxonomy to be kept.
- A species epithet must be a single alphabetic Latin word, so prose that
  merely follows a genus name ("Pseudomonas and", "Salmonella species",
  "Rhizobium medium", "Burkholderia two-component") is no longer reported as
  an organism. Strain suffixes are trimmed instead of discarded, so
  "Rhizobium meliloti-1021" is recorded as "Rhizobium meliloti".
- performance_method reports the experimental assay/method used to prove the
  phenotype — either matched against the curated per-phenotype vocabulary, or
  (when nothing in that vocabulary matches) extracted directly from generic
  reporting phrases like "using the X assay", tagged "(Reported)".
- Negation words ("undetected", "unable to", "no activity", etc.) only count
  as evidence about YOUR phenotype if they share a clause with an actual
  phenotype keyword mention in the sentence — a negation about an unrelated
  trait elsewhere in the same sentence no longer zeroes real positive evidence
  (positive mode) or gets misread as proof of a negative result (--negative_evidence).
- genome_source tells you where genome_accession came from: "Reported in
  Paper" (the article itself states a deposition accession for that organism)
  vs "NCBI Lookup" (no accession was found in the text, so MicroScope searched
  NCBI Assembly by organism name instead — this finds *some* genome for the
  species, not necessarily the exact strain the paper studied).
- All NCBI network calls (PubMed search/fetch, full-text, taxonomy, genome
  lookup) automatically retry on transient failures — dropped connections,
  incomplete/truncated responses, timeouts — with exponential backoff (up to
  3 attempts) before giving up on that request. Watch for "retrying in Ns"
  WARNING log lines during a run; if you instead see the request's own ERROR
  line, all retries were exhausted and that item was skipped.

Output:
-------
Saves results (CSV and/or JSON) containing:
PMID | pub_year | Organism | Confidence (Very High/High/Medium/Low) | Evidence Score
performance_value (float) | performance_unit | performance_range | performance_metric (raw)
performance_method | performance_evidence_sentence | Genome ID

Programmer: Khunanon Chanasongkhram
"""

from __future__ import annotations

# ==========================================
# IMPORTS
# ==========================================
import sys
import re
import argparse
import time
import logging
import json
import dataclasses
import urllib.error
import http.client
import xml.etree.ElementTree as ET
from math import floor, log10
from pathlib import Path
from typing import Optional

__version__ = "0.8.0"

# ==========================================
# LOGGING
# ==========================================

logger = logging.getLogger("microscope")


def setup_logging(level_str: str = "INFO") -> None:
    level = getattr(logging, level_str.upper(), logging.INFO)
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        level=level,
    )
    logger.setLevel(level)


# ==========================================
# VOCABULARY
# ==========================================

METHOD_CATEGORIES: dict = {
    "IAA": {
        "salkowski": "Salkowski Reagent Assay",
        "colorimetric": "Colorimetric Assay",
        "tryptophan": "L-Tryptophan Supplemented Assay",
        "hplc": "HPLC",
        "uplc": "UPLC",
        "lc-ms": "LC-MS",
        "gc-ms": "GC-MS",
        "tlc": "TLC",
    },
    "NITROGEN": {
        "acetylene reduction": "Acetylene Reduction Assay (ARA)",
        "ara": "Acetylene Reduction Assay (ARA)",
        "ethylene": "Acetylene Reduction Assay (ARA)",
        "nitrogenase": "Nitrogenase Assay",
        "gc": "Gas Chromatography",
    },
    "PHOSPHATE": {
        "molybdenum": "Molybdenum Blue Method",
        "pikovoskaya": "Pikovoskaya (PVK) Agar",
        "pvk": "Pikovoskaya (PVK) Agar",
        "nbrip": "NBRIP Medium",
        "clear zone": "Halo Zone Measurement",
        "halo": "Halo Zone Measurement",
        "solubilization index": "Solubilization Index (SI)",
        "hplc": "HPLC",
    },
    "SIDEROPHORE": {
        "chrome azurol": "CAS Assay",
        "cas": "CAS Assay",
        "siderophore unit": "CAS Assay (SU)",
        "orange halo": "CAS Agar Halo",
        "spectrophotomet": "Spectrophotometry",
    },
    "ACC_DEAMINASE": {
        "acc deaminase": "ACC Deaminase Assay",
        "df medium": "DF Salt Medium",
        "df-medium": "DF Salt Medium",
        "acrylonitrile": "ACC Degradation Assay",
        "ketobutyrate": "α-Ketobutyrate Assay",
    },
    "BIOFILM": {
        "crystal violet": "Crystal Violet Staining",
        "microtiter plate": "Microtiter Plate Assay",
        "biofilm index": "Biofilm Index",
        "od595": "OD595 Measurement",
        "od570": "OD570 Measurement",
        "confocal": "Confocal Microscopy",
    },
    "HYDROGEN": {
        "gas chromatography": "Gas Chromatography",
        "nitrogenase": "Nitrogenase Activity Assay",
        "hydrogenase": "Hydrogenase Assay",
        "methylene blue": "Methylene Blue Reduction",
        "sodium dithionite": "Dithionite Reduction Assay",
    },
}

ALL_METHODS: dict = {}
for _cat in METHOD_CATEGORIES:
    ALL_METHODS.update(METHOD_CATEGORIES[_cat])
ALL_METHODS.update({
    "spectrophotomet": "Spectrophotometry",
    "optical density": "OD Measurement",
    "od600": "OD600",
})

PATTERNS_DB: dict = {
    "NITROGEN": [
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:nmol|µmol|umol)\s*(?:C2H4|ethylene|C2H2)?(?:[\w\/\-\s]*)(?:h|hr|min|vial|protein)",
        r"(\d+(?:\.\d+)?)\s*(?:nmol|µmol|umol)\s*(?:C2H4|ethylene)?",
        # ethylene/acetylene per vial or flask
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:nmol|µmol)\s*(?:c2h4|ethylene|acetylene)\s*(?:per\s+)?(?:vial|tube|flask|culture)[-1\s]*h",
        # µg / mg N fixed per g or per plant
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:µg|mg|ng)\s*n\s*(?:fixed\s+)?(?:\/|per)\s*(?:g|plant|ml)",
    ],
    "IAA": [
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:µg|ug|mg)\s*(?:IAA|auxin)?\s*[\/·\s]?\s*(?:mL|ml|L|l)(?:[-1\s]*)\s*(?:IAA|auxin)?",
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*ppm",
        # ng/mL (HPLC precision) and explicit slash notation
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:ng|µg|ug)\s*/\s*(?:ml|l)\b",
    ],
    "PHOSPHATE": [
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:µg|ug|mg|g)\s*[\/·\s]?\s*(?:mL|ml|L|l|kg|g)(?:[-1\s]*)?",
        r"(?:SI|PSI|index)\s*(?:=|of|:)\s*(\d+(?:\.\d+)?(?:\s*[-–]\s*\d+(?:\.\d+)?)?)",
        r"(\d+(?:\.\d+)?(?:\s*[-–]\s*\d+(?:\.\d+)?)?)\s*(?:mm|cm)\s*(?:halo|zone|diameter)",
        # "mg P/L" or "mg P L-1" — phosphorus denominator between unit and volume
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:µg|ug|mg)\s+p\s*(?:\/\s*l|l[-1\s]|per\s+l)",
        # percentage solubilization (e.g., "72% phosphate solubilization")
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*%\s*(?:phosphate\s+|phosphorus\s+)?solubili",
        # ppm (parts per million) — common unit in phosphate solubilization assays
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*ppm\b",
    ],
    "SIDEROPHORE": [
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*%\s*(?:SU|units|efficiency|reduction)",
        r"(\d+(?:\.\d+)?(?:\s*[-–]\s*\d+(?:\.\d+)?)?)\s*(?:mm|cm)\s*(?:halo|zone|diameter)",
        # percentage CAS decolorization (common CAS shuttle assay reporting)
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*%\s*(?:cas\s+)?(?:decolori[sz]ation|decolouri[sz]ation)",
        # explicit siderophore units (SU)
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:su|siderophore\s+units?)\b",
    ],
    "ACC_DEAMINASE": [
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*nmol\s*(?:α-ketobutyrate|ketobutyrate)\s*(?:mg|per\s*mg)(?:[-\w\/\s]*)(?:h|hr)",
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:nmol|µmol)\s*(?:mg|protein)[-1\s]*h[-1\s]?",
        # slash notation: nmol/mg/h or nmol/mg protein/h (substrate name not required)
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:nmol|µmol)\s*/\s*(?:mg\s*(?:protein\s*)?\/\s*h|mg\s*h)",
        # α-KB abbreviation common in recent literature
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:nmol|µmol)\s*(?:α-?kb|alpha-?kb|α-?ketobutyrate)\s*",
    ],
    "BIOFILM": [
        r"OD\s*(?:595|570)\s*(?:=|of|:)?\s*(\d+(?:\.\d+)?)",
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*%\s*(?:biofilm|inhibition|reduction)",
        # "absorbance at 570/595" phrasing seen in microtiter assay descriptions
        r"(?:absorbance|od)\s*(?:at\s+)?(?:595|570)\s*(?:nm)?\s*(?:of|=|:)?\s*(\d+(?:\.\d+)?)",
        # fold-increase in biofilm
        r"(\d+(?:\.\d+)?)\s*[-–]?\s*fold\s+(?:higher|greater|more)\s+biofilm",
    ],
    "HYDROGEN": [
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:ml|mL)\s*(?:H2|hydrogen)\s*(?:\/|per)\s*(?:h|hr|hour)",
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*µmol\s*(?:H2|hydrogen)\s*(?:mg|per\s*mg)(?:[\w\/\s]*)(?:h|hr)",
        # mmol H2/L/h or per flask/culture
        r"(\d+(?:\.\d+)?(?:\s*[-–±]\s*\d+(?:\.\d+)?)?)\s*(?:mmol|µmol|nmol)\s*(?:h2|hydrogen)\s*(?:\/|per)\s*(?:l|liter|litre|culture|flask)",
    ],
}

INTERFERENCE_KEYWORDS: dict = {
    "IAA": ["phosphate", "phosphorus", "potassium", "zinc", "nitrogen", "siderophore", "ammonia"],
    "NITROGEN": ["phosphate", "iaa", "auxin", "potassium", "siderophore"],
    "PHOSPHATE": ["nitrogen", "iaa", "auxin", "potassium", "zinc"],
    "SIDEROPHORE": ["phosphate", "iaa", "nitrogen"],
    "ACC_DEAMINASE": ["phosphate", "iaa", "nitrogen", "siderophore"],
    "BIOFILM": ["phosphate", "iaa", "nitrogen"],
    "HYDROGEN": ["phosphate", "iaa", "nitrogen", "siderophore"],
}

# --- Evidence lexicons -------------------------------------------------------
#
# Entries are word *stems* matched with a leading \b, not whole words, because
# this literature overwhelmingly reports phenotypes as nominalisations rather
# than as finite verbs: "IAA production", "phosphate solubilisation" and
# "siderophore biosynthesis" are the normal phrasing, while "produces IAA" is
# the exception. Whole-word matching credited only the latter, so the strongest
# sentences in a typical abstract scored no action-verb points at all. British
# spellings (-ise/-isation) are covered by the [sz] alternations for the same
# reason — a large share of this corpus is published in European journals.
#
# Each entry is a regex fragment. Fragments that could prefix an unrelated word
# pin their own right edge (e.g. "potent" must not match "potential", which is
# speculative language and the opposite of a standout result).

STRONG_VERB_STEMS = [
    r"produc",              # produce/produced/production/producer
    r"(?:bio)?synthesi",    # synthesize/synthesis/biosynthesis
    r"degrad",              # degrade/degradation
    r"solubili[sz]",        # solubilize/solubilisation
    r"fix",                 # fix/fixed/fixation
    r"accumulat",
    r"yield",
    r"secret",              # secrete/secretion
    r"excret",
    r"minerali[sz]",
    r"mobili[sz]",          # phosphate/iron mobilisation
    r"chelat",              # siderophore chelation
    r"cataly",              # catalyze/catalysis/catalytic
]

# Compiled individually: analyze_abstract counts how many *distinct* action
# stems a sentence contains and awards a bonus for more than one.
STRONG_VERB_RES = [re.compile(r"\b" + s, re.IGNORECASE) for s in STRONG_VERB_STEMS]

WEAK_CONTEXTS = {
    "hypothesized", "speculated", "suggests", "future studies", "aim of",
    "goal of", "previously", "review", "unknown",
}

NEGATION_PATTERNS = [
    r"\bnot\s+(?:able\s+to\s+)?(?:produc|degrad|fix|solubili[sz]|synthesi[sz]|secret)\w*\b",
    r"\b(?:did|does|do|was|were|could|can)\s+not\s+\w+",
    r"\bfail(?:ed|s|ing)?\s+to\b",
    r"\bunable\s+to\b",
    r"\bincapable\s+of\b",
    # "no siderophore production", "no detectable nitrogenase activity" — the
    # negated noun is usually separated from "no" by one or two modifiers, so
    # requiring strict adjacency (the previous "\bno (activity|production)\b")
    # missed most real negative statements.
    r"\bno\s+(?:\w+\s+){0,2}(?:activit|production|producer|solubili|fixation|formation|synthesis|secretion|halo|zone)\w*\b",
    r"\bnot\s+detect(?:ed|able)\b",
    r"\bnot\s+observed\b",
    r"\bundetect(?:ed|able)\b",
    r"\bnegligible\b",
    r"\babsence\s+of\b",
    r"\bnegative\s+for\b",
    r"\black(?:ed|s|ing)?\b",
    r"\bnon-?(?:producing|producers?|diazotrophic)\b",
]
NEGATION_RES = [re.compile(p, re.IGNORECASE) for p in NEGATION_PATTERNS]

# Stems confirming an observation was actually measured/reported (+2 per sentence)
CONFIRMATION_STEMS = [
    r"confirm", r"demonstrat", r"show", r"exhibit", r"reveal", r"observ",
    r"report", r"detect", r"isolat", r"identif", r"characteri[sz]",
    r"measur", r"quantif", r"determin",
]
CONFIRMATION_RE = re.compile(r"\b(?:" + "|".join(CONFIRMATION_STEMS) + r")", re.IGNORECASE)

# Stems indicating a standout or comparative result (+1 per sentence)
QUANTITATIVE_STEMS = [
    r"maxim", r"highest", r"optim", r"significan", r"substantial",
    r"considerabl", r"superior", r"best\b", r"greatest", r"efficien",
    r"potent(?:ly)?\b",   # not "potential" — that is speculative, not a result
    r"remarkabl", r"notabl",
]
QUANTITATIVE_RE = re.compile(r"\b(?:" + "|".join(QUANTITATIVE_STEMS) + r")", re.IGNORECASE)

# ==========================================
# PHENOTYPE REGISTRY
# ==========================================


@dataclasses.dataclass
class PhenotypeConfig:
    display_name: str
    keywords: list
    methods: dict
    patterns: list
    interference: list


PHENOTYPE_REGISTRY: dict = {
    "IAA": PhenotypeConfig(
        display_name="IAA / Auxin Production",
        keywords=["iaa", "auxin", "indole"],
        methods=METHOD_CATEGORIES["IAA"],
        patterns=PATTERNS_DB["IAA"],
        interference=INTERFERENCE_KEYWORDS["IAA"],
    ),
    "NITROGEN": PhenotypeConfig(
        display_name="Nitrogen Fixation",
        keywords=["nitrogen", "fixation", "diazotroph"],
        methods=METHOD_CATEGORIES["NITROGEN"],
        patterns=PATTERNS_DB["NITROGEN"],
        interference=INTERFERENCE_KEYWORDS["NITROGEN"],
    ),
    "PHOSPHATE": PhenotypeConfig(
        display_name="Phosphate Solubilization",
        keywords=["phosphate", "solubil", "phosphorus"],
        methods=METHOD_CATEGORIES["PHOSPHATE"],
        patterns=PATTERNS_DB["PHOSPHATE"],
        interference=INTERFERENCE_KEYWORDS["PHOSPHATE"],
    ),
    "SIDEROPHORE": PhenotypeConfig(
        display_name="Siderophore Production",
        keywords=["siderophore", "iron"],
        methods=METHOD_CATEGORIES["SIDEROPHORE"],
        patterns=PATTERNS_DB["SIDEROPHORE"],
        interference=INTERFERENCE_KEYWORDS["SIDEROPHORE"],
    ),
    "ACC_DEAMINASE": PhenotypeConfig(
        display_name="ACC Deaminase Activity",
        keywords=["acc deaminase", "acc", "1-aminocyclopropane"],
        methods=METHOD_CATEGORIES["ACC_DEAMINASE"],
        patterns=PATTERNS_DB["ACC_DEAMINASE"],
        interference=INTERFERENCE_KEYWORDS["ACC_DEAMINASE"],
    ),
    "BIOFILM": PhenotypeConfig(
        display_name="Biofilm Formation",
        keywords=["biofilm"],
        methods=METHOD_CATEGORIES["BIOFILM"],
        patterns=PATTERNS_DB["BIOFILM"],
        interference=INTERFERENCE_KEYWORDS["BIOFILM"],
    ),
    "HYDROGEN": PhenotypeConfig(
        display_name="Hydrogen Production",
        keywords=["hydrogen", "h2 production", "biohydrogen"],
        methods=METHOD_CATEGORIES["HYDROGEN"],
        patterns=PATTERNS_DB["HYDROGEN"],
        interference=INTERFERENCE_KEYWORDS["HYDROGEN"],
    ),
}


def detect_phenotype_key(phenotype_query: str) -> Optional[str]:
    p = phenotype_query.lower()
    for key, cfg in PHENOTYPE_REGISTRY.items():
        if any(kw in p for kw in cfg.keywords):
            return key
    return None


def get_phenotype_config(
    phenotype_query: str, forced_key: Optional[str] = None
) -> Optional[PhenotypeConfig]:
    key = forced_key or detect_phenotype_key(phenotype_query)
    return PHENOTYPE_REGISTRY.get(key) if key else None


def get_allowed_methods(phenotype: str, forced_key: Optional[str] = None) -> dict:
    cfg = get_phenotype_config(phenotype, forced_key)
    return cfg.methods if cfg else ALL_METHODS


def get_patterns_for_phenotype(phenotype: str, forced_key: Optional[str] = None) -> list:
    cfg = get_phenotype_config(phenotype, forced_key)
    if cfg:
        return cfg.patterns
    all_pats: list = []
    for pats in PATTERNS_DB.values():
        all_pats.extend(pats)
    return all_pats


def get_interference_words(phenotype: str, forced_key: Optional[str] = None) -> list:
    cfg = get_phenotype_config(phenotype, forced_key)
    return cfg.interference if cfg else []


# ==========================================
# ORGANISM LISTS
# ==========================================

BACTERIAL_GENERA = {
    # ── Core / broad ──────────────────────────────────────────────────────────
    "Acetobacter", "Acinetobacter", "Actinomyces", "Agrobacterium", "Alcaligenes", "Arthrobacter",
    "Azospirillum", "Azotobacter", "Bacillus", "Bacteroides", "Burkholderia", "Campylobacter",
    "Citrobacter", "Clostridium", "Corynebacterium", "Enterobacter", "Enterococcus", "Erwinia",
    "Escherichia", "Klebsiella", "Lactobacillus", "Methylobacterium", "Micrococcus", "Mycobacterium",
    "Nitrobacter", "Nitrosomonas", "Paenibacillus", "Pantoea", "Pseudomonas", "Rhizobium",
    "Rhodococcus", "Salmonella", "Serratia", "Sinorhizobium", "Sphingomonas", "Staphylococcus",
    "Stenotrophomonas", "Streptococcus", "Streptomyces", "Thiobacillus", "Vibrio", "Xanthomonas",
    "Bradyrhizobium", "Mesorhizobium", "Gluconacetobacter", "Herbaspirillum", "Azospira", "Frankia",
    "Geobacter", "Shewanella", "Ralstonia", "Cupriavidus", "Comamonas", "Variovorax", "Sphingobium",
    "Ideonella", "Azoarcus", "Paraburkholderia", "Phyllobacterium", "Rahnella", "Raoultella",
    "Shinella", "Ochrobactrum", "Flavobacterium", "Microbacterium", "Lysobacter", "Novosphingobium",
    "Nocardia", "Cellulomonas", "Brevibacillus", "Lysinibacillus", "Priestia", "Nostoc", "Anabaena",
    # ── Rhizobia / N₂-fixing symbioses ───────────────────────────────────────
    "Neorhizobium",       # N₂ fixation in Galega, Hedysarum
    "Allorhizobium",      # N₂ fixation; recently split from Rhizobium
    "Ensifer",            # Sinorhizobium synonym; widely used in legume literature
    "Devosia",            # microsymbiont of Neptunia
    "Aminobacter",        # N₂-fixing endophyte
    "Methyloversatilis",  # methylotrophic, some N₂ fixation
    # ── Phosphate solubilization / common rhizosphere ─────────────────────────
    "Chryseobacterium",   # very common PSB in literature; ACC deaminase
    "Janthinobacterium",  # violacein-producing rhizosphere PSB
    "Massilia",           # PSB isolated from bulk and rhizosphere soil
    "Collimonas",         # fungal-feeding soil bacterium; some PSB activity
    "Brevundimonas",      # PSB; isolated from various rhizosphere soils
    "Acidovorax",         # rhizosphere PSB; plant-associated
    "Pedobacter",         # soil PSB; cold-adapted rhizosphere bacterium
    "Luteimonas",         # soil PSB; siderophore producer
    "Sphingopyxis",       # rhizosphere bacterium; some PSB reports
    "Leifsonia",          # endophytic PSB; IAA and siderophore producer
    "Curtobacterium",     # endophytic PSB; common in seed microbiome
    # ── IAA / ACC deaminase producers ────────────────────────────────────────
    "Kosakonia",          # reclassified from Enterobacter; key sugarcane PGPB
    "Gluconobacter",      # acetic acid bacterium; IAA and GA producer
    "Rhodobacter",        # phototrophic; IAA and cytokinin producer
    "Rhodospirillum",     # phototrophic α-proteobacterium; IAA and N₂ fixation
    "Rhodopseudomonas",   # PGPB; IAA, N₂ fixation, phosphate solubilization
    "Halomonas",          # halotolerant PGPB; ACC deaminase and IAA under stress
    "Exiguobacterium",    # extremophile; IAA and phosphate solubilization
    "Planococcus",        # soil bacterium; some IAA and PSB reports
    "Glutamicibacter",    # reclassified from Arthrobacter; IAA and PSB
    "Paenarthrobacter",   # reclassified from Arthrobacter; rhizosphere PGPB
    # ── Siderophore producers ────────────────────────────────────────────────
    "Luteibacter",        # siderophore producer; rhizosphere fluorescent bacterium
    "Cohnella",           # endospore-forming soil bacterium; siderophore producer
    # ── Actinomycetes (non-Streptomyces) PGPB ────────────────────────────────
    "Micromonospora",     # N₂ fixation in non-legume nodules (Casuarina, Alnus)
    "Amycolatopsis",      # soil actinomycete; phosphate and siderophore activity
    "Kitasatospora",      # antifungal; some PSB and IAA reports
    "Pseudonocardia",     # actinomycete endophyte; some PGPB traits
    # ── Cyanobacteria (free-living N₂ fixers) ────────────────────────────────
    "Fischerella",        # heterocystous; soil and rice-paddy N₂ fixation
    "Calothrix",          # heterocystous; biofertilizer cyanobacterium
    "Tolypothrix",        # heterocystous; rice-field N₂ fixer
    "Cylindrospermum",    # heterocystous; N₂ fixation and IAA production
    "Scytonema",          # heterocystous; drought-tolerant soil cyanobacterium
    # ── Recently reclassified from Burkholderia ───────────────────────────────
    "Caballeronia",       # N₂ fixation; plant-associated (formerly Burkholderia)
    "Trinickia",          # nodulates Mimosa; N₂ fixation (formerly Burkholderia)
}
INVALID_SPECIES_NAMES = {
    "sp", "spp", "strain", "strains", "isolate", "isolates",
    "mutant", "group", "bacterium", "cells", "gene",
    # Citation/reference abbreviations — e.g. "(A. et al., 2020)" reads as a
    # single-letter genus abbreviation followed by "et", which would otherwise
    # resolve against an earlier same-initial genus mention and fabricate a
    # bogus organism like "Azotobacter et".
    "et", "al", "vs", "etc",
}
CLONING_HOSTS = {
    "Escherichia coli", "E. coli", "Saccharomyces cerevisiae",
    "S. cerevisiae", "Bacillus subtilis", "B. subtilis",
}
PLANT_GENERA = {"Arabidopsis", "Oryza", "Zea", "Triticum", "Solanum", "Nicotiana"}

# Capitalised words that routinely start a sentence (or a title) and would
# otherwise be read as a genus by the "Genus species" regex once the
# BACTERIAL_GENERA whitelist is lifted by --discover_genera. Also used for the
# species slot, where the same words appear as ordinary lowercase text.
# NOTE: only add words here that cannot also be a Latin species epithet.
# "algae" belongs to Shewanella algae, "plantarum" to Lactobacillus plantarum —
# entries like those would silently delete real organisms from every run.
COMMON_WORD_STOPLIST = frozenset({
    "about", "after", "also", "although", "among", "analysis", "and", "another",
    "are", "assay", "based", "because", "before", "being", "below", "between",
    "both", "can", "cells", "compared", "conclusion", "control", "cultured",
    "cultures", "data", "did", "different", "does", "during", "each", "effect",
    "effects", "either", "exhibited", "experiment", "experiments", "figure",
    "finally", "for", "from", "further", "grown", "growth", "has", "have",
    "here", "higher", "highest", "however", "into", "isolate", "isolated",
    "isolates", "isolation", "its", "level", "levels", "material", "materials",
    "may", "medium", "method", "methods", "microbial", "moreover", "most",
    "much", "nitrogen", "none", "not", "novel", "only", "other", "over",
    "overall", "phosphate", "plants", "present", "previous", "produced",
    "production", "recent", "recently", "result", "results", "roots", "sample",
    "samples", "several", "showed", "significant", "significantly", "since",
    "soils", "some", "species", "strain", "strains", "studies", "study", "such",
    "supplementary", "table", "that", "the", "their", "then", "there", "these",
    "they", "this", "those", "three", "through", "thus", "total", "treatment",
    "treatments", "under", "using", "value", "values", "various", "was", "were",
    "when", "where", "which", "while", "will", "with", "within", "without",
    "would",
})

# ==========================================
# HELPERS
# ==========================================


class RateLimiter:
    """Sleeps only the remaining time needed between API calls."""

    def __init__(self, min_interval: float) -> None:
        self.min_interval = min_interval
        self._last_call: float = 0.0

    def wait(self) -> None:
        elapsed = time.monotonic() - self._last_call
        deficit = self.min_interval - elapsed
        if deficit > 0:
            time.sleep(deficit)
        self._last_call = time.monotonic()


class DiskCache:
    """Persists per-PMID results to a JSON file for --resume support.

    Writes are batched: the file is flushed every WRITE_INTERVAL entries
    to avoid rewriting the full JSON on every single PMID.
    """

    WRITE_INTERVAL = 20

    def __init__(self, path: str) -> None:
        self.path = Path(path)
        self._data: dict = {}
        self._dirty: int = 0
        if self.path.exists():
            try:
                with open(self.path) as f:
                    self._data = json.load(f)
                logger.info("Loaded cache with %d entries from %s", len(self._data), self.path)
            except (json.JSONDecodeError, OSError) as e:
                logger.warning("Could not load cache %s: %s", self.path, e)

    def has(self, pmid: str) -> bool:
        return pmid in self._data

    def get(self, pmid: str) -> list:
        return self._data.get(pmid, [])

    def set(self, pmid: str, data: list) -> None:
        self._data[pmid] = data
        self._dirty += 1
        if self._dirty >= self.WRITE_INTERVAL:
            self.flush()

    def flush(self) -> None:
        if self._dirty == 0:
            return
        try:
            # Write to a temp file and rename over the target — Path.replace()
            # is atomic (POSIX rename / Windows ReplaceFile), so a crash or
            # Ctrl+C mid-write can never leave a half-written, corrupted cache
            # file. Without this, __init__'s JSONDecodeError handler would
            # silently discard the *entire* cache (all previously cached
            # results) the next time --resume is used.
            tmp_path = self.path.with_name(self.path.name + ".tmp")
            with open(tmp_path, "w") as f:
                json.dump(self._data, f)
            tmp_path.replace(self.path)
            self._dirty = 0
        except OSError as e:
            logger.warning("Could not write cache: %s", e)


def build_abbreviation_map(text: str) -> dict:
    """Scan text for 'Genus species (G. species)' patterns and return an abbreviation map."""
    abbr_map: dict = {}
    pattern = re.compile(r"\b([A-Z][a-z]+\s+[a-z]+)\s+\(([A-Z]\.\s*[a-z]+)\)")
    for m in pattern.finditer(text):
        full_name = m.group(1)
        abbr = re.sub(r"\s+", " ", m.group(2)).strip()
        abbr_map[abbr] = full_name
    return abbr_map


def resolve_abbreviations(text: str, abbr_map: dict) -> str:
    """Replace abbreviated organism names with their full binomials."""
    for abbr, full in abbr_map.items():
        text = text.replace(abbr, full)
    return text


_ABBREV_SUFFIXES = {
    "sp", "spp", "et", "al", "fig", "vs", "approx", "cf",
    "no", "vol", "avg", "std", "i.e", "e.g",
}


def _ends_with_abbreviation(s: str, next_fragment: str = "") -> bool:
    last_token = s.rstrip().rsplit(None, 1)[-1] if s.strip() else ""
    # Single-letter genus abbreviation ("A.", "B.", "P.", ...) as used in
    # binomial species lists, e.g. "Azotobacter vinelandii, A. chroococcum,
    # and A. paspali" — never a genuine sentence end, but the base split
    # regex below only guards two-letter abbreviations like "Mr.", so this
    # merge step is what stops the list from being fragmented mid-sentence
    # (which would sever an organism's name from the evidence describing it).
    #
    # A trailing single capital is only an abbreviation when a lowercase word
    # follows it: that is the species epithet of a binomial. Element symbols
    # ending a sentence — "...the amount of N. Growth was measured...", very
    # common in this literature for N, P, K, C and Fe — are followed by a new
    # capitalised sentence instead, and merging there glued two unrelated
    # sentences together, letting evidence from one be attributed to an
    # organism named in the other.
    if re.fullmatch(r"[A-Z]\.", last_token):
        return bool(re.match(r"[a-z]", next_fragment.lstrip()))
    return last_token.rstrip(".").lower() in _ABBREV_SUFFIXES


def split_sentences(text: str) -> list:
    """Split text into sentences, merging fragments that follow known abbreviations."""
    raw = re.split(r"(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?|!)\s", text)
    merged: list = []
    for fragment in raw:
        if merged and _ends_with_abbreviation(merged[-1], fragment):
            merged[-1] += " " + fragment
        else:
            merged.append(fragment)
    return merged


def score_to_confidence(score: int) -> str:
    """Convert a numeric evidence score to a human-readable confidence tier."""
    if score >= 14:
        return "Very High"
    if score >= 9:
        return "High"
    if score >= 5:
        return "Medium"
    return "Low"


@dataclasses.dataclass
class ParsedMetric:
    """Structured representation of an extracted performance value."""
    raw: str
    value: float          # best single representative (midpoint if range)
    value_min: Optional[float]
    value_max: Optional[float]
    unit: str             # normalized unit string


# Unit normalization table (lowercase key → display string)
_UNIT_NORM: dict = {
    # IAA / general concentration
    "µg/ml": "µg/mL", "ug/ml": "µg/mL", "µg/l": "µg/L", "ug/l": "µg/L",
    "mg/l": "mg/L", "mg/ml": "mg/mL", "g/l": "g/L",
    "ppm": "ppm",
    # Halo zones
    "mm": "mm", "cm": "cm",
    # Enzymatic / gas production
    "nmol/mg/h": "nmol/mg/h", "µmol/mg/h": "µmol/mg/h", "umol/mg/h": "µmol/mg/h",
    "nmol": "nmol", "µmol": "µmol", "umol": "µmol",
    "ml/h": "mL/h", "ml/hr": "mL/h",
    # Phosphorus-annotated units (e.g. "mg P/L")
    "mg p/l": "mg P/L", "µg p/l": "µg P/L", "ug p/l": "µg P/L",
    # Nano-scale concentration
    "ng/ml": "ng/mL", "ng/l": "ng/L",
    # Rate units
    "nmol/h": "nmol/h", "µmol/h": "µmol/h", "umol/h": "µmol/h",
    "mmol/h": "mmol/h",
    # Dimensionless / ratios
    "%": "%", "si": "SI", "psi": "SI",
    "su": "SU", "fold": "fold",
}


_UNIT_CONTEXT_STRIP = re.compile(
    r"\b(?:halo|zone|diameter|production|activity|content|level|concentration"
    r"|su|units|efficiency|reduction|biofilm|inhibition"
    # assay qualifiers
    r"|cas|decolori[sz]ation|decolouri[sz]ation|solubili\w*|fixation"
    # phenotype substance qualifiers that can trail the unit (e.g. "µg/mL IAA")
    r"|iaa|auxin|indole|phosphate|phosphorus|nitrogen|hydrogen|siderophore"
    r"|ethylene|ketobutyrate|protein)\b.*$",
    flags=re.IGNORECASE,
)

# Normalise Unicode variants common in PDF-extracted scientific text.
# μ (U+03BC Greek mu) and µ (U+00B5 micro sign) are visually identical but
# byte-distinct — PDF extraction picks whichever the font encodes.
# − (U+2212 minus sign) often appears in superscripts like mL−1.
# × (U+00D7) is the multiplication sign used in scientific notation like 1.5×10³.
_UNICODE_TRANS = str.maketrans({
    '\u03bc': '\u00b5',  # μ → µ  (Greek mu → micro sign)
    '\u2212': '-',        # − → -  (minus sign → hyphen-minus)
    '\u00d7': 'x',        # × → x  (multiplication sign → plain x for sci notation)
})

# Superscript digit/sign characters → plain ASCII equivalents.
# Used to resolve exponents like 10³ → 10^3 before arithmetic expansion.
_SUPERSCRIPT_TRANS = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')

# Qualifiers that precede a value but carry no numeric meaning.
_APPROX_PREFIX = re.compile(
    r'^(?:~|≈|ca\.\s*|ca\s+|about\s+|approximately\s+|up\s+to\s+|at\s+least\s+|'
    r'more\s+than\s+|greater\s+than\s+|less\s+than\s+)',
    flags=re.IGNORECASE,
)


def _expand_scientific_notation(text: str) -> str:
    """Convert scientific notation in a metric string to a plain decimal.

    Handles:
      1.5×10³   (Unicode multiplication + superscript exponent)
      1.5x10^3  (plain x + caret)
      1.5x10-3  (plain x + signed exponent)
      1.5e3 / 1.5E-3  (standard e-notation)
    Returns the string with scientific notation replaced by its decimal value.
    """
    # Step 1 — resolve superscript exponents after "10" (e.g., 10³ → 103 context)
    # Translate superscript digits to normal digits first, but only in 10^sup patterns.
    def _sup_replace(m: re.Match) -> str:
        mantissa = float(m.group(1))
        exp = int(m.group(2).translate(_SUPERSCRIPT_TRANS))
        return f"{mantissa * (10.0 ** exp):.6g}"

    text = re.sub(
        r'(\d+(?:\.\d+)?)\s*x\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)',
        _sup_replace, text
    )

    # Step 2 — resolve x10^n or x10n with regular digit exponents
    def _sci_replace(m: re.Match) -> str:
        mantissa = float(m.group(1))
        exp = int(m.group(2))
        return f"{mantissa * (10.0 ** exp):.6g}"

    text = re.sub(
        r'(\d+(?:\.\d+)?)\s*x\s*10\s*\^?\s*([+-]?\d+)',
        _sci_replace, text
    )

    # Step 3 — standard e-notation (1.5e3, 2.4E-5)
    def _e_replace(m: re.Match) -> str:
        return f"{float(m.group(0)):.6g}"

    text = re.sub(r'\d+(?:\.\d+)?[eE][+-]?\d+', _e_replace, text)

    return text


def _sig_round(value: float, sig: int = 4) -> float:
    """Round to `sig` significant figures to suppress floating-point noise."""
    if value == 0.0:
        return 0.0
    magnitude = floor(log10(abs(value)))
    factor = 10 ** (sig - 1 - magnitude)
    return round(value * factor) / factor


def _normalize_text(text: str) -> str:
    """Normalise Unicode variants so unit patterns match regardless of PDF encoding."""
    return text.translate(_UNICODE_TRANS)


def _normalize_unit(raw_unit: str) -> str:
    # Strip trailing measurement-context words (e.g. "mm halo" → "mm")
    u = _UNIT_CONTEXT_STRIP.sub("", raw_unit).strip()
    # Strip superscript -1 notations
    u = re.sub(r"[-−]\s*1", "", u).strip().lower()
    return _UNIT_NORM.get(u, u if u else raw_unit.strip())


def parse_metric_value(raw: str) -> Optional[ParsedMetric]:
    """Parse a raw metric regex match into a structured ParsedMetric."""
    stripped = raw.strip()
    # Detect SI/index prefix before stripping so we can use it as the fallback unit
    is_si = bool(re.match(r"^(?:SI|PSI|index)\s*(?:=|of|:)\s*", stripped, flags=re.IGNORECASE))
    cleaned = re.sub(
        r"^(?:SI|PSI|index)\s*(?:=|of|:)\s*", "", stripped, flags=re.IGNORECASE
    )

    # Strip approximate qualifiers (~, ≈, "about", "up to", etc.) — keep the number.
    cleaned = _APPROX_PREFIX.sub("", cleaned).strip()

    # Resolve scientific notation (1.5×10³, 1.5e3, etc.) to plain decimals.
    cleaned = _expand_scientific_notation(cleaned)

    pm_pat = re.compile(r"(\d+(?:\.\d+)?)\s*±\s*(\d+(?:\.\d+)?)")
    range_pat = re.compile(r"(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)")
    single_pat = re.compile(r"(\d+(?:\.\d+)?)")

    # ± means "value ± uncertainty", not a range midpoint
    m_pm = pm_pat.search(cleaned)
    if m_pm:
        val_center = _sig_round(float(m_pm.group(1)))
        val_err = _sig_round(float(m_pm.group(2)))
        if 1900 <= val_center <= 2100:
            return None
        unit_raw = cleaned[m_pm.end():].strip().lstrip("/").strip()
        unit = _normalize_unit(unit_raw) or ("SI" if is_si else "")
        return ParsedMetric(raw=raw, value=val_center,
                            value_min=_sig_round(val_center - val_err),
                            value_max=_sig_round(val_center + val_err),
                            unit=unit)

    m_range = range_pat.search(cleaned)
    if m_range:
        vmin, vmax = _sig_round(float(m_range.group(1))), _sig_round(float(m_range.group(2)))
        if 1900 <= vmin <= 2100 or 1900 <= vmax <= 2100:
            return None
        value = _sig_round((vmin + vmax) / 2)
        unit_raw = cleaned[m_range.end():].strip().lstrip("/").strip()
        unit = _normalize_unit(unit_raw) or ("SI" if is_si else "")
        return ParsedMetric(raw=raw, value=value, value_min=vmin, value_max=vmax, unit=unit)

    m_single = single_pat.search(cleaned)
    if m_single:
        value = _sig_round(float(m_single.group(1)))
        if 1900 <= value <= 2100:
            return None
        unit_raw = cleaned[m_single.end():].strip().lstrip("/").strip()
        unit = _normalize_unit(unit_raw) or ("SI" if is_si else "")
        return ParsedMetric(raw=raw, value=value, value_min=None, value_max=None, unit=unit)

    return None


def save_results(df, base_path: str, fmt: str) -> None:
    """Save results DataFrame to csv, json, or both."""
    stem = base_path[: -len(Path(base_path).suffix)] if Path(base_path).suffix else base_path
    # Create the target directory rather than losing a completed run to a
    # missing folder (e.g. --out Result/results_siderophore.csv on a fresh clone).
    parent = Path(stem).parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)
    if fmt in ("csv", "both"):
        csv_path = f"{stem}.csv"
        # utf-8-sig writes a UTF-8 BOM so Excel on any locale (including Thai
        # Windows) opens the file correctly instead of misreading ± and µ as
        # Thai characters (Windows-874 encoding artefact).
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
        logger.info("Saved CSV → %s", csv_path)
    if fmt in ("json", "both"):
        json_path = f"{stem}.json"
        df.to_json(json_path, orient="records", indent=2)
        logger.info("Saved JSON → %s", json_path)


# "meliloti-1021", "subtilis-168": a species epithet carrying a strain
# designation. The epithet itself is still usable, so the suffix is trimmed
# rather than throwing the organism away. A purely alphabetic second half
# ("two-component") is NOT a strain number and stays rejected.
_STRAIN_SUFFIX_RE = re.compile(r"^([a-z]+)[.-]\d\w*$")


def canonical_organism_name(
    name: str, allow_unlisted_genera: bool = False
) -> Optional[str]:
    """Normalise `name` to "Genus species", or return None if it is not one.

    The species slot is where "Genus + next lowercase word" scanning goes
    wrong: prose such as "Salmonella in", "Pseudomonas that", "Rhizobium
    medium" and "Burkholderia two-component" all match the binomial shape and
    used to reach the output as organisms. A real epithet is a single
    alphabetic Latin word, so anything hyphenated, numeric or in the common
    English stoplist is rejected here.

    By default the genus must appear in the curated BACTERIAL_GENERA whitelist,
    which is precise but caps recall at the genera listed there. With
    `allow_unlisted_genera` (the --discover_genera flag) any plausible-looking
    binomial is accepted instead, and the caller is expected to confirm it
    against NCBI Taxonomy — the flag requires --verify_taxonomy for that reason.
    """
    clean_name = name.strip(".,;:()[]")
    parts = clean_name.split()
    if len(parts) < 2:
        return None
    genus = parts[0]
    species = parts[1].lower()

    strain_match = _STRAIN_SUFFIX_RE.match(species)
    if strain_match:
        species = strain_match.group(1)

    if species in INVALID_SPECIES_NAMES or species in COMMON_WORD_STOPLIST:
        return None
    if len(species) < 3 or not species.isalpha():
        return None

    if genus in BACTERIAL_GENERA:
        return f"{genus} {species}"
    if not allow_unlisted_genera:
        return None
    # Discovery mode: shape-based screen only. Anything surviving this still has
    # to resolve inside the Bacteria/Archaea subtree of NCBI Taxonomy before it
    # reaches the output, so the screen only needs to reject obvious prose
    # ("Results showed", "Soil samples") and known non-bacterial genera cheaply.
    if genus in PLANT_GENERA or genus.lower() in COMMON_WORD_STOPLIST:
        return None
    if len(genus) < 4 or not genus[1:].isalpha() or not genus[1:].islower():
        return None
    if len(species) < 4:
        return None
    return f"{genus} {species}"


def is_valid_name(name: str, allow_unlisted_genera: bool = False) -> bool:
    """True if `name` is a usable bacterial binomial (see canonical_organism_name)."""
    return canonical_organism_name(name, allow_unlisted_genera) is not None


_FULL_ORGANISM_RE = re.compile(r"\b([A-Z][a-z]+)\s+([a-z0-9]+(?:[\.-][a-z0-9]+)?)\b")
_ABBREV_ORGANISM_RE = re.compile(r"\b([A-Z])\.\s*([a-z]+(?:[\.-][a-z]+)?)\b")


def extract_organism_candidates(text: str, allow_unlisted_genera: bool = False) -> set:
    """Find every candidate organism name in `text`.

    A plain "Genus species" scan misses abbreviated binomials that appear in
    comma-separated species lists without an explicit "(G. species)" definition,
    e.g. "Bacillus subtilis, B. licheniformis, and B. amyloliquefaciens produced
    IAA" — only "Bacillus subtilis" has a capitalized full genus word, so
    "B. licheniformis" and "B. amyloliquefaciens" were previously dropped
    entirely rather than being extracted as their own organisms.

    This resolves each standalone "X. species" abbreviation against a preceding
    full genus mention sharing the same first letter — the standard convention
    once a genus has been spelled out once in scientific writing. If two
    *different* genera sharing that initial (e.g. Bacillus and Burkholderia,
    both "B.") were mentioned earlier in the same article, the abbreviation is
    ambiguous and is skipped rather than guessed — attributing a result to the
    wrong organism is worse than omitting it.
    """
    candidates: set = set()
    full_mentions: list = []  # (position, genus) for every full genus word seen

    for m in _FULL_ORGANISM_RE.finditer(text):
        full_mentions.append((m.start(), m.group(1)))
        canonical = canonical_organism_name(m.group(0), allow_unlisted_genera)
        if canonical:
            candidates.add(canonical)

    for m in _ABBREV_ORGANISM_RE.finditer(text):
        initial, species = m.group(1), m.group(2)
        pos = m.start()
        preceding_genera = {
            genus for gpos, genus in full_mentions if gpos < pos and genus[0] == initial
        }
        if len(preceding_genera) != 1:
            continue
        candidate = canonical_organism_name(
            f"{next(iter(preceding_genera))} {species}", allow_unlisted_genera
        )
        if candidate:
            candidates.add(candidate)

    return candidates


def extract_method_data(text: str, allowed_methods: dict) -> tuple:
    """Return (formatted_method_string, distinct_method_count)."""
    text_lower = text.lower()
    found = {
        method_name
        for keyword, method_name in allowed_methods.items()
        if keyword in text_lower
    }
    if not found:
        return "Not Explicitly Stated", 0
    return ", ".join(sorted(found)), len(found)


# Substrings that identify gene/molecular-identification methods.
# Any detected method whose display name contains one of these is excluded —
# the user wants only wet-lab assays that directly prove the phenotype.
_MOLECULAR_EXCLUSIONS: frozenset = frozenset({
    "pcr", "sequenc", "blast", "16s", "phylogen", "gene", "primer",
    "amplif", "cloning", "restriction", "rt-pcr", "qpcr", "genotyp",
    "rrna", "ribosom", "nucleotide", "dna extract", "southern", "northern blot",
})


def _is_molecular_method(display_name: str) -> bool:
    """Return True if the method name describes a gene/molecular technique, not a phenotypic assay."""
    name_lower = display_name.lower()
    return any(excl in name_lower for excl in _MOLECULAR_EXCLUSIONS)


def extract_method_with_proximity(
    sentence: str,
    allowed_methods: dict,
    phenotype_query: str,
    forced_key: Optional[str] = None,
) -> tuple:
    """Return (method_label, method_count) ranking methods by proximity to the phenotype keyword.

    When the phenotype keyword is present in the sentence, methods are sorted by
    their character distance to the nearest keyword occurrence — the closest method
    (most likely the one actually used to prove the phenotype) is listed first.
    When the phenotype keyword is absent, all methods are treated as equally relevant.
    Gene/molecular methods (PCR, sequencing, BLAST, etc.) are always excluded.
    """
    sent_lower = sentence.lower()

    _cfg = get_phenotype_config(phenotype_query, forced_key)
    _search_terms = _cfg.keywords if _cfg else [phenotype_query.lower().split()[0]]
    target_indices = [
        m.start()
        for term in _search_terms
        for m in re.finditer(re.escape(term.lower()), sent_lower)
    ]

    # Map each distinct method name to its minimum character distance to the phenotype keyword.
    # Skip any method that is gene/molecular-based, not a phenotypic assay.
    method_distances: dict = {}
    for keyword, method_name in allowed_methods.items():
        if _is_molecular_method(method_name):
            continue
        for m in re.finditer(re.escape(keyword), sent_lower):
            dist = (
                min(abs(m.start() - t) for t in target_indices)
                if target_indices else 0
            )
            if method_name not in method_distances or dist < method_distances[method_name]:
                method_distances[method_name] = dist

    if not method_distances:
        return "Not Explicitly Stated", 0

    # Sort by proximity to the phenotype keyword; closest (most relevant) first.
    ordered = sorted(method_distances, key=lambda n: method_distances[n])
    return ", ".join(ordered), len(ordered)


# Generic-phrasing fallback for methods/assays that aren't in the curated
# METHOD_CATEGORIES vocabulary (e.g. "disc diffusion assay", "seedling
# bioassay", "pot experiment") — pulls the assay name straight out of common
# reporting phrases instead of leaving it as "Not Explicitly Stated".
_METHOD_NOUN = r"assay|method|technique|protocol|procedure|bioassay"
_GENERIC_METHOD_PATTERNS = [
    re.compile(
        r"(?:using|via|by|through)\s+(?:the\s+|an?\s+)?"
        r"([A-Za-z][A-Za-z0-9\-]*(?:\s+[A-Za-z0-9\-]+){0,4}\s+(?:" + _METHOD_NOUN + r"))\b",
        re.IGNORECASE,
    ),
    re.compile(
        r"\b([A-Za-z][A-Za-z0-9\-]*(?:\s+[A-Za-z0-9\-]+){0,4}\s+(?:" + _METHOD_NOUN + r"))\s+"
        r"(?:was|were)\s+(?:used|performed|conducted|employed|applied)\b",
        re.IGNORECASE,
    ),
]


def extract_generic_method(sentence: str) -> str:
    """Best-effort extraction of an assay/method name from generic phrasing.

    Only used as a fallback when the curated dictionary lookup finds nothing —
    keeps the curated match as the authoritative source when both could apply.
    """
    for pattern in _GENERIC_METHOD_PATTERNS:
        m = pattern.search(sentence)
        if m:
            label = re.sub(r"\s+", " ", m.group(1)).strip(" -")
            if label and not _is_molecular_method(label):
                return " ".join(w.capitalize() if w.islower() else w for w in label.split())
    return ""


def extract_method_label(
    sentence: str,
    allowed_methods: dict,
    phenotype_query: str,
    forced_key: Optional[str] = None,
) -> tuple:
    """Method/assay lookup with a generic-phrasing fallback.

    Tries the curated per-phenotype vocabulary first (extract_method_with_proximity).
    If that finds nothing — e.g. the assay isn't one of the built-in categories,
    or the phenotype itself doesn't match a registered category — falls back to
    extract_generic_method so a method actually named in the text still surfaces
    in the output instead of "Not Explicitly Stated".
    """
    method_fmt, method_count = extract_method_with_proximity(
        sentence, allowed_methods, phenotype_query, forced_key
    )
    if method_fmt != "Not Explicitly Stated":
        return method_fmt, method_count
    generic = extract_generic_method(sentence)
    if generic:
        return f"{generic} (Reported)", 1
    return method_fmt, method_count


def extract_metrics_with_proximity(
    sentence: str, phenotype_query: str, forced_key: Optional[str] = None
) -> tuple:
    """
    Finds metrics filtered by distance to the phenotype keyword.
    Returns (raw_str, ParsedMetric|None) — the best match is also parsed into
    a structured numeric form for proper numeric sorting and unit normalization.
    Example: "Phosphate (99 mg/L), IAA (40 mg/L)".
    If searching for IAA, 40 mg/L is closer to IAA than 99 mg/L, so it wins.
    """
    # Normalize thousands-separator commas so "4,436.9" is matched as "4436.9"
    sentence = re.sub(r'(\d),(\d)', r'\1\2', sentence)
    sent_lower = sentence.lower()
    # Use all config keywords as proximity targets (e.g., "iaa", "auxin", "indole" for IAA),
    # not just the first word of the query.
    _cfg = get_phenotype_config(phenotype_query, forced_key)
    _search_terms = _cfg.keywords if _cfg else [phenotype_query.lower().split()[0]]
    target_indices: list = []
    for _term in _search_terms:
        target_indices.extend(m.start() for m in re.finditer(re.escape(_term.lower()), sent_lower))

    enemies = get_interference_words(phenotype_query, forced_key)
    enemy_indices: list = []
    for enemy in enemies:
        enemy_indices.extend([m.start() for m in re.finditer(re.escape(enemy), sent_lower)])

    patterns = get_patterns_for_phenotype(phenotype_query, forced_key)
    candidates: list = []

    for pattern in patterns:
        for m in re.finditer(pattern, sent_lower, re.IGNORECASE):
            # When group(0) starts with a non-digit, non-SI prefix (e.g. "OD595 = 1.23",
            # "absorbance at 595 ... 1.23") the captured numeric value is in group(1).
            # Using group(0) in those cases makes parse_metric_value pick up the
            # wavelength (595/570) instead of the actual reading.
            # Exception: "SI of 3.2" / "PSI = 1.8" — parse_metric_value strips that
            # prefix itself via its is_si logic, so group(0) must be passed.
            _g0 = m.group(0) or ""
            _use_g1 = (
                m.lastindex and m.group(1) is not None
                and _g0 and not _g0[0].isdigit()
                and not re.match(r'^(?:si|psi|index)\b', _g0, re.IGNORECASE)
            )
            val = m.group(1) if _use_g1 else _g0
            # Skip 4-digit year-like numbers
            if re.fullmatch(r"(?:19|20)\d{2}", val.strip()):
                continue
            metric_pos = m.start()
            dist_target = min(abs(metric_pos - t) for t in target_indices) if target_indices else 9999
            dist_enemy = min(abs(metric_pos - e) for e in enemy_indices) if enemy_indices else 9999
            # <= instead of < so that dist_target == dist_enemy == 9999 (neither keyword
            # present) still captures a candidate rather than silently dropping it.
            if dist_target <= dist_enemy:
                candidates.append((val, dist_target))

    if not candidates:
        return "N/A", None

    candidates.sort(key=lambda x: x[1])
    raw_str = "; ".join(list(dict.fromkeys([c[0] for c in candidates])))

    # Parse the closest (best) candidate into a structured metric
    best_parsed: Optional[ParsedMetric] = None
    for val, _ in candidates:
        parsed = parse_metric_value(val)
        if parsed is not None:
            best_parsed = parsed
            break

    return raw_str, best_parsed


def _find_organism_position(sent_lower: str, organism: str) -> Optional[int]:
    """Return the character offset of organism's most specific mention in the
    sentence — full binomial name, or abbreviated binomial (e.g. "b. subtilis")
    — or None if neither specific form appears.
    """
    parts = organism.lower().split()
    idx = sent_lower.find(" ".join(parts))
    if idx != -1:
        return idx
    if len(parts) >= 2:
        m = re.search(re.escape(parts[0][0]) + r'\.\s*' + re.escape(parts[1]), sent_lower)
        if m:
            return m.start()
    return None


def check_interference(sentence: str, target_organism: str, all_organisms: set) -> bool:
    """Return True if a *different genus* appears to own this sentence's evidence.

    Only cross-genus mentions count as interference (e.g. "Unlike E. coli, B.
    subtilis produced high IAA levels" — a different genus mentioned first
    suggests the evidence may actually belong to it). Other species of the
    *same* genus mentioned in the same breath (e.g. "Azotobacter vinelandii,
    A. chroococcum, and A. paspali were confirmed to fix nitrogen") are
    co-subjects of shared evidence, not competitors, so they're never treated
    as interference — and comparing them by genus-substring position alone
    would falsely flag every one of them anyway, since they'd all collide on
    the same genus offset.
    """
    sent_lower = sentence.lower()
    try:
        target_genus = target_organism.lower().split()[0]
        target_idx = _find_organism_position(sent_lower, target_organism)
        if target_idx is None:
            target_idx = sent_lower.find(target_genus)
        if target_idx == -1:
            return False
        for other in all_organisms:
            if other == target_organism:
                continue
            other_genus = other.lower().split()[0]
            if other_genus == target_genus:
                continue
            other_idx = _find_organism_position(sent_lower, other)
            if other_idx is None:
                other_idx = sent_lower.find(other_genus)
            if other_idx == -1:
                continue
            # Only flag interference if the competing organism is closer to the
            # start of the sentence than (or as close as) the target organism.
            if other_idx <= target_idx:
                return True
    except (AttributeError, ValueError):
        return False
    return False


_GENERIC_STRAIN_REFS = frozenset({
    "the strain", "this strain", "the isolate", "this isolate",
    "the bacterium", "this bacterium", "the organism", "this organism",
})


def _sentence_mentions_organism(sent_lower: str, organism: str) -> bool:
    """Return True if a sentence references the organism by any of:
    - full genus name  (e.g. "bacillus")
    - abbreviated form (e.g. "b. subtilis")
    - generic strain/isolate reference
    This is more permissive than a raw genus-in-sentence check, catching
    sentences that use pronouns or abbreviated binomials after the first mention.
    """
    parts = organism.lower().split()
    genus = parts[0]
    if genus in sent_lower:
        return True
    # Abbreviated form: first letter + "." + species (e.g. "b. subtilis")
    if len(parts) >= 2:
        if re.search(re.escape(parts[0][0]) + r'\.\s*' + re.escape(parts[1]), sent_lower):
            return True
    return any(ref in sent_lower for ref in _GENERIC_STRAIN_REFS)


# Clause-boundary markers: commas/semicolons/colons and contrast conjunctions.
# A single sentence often reports two different observations about two
# different traits ("X solubilized phosphate efficiently, but showed no
# siderophore production") — plain character distance isn't reliable enough
# to tell those apart in short sentences, so negation relevance is checked
# per-clause instead.
_CLAUSE_SPLIT_RE = re.compile(
    r"\b(?:though|but|however|while|whereas|although|yet|whilst)\b|[,;:]"
)


def _clause_index(pos: int, clause_bounds: list) -> int:
    for i in range(len(clause_bounds) - 1):
        if clause_bounds[i] <= pos < clause_bounds[i + 1]:
            return i
    return len(clause_bounds) - 2


def _find_relevant_negation(sent_lower: str, keyword_positions: list) -> Optional[re.Match]:
    """Return the first NEGATION_PATTERNS match that shares a clause with a
    phenotype keyword occurrence in this sentence, or None if there isn't one.

    A negation word appearing anywhere in a sentence isn't necessarily about
    the phenotype being scored: "no activity", "unable to", "undetected", and
    "negligible" are all generic enough to describe a completely different
    trait in the same (often compound) sentence — e.g. "...showed better
    antifungal activity ... indicating presence of undetected metabolite(s)"
    contains "undetected" but says nothing about phosphate solubilization.
    Requiring the negation to share a clause with an actual phenotype keyword
    mention keeps it tied to the phenotype actually under evaluation, rather
    than just being somewhere in the same sentence.
    """
    if not keyword_positions:
        return None

    # Deduplicated so a clause marker sitting at offset 0 (or two adjacent
    # markers, e.g. ", but") cannot create an empty clause that swallows
    # positions into the wrong index.
    clause_bounds = sorted(
        {0, len(sent_lower)} | {m.start() for m in _CLAUSE_SPLIT_RE.finditer(sent_lower)}
    )
    if len(clause_bounds) < 2:
        clause_bounds = [0, max(1, len(sent_lower))]
    keyword_clauses = {_clause_index(p, clause_bounds) for p in keyword_positions}

    # Every occurrence is checked, not just the first: a sentence can negate an
    # unrelated trait early ("no antifungal activity, but solubilised no
    # phosphate either") and the phenotype-relevant negation later. Stopping at
    # the first match of each pattern silently dropped those.
    for pattern in NEGATION_RES:
        for m in pattern.finditer(sent_lower):
            if _clause_index(m.start(), clause_bounds) in keyword_clauses:
                return m
    return None


def analyze_abstract(
    text: str,
    organism: str,
    phenotype: str,
    all_organisms: set,
    forced_key: Optional[str] = None,
    _pre_normalized: bool = False,
    _sentences: Optional[list] = None,
) -> tuple:
    if not _pre_normalized:
        # Expand abbreviations and normalize Unicode (done upstream in run() for efficiency)
        abbr_map = build_abbreviation_map(text)
        if abbr_map:
            text = resolve_abbreviations(text, abbr_map)
        text = _normalize_text(text)

    # Accept pre-split sentences to avoid re-splitting the same text per organism.
    sentences = _sentences if _sentences is not None else split_sentences(text)
    allowed_methods = get_allowed_methods(phenotype, forced_key)

    best_quant = {"val": "N/A", "parsed": None, "sent": "N/A", "method": "N/A", "score": -1, "method_src": ""}
    best_qual = {"sent": "", "method": "Not Explicitly Stated", "score": -1}

    term_root = phenotype.lower().split()[0]
    _cfg = get_phenotype_config(phenotype, forced_key)
    _search_terms = _cfg.keywords if _cfg else [term_root]

    for sent in sentences:
        sent_lower = sent.lower()
        if not _sentence_mentions_organism(sent_lower, organism):
            continue

        keyword_positions = [
            m.start()
            for term in _search_terms
            for m in re.finditer(re.escape(term.lower()), sent_lower)
        ]

        score = 0

        # Phenotype term present in sentence
        if keyword_positions:
            score += 2

        # Strong action stems — additive, bonus for multiple
        verb_hits = sum(1 for pattern in STRONG_VERB_RES if pattern.search(sent_lower))
        if verb_hits >= 1:
            score += 3
        if verb_hits >= 2:
            score += 1  # multi-verb bonus

        # Confirmation language (demonstrated, confirmed, showed, etc.)
        if CONFIRMATION_RE.search(sent_lower):
            score += 2

        # Standout/comparative signal words (highest, optimal, significant, etc.)
        if QUANTITATIVE_RE.search(sent_lower):
            score += 1

        # Method scoring: +2 base, +1 per additional distinct method (cap +2 extra).
        # Proximity-aware: methods nearest the phenotype keyword are ranked first,
        # so the label reflects the assay most likely used to prove the phenotype.
        method_fmt, method_count = extract_method_label(
            sent, allowed_methods, phenotype, forced_key
        )
        if method_count > 0:
            score += 2 + min(method_count - 1, 2)

        # Penalties
        if any(wk in sent_lower for wk in WEAK_CONTEXTS):
            score -= 11
        if _find_relevant_negation(sent_lower, keyword_positions):
            score = 0
        if check_interference(sent, organism, all_organisms):
            score -= 10

        if score <= 0:
            continue

        if score > best_qual["score"]:
            best_qual = {"sent": sent, "method": method_fmt, "score": score}

        metrics, parsed_metric = extract_metrics_with_proximity(sent, phenotype, forced_key)
        if metrics != "N/A":
            quant_score = score + 5
            if quant_score > best_quant["score"]:
                best_quant = {
                    "val": metrics,
                    "parsed": parsed_metric,
                    "sent": sent,
                    "method": method_fmt,
                    "score": quant_score,
                    "method_src": sent,
                }

    if best_quant["val"] != "N/A" and best_quant["method"] == "Not Explicitly Stated":
        for sent in sentences:
            sent_lower_fb = sent.lower()
            # Only infer from sentences that are about the phenotype and mention the organism,
            # to avoid picking up methods used for unrelated purposes in the abstract.
            if term_root not in sent_lower_fb:
                continue
            if not _sentence_mentions_organism(sent_lower_fb, organism):
                continue
            method_in_other, _ = extract_method_label(
                sent, allowed_methods, phenotype, forced_key
            )
            if method_in_other != "Not Explicitly Stated":
                # Strip a generic-fallback "(Reported)" tag before adding
                # "(Inferred)" so labels don't stack into "X (Reported) (Inferred)".
                if method_in_other.endswith(" (Reported)"):
                    method_in_other = method_in_other[: -len(" (Reported)")]
                best_quant["method"] = method_in_other + " (Inferred)"
                best_quant["method_src"] = sent
                break

    return best_qual, best_quant


def analyze_abstract_negative(
    text: str,
    organism: str,
    phenotype: str,
    all_organisms: set,
    forced_key: Optional[str] = None,
    _sentences: Optional[list] = None,
) -> dict:
    """Score sentences reporting that `organism` was tested and did NOT show the phenotype.

    This is the mirror image of analyze_abstract: instead of rewarding sentences
    for strong action verbs, it requires an explicit NEGATION_PATTERNS match
    (e.g. "failed to fix", "no nitrogenase activity", "undetected") that sits
    near an actual phenotype keyword mention in the sentence — via
    _find_relevant_negation, same as the positive-evidence path. This matters
    because most NEGATION_PATTERNS entries ("undetected", "negligible", "unable
    to", "no activity") are generic negation words with no phenotype-specific
    verb, so without this check a negation about a completely unrelated trait
    (e.g. "...showed better antifungal activity ... indicating presence of
    undetected metabolite(s)...") would get misread as proof this organism
    can't do the phenotype actually being searched for. A relevant negation
    alone still isn't enough — WEAK_CONTEXTS ("hypothesized", "future studies")
    and interference from a competing organism in the same sentence are still
    disqualifying, since the point is a confirmed wet-lab negative result for
    *this* organism specifically, not speculation or a result that actually
    belongs to a different organism mentioned nearby.
    """
    sentences = _sentences if _sentences is not None else split_sentences(text)
    allowed_methods = get_allowed_methods(phenotype, forced_key)
    term_root = phenotype.lower().split()[0]
    _cfg = get_phenotype_config(phenotype, forced_key)
    _search_terms = _cfg.keywords if _cfg else [term_root]

    best_neg = {"sent": "", "method": "Not Explicitly Stated", "score": -1}

    for sent in sentences:
        sent_lower = sent.lower()
        if not _sentence_mentions_organism(sent_lower, organism):
            continue

        keyword_positions = [
            m.start()
            for term in _search_terms
            for m in re.finditer(re.escape(term.lower()), sent_lower)
        ]
        if not _find_relevant_negation(sent_lower, keyword_positions):
            continue
        if any(wk in sent_lower for wk in WEAK_CONTEXTS):
            continue
        if check_interference(sent, organism, all_organisms):
            continue

        # base (negation confirmed) + phenotype term confirmed present
        # (guaranteed at this point, since _find_relevant_negation requires it)
        score = 3 + 2
        if CONFIRMATION_RE.search(sent_lower):
            score += 2

        method_fmt, method_count = extract_method_label(
            sent, allowed_methods, phenotype, forced_key
        )
        if method_count > 0:
            score += 2

        if score > best_neg["score"]:
            best_neg = {"sent": sent, "method": method_fmt, "score": score}

    return best_neg


# ==========================================
# MINER CLASS
# ==========================================


def _extract_pub_year(article_node) -> str:
    """Pull a 4-digit publication year out of a PubmedArticle XML node.

    PubMed dates aren't uniform: most articles have a clean <Year> element,
    but older/some journal records only provide a free-text <MedlineDate>
    (e.g. "2019 Winter" or "2020 Jan-Feb"), so that's regex-scanned as a fallback.
    """
    year = article_node.findtext(".//Article/Journal/JournalIssue/PubDate/Year")
    if year:
        return year
    medline_date = article_node.findtext(".//Article/Journal/JournalIssue/PubDate/MedlineDate")
    if medline_date:
        m = re.search(r"(19|20)\d{2}", medline_date)
        if m:
            return m.group(0)
    year = article_node.findtext(".//Article/ArticleDate/Year")
    if year:
        return year
    return ""


# Context words that indicate a nearby accession-looking token is actually a
# genome/sequence deposition, not an unrelated code (a catalog number, a dose,
# a plot ID, etc.) — the accession regexes below are loose enough on their
# own to false-positive without this guard.
_GENOME_CONTEXT_WORDS = (
    "genome", "assembly", "deposited", "accession", "genbank", "refseq",
    "bioproject", "biosample", "whole genome shotgun", "wgs", "ena", "ddbj",
)

_GENOME_ACCESSION_PATTERNS = [
    re.compile(r"\bGC[AF]_\d{9}\.\d+\b"),                # NCBI Assembly (GCF_/GCA_)
    re.compile(r"\b[A-Z]{4,6}\d{8,10}\b"),                # WGS master record
    re.compile(r"\bPRJ[EDN][A-Z]\d+\b"),                  # BioProject
    re.compile(r"\bSAM[EDN][A-Z]?\d+\b"),                 # BioSample
    re.compile(r"\b[A-Z]{1,2}\d{5,8}\.\d+\b"),            # Single GenBank nucleotide accession
]


def _all_organism_positions(text_lower: str, organism: str) -> list:
    """Every character offset where `organism` is mentioned — full binomial or
    abbreviated ("b. subtilis") form — used to find the mention *nearest* a
    given position rather than just the first one in the article.
    """
    parts = organism.lower().split()
    positions = [m.start() for m in re.finditer(re.escape(" ".join(parts)), text_lower)]
    if len(parts) >= 2:
        pattern = re.escape(parts[0][0]) + r"\.\s*" + re.escape(parts[1])
        positions.extend(m.start() for m in re.finditer(pattern, text_lower))
    return positions


def find_reported_genome_accessions(text: str, all_organisms: set) -> dict:
    """Resolve genome/assembly accessions the *paper itself* reports, per organism.

    Genome-announcement papers typically name the organism once early on and
    report the deposition accession in a later sentence without repeating the
    name there, so this doesn't require them to co-occur in the same sentence.
    Instead — mirroring the abbreviation-resolution convention used elsewhere
    in this file — each accession is attributed to whichever candidate
    organism was *most recently mentioned* (by any of its forms) immediately
    before it. A generic NCBI Assembly name-lookup finds *some* genome for a
    species, not necessarily the exact strain the article studied, so a
    paper-reported accession is meaningfully stronger evidence and is always
    preferred over that fallback (see genome_source in the output).

    Returns {organism_name: accession} — one pass per article, reused across
    every organism in it rather than recomputed per organism.
    """
    text_lower = text.lower()
    org_positions = {org: _all_organism_positions(text_lower, org) for org in all_organisms}

    result: dict = {}
    for pattern in _GENOME_ACCESSION_PATTERNS:
        for m in pattern.finditer(text):
            window_start = max(0, m.start() - 200)
            if not any(kw in text_lower[window_start:m.start()] for kw in _GENOME_CONTEXT_WORDS):
                continue

            nearest_org, nearest_pos = None, -1
            for org, positions in org_positions.items():
                preceding = [p for p in positions if p < m.start()]
                if preceding and max(preceding) > nearest_pos:
                    nearest_pos = max(preceding)
                    nearest_org = org

            if nearest_org is not None and nearest_org not in result:
                result[nearest_org] = m.group(0)

    return result


# Transient failure modes from a bulk NCBI E-utilities pull: dropped
# connections, incomplete/truncated HTTP responses, timeouts, and the
# malformed XML a truncated response often parses into. None of these mean
# the request itself was invalid — retrying usually succeeds.
_TRANSIENT_EXCEPTIONS = (
    urllib.error.URLError,
    http.client.HTTPException,  # includes http.client.IncompleteRead
    ConnectionError,
    TimeoutError,
    ET.ParseError,
)


def _retry_call(func, *, max_attempts: int = 3, base_delay: float = 2.0, description: str = "NCBI request"):
    """Call `func()`, retrying with exponential backoff on transient network failures.

    Without this, a single dropped connection mid-transfer (e.g. the
    "IncompleteRead" errors seen on large --max_articles pulls) silently costs
    an entire batch of articles with no retry — the caller's own try/except
    just logs it and moves on. This re-issues the whole request (a stale,
    partially-read handle can't just be resumed) up to `max_attempts` times
    before letting the exception propagate to the caller as before.
    """
    last_exc: BaseException = RuntimeError("_retry_call invoked with max_attempts < 1")
    for attempt in range(1, max_attempts + 1):
        try:
            return func()
        except _TRANSIENT_EXCEPTIONS as e:
            last_exc = e
            if attempt < max_attempts:
                delay = base_delay * (2 ** (attempt - 1))
                logger.warning(
                    "%s failed (attempt %d/%d): %s — retrying in %.1fs",
                    description, attempt, max_attempts, e, delay,
                )
                time.sleep(delay)
    raise last_exc


class MicroScopeMiner:
    def __init__(self, email: str, api_key: Optional[str] = None) -> None:
        from Bio import Entrez as _Entrez

        self._entrez = _Entrez
        self._entrez.email = email
        if api_key:
            self._entrez.api_key = api_key
        self.verified_cache: dict = {}
        self.genome_cache: dict = {}

    # NCBI returns at most 10,000 UIDs per esearch call.
    ESEARCH_PAGE_SIZE = 10000

    def search_pubmed(
        self, query: str, max_results: int, rate_limiter: Optional[RateLimiter] = None
    ) -> list:
        """Return up to `max_results` PMIDs, newest-first, paging as needed.

        sort="pub_date" returns PMIDs in NCBI's descending publication-date
        order instead of the "Best Match" relevance ranking, which otherwise
        surfaces old, highly-cited papers ahead of recent ones within the
        --max_articles cap.

        A single esearch caps out at 10,000 UIDs, and plain retstart paging is
        itself capped just below that — so anything larger is paged through the
        Entrez history server (usehistory + WebEnv/query_key). Before this,
        --max_articles above 10,000 silently returned only the first page.
        """
        collected: list = []
        webenv: Optional[str] = None
        query_key: Optional[str] = None

        while len(collected) < max_results:
            retstart = len(collected)
            retmax = min(self.ESEARCH_PAGE_SIZE, max_results - retstart)

            def _do():
                kwargs = {
                    "db": "pubmed",
                    "term": query,
                    "retmax": retmax,
                    "retstart": retstart,
                    "sort": "pub_date",
                }
                if webenv and query_key:
                    kwargs["WebEnv"] = webenv
                    kwargs["query_key"] = query_key
                else:
                    kwargs["usehistory"] = "y"
                handle = self._entrez.esearch(**kwargs)
                return self._entrez.read(handle)

            try:
                if rate_limiter is not None and collected:
                    rate_limiter.wait()
                record = _retry_call(_do, description="PubMed search")
            except urllib.error.URLError as e:
                logger.error("Network error searching PubMed: %s", e)
                break
            except Exception as e:
                logger.error("PubMed search failed: %s", e)
                break

            webenv = record.get("WebEnv", webenv)
            query_key = record.get("QueryKey", query_key)
            ids = record.get("IdList", [])
            collected.extend(ids)

            # Short page ⇒ the result set is exhausted; stop before asking NCBI
            # for a range that does not exist.
            if len(ids) < retmax:
                break

        return collected[:max_results]

    def fetch_details(self, pmids: list) -> list:
        if not pmids:
            return []

        def _do():
            handle = self._entrez.efetch(db="pubmed", id=",".join(pmids), retmode="xml")
            return ET.fromstring(handle.read())

        try:
            root = _retry_call(_do, description=f"fetch_details ({len(pmids)} PMIDs)")
            articles = []
            for art in root.findall(".//PubmedArticle"):
                pmid = art.findtext(".//PMID")
                title = art.findtext(".//ArticleTitle", default="")
                abstract = " ".join(
                    [e.text for e in art.findall(".//AbstractText") if e.text]
                )
                pmc_id = None
                for id_node in art.findall(".//ArticleId"):
                    if id_node.get("IdType") == "pmc":
                        pmc_id = id_node.text
                        break
                articles.append({
                    "pmid": pmid,
                    "pmc_id": pmc_id,
                    "title": title,
                    "text": f"{title}. {abstract}",
                    "pub_year": _extract_pub_year(art),
                })
            return articles
        except ET.ParseError as e:
            logger.error("XML parse error in fetch_details: %s", e)
            return []
        except urllib.error.URLError as e:
            logger.error("Network error in fetch_details: %s", e)
            return []
        except Exception as e:
            logger.error("fetch_details failed: %s", e)
            return []

    def fetch_full_text(self, pmc_id: str) -> str:
        def _do():
            handle = self._entrez.efetch(db="pmc", id=pmc_id, retmode="xml")
            return ET.parse(handle)

        try:
            tree = _retry_call(_do, description=f"fetch_full_text ({pmc_id})")
            root = tree.getroot()
            full_text = []
            for sec in root.findall(".//sec"):
                title_node = sec.find("title")
                if title_node is None or title_node.text is None:
                    continue
                title_lower = title_node.text.lower()
                if (
                    "result" in title_lower or "method" in title_lower
                    or "experimental" in title_lower
                    # "Data Availability" / "Accession numbers" sections are where
                    # genome/assembly deposition accessions are most often reported.
                    or "availability" in title_lower or "accession" in title_lower
                ):
                    paragraphs = [p.text for p in sec.findall(".//p") if p.text]
                    if paragraphs:
                        full_text.append(f"--- SECTION: {title_node.text} ---")
                        full_text.append(" ".join(paragraphs))
            return " ".join(full_text)
        except ET.ParseError as e:
            logger.debug("XML parse error for PMC %s: %s", pmc_id, e)
            return ""
        except urllib.error.URLError as e:
            logger.debug("Network error fetching PMC %s: %s", pmc_id, e)
            return ""
        except Exception as e:
            logger.debug("fetch_full_text failed for PMC %s: %s", pmc_id, e)
            return ""

    def verify_taxonomy(
        self, name: str, rate_limiter: RateLimiter, require_prokaryote: bool = False
    ) -> Optional[str]:
        """Resolve `name` to an NCBI Taxonomy ID, or None if it does not exist.

        With `require_prokaryote` (used by --discover_genera, where the genus
        whitelist no longer vouches for the name) the lookup is restricted to
        the Bacteria/Archaea subtree, so plants, fungi and insects picked up by
        the shape-based binomial screen are rejected instead of being reported
        as microbes.
        """
        if name in self.verified_cache:
            return self.verified_cache[name]

        def _do():
            if require_prokaryote:
                term = (
                    f'"{name}"[Scientific Name] '
                    f"AND (Bacteria[Subtree] OR Archaea[Subtree])"
                )
            else:
                term = f'"{name}"'
            handle = self._entrez.esearch(db="taxonomy", term=term, retmax=1)
            return self._entrez.read(handle)

        try:
            rate_limiter.wait()
            record = _retry_call(_do, description=f"verify_taxonomy ({name})")
            if record["IdList"]:
                self.verified_cache[name] = record["IdList"][0]
                return record["IdList"][0]
        except urllib.error.URLError as e:
            logger.warning("Network error verifying taxonomy for %s: %s", name, e)
        except Exception as e:
            logger.warning("Taxonomy lookup failed for %s: %s", name, e)
        self.verified_cache[name] = None
        return None

    def find_genome(self, organism: str, rate_limiter: RateLimiter) -> str:
        if organism in self.genome_cache:
            return self.genome_cache[organism]
        result = ""

        def _do():
            term = f'"{organism}"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])'
            handle = self._entrez.esearch(db="assembly", term=term, retmax=1)
            ids = self._entrez.read(handle)["IdList"]
            if not ids:
                return ""
            summ = self._entrez.esummary(db="assembly", id=ids[0], report="full")
            return self._entrez.read(summ)["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]

        try:
            rate_limiter.wait()
            result = _retry_call(_do, description=f"find_genome ({organism})")
        except urllib.error.URLError as e:
            logger.debug("Network error in find_genome for %s: %s", organism, e)
        except Exception as e:
            logger.debug("find_genome failed for %s: %s", organism, e)
        self.genome_cache[organism] = result
        return result

    def run(self, args) -> None:
        import pandas as pd

        term = args.phenotype
        forced_key: Optional[str] = getattr(args, "phenotype_category", None)
        if forced_key:
            forced_key = forced_key.upper()
            if forced_key not in PHENOTYPE_REGISTRY:
                logger.warning(
                    "Unknown --phenotype_category '%s'. Valid options: %s",
                    forced_key, ", ".join(PHENOTYPE_REGISTRY.keys()),
                )
                forced_key = None
            else:
                logger.info(
                    "Forcing phenotype category: %s (%s)",
                    forced_key, PHENOTYPE_REGISTRY[forced_key].display_name,
                )

        query = (
            f'("{term}"[Title/Abstract]) '
            f"AND (bacteria[MeSH Terms] OR bacterium[Title/Abstract] OR bacteria[Title/Abstract] OR bacterial[Title/Abstract])"
        )

        year_from = getattr(args, "year_from", None)
        year_to = getattr(args, "year_to", None)
        if year_from or year_to:
            lo = year_from or 1800
            hi = year_to or 3000
            query += f' AND ("{lo}"[Date - Publication] : "{hi}"[Date - Publication])'
            logger.info("Restricting to publication years %s–%s", lo, hi)

        rate_limiter = RateLimiter(0.11 if getattr(args, "api_key", None) else 0.34)

        logger.info("Searching PubMed for: %s", term)
        pmids = self.search_pubmed(query, args.max_articles, rate_limiter)
        logger.info("Found %d articles.", len(pmids))

        use_resume = getattr(args, "resume", False)
        cache_file = getattr(args, "cache_file", "microscope_cache.json")
        cache = DiskCache(cache_file) if use_resume else None

        results: list = []

        # Load cached results for PMIDs in this search
        if cache:
            for pmid in pmids:
                if cache.has(pmid):
                    results.extend(cache.get(pmid))
            if results:
                logger.info("Resumed %d cached results.", len(results))

        pmids_to_process = [p for p in pmids if not (cache and cache.has(p))]
        logger.info("Processing %d uncached articles...", len(pmids_to_process))

        # Fetch all article metadata upfront (rate-limited between batches)
        all_articles: list = []
        batch_size = 20
        for i in range(0, len(pmids_to_process), batch_size):
            batch = pmids_to_process[i : i + batch_size]
            rate_limiter.wait()
            all_articles.extend(self.fetch_details(batch))

        try:
            self._mine_articles(all_articles, args, term, forced_key, cache, results, rate_limiter)
        except KeyboardInterrupt:
            logger.warning(
                "Interrupted — saving %d result(s) collected so far before exiting.",
                len(results),
            )

        if cache:
            cache.flush()

        if results:
            cols = [
                "pmid", "pmc_id", "pub_year", "organism", "taxonomy_id", "phenotype",
                "phenotype_status", "confidence",
                "performance_value", "performance_unit", "performance_range",
                "performance_metric", "performance_method", "performance_evidence_sentence",
                "genome_accession", "genome_source", "evidence_score", "title",
            ]
            df = pd.DataFrame(results)
            # A cache written by an older MicroScope version can be missing
            # columns added since; fill them in rather than aborting a finished
            # run with a KeyError on the very last step.
            for missing in [c for c in cols if c not in df.columns]:
                logger.warning("Column '%s' missing from cached results — filling blank.", missing)
                df[missing] = None
            # Sort by numeric value (float) so 150 mg/L ranks above 10 mg/L correctly,
            # then by evidence score, then by publication year (most recent research
            # first among otherwise-tied results), then alphabetically by organism.
            df["_pub_year_sort"] = pd.to_numeric(df["pub_year"], errors="coerce")
            df = df.sort_values(
                by=["performance_value", "evidence_score", "_pub_year_sort", "organism"],
                ascending=[False, False, False, True],
                na_position="last",
            )
            df = df[cols]
            out_fmt = getattr(args, "output_format", "csv")
            save_results(df, args.out, out_fmt)
            logger.info(
                "%d result rows covering %d organism(s) from %d article(s).",
                len(df), df["organism"].nunique(), df["pmid"].nunique(),
            )
        else:
            logger.warning("No results found.")

    def _mine_articles(
        self, all_articles, args, term, forced_key, cache, results: list, rate_limiter: RateLimiter
    ) -> None:
        from tqdm import tqdm as _tqdm

        for art in _tqdm(all_articles, total=len(all_articles), desc="Mining", dynamic_ncols=True):
            text_to_analyze = art["text"]
            if art["pmc_id"]:
                rate_limiter.wait()
                full_text_content = self.fetch_full_text(art["pmc_id"])
                if full_text_content:
                    text_to_analyze += " " + full_text_content

            # Expand abbreviations and normalize Unicode once per article,
            # not once per organism (was repeated N times for the same text).
            _abbr_map = build_abbreviation_map(text_to_analyze)
            if _abbr_map:
                text_to_analyze = resolve_abbreviations(text_to_analyze, _abbr_map)
            text_to_analyze = _normalize_text(text_to_analyze)

            # Split sentences once per article — shared across all organisms.
            article_sentences = split_sentences(text_to_analyze)

            discover = getattr(args, "discover_genera", False)
            all_candidates: set = extract_organism_candidates(text_to_analyze, discover)
            reported_genomes = find_reported_genome_accessions(text_to_analyze, all_candidates)

            art_results: list = []
            # Sort for deterministic output order when scores are equal.
            for organism in sorted(all_candidates):
                if organism in CLONING_HOSTS and "isolated" not in text_to_analyze.lower():
                    continue

                qual, quant = analyze_abstract(
                    text_to_analyze, organism, term, all_candidates, forced_key,
                    _pre_normalized=True,
                    _sentences=article_sentences,
                )

                if args.negative_evidence:
                    # Inverted mode: keep only organisms wet-lab-proven to LACK the
                    # phenotype. Any organism with positive evidence anywhere in the
                    # article "can do it", so it's excluded here even if a different
                    # sentence also reads as negative (e.g. a different strain/context).
                    if qual["score"] > 0 or quant["val"] != "N/A":
                        continue
                    neg = analyze_abstract_negative(
                        text_to_analyze, organism, term, all_candidates, forced_key,
                        _sentences=article_sentences,
                    )
                    if neg["score"] < args.min_score:
                        continue
                    final_data = neg
                    score = neg["score"]
                    quant = {"val": "N/A", "parsed": None}
                else:
                    if args.strict_metric and quant["val"] == "N/A":
                        continue

                    final_data = quant if quant["val"] != "N/A" else qual
                    score = final_data["score"]
                    if score < args.min_score:
                        continue

                # Taxonomy verification runs only on organisms that already
                # cleared the evidence threshold. Verifying every candidate
                # binomial first meant one NCBI round-trip per name mentioned
                # anywhere in the article, the overwhelming majority of which
                # were then discarded for lack of evidence — the same rows come
                # out either way, just far fewer requests to get there.
                tax_id = ""
                if args.verify_taxonomy:
                    # In discovery mode the genus was never whitelisted, so the
                    # taxonomy check must also confirm the name is a prokaryote.
                    tax_id = self.verify_taxonomy(
                        organism, rate_limiter,
                        require_prokaryote=discover and organism.split()[0] not in BACTERIAL_GENERA,
                    )
                    if not tax_id:
                        continue

                final_sent = final_data["sent"]
                if (
                    "method_src" in final_data
                    and final_data["method_src"] not in ("", final_data["sent"])
                ):
                    final_sent = (
                        f"[RESULT]: {final_data['sent']} "
                        f"... [METHOD CONTEXT]: {final_data['method_src']}"
                    )

                # Prefer an accession the paper itself reports for this organism
                # over a generic NCBI Assembly name-lookup — a name-lookup finds
                # *some* genome for the species, not necessarily the exact strain
                # the article studied, so the two are not equivalent evidence.
                genome = reported_genomes.get(organism, "")
                if genome:
                    genome_source = "Reported in Paper"
                elif not args.skip_genome:
                    genome = self.find_genome(organism, rate_limiter)
                    genome_source = "NCBI Lookup" if genome else ""
                else:
                    genome_source = ""

                # Extract structured numeric value and unit from parsed metric
                parsed: Optional[ParsedMetric] = quant.get("parsed")
                perf_value: Optional[float] = parsed.value if parsed else None
                perf_unit: str = parsed.unit if parsed else ""
                perf_range: str = (
                    f"{parsed.value_min:g} – {parsed.value_max:g}"
                    if parsed and parsed.value_min is not None
                    else ""
                )

                art_results.append({
                    "pmid": art["pmid"],
                    "pmc_id": art["pmc_id"] if art["pmc_id"] else "N/A",
                    "pub_year": art.get("pub_year", ""),
                    "organism": organism,
                    "taxonomy_id": tax_id,
                    "phenotype": term,
                    "phenotype_status": "Negative" if args.negative_evidence else "Positive",
                    "confidence": score_to_confidence(score),
                    "performance_metric": quant["val"],
                    "performance_value": perf_value,
                    "performance_unit": perf_unit,
                    "performance_range": perf_range,
                    "performance_method": final_data["method"],
                    "performance_evidence_sentence": final_sent,
                    "evidence_score": score,
                    "genome_accession": genome,
                    "genome_source": genome_source,
                    "title": art["title"],
                })

            # Cache the empty list too: articles that yield nothing are the
            # majority of any run, and skipping them here meant --resume
            # re-downloaded and re-analysed every one of them on the next pass.
            if cache:
                cache.set(art["pmid"], art_results)
            results.extend(art_results)


# ==========================================
# CLI
# ==========================================


def install_dependencies() -> None:
    import importlib.util
    import subprocess

    in_venv = sys.prefix != sys.base_prefix

    def _check_and_install(package_name: str, import_name: Optional[str] = None) -> None:
        if import_name is None:
            import_name = package_name
        if importlib.util.find_spec(import_name) is not None:
            return
        logger.info("Installing %s...", package_name)
        base = [sys.executable, "-m", "pip", "install", package_name]
        # Inside a virtualenv pip installs cleanly; --break-system-packages is
        # only needed on externally-managed system interpreters (PEP 668), and
        # forcing it there unconditionally risks writing into the system
        # site-packages of a distro-managed Python.
        attempts = [base] if in_venv else [base, base[:-1] + ["--break-system-packages", package_name]]
        for cmd in attempts:
            try:
                subprocess.check_call(cmd)
                return
            except subprocess.CalledProcessError:
                continue
        logger.error("Failed to install %s. Please install it manually.", package_name)
        sys.exit(1)

    _check_and_install("biopython", "Bio")
    _check_and_install("pandas")
    _check_and_install("tqdm")
    _check_and_install("certifi")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=f"MicroScope v{__version__} — PubMed text-mining for microbial phenotypes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"MicroScope {__version__}")
    parser.add_argument("--phenotype", required=True, help="Phenotype to search (e.g. 'nitrogen fixation')")
    parser.add_argument("--email", required=True, help="Your email for NCBI Entrez")
    parser.add_argument("--api_key", help="NCBI API Key (enables 10 req/s rate limit)")
    parser.add_argument("--max_articles", type=int, default=100, help="Max PubMed articles to retrieve")
    parser.add_argument("--min_score", type=int, default=5, help="Minimum evidence score to include")
    parser.add_argument("--strict_metric", action="store_true", help="Only save results with numerical data")
    parser.add_argument(
        "--negative_evidence",
        action="store_true",
        help=(
            "Invert the search: find organisms wet-lab-proven to NOT have the "
            "phenotype (e.g. 'failed to fix nitrogen', 'no siderophore activity "
            "detected'), and exclude any organism that shows positive evidence "
            "for the phenotype elsewhere in the same article."
        ),
    )
    parser.add_argument("--verify_taxonomy", action="store_true", help="Verify organisms against NCBI Taxonomy")
    parser.add_argument(
        "--discover_genera",
        action="store_true",
        help=(
            "Accept bacterial genera that are not in the built-in whitelist. Each "
            "such name must resolve inside the Bacteria/Archaea subtree of NCBI "
            "Taxonomy, so this requires --verify_taxonomy and costs one extra "
            "lookup per newly seen name."
        ),
    )
    parser.add_argument("--year_from", type=int, help="Earliest publication year to include")
    parser.add_argument("--year_to", type=int, help="Latest publication year to include")
    parser.add_argument("--skip_genome", action="store_true", help="Skip genome assembly lookup")
    parser.add_argument("--out", default="results_MicroScope.csv", help="Output file base path")
    parser.add_argument("--resume", action="store_true", help="Skip PMIDs already in cache")
    parser.add_argument("--cache_file", default="microscope_cache.json", help="Path to disk cache file")
    parser.add_argument(
        "--output_format",
        choices=["csv", "json", "both"],
        default="csv",
        help="Output format (default: csv)",
    )
    parser.add_argument(
        "--log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity (default: INFO)",
    )
    parser.add_argument(
        "--phenotype_category",
        choices=list(PHENOTYPE_REGISTRY.keys()),
        help="Force explicit phenotype category instead of auto-detect",
    )
    args = parser.parse_args()

    if args.discover_genera and not args.verify_taxonomy:
        parser.error(
            "--discover_genera requires --verify_taxonomy: without an NCBI Taxonomy "
            "check, unlisted genus names are accepted on shape alone and the output "
            "fills with non-bacterial false positives."
        )
    if args.year_from and args.year_to and args.year_from > args.year_to:
        parser.error("--year_from cannot be later than --year_to")
    if args.max_articles < 1:
        parser.error("--max_articles must be at least 1")

    setup_logging(args.log_level)
    install_dependencies()

    import ssl
    import certifi

    try:
        ssl._create_default_https_context = lambda: ssl.create_default_context(
            cafile=certifi.where()
        )
    except AttributeError:
        pass

    logger.info("MicroScope v%s starting...", __version__)
    miner = MicroScopeMiner(args.email, args.api_key)
    try:
        miner.run(args)
    except KeyboardInterrupt:
        logger.info("Interrupted by user.")


if __name__ == "__main__":
    main()
