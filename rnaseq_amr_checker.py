# rnaseq_amr_checker.py
# =========================
# Multi-species Conservative RNA-seq + AMR checker (LOCAL)
# LAB-ONLY AMR (strict + evidence/source backup, blocks computational)
# Runs species ONE BY ONE (sequential)
#
# UPDATED (SRA + ENA, STRONG BioSample-only):
# - RNA-seq run discovery uses BOTH:
#     (A) NCBI SRA runinfo (targeted BioSample queries, not species-wide)
#     (B) ENA filereport by BioSample
# - STRONG ONLY: join is ONLY via the SAME BioSample accession (no BioProject fallback)
# - Accepts SRR/ERR/DRR (not SRR-only)
#
# OUTPUT RULE:
# - If exactly 1 species succeeds: CSV outputs only (NO ZIP)
# - If 2+ species succeed: ONE ZIP ONLY in outputs folder.
#     - ZIP contains: all per-species outputs + combined rnaseq_true CSV
#     - Per-species CSVs are cleaned up after zipping (not left loose in outputs/)
# =========================

import os
import time
import random
import json
import shutil
import zipfile
import re
import argparse
from io import StringIO
from urllib.parse import quote
from typing import Optional, List, Tuple

import pandas as pd
import requests
from tqdm import tqdm

# -------------------------
# DEFAULT CONFIG (override via CLI)
# -------------------------
SPECIES_INPUTS_DEFAULT = (
    "Enterococcus faecium, Staphylococcus aureus, Klebsiella pneumoniae, "
    "Acinetobacter baumannii, Pseudomonas aeruginosa, Enterobacter"
)

RESET_CACHE_FOR_EACH_SPECIES = False
ONLY_GENOMES_WITH_AMR = True

# Console verbosity controls (noise control)
PRINT_STEP_STATS = True
PRINT_EVERY_N_TAXON_IDS = 25
SHOW_TQDM_BARS = True
BVBRC_VERBOSE = False

# Optional caps (None = fetch everything)
MAX_BVBRC_GENOMES = None
MAX_BVBRC_AMR_RECORDS_PER_TAXON = None

# Caps for RNA side (optional safety)
MAX_BIOSAMPLES_PER_SPECIES = None
MAX_SRA_RECORDS = None

# Paging sizes
BVBRC_PAGE = 25000
SRA_EFETCH_PAGE = 500

# NCBI env
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "").strip()
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
RPS = 2 if not NCBI_API_KEY else 6
MIN_INTERVAL = 1.0 / RPS

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BVBRC_API = "https://www.bv-brc.org/api/"

# ENA
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Regex
RUN_RE = re.compile(r"(?:SRR|ERR|DRR)\d+", re.IGNORECASE)
BS_TOKEN_RE = re.compile(r"(?:SAMN|SAMEA|SAMD)\d+", re.IGNORECASE)
ASM_RE = re.compile(r"(?:GCA|GCF)_\d+\.\d+", re.IGNORECASE)

sess = requests.Session()
sess.headers.update({"User-Agent": "rnaseq_amr_checker_local/6.4"})
_last_call = 0.0


def norm_species_name(name: str) -> str:
    name = " ".join(str(name).strip().split())
    if not name:
        raise ValueError("Empty species name.")
    words = name.split()
    words[0] = words[0].capitalize()
    for i in range(1, len(words)):
        words[i] = words[i].lower()
    return " ".join(words)


def rql_val(s: str) -> str:
    return quote(str(s), safe="")


def _params(extra: dict):
    p = dict(extra)
    p["tool"] = "rnaseq_amr_checker"
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    return p


def eutils_get(endpoint: str, params: dict, timeout: int = 60, max_retries: int = 8) -> requests.Response:
    global _last_call
    url = EUTILS + endpoint
    for attempt in range(max_retries):
        now = time.monotonic()
        wait = (_last_call + MIN_INTERVAL) - now
        if wait > 0:
            time.sleep(wait)

        r = sess.get(url, params=_params(params), timeout=timeout)
        _last_call = time.monotonic()

        if r.status_code == 429 or (500 <= r.status_code < 600):
            ra = r.headers.get("Retry-After")
            if ra and ra.isdigit():
                sleep_s = int(ra)
            else:
                sleep_s = min(60, 2 ** attempt) + random.uniform(0, 0.5)
            time.sleep(sleep_s)
            continue

        r.raise_for_status()
        return r

    r.raise_for_status()
    return r


def http_get_retry(url: str, timeout: int = 120, max_retries: int = 10, backoff_base: float = 1.6) -> requests.Response:
    last = None
    for attempt in range(max_retries):
        try:
            r = sess.get(url, timeout=timeout)
            if r.status_code in (429, 500, 502, 503, 504):
                ra = r.headers.get("Retry-After")
                if ra and ra.isdigit():
                    sleep_s = int(ra)
                else:
                    sleep_s = min(90, (backoff_base ** attempt)) + random.uniform(0, 0.8)
                time.sleep(sleep_s)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            last = e
            sleep_s = min(90, (backoff_base ** attempt)) + random.uniform(0, 0.8)
            time.sleep(sleep_s)
    raise last


def bvbrc_count(endpoint: str, rql_core: str, timeout: int = 120, label: Optional[str] = None) -> int:
    url = f"{BVBRC_API}{endpoint}/?{rql_core}&limit(0,0)&http_accept=application/solr+json"
    r = http_get_retry(url, timeout=timeout)
    j = r.json()
    n = int(j.get("response", {}).get("numFound", 0))
    if BVBRC_VERBOSE:
        tag = f" | {label}" if label else ""
        print(f"[BV-BRC]{tag} {endpoint} expected total: {n}")
    return n


def bvbrc_fetch_all(
    endpoint: str,
    rql_core: str,
    select_fields: list,
    page: int = BVBRC_PAGE,
    max_rows: Optional[int] = None,
    timeout: int = 120,
    label: Optional[str] = None
) -> Tuple[pd.DataFrame, int]:
    expected = bvbrc_count(endpoint, rql_core, timeout=timeout, label=label)

    out = []
    start = 0
    while True:
        if max_rows is not None and start >= max_rows:
            break

        this_limit = page if max_rows is None else min(page, max_rows - start)

        url = (
            f"{BVBRC_API}{endpoint}/?"
            f"{rql_core}"
            f"&select({','.join(select_fields)})"
            f"&limit({this_limit},{start})"
            f"&http_accept=application/json"
        )

        r = http_get_retry(url, timeout=timeout)
        batch = r.json()
        if not batch:
            break

        out.extend(batch)
        got = len(batch)
        start += got

        if got < this_limit:
            break

    df = pd.DataFrame(out)

    for col in select_fields:
        if col not in df.columns:
            df[col] = pd.Series(dtype="object")

    df = df[select_fields + [c for c in df.columns if c not in select_fields]]
    return df, expected


def norm_pheno(x) -> str:
    x = str(x).strip().lower()
    if x in {"resistant", "r"}:
        return "R"
    if x in {"susceptible", "s"}:
        return "S"
    if x in {"intermediate", "i"}:
        return "I"
    return ""


def summarise_amr_per_genome(ab_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gid, g in ab_df.groupby("genome_id"):
        mp = dict(zip(g["antibiotic"], g["pheno"]))
        conflicts = sorted([a for a, p in mp.items() if len(p) > 1])
        resistant_any = sorted([a for a, p in mp.items() if "R" in p])
        susceptible_any = sorted([a for a, p in mp.items() if "S" in p])
        intermediate_any = sorted([a for a, p in mp.items() if "I" in p])

        rows.append({
            "genome_id": str(gid),
            "amr_antibiotics_count": len(mp),
            "amr_resistant_count": len(resistant_any),
            "amr_susceptible_count": len(susceptible_any),
            "amr_intermediate_count": len(intermediate_any),
            "amr_conflict_count": len(conflicts),
            "amr_conflict_antibiotics": ";".join(conflicts),
            "antibiotics_resistant": ";".join(resistant_any),
            "antibiotics_susceptible": ";".join(susceptible_any),
            "antibiotics_intermediate": ";".join(intermediate_any),
        })
    return pd.DataFrame(rows)


def pick_first_token(x: str, token_re: re.Pattern) -> str:
    x = str(x or "")
    m = token_re.search(x)
    return m.group(0).upper() if m else ""


def extract_run_ids(x: str) -> list:
    return sorted(set([m.group(0).upper() for m in RUN_RE.finditer(str(x or ""))]))


def union_run_strings(*vals) -> str:
    s = set()
    for v in vals:
        for rid in extract_run_ids(v):
            s.add(rid)
    return ";".join(sorted(s))


def ena_filereport(accessions: List[str], fields: List[str], chunk_size: int = 25, timeout: int = 120) -> pd.DataFrame:
    if not accessions:
        return pd.DataFrame()

    frames = []
    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i:i + chunk_size]
        params = {
            "accession": ",".join(chunk),
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
            "download": "false",
        }
        r = sess.get(ENA_PORTAL, params=params, timeout=timeout)
        if r.status_code in (429, 500, 502, 503, 504):
            r = http_get_retry(r.url, timeout=timeout)

        r.raise_for_status()
        txt = r.text.strip()
        if not txt:
            continue

        df = pd.read_csv(StringIO(txt), sep="\t")
        if not df.empty:
            frames.append(df)

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def build_sra_runinfo_by_biosamples(biosamples: List[str], chunk_size: int = 60, max_records: Optional[int] = None) -> pd.DataFrame:
    if not biosamples:
        return pd.DataFrame()

    dfs = []
    fetched_total = 0

    itr = range(0, len(biosamples), chunk_size)
    itr = tqdm(itr, desc="SRA BioSample chunks", disable=(not SHOW_TQDM_BARS))

    for i in itr:
        chunk = biosamples[i:i + chunk_size]
        term = "(" + " OR ".join([f"{bs}[BioSample]" for bs in chunk]) + ") AND rna seq[Strategy]"

        j0 = eutils_get("esearch.fcgi", {
            "db": "sra", "term": term, "retmode": "json", "retmax": 0, "usehistory": "y"
        }).json()

        es = j0.get("esearchresult", {})
        total = int(es.get("count", 0))
        webenv = es.get("webenv", "")
        query_key = es.get("querykey", "")

        if total == 0:
            continue

        to_fetch = total
        if max_records is not None:
            remaining = max_records - fetched_total
            if remaining <= 0:
                break
            to_fetch = min(to_fetch, remaining)

        for start in range(0, to_fetch, SRA_EFETCH_PAGE):
            this_max = min(SRA_EFETCH_PAGE, to_fetch - start)
            txt = eutils_get("efetch.fcgi", {
                "db": "sra",
                "query_key": query_key,
                "WebEnv": webenv,
                "rettype": "runinfo",
                "retmode": "text",
                "retstart": start,
                "retmax": this_max
            }, timeout=120).text.strip()

            if txt:
                df = pd.read_csv(StringIO(txt))
                dfs.append(df)
                fetched_total += len(df)

            if max_records is not None and fetched_total >= max_records:
                break

        if max_records is not None and fetched_total >= max_records:
            break

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def run_map_biosample_only_sra_ena(
    genomes_amr_filtered: pd.DataFrame,
    cache_dir: str,
    max_biosamples: Optional[int] = None,
    max_sra_records: Optional[int] = None
) -> pd.DataFrame:
    g = genomes_amr_filtered.copy()
    g["biosample_accession"] = g["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))

    biosamples = sorted([x for x in g["biosample_accession"].unique().tolist() if x])
    if max_biosamples is not None:
        biosamples = biosamples[:max_biosamples]

    # --- SRA ---
    runinfo = build_sra_runinfo_by_biosamples(biosamples, chunk_size=60, max_records=max_sra_records)

    sra_bs_map = pd.DataFrame(columns=["biosample_accession", "rnaseq_run_ids_sra"])
    if not runinfo.empty:
        for col in ["Run", "BioSample", "LibraryStrategy"]:
            if col not in runinfo.columns:
                runinfo[col] = ""

        runinfo["Run"] = runinfo["Run"].fillna("").astype(str).str.strip().str.upper()
        runinfo["BioSample"] = runinfo["BioSample"].fillna("").astype(str).str.strip().str.upper()
        runinfo["LibraryStrategy"] = runinfo["LibraryStrategy"].fillna("").astype(str).str.strip().str.lower()

        ri = runinfo[
            (runinfo["BioSample"] != "") &
            (runinfo["Run"] != "") &
            (runinfo["LibraryStrategy"].str.contains("rna"))
        ].copy()

        sra_bs_map = (
            ri.groupby("BioSample")["Run"]
            .apply(lambda s: ";".join(sorted(set([x for x in s if RUN_RE.search(str(x))]))[:500]))
            .reset_index()
            .rename(columns={"BioSample": "biosample_accession", "Run": "rnaseq_run_ids_sra"})
        )

    # --- ENA ---
    ena_fields = ["run_accession", "sample_accession", "secondary_sample_accession", "library_strategy", "fastq_ftp", "fastq_md5"]
    ena = ena_filereport(biosamples, fields=ena_fields, chunk_size=25, timeout=120)

    ena_bs_map = pd.DataFrame(columns=["biosample_accession", "rnaseq_run_ids_ena", "ena_fastq_ftp", "ena_fastq_md5"])
    if not ena.empty:
        for c in ena_fields:
            if c not in ena.columns:
                ena[c] = ""

        ena["library_strategy"] = ena["library_strategy"].fillna("").astype(str).str.lower()
        ena = ena[ena["library_strategy"].str.contains("rna")].copy()

        ena["biosample_accession"] = ena["secondary_sample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
        mask_empty = ena["biosample_accession"].eq("")
        ena.loc[mask_empty, "biosample_accession"] = ena.loc[mask_empty, "sample_accession"].apply(
            lambda x: pick_first_token(x, BS_TOKEN_RE)
        )

        ena["run_accession"] = ena["run_accession"].fillna("").astype(str).str.strip().str.upper()
        ena["fastq_ftp"] = ena["fastq_ftp"].fillna("").astype(str).str.strip()
        ena["fastq_md5"] = ena["fastq_md5"].fillna("").astype(str).str.strip()

        ena = ena[(ena["biosample_accession"] != "") & (ena["run_accession"] != "")]
        ena_bs_map = (
            ena.groupby("biosample_accession")
            .agg(
                rnaseq_run_ids_ena=("run_accession", lambda s: ";".join(sorted(set([x for x in s if RUN_RE.search(str(x))]))[:500])),
                ena_fastq_ftp=("fastq_ftp", lambda s: ";".join(sorted(set([x for x in s if x]))[:50])),
                ena_fastq_md5=("fastq_md5", lambda s: ";".join(sorted(set([x for x in s if x]))[:50])),
            )
            .reset_index()
        )

    bs_map = sra_bs_map.merge(ena_bs_map, on="biosample_accession", how="outer")
    bs_map["rnaseq_run_ids"] = bs_map.apply(
        lambda r: union_run_strings(r.get("rnaseq_run_ids_sra", ""), r.get("rnaseq_run_ids_ena", "")),
        axis=1
    )

    bs_map = bs_map[[
        "biosample_accession",
        "rnaseq_run_ids",
        "rnaseq_run_ids_sra",
        "rnaseq_run_ids_ena",
        "ena_fastq_ftp",
        "ena_fastq_md5"
    ]].copy()

    bs_map.to_parquet(os.path.join(cache_dir, "biosample_run_map_SRA_ENA_STRONG.parquet"), index=False)
    return bs_map


# -------------------------
# LAB-ONLY AMR filter
# -------------------------
LAB_POS = re.compile(
    r"(?:disk|disc)\s*diffusion|kirby\s*[- ]?bauer|broth\s*dilution|agar\s*dilution|"
    r"\bmic\b|minimum\s*inhibitory\s*concentration|e-?\s*test|etest|"
    r"\bvitek\b|\bphoenix\b|sensititre|"
    r"\bclsi\b|\beucast\b|"
    r"\bast\b|antimicrobial\s*susceptibility\s*test|phenotypic|phenotype|"
    r"zone\s*diameter|inhibition\s*zone",
    re.IGNORECASE
)

NONLAB_NEG = re.compile(
    r"predicted|prediction|in\s*silico|computational|genotype|genomic|"
    r"resistance\s*gene|amr\s*gene|gene\s*presence|"
    r"amrfinder|resfinder|card|rgi|abricate|"
    r"k-?mer|machine\s*learning|\bml\b|model|"
    r"homology|blast|annotation|assembly|variant|mutation|snp|"
    r"rule[- ]based|inferred",
    re.IGNORECASE
)


def lab_only_mask(amr: pd.DataFrame) -> pd.Series:
    n = len(amr)
    if n == 0:
        return pd.Series([], dtype=bool)

    def col_text(c):
        if c not in amr.columns:
            return pd.Series([""] * n)
        return amr[c].fillna("").astype(str).str.strip()

    lab_method = col_text("laboratory_typing_method")
    standard = col_text("testing_standard")
    meas = col_text("measurement_value")
    evidence = col_text("evidence")
    source = col_text("source")

    strong_lab = (lab_method.ne("") | standard.ne("") | meas.ne(""))
    txt = (lab_method + " " + standard + " " + meas + " " + evidence + " " + source).str.lower()

    has_nonlab = txt.str.contains(NONLAB_NEG)
    has_labpos = txt.str.contains(LAB_POS)

    keep = (~has_nonlab) & (strong_lab | has_labpos)
    return keep


def make_species_cache_dir(cache_root: str, species: str) -> str:
    species_key = species.lower().replace(" ", "_")
    d = os.path.join(cache_root, species_key)
    os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join(d, "amr_by_taxon"), exist_ok=True)
    return d


def run_one_species(species_input: str, cache_root: str, out_root: str) -> Tuple[str, str, str]:
    stats = []
    step = 0

    def stat(msg: str, force_print: bool = True, **kv):
        nonlocal step
        step += 1
        stats.append({"step": step, "msg": msg, **kv})
        if PRINT_STEP_STATS and force_print:
            line = f"[{step:04d}] {msg}"
            if kv:
                line += " | " + " | ".join([f"{k}={v}" for k, v in kv.items()])
            print(line)

    def warn(msg: str, **kv):
        stats.append({"step": "WARN", "msg": msg, **kv})
        line = f"[!!] {msg}"
        if kv:
            line += " | " + " | ".join([f"{k}={v}" for k, v in kv.items()])
        print(line)

    species = norm_species_name(species_input)
    cache_dir = make_species_cache_dir(cache_root, species)

    if RESET_CACHE_FOR_EACH_SPECIES and os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
        cache_dir = make_species_cache_dir(cache_root, species)

    print(f"\n[Species] {species} | cache={cache_dir}")

    def cache_path(name: str) -> str:
        return os.path.join(cache_dir, name)

    def amr_taxon_cache_paths(tid: int) -> Tuple[str, str]:
        base = os.path.join(cache_dir, "amr_by_taxon", f"taxon_{tid}")
        return base + ".parquet", base + ".meta.json"

    # 1) Genomes (cached)
    genomes_file = cache_path("genomes.parquet")
    if os.path.exists(genomes_file):
        genomes = pd.read_parquet(genomes_file)
        stat("Loaded genomes from cache", True, rows=len(genomes))
    else:
        stat("Fetching BV-BRC genomes by species string (cap-safe)", True, species=species)
        rql_species = f"eq(species,{rql_val(species)})"

        required_fields = ["genome_id", "genome_name", "strain", "taxon_id", "biosample_accession", "bioproject_accession"]
        accession_fields = ["assembly_accession", "genbank_accessions", "refseq_accessions", "sra_accession"]
        optional_fields = ["genome_status", "genome_quality", "genome_quality_flags"]

        try:
            genomes, exp_g = bvbrc_fetch_all(
                endpoint="genome",
                rql_core=rql_species,
                select_fields=required_fields + optional_fields + accession_fields,
                max_rows=MAX_BVBRC_GENOMES,
                label=f"species={species}"
            )
        except Exception as e:
            warn("Genome fetch with accession fields failed; retrying without them", error=str(e)[:200])
            try:
                genomes, exp_g = bvbrc_fetch_all(
                    endpoint="genome",
                    rql_core=rql_species,
                    select_fields=required_fields + optional_fields,
                    max_rows=MAX_BVBRC_GENOMES,
                    label=f"species={species}"
                )
            except Exception as e2:
                warn("Optional genome fields failed; retrying required fields only", error=str(e2)[:200])
                genomes, exp_g = bvbrc_fetch_all(
                    endpoint="genome",
                    rql_core=rql_species,
                    select_fields=required_fields,
                    max_rows=MAX_BVBRC_GENOMES,
                    label=f"species={species}"
                )

        if genomes.empty:
            raise ValueError(f"No genomes returned for species='{species}'. (Genus-only like 'Enterobacter' may return 0.)")

        # Ensure/clean string columns
        str_cols = required_fields + accession_fields
        for c in str_cols:
            if c in genomes.columns:
                genomes[c] = genomes[c].fillna("").astype(str).str.strip()
            else:
                genomes[c] = ""

        genomes["taxon_id"] = pd.to_numeric(genomes["taxon_id"], errors="coerce")

        # Convenience columns: extract assembly accessions if present in genbank/refseq accession strings
        genomes["genbank_assembly_accession"] = genomes["genbank_accessions"].apply(lambda x: pick_first_token(x, ASM_RE))
        genomes["refseq_assembly_accession"] = genomes["refseq_accessions"].apply(lambda x: pick_first_token(x, ASM_RE))

        # Prefer explicit assembly_accession first, then RefSeq (GCF), then GenBank (GCA)
        genomes["assembly_accession_best"] = genomes["assembly_accession"].fillna("").astype(str).str.strip()
        mask_empty = genomes["assembly_accession_best"].eq("")
        genomes.loc[mask_empty, "assembly_accession_best"] = genomes.loc[mask_empty, "refseq_assembly_accession"]
        mask_empty = genomes["assembly_accession_best"].eq("")
        genomes.loc[mask_empty, "assembly_accession_best"] = genomes.loc[mask_empty, "genbank_assembly_accession"]

        genomes.to_parquet(genomes_file, index=False)
        stat("Saved genomes to cache", True, rows=len(genomes), expected=exp_g)

    stat(
        "Genomes stats", True,
        total=len(genomes),
        biosample_nonempty=int((genomes["biosample_accession"].astype(str).str.len() > 0).sum()),
        taxon_ids=int(genomes["taxon_id"].dropna().nunique())
    )

    # 2) Genome filters
    before = len(genomes)

    if "genome_status" in genomes.columns:
        gs_before = len(genomes)
        genomes["genome_status_norm"] = genomes["genome_status"].fillna("").astype(str).str.strip().str.lower()
        genomes = genomes[genomes["genome_status_norm"].isin({"complete", "wgs"})].copy()
        stat("Filter genome_status in {Complete,WGS}", True, before=gs_before, after=len(genomes), removed=gs_before - len(genomes))
    else:
        warn("genome_status missing -> skipping Complete/WGS filter")

    if "genome_quality" in genomes.columns:
        gq_before = len(genomes)
        genomes["genome_quality_norm"] = genomes["genome_quality"].fillna("").astype(str).str.strip().str.lower()
        genomes = genomes[genomes["genome_quality_norm"].eq("good")].copy()
        stat("Filter genome_quality == Good", True, before=gq_before, after=len(genomes), removed=gq_before - len(genomes))
    elif "genome_quality_flags" in genomes.columns:
        gqf_before = len(genomes)
        genomes["genome_quality_flags_norm"] = genomes["genome_quality_flags"].fillna("").astype(str).str.strip().str.lower()
        genomes = genomes[(genomes["genome_quality_flags_norm"].eq("")) | (genomes["genome_quality_flags_norm"].str.contains("good"))].copy()
        stat("Filter genome_quality_flags ~ good/empty", True, before=gqf_before, after=len(genomes), removed=gqf_before - len(genomes))
    else:
        warn("genome_quality + flags missing -> skipping Good filter")

    stat("Genomes after filters", True, before=before, after=len(genomes), removed=before - len(genomes))

    taxon_ids = sorted([int(x) for x in genomes["taxon_id"].dropna().unique().tolist()])
    stat("Taxon IDs to process", True, n=len(taxon_ids), sample=taxon_ids[:10])

    # 3) AMR stage (lab-only)
    stat("AMR stage start", True, taxon_ids=len(taxon_ids), print_every_n=PRINT_EVERY_N_TAXON_IDS)

    AMR_FIELDS = [
        "genome_id", "antibiotic", "resistant_phenotype",
        "laboratory_typing_method", "testing_standard", "measurement_value",
        "evidence", "source"
    ]

    amr_sum_parts = []
    taxon_summary_rows = []
    total_taxa = len(taxon_ids)

    taxa_with_any_lab = 0
    taxa_all_dropped = 0
    total_lab_rows_kept = 0

    for idx, tid in enumerate(taxon_ids, start=1):
        do_print_taxon = (idx == 1) or (idx == total_taxa) or (PRINT_EVERY_N_TAXON_IDS and idx % PRINT_EVERY_N_TAXON_IDS == 0)

        pq, meta = amr_taxon_cache_paths(tid)

        if os.path.exists(pq) and os.path.exists(meta):
            with open(meta, "r", encoding="utf-8") as f:
                m = json.load(f)
            if m.get("complete") is True:
                df = pd.read_parquet(pq)
                amr_sum_parts.append(df)

                kept_lab = int(m.get("lab_only_kept", 0))
                total_lab_rows_kept += kept_lab
                if kept_lab > 0:
                    taxa_with_any_lab += 1
                else:
                    taxa_all_dropped += 1

                taxon_summary_rows.append({
                    "taxon_id": tid,
                    "taxon_idx": idx,
                    "expected_raw": m.get("expected_raw"),
                    "fetched_raw": m.get("fetched_raw"),
                    "lab_only_kept": kept_lab,
                    "from_cache": True
                })

                stat("Taxon AMR (cache)", force_print=do_print_taxon,
                     taxon_id=tid, taxon_idx=idx, taxon_total=total_taxa,
                     lab_kept=kept_lab, amr_summary_rows=len(df))
                continue

        stat("Taxon AMR (fetch)", force_print=do_print_taxon, taxon_id=tid, taxon_idx=idx, taxon_total=total_taxa)

        rql_amr = f"eq(taxon_id,{tid})"
        amr, expected = bvbrc_fetch_all(
            endpoint="genome_amr",
            rql_core=rql_amr,
            select_fields=AMR_FIELDS,
            max_rows=MAX_BVBRC_AMR_RECORDS_PER_TAXON,
            label=f"taxon_id={tid}"
        )

        fetched = len(amr)

        for c in AMR_FIELDS:
            if c in amr.columns:
                amr[c] = amr[c].fillna("").astype(str).str.strip()

        if len(amr) > 0:
            m_lab = lab_only_mask(amr)
            amr_lab = amr[m_lab].copy()
        else:
            amr_lab = amr.copy()

        kept_lab = len(amr_lab)
        total_lab_rows_kept += kept_lab
        if kept_lab > 0:
            taxa_with_any_lab += 1
        else:
            taxa_all_dropped += 1

        if not amr_lab.empty:
            amr_lab["pheno"] = amr_lab["resistant_phenotype"].apply(norm_pheno)
            amr_lab = amr_lab[amr_lab["pheno"] != ""].copy()

        if amr_lab.empty:
            amr_sum = pd.DataFrame(columns=[
                "genome_id", "amr_antibiotics_count", "amr_resistant_count", "amr_susceptible_count", "amr_intermediate_count",
                "amr_conflict_count", "amr_conflict_antibiotics", "antibiotics_resistant", "antibiotics_susceptible", "antibiotics_intermediate"
            ])
        else:
            ab = (
                amr_lab.groupby(["genome_id", "antibiotic"])["pheno"]
                .apply(lambda s: "".join(sorted(set(s))))
                .reset_index()
            )
            amr_sum = summarise_amr_per_genome(ab)

        amr_sum.to_parquet(pq, index=False)
        with open(meta, "w", encoding="utf-8") as f:
            json.dump({
                "taxon_id": tid,
                "expected_raw": expected,
                "fetched_raw": fetched,
                "lab_only_kept": int(kept_lab),
                "complete": True
            }, f)

        amr_sum_parts.append(amr_sum)

        taxon_summary_rows.append({
            "taxon_id": tid,
            "taxon_idx": idx,
            "expected_raw": expected,
            "fetched_raw": fetched,
            "lab_only_kept": kept_lab,
            "from_cache": False
        })

        stat("Taxon AMR done", force_print=do_print_taxon,
             taxon_id=tid, taxon_idx=idx, expected=expected, fetched=fetched,
             lab_kept=kept_lab, amr_summary_rows=len(amr_sum))

    # Save taxon summary inside cache
    taxon_summary_csv_cache = cache_path(f"{species.replace(' ', '_')}_amr_taxon_summary.csv")
    pd.DataFrame(taxon_summary_rows).to_csv(taxon_summary_csv_cache, index=False)

    amr_sum = pd.concat(amr_sum_parts, ignore_index=True) if amr_sum_parts else pd.DataFrame(columns=["genome_id"])
    dup = int(amr_sum["genome_id"].duplicated().sum()) if not amr_sum.empty else 0

    stat("AMR stage summary", True,
         taxon_ids=total_taxa,
         taxa_with_any_lab=taxa_with_any_lab,
         taxa_all_dropped=taxa_all_dropped,
         total_lab_rows_kept=total_lab_rows_kept,
         amr_summary_rows=len(amr_sum),
         duplicate_genome_ids=dup)

    if dup > 0:
        warn("Duplicate genome_id in AMR summaries. Keeping first.")
        amr_sum = amr_sum.drop_duplicates("genome_id", keep="first").copy()

    if ONLY_GENOMES_WITH_AMR and not amr_sum.empty:
        before_g = len(genomes)
        genomes = genomes[genomes["genome_id"].isin(set(amr_sum["genome_id"]))].copy()
        stat("Filtered genomes to ONLY_GENOMES_WITH_AMR (lab-only)", True,
             before=before_g, after=len(genomes), removed=before_g - len(genomes))

    # 4) STRONG SRA + ENA BioSample map (cached)
    genomes["biosample_accession"] = genomes["biosample_accession"].apply(lambda x: pick_first_token(x, BS_TOKEN_RE))
    nonempty_bs = int((genomes["biosample_accession"].astype(str).str.len() > 0).sum())
    stat("BioSample diagnostic (post-AMR filter)", True, genomes=len(genomes), biosample_nonempty=nonempty_bs)

    bs_map_file = cache_path("biosample_run_map_SRA_ENA_STRONG.parquet")
    if os.path.exists(bs_map_file):
        bs_map = pd.read_parquet(bs_map_file)
        stat("Loaded STRONG SRA+ENA BioSample map from cache", True, rows=len(bs_map))
    else:
        stat("Building STRONG SRA+ENA BioSample map (AMR-targeted)", True)
        bs_map = run_map_biosample_only_sra_ena(
            genomes_amr_filtered=genomes,
            cache_dir=cache_dir,
            max_biosamples=MAX_BIOSAMPLES_PER_SPECIES,
            max_sra_records=MAX_SRA_RECORDS
        )
        stat("Saved STRONG SRA+ENA BioSample map", True, rows=len(bs_map))

    # 5) Join + export
    out = genomes.merge(bs_map, on="biosample_accession", how="left")

    out["rnaseq_available"] = out["rnaseq_run_ids"].fillna("").str.contains(RUN_RE)
    out["evidence_tier"] = out["rnaseq_available"].map(lambda v: "Strong (same BioSample: SRA+ENA)" if v else "None")
    out["access_method"] = out["rnaseq_available"].map(
        lambda v: "Direct BioSample match → SRA Toolkit (prefetch/fasterq-dump) or ENA FTP FASTQ" if v else ""
    )
    out["rnaseq_run_count"] = out["rnaseq_run_ids"].fillna("").apply(lambda s: len(extract_run_ids(s)))

    stat("Post-join RNA-seq availability (STRONG ONLY)", True, rnaseq_true=int(out["rnaseq_available"].sum()), total=len(out))

    final = out.merge(amr_sum, on="genome_id", how="left")

    final_cols = [
        "genome_id", "genome_name", "strain", "taxon_id",
        "biosample_accession", "bioproject_accession",
        "assembly_accession", "genbank_accessions", "refseq_accessions", "sra_accession",
        "genbank_assembly_accession", "refseq_assembly_accession", "assembly_accession_best",
        "rnaseq_available", "evidence_tier", "rnaseq_run_count",
        "rnaseq_run_ids", "rnaseq_run_ids_sra", "rnaseq_run_ids_ena", "ena_fastq_ftp", "ena_fastq_md5",
        "access_method",
        "amr_antibiotics_count", "amr_resistant_count", "amr_susceptible_count", "amr_intermediate_count",
        "amr_conflict_count", "amr_conflict_antibiotics",
        "antibiotics_resistant", "antibiotics_susceptible", "antibiotics_intermediate"
    ]

    for c in final_cols:
        if c not in final.columns:
            final[c] = ""

    final_table = final[final_cols].copy()
    stat("Final table ready", True, rows=len(final_table), rnaseq_true=int(final_table["rnaseq_available"].sum()))

    os.makedirs(out_root, exist_ok=True)
    out_csv = os.path.join(out_root, f"{species.replace(' ', '_')}_LABONLY_conservative_rnaseq_amr.csv")
    dbg_csv = os.path.join(out_root, f"{species.replace(' ', '_')}_debug_stats.csv")
    tax_csv = os.path.join(out_root, f"{species.replace(' ', '_')}_amr_taxon_summary.csv")

    final_table.to_csv(out_csv, index=False)
    pd.DataFrame(stats).to_csv(dbg_csv, index=False)
    shutil.copyfile(taxon_summary_csv_cache, tax_csv)

    stat("Saved outputs", True, out_csv=out_csv, dbg_csv=dbg_csv, taxon_summary_csv=tax_csv)
    return out_csv, dbg_csv, tax_csv


def parse_species_list(species_arg: str, species_file: Optional[str]) -> List[str]:
    if species_file:
        with open(species_file, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f.readlines()]
        return [ln for ln in lines if ln and not ln.startswith("#")]

    if species_arg:
        return [s.strip() for s in species_arg.split(",") if s.strip()]

    return [s.strip() for s in SPECIES_INPUTS_DEFAULT.split(",") if s.strip()]


def zip_outputs(out_root: str, paths: List[str], zip_path: str) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for f in paths:
            if os.path.exists(f):
                z.write(f, arcname=os.path.relpath(f, out_root))


def cleanup_files(paths: List[str]) -> None:
    seen = set()
    for p in paths:
        if not p or p in seen:
            continue
        seen.add(p)
        try:
            if os.path.exists(p):
                os.remove(p)
        except Exception as e:
            print(f"[WARN] Could not delete {p}: {e}")


def build_combined_rnaseq_true_csv(species_csvs: List[tuple], out_root: str) -> str:
    # Always create combined CSV for multi-species runs, even if 0 rnaseq_true rows.
    keep_frames = []
    template_cols = None

    for sp, path in species_csvs:
        df = pd.read_csv(path)
        if template_cols is None:
            template_cols = ["species"] + df.columns.tolist()

        rn = df["rnaseq_available"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
        df_true = df[rn].copy()

        if not df_true.empty:
            df_true.insert(0, "species", sp)
            keep_frames.append(df_true)

        print(f"[Combined] {sp} | total={len(df)} | rnaseq_true={int(rn.sum())}")

    if keep_frames:
        combined = pd.concat(keep_frames, ignore_index=True)
    else:
        combined = pd.DataFrame(columns=template_cols if template_cols else ["species"])

    combined_true_csv = os.path.join(out_root, f"rnaseq_true_combined_{int(time.time())}.csv")
    combined.to_csv(combined_true_csv, index=False)
    print("Saved combined rnaseq_true CSV:", combined_true_csv, "| rows:", len(combined))
    return combined_true_csv


def main() -> None:
    ap = argparse.ArgumentParser(description="Conservative RNA-seq + AMR checker (local)")
    ap.add_argument("--species", type=str, default="", help="Comma-separated species list")
    ap.add_argument("--species-file", type=str, default="", help="Text file with one species per line")
    ap.add_argument("--cache-root", type=str, default="./cache", help="Cache directory")
    ap.add_argument("--out-root", type=str, default="./outputs", help="Output directory")
    args = ap.parse_args()

    species_file = args.species_file.strip() or None
    species_list = parse_species_list(args.species, species_file)
    print("Species list:", species_list)

    os.makedirs(args.cache_root, exist_ok=True)
    os.makedirs(args.out_root, exist_ok=True)

    all_outputs: List[str] = []
    species_csvs: List[tuple] = []

    for s in species_list:
        print("\n" + "=" * 90)
        print("RUNNING SPECIES:", s)
        print("=" * 90)

        try:
            out_csv, dbg_csv, tax_csv = run_one_species(s, cache_root=args.cache_root, out_root=args.out_root)
        except Exception as e:
            print(f"[!!] FAILED species='{s}' | error={type(e).__name__}: {str(e)[:300]}")
            continue

        all_outputs.extend([out_csv, dbg_csv, tax_csv])
        species_csvs.append((norm_species_name(s), out_csv))

    # Multi-species behaviour only if 2+ species succeeded
    if len(species_csvs) > 1:
        print("\n" + "=" * 90)
        print("MULTI-SPECIES RUN → ZIP ONLY (no loose CSVs kept)")
        print("=" * 90)

        combined_true_csv = build_combined_rnaseq_true_csv(species_csvs, args.out_root)
        all_outputs.append(combined_true_csv)
        
        ts = time.strftime("%Y%m%d_%H%M%S")
        zip_path = os.path.join(args.out_root, f"rnaseq_amr_outputs_{ts}.zip")

        zip_outputs(args.out_root, all_outputs, zip_path)
        print("Zipped outputs:", zip_path)

        # Cleanup after zipping: remove per-species outputs + combined csv
        cleanup_files(all_outputs)
        print("Cleaned up temporary outputs. Only ZIP remains for this run.")

    elif len(species_csvs) == 1:
        print("\nSingle species run → CSV outputs only (no zip, no combined file).")
    else:
        print("\nNo species succeeded → no outputs created.")

    print("\nDONE.")
    print("Outputs folder:", os.path.abspath(args.out_root))
    print("Cache folder:", os.path.abspath(args.cache_root))


if __name__ == "__main__":
    main()
