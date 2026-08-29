from __future__ import annotations

import re
from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Literal, Tuple
from urllib.request import urlopen
import sys
import pathlib as ph
import shutil
import json

_MODULE_DIR = ph.Path(__file__).resolve().parent
OUTPUT_DIR = _MODULE_DIR / "dsir_results"

DSIR_MODEL_URLS = {
    "19nt": "https://biodev.cea.fr/DSIR/data/dsir-model-19.txt",
    "21nt": "https://biodev.cea.fr/DSIR/data/dsir-model-21.txt",
}

RNA_COMP = str.maketrans({
    "A": "U",
    "U": "A",
    "T": "A",  # allow DNA input
    "C": "G",
    "G": "C",
})


@dataclass
class DsirCandidate:
    pos: int               # 1-based start position in target
    target_window: str     # sequence in target (RNA alphabet)
    guide: str             # antisense/guide strand, 5'->3'
    passenger: str         # sense/passenger strand, 5'->3'
    score: float

    def string(self): #método que retorna candidato em string
        return str(f"Pos {self.pos} | score   {self.score} | guide {self.guide} | target {self.target_window}")


def normalize_rna(seq: str) -> str:
    seq = seq.upper().replace(" ", "").replace("\n", "").replace("\r", "")
    seq = seq.replace("T", "U")
    bad = set(seq) - set("ACGU")
    if bad:
        raise ValueError(f"Sequence contains invalid characters: {sorted(bad)}")
    return seq


def revcomp_rna(seq: str) -> str:
    seq = normalize_rna(seq)
    return seq.translate(RNA_COMP)[::-1]


DSIR_WEIGHTS_CACHE_DIR = _MODULE_DIR / "dsir_weights_cache"


def load_dsir_weights(model: Literal["19nt", "21nt"] = "21nt") -> Dict[str, float]:
    """
    Busca e faz o parse do arquivo oficial de pesos do DSIR no servidor do CEA.
    Usa um cache local em disco: se o arquivo já foi baixado uma vez, não tenta
    baixar de novo (evita timeouts de rede em runs futuras).
    """
    if model not in DSIR_MODEL_URLS:
        raise ValueError("model must be '19nt' or '21nt'")

    DSIR_WEIGHTS_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = DSIR_WEIGHTS_CACHE_DIR / f"dsir-model-{model.replace('nt', '')}.txt"

    if cache_file.exists():
        text = cache_file.read_text(encoding="utf-8", errors="replace")
    else:
        url = DSIR_MODEL_URLS[model]
        with urlopen(url, timeout=30) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        cache_file.write_text(text, encoding="utf-8")

    # The file is a single line of "Feature : weight" pairs.
    pairs = re.findall(r"([A-Za-z0-9]+)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", text)
    if not pairs:
        raise RuntimeError(f"Could not parse DSIR weights from {'cache file' if cache_file.exists() else url}")

    weights = {k: float(v) for k, v in pairs}
    return weights


def sparse_features(seq: str) -> Dict[str, int]:
    """
    Position-specific one-hot features.
    Example for position 1: A1, C1, G1, U1
    """
    seq = normalize_rna(seq)
    feats: Dict[str, int] = {}
    for i, nt in enumerate(seq, start=1):
        for base in "ACGU":
            feats[f"{base}{i}"] = 1 if nt == base else 0
    return feats


def spectral_features(seq: str) -> Dict[str, int]:
    """
    Counts of 1-mer, 2-mer, and 3-mer motifs across the sequence.
    Example features: A, C, G, U, AA, AC, ..., AAA, AAC, ...
    """
    seq = normalize_rna(seq)
    feats: Dict[str, int] = {}

    # 1-mers
    for base in "ACGU":
        feats[base] = seq.count(base)

    # 2-mers
    for a, b in product("ACGU", repeat=2):
        motif = a + b
        feats[motif] = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == motif)

    # 3-mers
    for a, b, c in product("ACGU", repeat=3):
        motif = a + b + c
        feats[motif] = sum(1 for i in range(len(seq) - 2) if seq[i:i+3] == motif)

    return feats


def dsir_score(seq: str, model: Literal["19nt", "21nt"] = "21nt",
               weights: Dict[str, float] | None = None) -> float:
    """
    Compute the DSIR score for one guide strand sequence.

    This matches the official server's linear model:
    score = Offset + sum(feature_value * feature_weight)

    The score is already in the same units used by DSIR (percent-style efficacy score).
    """
    seq = normalize_rna(seq)
    if model == "19nt" and len(seq) != 19:
        raise ValueError(f"19nt model expects length 19, got {len(seq)}")
    if model == "21nt" and len(seq) != 21:
        raise ValueError(f"21nt model expects length 21, got {len(seq)}")

    if weights is None:
        weights = load_dsir_weights(model)

    feats = {}
    feats.update(sparse_features(seq))
    feats.update(spectral_features(seq))

    score = weights.get("Offset", 0.0)
    for name, value in feats.items():
        score += weights.get(name, 0.0) * value

    return score


def has_poly_run(seq: str, n: int = 4) -> bool:
    seq = normalize_rna(seq)
    return bool(re.search(rf"(A{{{n},}}|C{{{n},}}|G{{{n},}}|U{{{n},}})", seq))


def has_immunostimulatory_motif(seq: str) -> bool:
    seq = normalize_rna(seq)
    motifs = ("UGUGU", "GUCCUUCAA")
    return any(m in seq for m in motifs)


def design_dsir_from_target(
    target_seq: str,
    model: Literal["19nt", "21nt"] = "21nt",
    threshold: float = 90.0,
    filter_poly_runs: bool = True,
    filter_motifs: bool = True,
    weights: Dict[str, float] | None = None,
) -> List[DsirCandidate]:
    """
    Generate and score all DSIR candidates from a target RNA/DNA sequence.

    target_seq:
        Input mRNA/cDNA sequence in 5'->3' orientation.
    model:
        "19nt" or "21nt"
    threshold:
        Keep only candidates with score >= threshold, matching the DSIR default style.
    """
    target_seq = normalize_rna(target_seq)
    n = 19 if model == "19nt" else 21

    if len(target_seq) < n:
        raise ValueError(f"Target sequence must be at least {n} nt long.")

    if weights is None:
        weights = load_dsir_weights(model)

    results: List[DsirCandidate] = []

    for start in range(0, len(target_seq) - n + 1):
        window = target_seq[start:start + n]
        guide = revcomp_rna(window)     # antisense / guide strand, 5'->3'
        passenger = revcomp_rna(guide)  # sense / passenger strand, 5'->3'

        if filter_poly_runs and has_poly_run(guide, 4):
            continue
        if filter_motifs and has_immunostimulatory_motif(guide):
            continue

        score = dsir_score(guide, model=model, weights=weights)
        if score >= threshold:
            results.append(
                DsirCandidate(
                    pos=start + 1,
                    target_window=window,
                    guide=guide,
                    passenger=passenger,
                    score=score,
                )
            )

    results.sort(key=lambda x: x.score, reverse=True)
    return results

def get_targets(fasta_file) -> dict:
    if(type(fasta_file) == str):
        fasta_file = ph.Path(fasta_file).resolve()
    with fasta_file.open() as f:
            l = f.readline()
            target_seq = []
            target_seqs = {}  
            current_header = None  # guarda o header atual

            if('>' not in l):
                raise ValueError("Error: file given is not in fasta format")

            current_header = l.strip()  # salva o primeiro header

            for l in f:
                if('>' in l):
                    # salva a sequência atual com o header correspondente
                    target_seqs[current_header] = ''.join(target_seq)
                    target_seq = []
                    current_header = l.strip()  # novo header
                    l = ""
                target_seq.append(l)

            # salva a última sequência
            if(current_header is not None):
                target_seqs[current_header] = ''.join(target_seq)

            return target_seqs  # sempre dict
    
def run_targets(target, mode : Literal["19nt", "21nt"] = "19nt", threshold : float = 0.90, silent = True) -> dict:
        cds = {}
        for h, seq in target.items():
            cds[h] = design_dsir_from_target(target_seq=seq, model=mode, threshold=threshold)
        if(not silent):
            x = 0
            print("#"*100)
            for h, cd in cds.items():
                
                print(f"Found {len(cd)} candidates for seq {h} (seq {x})")
                for c in cd:
                    print(f"Pos {c.pos:>3} | score {c.score:>6.2f} | guide {c.guide} | target {c.target_window}")
                print("#"*100)
                x += 1
        return cds


def run(fasta_path, mode : Literal["19nt", "21nt"] = "19nt", threshold : float = 0.90, silent = True, override : bool= True):
    if(type(fasta_path) == str):
        fasta_path = ph.Path(fasta_path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    folder = ph.Path(OUTPUT_DIR / fasta_path.name.replace(".fasta", ""))
    if(override and folder.exists() and folder.is_dir()):
        shutil.rmtree(folder) #destrói run passada
    elif(folder.exists() and not folder.is_dir()):
        raise RuntimeError("Fasta name is already present in output dir, but isnt a folder, please clean up the output directory")
    elif(not override and folder.exists() and folder.is_dir()):
        if not silent: print("Folder already exists and override is disalbled.")
        return
    folder.mkdir(parents=True)


    results = run_targets(target=get_targets(fasta_path), mode=mode, threshold=threshold ,silent=silent)
    
    for h, seq in results.items():
        file_dir = ph.Path(folder / h.replace(">", "", 1))
        file_dir.mkdir(parents=True)
        file = file_dir / f"{file_dir.name}.txt"
        if not silent: print(file)
        with file.open(mode="w+") as f:
            f.write(f"Found {len(seq)} candidates for seq {h} \n")
            for s in seq:
                f.write(s.string() + "\n")
        dsir_to_json(file)


def dsir_to_json(input_file, output_file=None):
    input_file = ph.Path(input_file)

    results = []

    with input_file.open() as f:
        lines = f.readlines()

    # --- parse entries ---
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        # split by "|"
        parts = [p.strip() for p in line.split("|")]

        try:
            pos = int(parts[0].replace("Pos", "").strip())
            score = float(parts[1].replace("score", "").strip())
            guide = parts[2].replace("guide", "").strip()
            target = parts[3].replace("target", "").strip()
        except Exception as e:
            # skip malformed lines safely
            continue

        entry = {
            "sirna": guide,
            "mrna_segment": target,
            "position": pos,
            "efficacy": score,
            "source": "dsir",
            "origin_seq": input_file.name
        }

        results.append(entry)

    # --- decide output path ---
    if output_file:
        output_file = ph.Path(output_file)
    else:
        output_file = input_file.with_suffix(".json")

    # --- write file ---
    with output_file.open("w") as f:
        json.dump(results, f, indent=2)

    return results
    

def show_results(k : str = ""):
    if(not(OUTPUT_DIR.exists() and OUTPUT_DIR.is_dir())):
        raise ValueError ("No directory for output detected")
    def print_table(items, title):
        print(f"\n📂 {title}")
        print("#" * 60)
        print(f"{'Index':<6} | {'Name'}")
        print("-" * 60)
        for i, r in enumerate(items, 1):
            print(f"{i:<6} | {r.name}")
        print("#" * 60)

    if(k.endswith(".fasta")):
        k = k.replace(".fasta", "")
    fastas = list(OUTPUT_DIR.iterdir())
    fnames = [f.name for f in fastas]
    
    if(k == ""):
        print_table(fastas, "Dsir output directory")
        return fastas
    if(k in fnames):
        kp = OUTPUT_DIR / k
        l = list(kp.iterdir())
        print_table(l, f"Dsir {k} results")
        return l
    

def load_from_json(paths: tuple | str) -> list[dict]: 
    if(not(OUTPUT_DIR.exists() and OUTPUT_DIR.is_dir())):
        raise ValueError ("No directory for output detected")
    
    json_files = []

    if(type(paths) == str):
        if(paths.endswith(".fasta")):
            paths = paths.replace(".fasta", "")
        l = list(OUTPUT_DIR.iterdir())
        paths = OUTPUT_DIR / paths
        if(paths in l):
            json_files = list(paths.rglob("*.json*"))
        else:
            raise ValueError ("path not found in output directory")
    else:
        fasta = paths[0]
        seq = paths[1]
        if(fasta.endswith(".fasta")):
            fasta = fasta.replace(".fasta", "")
        fasta = OUTPUT_DIR / fasta
        seq = fasta / seq
        if(fasta in list(OUTPUT_DIR.iterdir()) and seq in list(fasta.iterdir())):
            json_files = list(seq.glob("*.json*"))
        else:
            raise ValueError("path not found in output directory")
    results = []
    for file in json_files:
        with file.open("r", encoding="utf-8") as f:
            results.append(json.load(f))
    return results

    
    





# ---- Example usage ----
if __name__ == "__main__":
    target = "AUGCGCGAUCUCGAUGCAUGUGCGAUCGAUGCGUAUCGAUUGCUAGCUAGCUAGCUAGCUA"
    if(len(sys.argv) > 1):
        #assume que o arg extra é um arquivo fasta.
        target = get_targets(sys.argv[1])

    candidates = run_targets(target, silent=False)