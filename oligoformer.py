import subprocess
import pathlib
import sys
import shutil
import time
import warnings
import json
import os


# --- Configuração ---
#
# OUTPUT_DIR fica dentro do próprio repositório (pasta onde este arquivo está), então
# nenhuma configuração é necessária: a pasta é criada automaticamente na primeira execução.
#
# OLIGOFORMER_DIR e OLIGOFORMER_CONDA_ENV dependem de uma instalação externa do OligoFormer
# (https://github.com/) que não faz parte deste repositório, e por isso são lidas de variáveis
# de ambiente em vez de hard-coded. Veja o README para instruções de configuração.

_MODULE_DIR = pathlib.Path(__file__).resolve().parent
OUTPUT_DIR = _MODULE_DIR / "oligoformer_results"


def _get_oligoformer_dir() -> pathlib.Path:
    """Lê o diretório de instalação do OligoFormer a partir da variável de ambiente
    OLIGOFORMER_DIR. Só é chamada pelas funções que efetivamente precisam rodar o modelo
    (run/run_specific) — funções que só leem resultados já salvos não exigem essa variável."""
    val = os.environ.get("OLIGOFORMER_DIR")
    if not val:
        raise EnvironmentError(
            "A variável de ambiente OLIGOFORMER_DIR não está definida. Ela deve apontar para "
            "o diretório onde o OligoFormer foi clonado/instalado. Veja o README para instruções."
        )
    return pathlib.Path(val)


def _get_conda_env_name() -> str:
    """Lê o nome do ambiente conda do OligoFormer a partir da variável de ambiente
    OLIGOFORMER_CONDA_ENV."""
    val = os.environ.get("OLIGOFORMER_CONDA_ENV")
    if not val:
        raise EnvironmentError(
            "A variável de ambiente OLIGOFORMER_CONDA_ENV não está definida. Ela deve conter "
            "o nome do ambiente conda com as dependências do OligoFormer instaladas. Veja o "
            "README para instruções."
        )
    return val


def run(fasta_path, silent=False, override=True, model_silent=False):
    oligoformer_dir = _get_oligoformer_dir()
    fasta_path = pathlib.Path(fasta_path).resolve()
    #path final onde colocará os resultados
    final_dir = pathlib.Path(OUTPUT_DIR) / pathlib.Path(fasta_path).stem #o .stem tira apenas o último ponto e o que vem depois 
    result_dir = oligoformer_dir / "result"

    with fasta_path.open() as f: #abre o fasta e verifica rapidamente
        line = f.readline().strip()
        if(line == ""):
            raise ValueError ("Error: Input file empty")
        if (line[0] != '>'):
            raise ValueError ("Error: Input file is not fasta")

    old_dir =  set(result_dir.glob("*")) #o que tinha antes da run
    start = time.time()
    


    #se na pasta resultados já tem um diretório como esse, é porque foi gerado pelo mesmo arquivo fasta, logo overrite é recomendado (deleta e refaz)
    if final_dir.exists() and final_dir.is_dir():
        if override:
            shutil.rmtree(final_dir)
        else:
            if(not silent):
                print("Results already exist in the directory, overriding has been disabled.")
                print(":::Showing stored results:::")
                show_results(final_dir.name)
                print("\n:::Result End:::")
            return


    if(not silent):
        print("Text feedback and outputs from the model are buffered! \nplease wait for inference to end or an error message to appear, the code is not looping.\n")


    #o comando do oligoformer (só deixo esses args mesmo)
    cmd = [
        "conda", "run", "-n", _get_conda_env_name(),
        "python", "scripts/main.py", "--infer", "1", "-i1", str(fasta_path), 
    ]

    result = subprocess.run(
            cmd,
            cwd=oligoformer_dir,   # precisa ser o diretório do OligoFormer para que "scripts/main.py" seja encontrado
            text=True,
            capture_output=(silent or model_silent)
        )

    if result.returncode != 0: #em caso de erro no subprocesso 
            print(result.stderr)
            raise RuntimeError("OligoFormer failed")
    if(not silent):
        print("Inference ended successfully.")
    final_dir.mkdir(parents=True, exist_ok=True)

    #pega os arquivos de resultados no path oligoformer/result e coloca numa lista os que tem id correspondente
    after = set(result_dir.glob("*"))
    files = [f for f in (after - old_dir) if f.is_file()]

    if (len(files) == 0):
        files = [f for f in result_dir.glob("*") if f.is_file() and f.stat().st_mtime > start]
        warnings.warn("No new files found after run, falling back to timestamps for filtering, please delete files on the oligoformer result directory that correlate with this run to avoid potential errors")
        if(len(files) == 0):
            raise ValueError("Could not find new files after falling back to timestamps")

    folders = [f.name.replace("_ranked_filtered.txt", "") for f in files if f.name.endswith("_ranked_filtered.txt")] #pega um de cada header file e extrai o nome
    for f in folders:
        x = final_dir / pathlib.Path(f) #dentro do path final cria uma pasta para cada subsequência
        x.mkdir()
        mv = [v for v in files if v.name.startswith(f)]
        for v in mv:
            new_path = shutil.move(v, x)
            oligoformer_to_json(new_path)

    files = list(final_dir.glob("*"))
    if (not silent):
        print("results saved in the following files:\n")
        for f in files:
            print(f)
    return files


def run_specific(fasta_path, sirna_fasta_path, silent=False, override=True, model_silent=False):
    #Option 2 do OligoFormer: -i1 (mRNA) + -i2 (siRNAs específicas a prever)
    #Estrutura idêntica ao run() normal, só muda o cmd (adiciona -i2) e o nome da pasta final,
    #que combina os stems dos dois fastas para não colidir com resultados do run normal
    fasta_path = pathlib.Path(fasta_path).resolve()
    sirna_fasta_path = pathlib.Path(sirna_fasta_path).resolve()
    #path final onde colocará os resultados (mRNA_stem__siRNA_stem, combinado para evitar colisão com run() normal)
    oligoformer_dir = _get_oligoformer_dir()
    final_dir = pathlib.Path(OUTPUT_DIR) / f"{fasta_path.stem}__{sirna_fasta_path.stem}"
    result_dir = oligoformer_dir / "result"

    with fasta_path.open() as f: #abre o fasta do mRNA e verifica rapidamente
        line = f.readline().strip()
        if(line == ""):
            raise ValueError ("Error: Input file empty")
        if (line[0] != '>'):
            raise ValueError ("Error: Input file is not fasta")

    with sirna_fasta_path.open() as f: #abre o fasta das siRNAs e verifica rapidamente
        line = f.readline().strip()
        if(line == ""):
            raise ValueError ("Error: siRNA input file empty")
        if (line[0] != '>'):
            raise ValueError ("Error: siRNA input file is not fasta")

    old_dir = set(result_dir.glob("*")) #o que tinha antes da run
    start = time.time()

    #se na pasta resultados já tem um diretório como esse, é porque foi gerado pelo mesmo par de arquivos, logo overrite é recomendado (deleta e refaz)
    if final_dir.exists() and final_dir.is_dir():
        if override:
            shutil.rmtree(final_dir)
        else:
            if(not silent):
                print("Results already exist in the directory, overriding has been disabled.")
                print(":::Showing stored results:::")
                show_results(final_dir.name)
                print("\n:::Result End:::")
            return

    if(not silent):
        print("Text feedback and outputs from the model are buffered! \nplease wait for inference to end or an error message to appear, the code is not looping.\n")

    #o comando do oligoformer, Option 2: mRNA + siRNAs específicas
    cmd = [
        "conda", "run", "-n", _get_conda_env_name(),
        "python", "scripts/main.py", "--infer", "1", "-i1", str(fasta_path), "-i2", str(sirna_fasta_path),
    ]

    result = subprocess.run(
            cmd,
            cwd=oligoformer_dir,   # precisa ser o diretório do OligoFormer para que "scripts/main.py" seja encontrado
            text=True,
            capture_output=(silent or model_silent)
        )

    if result.returncode != 0: #em caso de erro no subprocesso
            print(result.stderr)
            raise RuntimeError("OligoFormer failed")
    if(not silent):
        print("Inference ended successfully.")
    final_dir.mkdir(parents=True, exist_ok=True)

    #pega os arquivos de resultados no path oligoformer/result e coloca numa lista os que tem id correspondente
    after = set(result_dir.glob("*"))
    files = [f for f in (after - old_dir) if f.is_file()]

    if (len(files) == 0):
        files = [f for f in result_dir.glob("*") if f.is_file() and f.stat().st_mtime > start]
        warnings.warn("No new files found after run, falling back to timestamps for filtering, please delete files on the oligoformer result directory that correlate with this run to avoid potential errors")
        if(len(files) == 0):
            raise ValueError("Could not find new files after falling back to timestamps")

    folders = [f.name.replace("_ranked_filtered.txt", "") for f in files if f.name.endswith("_ranked_filtered.txt")] #pega um de cada header file e extrai o nome
    for f in folders:
        x = final_dir / pathlib.Path(f) #dentro do path final cria uma pasta para cada subsequência
        x.mkdir()
        mv = [v for v in files if v.name.startswith(f)]
        for v in mv:
            new_path = shutil.move(v, x)
            oligoformer_to_json(new_path)

    files = list(final_dir.glob("*"))
    if (not silent):
        print("results saved in the following files:\n")
        for f in files:
            print(f)
    return files


def oligoformer_to_json(input_file, output_file=None):
    input_file = pathlib.Path(input_file)

    results = []

    with input_file.open() as f:
        header = f.readline().strip().split("\t")

        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split("\t")
            row = dict(zip(header, fields))

            entry = {
                "sirna": row["siRNA"],
                "mrna_segment": row["sense"],
                "position": int(row["pos"]),
                "efficacy": float(row["efficacy"]),
                "filters": {
                    "func_filter": int(row["func_filter"]),
                    "filter": int(row["filter"])
                },
                "source": "oligoformer",
                "origin_file": input_file.name
            }

            results.append(entry)

    # write to file if requested
    if output_file:
        output_file = pathlib.Path(output_file)
    else:
        output_file = input_file.with_suffix(".json")

    with output_file.open("w") as f:
        json.dump(results, f, indent=2)

    return results


def load_oligoformer_json(paths: tuple | str, ranked :bool = True, filtered :bool = True) -> list[dict]:
    if(not(OUTPUT_DIR.exists() and OUTPUT_DIR.is_dir())):
        raise ValueError ("No directory for output detected")
    if(not ranked and filtered):
        raise KeyError ("No result files are filtered but not ranked, please use neither or ranked only or both")
    json_files = []

    if(type(paths) == str):
        if(paths.endswith(".fasta")):
            paths = paths.replace(".fasta", "")
        l = list(OUTPUT_DIR.iterdir())
        paths = OUTPUT_DIR / paths
        if(paths in l):
            temp = list(paths.rglob("*.json*"))
            for j in temp:
                if(ranked and filtered):
                    if(j.name.endswith("_ranked_filtered.json")): 
                        json_files.append(j)
                        continue
                if(ranked and not filtered):
                    if(j.name.endswith("_ranked.json")): 
                        json_files.append(j)
                        continue
                if(not ranked and not filtered):
                    if(j.name.endswith(".json") and (not("ranked" in j.name))): 
                        json_files.append(j)
                        continue
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
            temp = list(seq.glob("*.json*"))
            for j in temp:
                if(ranked and filtered):
                    if(j.name.endswith("_ranked_filtered.json")): 
                        json_files.append(j)
                        continue
                if(ranked and not filtered):
                    if(j.name.endswith("_ranked.json")): 
                        json_files.append(j)
                        continue
                if(not ranked and not filtered):
                    if(j.name.endswith(".json") and (not("ranked" in j.name))): 
                        json_files.append(j)
                        continue
        else:
            raise ValueError("path not found in output directory")
    results = []
    for file in json_files:
        with file.open("r", encoding="utf-8") as f:
            results.append(json.load(f))
    return results


def show_results(name: str | int | tuple = ""):
    def print_table(items, title):
        print(f"\n📂 {title}")
        print("#" * 60)
        print(f"{'Index':<6} | {'Name'}")
        print("-" * 60)
        for i, r in enumerate(items, 1):
            print(f"{i:<6} | {r.name}")
        print("#" * 60)

    # 🔹 nível 1 → FASTAS
    if name == "":
        results = sorted([p for p in OUTPUT_DIR.iterdir() if p.is_dir()])
        print_table(results, "Oligoformer Results")
        return results

    # 🔹 navegação por tupla (multi-nível)
    if isinstance(name, tuple):
        path = OUTPUT_DIR

        try:
            for idx in name:
                items = sorted([p for p in path.iterdir() if p.exists()])
                path = items[idx - 1]

            if path.is_dir():
                items = sorted(path.iterdir())
                print_table(items, path.name)
                return items
            else:
                print(f"\n📄 File: {path.name}")
                return [path]

        except Exception:
            print("Error: invalid navigation path")
            return

    # 🔹 nível 2 → FASTA → HEADERS
    if isinstance(name, int):
        results = sorted([p for p in OUTPUT_DIR.iterdir() if p.is_dir()])

        if name <= 0 or name > len(results):
            print("Error: type a valid numerical index")
            return

        selected = results[name - 1]

        headers = sorted([p for p in selected.iterdir() if p.is_dir()])
        print_table(headers, selected.name)
        return headers

    # 🔹 por nome (mantém compatibilidade)
    else:
        name = name.strip()

        # caso seja fasta
        if name.endswith(".fasta") or name.endswith(".fa"):
            try:
                fasta = pathlib.Path(name).resolve()
            except:
                print("Error: original file not found")
                return

            path = OUTPUT_DIR / fasta.stem

        else:
            matches = list(OUTPUT_DIR.glob(name))
            if len(matches) == 1:
                path = matches[0]
            elif len(matches) == 0:
                print("Error: no results found")
                return
            else:
                print("Error: multiple matches found")
                return

        if path.is_dir():
            items = sorted(path.iterdir())
            print_table(items, path.name)
            return items
        else:
            print("Error: not a directory")
            return

def show_file(filename: str, n: int = 10):
    matches = list(OUTPUT_DIR.rglob(filename))

    if len(matches) == 0:
        print("Error: file not found")
        return
    if len(matches) > 1:
        print("Error: multiple files found:")
        for m in matches:
            print(m)
        return

    path = matches[0]

    print(f"\nFirst {n} lines of {path.name} 📄:")
    print("-" * 60)

    with path.open() as f:
        for i, line in enumerate(f):
            if i >= n:
                print("...")
                break
            print(line.rstrip())

    print("-" * 60)
    return path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: oligoformer <fasta_file>")
        sys.exit(1)

    fasta = sys.argv[1]
    results = run(fasta)