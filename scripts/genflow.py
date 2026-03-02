"""
GenFlow - Phylogenomic analysis pipeline launcher
Genome analysis for bacteria and fungi with Snakemake
"""
import argparse
import os
import re
import shlex
import subprocess
import sys
import unicodedata
from pathlib import Path

def sanitize(text: str) -> str:
    """Convert text to ASCII alphanumeric, replacing special characters with underscores."""
    normalized = unicodedata.normalize('NFKD', text)
    ascii_text = normalized.encode('ascii', 'ignore').decode('ascii')
    ascii_text = re.sub(r'[^A-Za-z0-9]+', '_', ascii_text)
    ascii_text = re.sub(r'_+', '_', ascii_text).strip('_')
    return ascii_text or 'Unknown'

def edirect_organism_for_accession(acc: str, edirect_dir: str) -> str:
    """Query NCBI assembly/nucleotide database via EDirect for organism name from accession."""
    env = os.environ.copy()
    env['PATH'] = f"{edirect_dir}:{env.get('PATH', '')}"
    try:
        # Try assembly database first (faster)
        print(f"  [EDirect] Querying assembly database for {acc}...", file=sys.stderr)
        cmd = f"esearch -db assembly -query {acc} | esummary | xtract -pattern DocumentSummary -element Organism"
        org = subprocess.check_output(cmd, shell=True, env=env, stderr=subprocess.DEVNULL, text=True).strip()
        if org:
            print(f"  [EDirect] Found organism: {org}", file=sys.stderr)
            return org
        
        # Fallback to nucleotide database
        print(f"  [EDirect] Assembly not found, querying nucleotide database for {acc}...", file=sys.stderr)
        cmd2 = f"esearch -db nucleotide -query {acc} | esummary | xtract -pattern DocumentSummary -element Organism"
        org2 = subprocess.check_output(cmd2, shell=True, env=env, stderr=subprocess.DEVNULL, text=True).strip()
        if org2:
            print(f"  [EDirect] Found: {org2}", file=sys.stderr)
            return org2
        print(f"  [EDirect] No organism found for {acc}", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"  [EDirect] Error querying {acc}: {e}", file=sys.stderr)
    return ""

def derive_accession_from_filename(name: str) -> str:
    """Extract NCBI accession code (GCF_/GCA_ format) from filename if present."""
    m = re.search(r"(GC[AF]_[0-9]+\.[0-9]+)", name)
    if m:
        return m.group(1)
    m2 = re.search(r"(GC[AF]_[0-9]+)", name)
    if m2:
        return m2.group(1)
    return ""

def build_name_map(intermediate: Path, edirect_dir: str) -> dict:
    """Map original FASTA filenames to sanitized organism names using EDirect lookups."""
    files = []
    for ext in ('*.fasta', '*.fa', '*.fna'):
        files.extend(intermediate.glob(ext))
    files = sorted({f.resolve() for f in files if f.exists()})

    if not files:
        print("Warning: No sequence files found in Intermediate", file=sys.stderr)

    print(f"\n[GenFlow] Building name map with EDirect", file=sys.stderr)
    print(f"[GenFlow] Found {len(files)} sequence file(s)", file=sys.stderr)
    
    name_map = {}
    for i, f in enumerate(files, 1):
        stem = f.stem
        acc = derive_accession_from_filename(f.name)
        org = ''
        if acc:
            print(f"[GenFlow] [{i}/{len(files)}] Processing {f.name}", file=sys.stderr)
            org = edirect_organism_for_accession(acc, edirect_dir)
        else:
            print(f"[GenFlow] [{i}/{len(files)}] No accession in {f.name}, using filename", file=sys.stderr)
        
        base_name = org if org else stem
        clean = sanitize(base_name)
        name_map[stem] = clean
        print(f"           → Mapped: {stem} → {clean}", file=sys.stderr)

    print(f"[GenFlow] Name map complete ({len(name_map)} entries)\n", file=sys.stderr)
    return name_map

def rename_fasta_headers(in_fa: Path, mapping: dict, out_fa: Path) -> None:
    """Rename FASTA headers using provided mapping dictionary."""
    with open(in_fa) as inp, open(out_fa, 'w') as out:
        for line in inp:
            if line.startswith('>'):
                hdr = line[1:].strip()
                key = hdr.split()[0]
                new = mapping.get(key, hdr)
                out.write(f'>{new}\n')
            else:
                out.write(line)

def generate_individual_sequence_names(combined_fasta: Path, output_map: Path) -> None:
    """Extract individual sequence identifiers from combined FASTA file (Prokka format)."""
    name_map = {}
    
    with open(combined_fasta) as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                # Parse format: NODE_N_length_XXXX_cov_YYY.YYY
                match = re.search(r'NODE_(\d+)_length_(\d+)_cov_([\d.]+)', header)
                if match:
                    node_num = match.group(1)
                    length = int(match.group(2))
                    cov = float(match.group(3))
                    
                    # Extract sample name (prefix before NODE)
                    sample_match = re.match(r'(\w+)_NODE_', header)
                    sample_name = sample_match.group(1) if sample_match else "UNKNOWN"
                    
                    # Build identifier: SAMPLE_NODENUMBER_LXXXK_CYYY
                    length_kb = length // 1000
                    cov_int = int(cov)
                    new_id = f"{sample_name}_NODE{node_num}_L{length_kb}k_C{cov_int}"
                    name_map[header] = new_id
    
    # Write TSV: original_header → new_name
    with open(output_map, 'w') as f:
        for orig, new in name_map.items():
            f.write(f"{orig}\t{new}\n")
    
    print(f"[GenFlow] Generated {len(name_map)} individual sequence identifiers", file=sys.stderr)

def run_snakemake(cwd: Path, args: argparse.Namespace) -> None:
    """Configure and execute Snakemake phylogenomic pipeline."""
    # Create required directories
    (cwd / 'results').mkdir(exist_ok=True)
    (cwd / 'Intermediate').mkdir(exist_ok=True)
    (cwd / 'logs').mkdir(exist_ok=True)

    print(f"\n{'='*70}")
    print(f"[GenFlow] Preparing pipeline configuration")
    print(f"{'='*70}")

    # Assemble Snakemake configuration
    cfg_parts = {
        'mode': args.mode,
        'use_ncbi_download': str(args.ncbi).lower(),
        'run_ani': str(args.run_ani).lower(),
        'genomes': args.genomes,
        'threads': str(args.threads),
        'G': str(args.G),
        'F': str(args.F),
        'DNA_mode': str(args.dna).lower(),
        'mcl_inflation': str(args.mcl_inflation),
        'edirect_dir': str(cwd / 'scripts' / 'edirect'),
    }
    if getattr(args, 'fasta', None):
        cfg_parts['fasta'] = args.fasta
    cfg_str = ' '.join(f"{k}={shlex.quote(str(v))}" for k, v in cfg_parts.items())

    # Build Snakemake command
    cmd = (
        "snakemake -s scripts/Snakefile "
        "--use-conda "
        f"--cores {args.threads} "
        "--quiet --rerun-incomplete --latency-wait 60 "
        f"--config {cfg_str}"
    )

    # Display configuration
    print(f"[GenFlow] Mode: {args.mode}")
    print(f"[GenFlow] Threads: {args.threads}")
    print(f"[GenFlow] NCBI Download: {args.ncbi}")
    print(f"[GenFlow] DNA Mode: {args.dna}")
    if getattr(args, 'fasta', None):
        print(f"[GenFlow] FASTA files: {args.fasta}")
    print(f"{'='*70}\n")
    
    print(f"\n{'='*70}")
    print(f"[GenFlow] Executing phylogenomic pipeline")
    print(f"{'='*70}\n")
    
    # Execute Snakemake
    subprocess.check_call(cmd, shell=True, cwd=str(cwd))
    print(f"\n{'='*70}")
    print(f"[GenFlow] Pipeline completed successfully!")
    print(f"[GenFlow] Output tree: results/phylogenomic-tree.txt")
    print(f"[GenFlow] Output alignment: results/aligned.named.fasta")
    print(f"{'='*70}\n")

def main():
    """Parse command-line arguments and execute GenFlow pipeline."""
    ap = argparse.ArgumentParser(
        description='GenFlow - Phylogenomic analysis pipeline (bacteria/fungi)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 scripts/genflow.py -m bacteria -t 4 --fasta "Data/f1.fasta,Data/f2.fasta"
  python3 scripts/genflow.py -m fungi -t 4 -d --genomes Data/genomes.txt
  python3 scripts/genflow.py -m bacteria -t 8 --run-ani

Notes:
  - Fungi mode requires minimum 4 genomes
  - Local FASTA: use --fasta with comma-separated paths
  - NCBI download (-d): populate file with GCF_/GCA_ accessions
        """
    )
    
    ap.add_argument(
        '-m', '--mode',
        choices=['bacteria', 'fungi'],
        default='bacteria',
        help='Analysis mode: bacteria (default) or fungi'
    )
    ap.add_argument(
        '-t', '--threads',
        type=int,
        default=8,
        help='Number of threads to use (default: 8)'
    )
    ap.add_argument(
        '-d', '--ncbi',
        action='store_true',
        help='Enable NCBI genome download using Data/genomes.txt'
    )
    ap.add_argument(
        '-N', '--dna',
        action='store_true',
        help='Use DNA mode for tree building (default: protein)'
    )
    ap.add_argument(
        '--run-ani',
        action='store_true',
        help='Run ANI analysis (default: False)'
    )
    ap.add_argument(
        '-I', '--mcl-inflation',
        type=int,
        default=10,
        help='MCL inflation parameter for pangenome clustering (default: 10)'
    )
    ap.add_argument(
        '-G',
        type=float,
        default=0.8,
        dest='G',
        help='Geometric index threshold for core genes (default: 0.8)'
    )
    ap.add_argument(
        '-F',
        type=float,
        default=0.8,
        dest='F',
        help='Functional index threshold for core genes (default: 0.8)'
    )
    ap.add_argument(
        '--fasta',
        type=str,
        default=None,
        help='Comma-separated FASTA files (overrides automatic discovery)'
    )
    ap.add_argument(
        '--genomes',
        type=str,
        default='Data/genomes.txt',
        help='Path to NCBI accessions file for download (default: Data/genomes.txt)'
    )
    
    args = ap.parse_args()

    # Navigate to GenFlow root directory
    cwd = Path(__file__).resolve().parent.parent
    os.chdir(cwd)

    # Validate NCBI download file
    if args.ncbi:
        genomes = Path(args.genomes)
        if not genomes.exists():
            raise SystemExit(f"Error: {args.genomes} not found.")
    
    # Fungi mode requires 4+ genomes
    if args.mode == 'fungi':
        if args.fasta:
            fasta_count = len([f.strip() for f in args.fasta.split(',') if f.strip()])
            if fasta_count < 4:
                raise SystemExit(f"Error: Fungi mode requires 4+ genomes. Found {fasta_count}.")
        elif args.ncbi:
            with open(Path(args.genomes)) as f:
                accessions = [line.strip() for line in f if line.strip()]
            if len(accessions) < 4:
                raise SystemExit(f"Error: Fungi mode requires 4+ genomes. Found {len(accessions)}.")

    # Check EDirect availability
    edirect_dir = Path('scripts/edirect')
    if not edirect_dir.exists():
        print(f"Warning: {edirect_dir} not found. EDirect name resolution will be skipped.", file=sys.stderr)

    # Execute pipeline
    run_snakemake(cwd, args)

if __name__ == '__main__':
    main()