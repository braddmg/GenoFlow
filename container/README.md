GenFlow - Phylogenomic Pipeline

Docker
  docker-compose build
  docker-compose run --rm genflow -c "snakemake --cores 8 --use-conda --config fasta='../Data/*.fasta'"
  docker-compose run --rm genflow -i

Singularity
  ./docker-to-singularity.sh convert genflow:latest genflow.sif
  singularity run genflow.sif snakemake --cores 8 --use-conda --config fasta='../Data/*.fasta'
  singularity shell genflow.sif

Local
  conda create -n genflow snakemake=7 python=3.10
  conda activate genflow
  cd ../envs && for f in *.yml; do conda env create -f "$f"; done

Common Workflows
  Bacteria local: --config fasta='../Data/*.fasta'
  Bacteria NCBI: --config use_ncbi_download=true genomes='../Data/genomes.txt'
  Bacteria RAxML: --config fasta='../Data/*.fasta' use_raxml=true
  Bacteria ANI: --config fasta='../Data/*.fasta' run_ani=true
  Fungi: --config mode=fungi genomes='../Data/genomes.txt'

Parameters
  fasta - Input FASTA files
  mode - 'bacteria' or 'fungi'
  use_ncbi_download - Download from NCBI (true/false)
  use_raxml - Use RAxML instead of FastTree
  run_ani - Run ANI analysis

Installed Environments
  genflow-base, genflow-download-ani, genflow-prokka, genflow-roary, genflow-raxml, ufcg

Outputs
  results/phylogenomic-tree.txt
  results/aligned.fasta
  results/ani_heatmap.pdf
  prokka/, roary/, ufcg/, pyANI/

Troubleshooting
  Docker daemon: sudo systemctl start docker (Linux) or open Docker app (macOS)
  Permission denied: sudo usermod -aG docker $USER
  Conda error: docker-compose build --no-cache
