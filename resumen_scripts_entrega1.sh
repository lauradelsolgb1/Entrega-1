
cd /home/lauradelsolgb12/curso_bioinformatica
mkdir analysis
cd analysis
mkdir -p data/{ref,fastq,trimmed} qc scripts
cd data
cd fastq
wget -O data.tar.gz https://osf.io/2jc4a/download
ls
conda activate qc-reads
cd /home/lauradelsolgb12/curso_bioinformatica/analysis/data/fastq
tar -tzf data.tar.gz
tar -xvzf data.tar.gz
cd /home/lauradelsolgb12/curso_bioinformatica/analysis/scripts
conda activate bioinfo
nano /qc_pre.sh
------------------------------------------------------------------------------------------------- #Aquí están los nanos 
#!/usr/bin/env bash
set -euo pipefail

IN="/home/lauradelsolgb12/curso_bioinformatica/analysis/scripts/data"
OUTD="trimmed"
QCD="qc/fastp"
mkdir -p "${OUTD}" "${QCD}" #ya estaban creados los directorios , los voli a crear y despues los eliminé :/ 

# al anc le hice un trimming más suave

fastp \
  -i "${IN}/anc_R1.fastq.gz" -I "${IN}/anc_R2.fastq.gz" \
  -o "${OUTD}/anc_R1.trim.fastq.gz" -O "${OUTD}/anc_R2.trim.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 30 \
  --length_required 50 \
  --trim_front1 0 --trim_front2 0 \
  --cut_front --cut_tail --cut_mean_quality 20 \
  --thread 4 \
  --html "${QCD}/anc_fastp.html" \
  --json "${QCD}/anc_fastp.json"

# al evol1 le hice un trimming más agresivo
fastp \
  -i "${IN}/evol1_R1.fastq.gz" -I "${IN}/evol1_R2.fastq.gz" \
  -o "${OUTD}/evol1_R1.trim.fastq.gz" -O "${OUTD}/evol1_R2.trim.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 35 \
  --length_required 75 \
  --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --thread 4 \
  --html "${QCD}/evol1_fastp.html" \
  --json "${QCD}/evol1_fastp.json"

# evol2 un trimming masomenos 

fastp \
  -i "${IN}/evol2_R1.fastq.gz" -I "${IN}/evol2_R2.fastq.gz" \
  -o "${OUTD}/evol2_R1.trim.fastq.gz" -O "${OUTD}/evol2_R2.trim.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 30 \
  --length_required 75 \
  --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --thread 4 \
  --html "${QCD}/evol2_fastp.html" \
  --json "${QCD}/evol2_fastp.json"

echo "[OK] fastp completado. Salidas en ${OUTD}, reportes en ${QCD}"
-------------------------------------------------------------------------------------------------
chmod qc_pre.sh
./qc_pre.sh
cd /home/lauradelsolgb12/curso_bioinformatica/analysis/scripts
nano 04_qc_post.sh
-------------------------------------------------------------------------------------------------
#!/usr/bin/env bash
set -euo pipefail

IN="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed"
OUT="/home/lauradelsolgb12/curso_bioinformatica/analysis/qc/post"
mkdir -p "${OUT}"

# entarda 
FILES=(
  "${IN}/anc_R1.trim.fastq.gz"
  "${IN}/anc_R2.trim.fastq.gz"
  "${IN}/evol1_R1.trim.fastq.gz"
  "${IN}/evol1_R2.trim.fastq.gz"
  "${IN}/evol2_R1.trim.fastq.gz"
  "${IN}/evol2_R2.trim.fastq.gz"
)

# fastqc en todas las lecturas recortadas
fastqc -o "${OUT}" -t 4 "${FILES[@]}"

# multiqc para integrar reportes (fastp y post-trimming)
multiqc -o "${OUT}" "${OUT}" "qc/fastp" 2>/dev/null || true

echo "[OK] QC post-trimming completado en ${OUT}"
-------------------------------------------------------------------------------------------------
chmod +x 04_qc_post.sh
./04_qc_post.sh
conda activate bioinfo
multiqc -o "/home/lauradelsolgb12/curso_bioinformatica/analysis/qc/post"   "/home/lauradelsolgb12/curso_bioinformatica/analysis/qc/fastp"   --force |& tee /tmp/multiqc_run.log
ls -la /home/lauradelsolgb12/curso_bioinformatica/analysis/qc/post

##Creacion de assemby_env
mamba create -n assembly_env python=3.10 -c conda-forge -c bioconda spades quast -y
conda activate assembly_env
which python
python --version
spades.py --version
quast --version

## ensamblaje 
cd /home/lauradelsolgb12/curso_bioinformatica/analysis
cd /home/lauradelsolgb12/curso_bioinformatica/analysis/scripts
nano ensamblaje01.sh ---que puse en el nano 
-------------------------------------------------------------------------------------------------
#!/usr/bin/env bash
set -euo pipefail

R1="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed/anc_R1.trim.fastq.gz"
R2="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed/anc_R2.trim.fastq.gz"
OUTDIR="/home/lauradelsolgb12/curso_bioinformatica/analysis/assembly_01/anc_spades"

# CORRECCIÓN: sin espacios alrededor del =
THREADS=40         # número de hilos
MEM_GB=160         # memoria en GB para SPAdes

# Ejecutables (se asume que están en PATH)
SPADES_BIN="spades.py"
QUAST_BIN="quast"

mkdir -p "${OUTDIR}"

# 1) Ensamblaje con SPAdes
echo "[INFO] Ejecutando SPAdes..."
"${SPADES_BIN}" -1 "$R1" -2 "$R2" -o "$OUTDIR" -t "$THREADS" -m "$MEM_GB" --careful

# 2) Verificar que se generó contigs.fasta
CONTIGS="${OUTDIR}/contigs.fasta"
if [[ -f "$CONTIGS" ]]; then
  echo "[OK] Ensamblaje terminado: $CONTIGS"
else
  echo "[ERROR] No se generó contigs.fasta en $OUTDIR"
  exit 1
fi

# 3) Evaluación con QUAST (N50 y total length)
if command -v "${QUAST_BIN}" >/dev/null 2>&1; then
  QOUT="${OUTDIR}/quast_report"
  mkdir -p "$QOUT"
  echo "[INFO] Ejecutando QUAST..."
  "${QUAST_BIN}" "$CONTIGS" -o "$QOUT" -t "$THREADS"

  REPORT="${QOUT}/report.tsv"
  if [[ -f "$REPORT" ]]; then
    echo
    echo "=== Métricas extraídas de ${REPORT} ==="
    awk -F'\t' '
      $1 ~ /Total length/ { print "Total length (bp):", $2 }
      $1 == "N50"          { print "N50 (bp):", $2 }
    ' "$REPORT"
  else
    echo "[WARN] QUAST no generó report.tsv en $QOUT"
  fi
else
  echo "[WARN] quast no está en PATH. Instálalo o activa el env correcto."
fi

echo "[DONE] Ensamblaje y evaluación finalizados. Resultados en: $OUTDIR"

-------------------------------------------------------------------------------------------------
chmod +x ensamblaje01.sh
./ensamblaje01.sh
nano ensamblaje_secuencias_crudas.sh
-------------------------------------------------------------------------------------------------
#!/usr/bin/env bash
set -euo pipefail

R1="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/fastq/anc_R1.fastq.gz"
R2="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/fastq/anc_R2.fastq.gz"
OUTDIR="/home/lauradelsolgb12/curso_bioinformatica/analysis/assembly_datos_crudos/anc_spades"

# CORRECCIÓN: sin espacios alrededor del =
THREADS=40         # número de hilos
MEM_GB=160         # memoria en GB para SPAdes

# Ejecutables (se asume que están en PATH)
SPADES_BIN="spades.py"
QUAST_BIN="quast"

mkdir -p "${OUTDIR}"

# 1) Ensamblaje con SPAdes
echo "[INFO] Ejecutando SPAdes..."
"${SPADES_BIN}" -1 "$R1" -2 "$R2" -o "$OUTDIR" -t "$THREADS" -m "$MEM_GB" --careful

# 2) Verificar que se generó contigs.fasta
CONTIGS="${OUTDIR}/contigs.fasta"
if [[ -f "$CONTIGS" ]]; then
  echo "[OK] Ensamblaje terminado: $CONTIGS"
else
  echo "[ERROR] No se generó contigs.fasta en $OUTDIR"
  exit 1
fi

# 3) Evaluación con QUAST (N50 y total length)
if command -v "${QUAST_BIN}" >/dev/null 2>&1; then
  QOUT="${OUTDIR}/quast_report"
  mkdir -p "$QOUT"
  echo "[INFO] Ejecutando QUAST..."
  "${QUAST_BIN}" "$CONTIGS" -o "$QOUT" -t "$THREADS"

  REPORT="${QOUT}/report.tsv"
  if [[ -f "$REPORT" ]]; then
    echo
    echo "=== Métricas extraídas de ${REPORT} ==="
    awk -F'\t' '
      $1 ~ /Total length/ { print "Total length (bp):", $2 }
      $1 == "N50"          { print "N50 (bp):", $2 }
    ' "$REPORT"
  else
    echo "[WARN] QUAST no generó report.tsv en $QOUT"
  fi
else
  echo "[WARN] quast no está en PATH. Instálalo o activa el env correcto."
fi

echo "[DONE] Ensamblaje y evaluación finalizados. Resultados en: $OUTDIR"

-------------------------------------------------------------------------------------------------

chmod +x ensamblaje_secuencias_crudas.sh
./ensamblaje_secuencias_crudas.sh
conda activate bioinfo
conda install -c conda-forge bwa samtools qualimap -y
cd /home/lauradelsolgb12/curso_bioinformatica/analysis/scripts
nano 05_mapping_evol1.sh

-------------------------------------------------------------------------------------------------
#!/usr/bin/env bash
set -euo pipefail

REF="/home/lauradelsolgb12/curso_bioinformatica/analysis/assembly_01/anc_spades/contigs.fasta"
R1="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed/evol1_R1.trim.fastq.gz"
R2="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed/evol1_R2.trim.fastq.gz"
OUTDIR="/home/lauradelsolgb12/curso_bioinformatica/analysis/mapping_evol1"
mkdir -p "${OUTDIR}"

bwa index "${REF}"
bwa mem "${REF}" "${R1}" "${R2}" > "${OUTDIR}/evol1.sam"

samtools view -b "${OUTDIR}/evol1.sam" | samtools sort -o "${OUTDIR}/evol1.sorted.bam" #pata convertir

samtools index "${OUTDIR}/evol1.sorted.bam" # para indezar a bam
samtools flagstat "${OUTDIR}/evol1.sorted.bam" > "${OUTDIR}/evol1.flagstat.txt"

qualimap bamqc -bam "${OUTDIR}/evol1.sorted.bam" -outdir "${OUTDIR}/qualimap_report" -nt 4 #el qualimap
echo "[OK] Mapeo de evol1 completado. Resultados en ${OUTDIR}"

-------------------------------------------------------------------------------------------------
chmod +x 05_mapping_evol1.sh
./05_mapping_evol1.sh 
nano 05_mapping_evol2.sh
-------------------------------------------------------------------------------------------------
#!/usr/bin/env bash
set -euo pipefail

REF="/home/lauradelsolgb12/curso_bioinformatica/analysis/assembly_01/anc_spades/contigs.fasta"
R1="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed/evol2_R1.trim.fastq.gz"
R2="/home/lauradelsolgb12/curso_bioinformatica/analysis/data/trimmed/evol2_R2.trim.fastq.gz"
OUTDIR="/home/lauradelsolgb12/curso_bioinformatica/analysis/mapping_evol2"
mkdir -p "${OUTDIR}"

bwa index "${REF}"
bwa mem "${REF}" "${R1}" "${R2}" > "${OUTDIR}/evol2.sam"

samtools view -b "${OUTDIR}/evol2.sam" | samtools sort -o "${OUTDIR}/evol2.sorted.bam" #pata convertir

samtools index "${OUTDIR}/evol2.sorted.bam" # para indexar a bam
samtools flagstat "${OUTDIR}/evol2.sorted.bam" > "${OUTDIR}/evol2.flagstat.txt"

qualimap bamqc -bam "${OUTDIR}/evol2.sorted.bam" -outdir "${OUTDIR}/qualimap_report" -nt 4 #el qualimap
echo "[OK] Mapeo de evol2 completado. Resultados en ${OUTDIR}"
-------------------------------------------------------------------------------------------------
chmod +x 05_mapping_evol2.sh
./05_mapping_evol2.sh 
conda activate bioinfo 
igv 