#!/usr/bin/env Rscript

# ===============================
# 📌 PARÁMETROS GENERALES
# ===============================
library(glue)

# Ruta al archivo con SRRs (uno por línea)
srr_list <- "REEMPLAZAR_CON/RUTA/SRR_ids.txt"

# Directorios
sra_dir       <- "data/sra"
fastq_dir     <- "data/fastq"
ref_genome    <- "REEMPLAZAR_CON/RUTA/genome.fa"
genome_index  <- "ref/index_bwa"
bam_dir       <- "alignments"
vcf_dir       <- "variants"
n_threads     <- 4

# Crear carpetas necesarias
dir.create("bin", showWarnings = FALSE)
dir.create(bam_dir, showWarnings = FALSE)
dir.create(vcf_dir, showWarnings = FALSE)

# ===============================
# 🔧 PASO 1: CHEQUEO DE HERRAMIENTAS
# ===============================
check_script <- '
#!/bin/bash
tools=("fastqc" "fasterq-dump" "prefetch" "bwa" "samtools" "bcftools" "tabix" "python3")

MISSING=0
for tool in "${tools[@]}"; do
  if ! command -v "$tool" &> /dev/null; then
    echo "❌ $tool no está instalado."
    MISSING=1
  else
    echo "✅ $tool está instalado."
  fi
done

if [ "$MISSING" -ne 0 ]; then
  echo "🚫 Faltan herramientas. Abortando."
  exit 1
fi
'
writeLines(check_script, "bin/check_dependencies.sh")
Sys.chmod("bin/check_dependencies.sh", "755")
stopifnot(system("bash bin/check_dependencies.sh") == 0)

# ===============================
# ⬇️ PASO 2: DESCARGA Y CONVERSIÓN
# ===============================
download_script <- glue('
#!/bin/bash
set -euo pipefail
mkdir -p "{sra_dir}" "{fastq_dir}"

while read -r SRR; do
  if [ -z "$SRR" ]; then continue; fi
  if [ ! -f "{sra_dir}/${{SRR}}.sra" ]; then
    echo "📥 Descargando $SRR..."
    prefetch "$SRR" --output-directory "{sra_dir}"
  fi
  echo "🔄 Convirtiendo $SRR a FASTQ..."
  fasterq-dump "{sra_dir}/${{SRR}}.sra" -O "{fastq_dir}" --split-files --threads {n_threads}
  echo "📦 Comprimiendo $SRR..."
  gzip -f "{fastq_dir}/${{SRR}}_1.fastq"
  gzip -f "{fastq_dir}/${{SRR}}_2.fastq"
done < "{srr_list}"
')
writeLines(download_script, "bin/convert_srrs.sh")
Sys.chmod("bin/convert_srrs.sh", "755")
stopifnot(system("bash bin/convert_srrs.sh") == 0)

# ===============================
# 🧬 PASO 3: INDEXADO DEL GENOMA
# ===============================
index_script <- glue('
#!/bin/bash
set -euo pipefail
mkdir -p "{genome_index}"
bwa index -p "{genome_index}/genome" "{ref_genome}"
')
writeLines(index_script, "bin/index_ref.sh")
Sys.chmod("bin/index_ref.sh", "755")
stopifnot(system("bash bin/index_ref.sh") == 0)

# ===============================
# 🔗 PASO 4: ALINEAMIENTO
# ===============================
align_script <- glue('
#!/bin/bash
set -euo pipefail
mkdir -p "{bam_dir}"
shopt -s nullglob
for fq1 in "{fastq_dir}"/*_1.fastq.gz; do
  fq2="${{fq1/_1.fastq.gz/_2.fastq.gz}}"
  sample=$(basename "$fq1" _1.fastq.gz)
  echo "🧬 Alineando $sample"
  bwa mem -t {n_threads} "{genome_index}/genome" "$fq1" "$fq2" | \
  samtools sort -@ {n_threads} -o "{bam_dir}/${{sample}}.sorted.bam"
  samtools index "{bam_dir}/${{sample}}.sorted.bam"
done
')
writeLines(align_script, "bin/align.sh")
Sys.chmod("bin/align.sh", "755")
stopifnot(system("bash bin/align.sh") == 0)

# ===============================
# 🧬 PASO 5: LLAMADO DE VARIANTES (VCF comprimido + index)
# ===============================
vcf_script <- glue('
#!/bin/bash
set -euo pipefail
mkdir -p "{vcf_dir}"
shopt -s nullglob
for bam in "{bam_dir}"/*.sorted.bam; do
  sample=$(basename "$bam" .sorted.bam)
  echo "🧪 Llamando variantes en $sample"
  bcftools mpileup -f "{ref_genome}" "$bam" | \
  bcftools call -mv -Oz -o "{vcf_dir}/${{sample}}.vcf.gz"
  tabix -p vcf "{vcf_dir}/${{sample}}.vcf.gz"
done
')
writeLines(vcf_script, "bin/call_variants.sh")
Sys.chmod("bin/call_variants.sh", "755")
stopifnot(system("bash bin/call_variants.sh") == 0)

# ===============================
# 🧬 PASO 6: Reconstrucción de secuencia mutante por muestra (Python)
# ===============================
mut_script <- glue('
#!/usr/bin/env python3
import sys
from pathlib import Path
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from cyvcf2 import VCF
except Exception as e:
    sys.stderr.write(f"ERROR: faltan dependencias de Python (biopython, cyvcf2). Detalle: {e}\\n")
    sys.exit(1)

GENOME_PATH = "<<ref_genome>>"
VCF_DIR = Path("<<vcf_dir>>")
OUTPUT_DIR = VCF_DIR / "GEN_saxiphilin_mutantes"
OUTPUT_DIR.mkdir(exist_ok=True)

# ---- Parámetros del gen (editar si cambian) ----
GENE_ID   = "MCH030139.1"
CHROM     = "PGA_scaffold8"
GEN_START = 15590383
GEN_END   = 15592441
STRAND    = "-"   # "+" o "-"

# Cargar genoma y subsecuencia (0-based slicing)
genome = {rec.id: rec.seq for rec in SeqIO.parse(GENOME_PATH, "fasta")}
if CHROM not in genome:
    sys.stderr.write(f"ERROR: cromosoma {CHROM} no encontrado en {GENOME_PATH}\\n")
    sys.exit(1)
ref_seq = list(str(genome[CHROM][GEN_START-1:GEN_END]))
pos_map = {f"{CHROM}:{pos}": i for i, pos in enumerate(range(GEN_START, GEN_END+1))}

# Tomar VCFs (acepta *_annot.vcf.gz o *.vcf.gz)
vcfs = sorted(set(list(VCF_DIR.glob("*_annot.vcf.gz")) + list(VCF_DIR.glob("*.vcf.gz"))))
if not vcfs:
    sys.stderr.write(f"ADVERTENCIA: no se encontraron VCFs en {VCF_DIR}\\n")

for vcf_path in vcfs:
    sample = vcf_path.name.replace("_annot.vcf.gz","").replace(".vcf.gz","")
    vcf = VCF(str(vcf_path))
    mut_seq = ref_seq.copy()
    applied = 0
    for var in vcf:
        if not var.ALT:
            continue
        key = f"{var.CHROM}:{var.POS}"
        # Solo SNPs simples dentro de la región
        if key in pos_map and len(var.REF)==1 and len(var.ALT[0])==1:
            idx = pos_map[key]
            if mut_seq[idx]==var.REF:
                mut_seq[idx]=var.ALT[0]
                applied+=1

    final_seq = Seq("".join(mut_seq))
    if STRAND=="-":
        final_seq = final_seq.reverse_complement()

    out = OUTPUT_DIR / f"{sample}_GEN_saxi_mut.fasta"
    rec = SeqRecord(final_seq, id=f"{sample}_mutGEN",
                    description=f"{GENE_ID} | {applied} SNPs en región completa")
    SeqIO.write(rec, out, "fasta")
    print(f"✅ {sample}: {applied} SNPs aplicados → {out}")
', .open='<<', .close='>>')
writeLines(mut_script, "bin/generate_mutants.py")
Sys.chmod("bin/generate_mutants.py", "755")
stopifnot(system("python3 bin/generate_mutants.py") == 0)

# ===============================
# ✅ FINALIZADO
# ===============================
cat("🎉 Pipeline DNA-seq completado: VCFs (.vcf.gz) y secuencias mutantes por muestra generadas.\n")
