#!/bin/bash
# ===============================================
# 🧬 DNA-seq Pipeline - Flujo de llamado de variantes
# Autor: Álvaro Gómez (TACO)
# -----------------------------------------------
# Este pipeline automatiza un flujo completo de DNA-seq,
# incluyendo descarga de SRR, indexado del genoma,
# alineamiento, llamado de variantes y reconstrucción
# de secuencias génicas mutantes.
# ===============================================

# ESTRUCTURA DEL REPOSITORIO:
# DNA-seq-pipeline/
# ├── DNA-seqVariantes.R            # Script maestro en R
# ├── bin/                          # Scripts intermedios
# │   ├── check_dependencies.sh     # Verifica herramientas necesarias
# │   ├── convert_srrs.sh           # Descarga y conversión SRA → FASTQ
# │   ├── index_ref.sh              # Indexado del genoma con BWA
# │   ├── align.sh                  # Alineamiento con BWA + samtools
# │   ├── call_variants.sh          # Llamado de variantes con bcftools
# │   └── generate_mutants.py       # Reconstrucción de secuencias mutantes
# ├── data/                         # Carpeta sugerida para datos crudos
# ├── alignments/                   # Archivos BAM alineados
# ├── variants/                     # VCFs y FASTA mutantes
# └── README.sh                     # Este archivo

# -----------------------------------------------
# 📦 REQUISITOS
# -----------------------------------------------
# Herramientas requeridas:
#   - fastqc
#   - prefetch + fasterq-dump (SRA Toolkit)
#   - bwa
#   - samtools
#   - bcftools
#   - bedtools
#   - Python3 + biopython + cyvcf2
#   - R + paquete "glue"
#
# Instalar "glue" desde R:
# Rscript -e 'install.packages("glue")'
#
# Instalar dependencias de Python:
# pip install biopython cyvcf2

# -----------------------------------------------
# ▶️ USO
# -----------------------------------------------
# 1. Edita las rutas en DNA-seqVariantes.R:
#    - Lista de SRRs
#    - Genoma de referencia (FASTA)
#    - Directorios de salida
#
# 2. Ejecuta el pipeline completo:
#    Rscript DNA-seqVariantes.R
#
# 3. El script generará y ejecutará los pasos:
#    🔹 Verificación de herramientas
#    🔹 Descarga + conversión a FASTQ
#    🔹 Indexado del genoma
#    🔹 Alineamiento (BWA + samtools)
#    🔹 Llamado de variantes (bcftools)
#    🔹 Reconstrucción de secuencias mutantes (Python)

# -----------------------------------------------
# 📁 SALIDAS
# -----------------------------------------------
# - data/fastq/: Archivos FASTQ limpios
# - ref/index_bwa/: Índices del genoma
# - alignments/: Archivos BAM ordenados + indexados
# - variants/: Archivos VCF comprimidos + indexados
# - variants/GEN_saxiphilin_mutantes/: FASTA mutantes por muestra

# -----------------------------------------------
# 🔬 DETALLES TÉCNICOS
# -----------------------------------------------
# - Optimizado para lecturas paired-end (PE)
# - Por defecto usa 4 hilos (configurable en el script)
# - SNP calling limitado a variantes bialélicas simples
# - Reconstrucción génica definida en el script Python

# -----------------------------------------------
# 👨‍🔬 AUTOR
# -----------------------------------------------
# Álvaro Gómez (TACO)
# Bioingeniería · Universidad de Concepción
# Tesis en Bioinformática y Genómica
# GitHub: https://github.com/CDMTaco
