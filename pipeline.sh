#!/bin/bash
SAMPLE="ERR11866391"
THREADS=8

echo "=== Iniciando Pipeline TFG para $SAMPLE ==="

# Crear carpetas
mkdir -p ~/tfg_ckd/raw_data
mkdir -p ~/tfg_ckd/cleaned_data/${SAMPLE}_knead
mkdir -p ~/tfg_ckd/humann_output_${SAMPLE}

# --- A. DESCARGA ---
echo ">> [A] Descargando reads crudos..."
source $HOME/miniforge3/bin/activate env_descargas
cd ~/tfg_ckd/raw_data
prefetch $SAMPLE
fasterq-dump $SAMPLE --split-files --threads $THREADS
pigz -p $THREADS ${SAMPLE}_1.fastq
pigz -p $THREADS ${SAMPLE}_2.fastq
conda deactivate

# --- B. CONTROL DE CALIDAD Y FILTRADO ---
echo ">> [B] Ejecutando KneadData..."
source $HOME/miniforge3/bin/activate env_perfect
kneaddata \
  -i1 ~/tfg_ckd/raw_data/${SAMPLE}_1.fastq.gz \
  -i2 ~/tfg_ckd/raw_data/${SAMPLE}_2.fastq.gz \
  -o ~/tfg_ckd/cleaned_data/${SAMPLE}_knead \
  -db ~/db/kneaddata \
  --threads $THREADS \
  --remove-intermediate-output

echo ">> Concatenando archivos PAIRED..."
cat ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_1.fastq \
    ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_2.fastq \
    > ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_humann_input.fastq

# --- C. PERFIL TAXONÓMICO ---
echo ">> [C] Ejecutando MetaPhlAn 4.1.1..."
metaphlan \
  ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_1.fastq,~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_2.fastq \
  --input_type fastq \
  --nproc $THREADS \
  --bt2_ps sensitive-local \
  -o ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_metaphlan_profile.txt

# --- D. PERFIL FUNCIONAL ---
echo ">> [D] Ejecutando HUMAnN 3.9..."
humann \
  --input ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_humann_input.fastq \
  --output ~/tfg_ckd/humann_output_${SAMPLE} \
  --threads $THREADS \
  --nucleotide-database ~/db/humann/chocophlan \
  --protein-database ~/db/humann/uniref \
  --taxonomic-profile ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_metaphlan_profile.txt \
  --memory-use maximum

# --- E. EXPORTAR A KEGG Y NORMALIZAR ---
echo ">> [E] Reagrupando a KOs y normalizando..."
humann_regroup_table \
  --input ~/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_humann_input_genefamilies.tsv \
  --groups uniref90_ko \
  --output ~/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko.tsv

humann_renorm_table \
  --input ~/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko.tsv \
  --output ~/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko_relab.tsv \
  --units relab

humann_renorm_table \
  --input ~/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_humann_input_pathabundance.tsv \
  --output ~/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_pathabundance_relab.tsv \
  --units relab

echo "=== ¡PIPELINE COMPLETADO CON ÉXITO! ==="