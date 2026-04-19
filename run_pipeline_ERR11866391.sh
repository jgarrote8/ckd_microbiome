#!/bin/bash
set -e
SAMPLE="ERR11866391"
THREADS=8

echo "=== Iniciando Pipeline TFG para $SAMPLE ==="

# Crear carpetas
mkdir -p ~/tfg_ckd/raw_data
mkdir -p ~/tfg_ckd/cleaned_data/${SAMPLE}_knead
mkdir -p ~/tfg_ckd/humann_output_${SAMPLE}

# --- A. DESCARGA DE ARCHIVOS FASTQ (VÍA ENA API) ---
echo ">> [A] Consultando la API de ENA para $SAMPLE..."
cd ~/tfg_ckd/raw_data

if [ ! -f "${SAMPLE}_1.fastq.gz" ] || [ ! -f "${SAMPLE}_2.fastq.gz" ]; then
    LINKS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SAMPLE}&result=read_run&fields=fastq_ftp&format=tsv" | awk 'NR==2 {print $2}')
    LINK1=$(echo $LINKS | cut -d ';' -f 1)
    LINK2=$(echo $LINKS | cut -d ';' -f 2)
    echo "Descargando archivo 1..."
    wget -q --show-progress "https://$LINK1" -O ${SAMPLE}_1.fastq.gz
    if [ "$LINK1" != "$LINK2" ] && [ -n "$LINK2" ]; then
        echo "Muestra Paired-end detectada. Descargando archivo 2..."
        wget -q --show-progress "https://$LINK2" -O ${SAMPLE}_2.fastq.gz
    else
        echo "Muestra Single-end detectada. Solo hay un archivo."
    fi
else
    echo ">>> Archivos fastq.gz ya existen, saltando descarga."
fi

# Activar el entorno de bioinformática para el resto del proceso para compatibilidad de versiones
source $HOME/miniforge3/bin/activate env_perfect

# --- B. CONTROL DE CALIDAD Y FILTRADO ---
echo ">> [B] Ejecutando KneadData..."
if [ ! -f $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_1.fastq ]; then
    kneaddata \
      -i1 $HOME/tfg_ckd/raw_data/${SAMPLE}_1.fastq.gz \
      -i2 $HOME/tfg_ckd/raw_data/${SAMPLE}_2.fastq.gz \
      -o $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead \
      -db $HOME/db/kneaddata \
      --threads $THREADS \
      --remove-intermediate-output
else
    echo ">>> KneadData ya completado, saltando."
fi

if [ ! -f $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_humann_input.fastq ]; then
    echo ">> Concatenando archivos PAIRED..."
    cat $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_1.fastq \
        $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_2.fastq \
        > $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_humann_input.fastq
else
    echo ">>> Concatenación ya hecha, saltando."
fi

# --- C. PERFIL TAXONÓMICO ---
echo ">> [C] Ejecutando MetaPhlAn 4.1.1..."
if [ ! -f $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_metaphlan_profile.txt ]; then
    metaphlan \
      $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_1.fastq,$HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_1_kneaddata_paired_2.fastq \
      --input_type fastq \
      --nproc $THREADS \
      --bt2_ps sensitive-local \
      --bowtie2out $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_bowtie2.bz2 \
      -o $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_metaphlan_profile.txt
else
    echo ">>> MetaPhlAn ya completado, saltando."
fi

# --- D. PERFIL FUNCIONAL ---
echo ">> [D] Ejecutando HUMAnN 3.9..."

# --- Forzamos a HUMAnN para que acepte la base de datos de 2025 porque sí es compatible con estas versiones de los paquetes pero en el prescreen la rechaza por fecha ---
sed -i '1s/.*/#mpa_vJun23_CHOCOPhlAnSGB_202306/' ~/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_metaphlan_profile.txt

if [ ! -f $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_humann_input_genefamilies.tsv ]; then
    humann \
      --input $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_humann_input.fastq \
      --output $HOME/tfg_ckd/humann_output_${SAMPLE} \
      --threads $THREADS \
      --nucleotide-database $HOME/db/humann/chocophlan \
      --protein-database $HOME/db/humann/uniref \
      --taxonomic-profile $HOME/tfg_ckd/cleaned_data/${SAMPLE}_knead/${SAMPLE}_metaphlan_profile.txt \
      --memory-use maximum
else
    echo ">>> HUMAnN ya completado, saltando."
fi

# --- E. EXPORTAR A KEGG Y NORMALIZAR ---
echo ">> [E] Reagrupando a KOs y normalizando..."
if [ ! -f $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko.tsv ]; then
    humann_regroup_table \
      --input $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_humann_input_genefamilies.tsv \
      --groups uniref90_ko \
      --output $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko.tsv
else
    echo ">>> Reagrupamiento KO ya hecho, saltando."
fi

if [ ! -f $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko_relab.tsv ]; then
    humann_renorm_table \
      --input $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko.tsv \
      --output $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_genefamilies_ko_relab.tsv \
      --units relab
else
    echo ">>> Normalización KO ya hecha, saltando."
fi

if [ ! -f $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_pathabundance_relab.tsv ]; then
    humann_renorm_table \
      --input $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_humann_input_pathabundance.tsv \
      --output $HOME/tfg_ckd/humann_output_${SAMPLE}/${SAMPLE}_pathabundance_relab.tsv \
      --units relab
else
    echo ">>> Normalización pathabundance ya hecha, saltando."
fi

echo "=== ¡PIPELINE COMPLETADO CON ÉXITO! ==="
