#!/bin/bash
echo "=== Iniciando instalación ==="

# 1. Entorno de descargas
mamba create -n env_descargas sra-tools -y

# 2. Entorno unificado para evitar conflictos de paths y versiones
mamba create -n env_perfect python=3.9 metaphlan=4.1.1 kneaddata bowtie2 diamond glpk=5.0 tbb=2020.2 -c bioconda -c conda-forge -y

# 3. HUMAnN 3.9 con instalación vía pip
source $HOME/miniforge3/bin/activate env_perfect
pip install --no-cache-dir humann==3.9

# 4. Crear carpetas
mkdir -p ~/db/kneaddata
mkdir -p ~/db/humann

echo "=== Descargando bases de datos ==="

# 5. Descargas
kneaddata_database --download human_genome bowtie2 ~/db/kneaddata
metaphlan --install
humann_databases --download chocophlan full ~/db/humann
humann_databases --download uniref uniref90_diamond ~/db/humann
humann_databases --download utility_mapping full ~/db/humann

echo "=== ¡INSTALACIÓN Y DESCARGAS COMPLETADAS! ==="