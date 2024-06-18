#!/bin/bash

#This package dpends on pplacer, prodigal, hmmer

genomes="./"
cpus=35
extension="fasta"
out_dir="./gtdbtk_out"
db_dir=$(echo $GTDBTK_DATA_PATH)

display_help(){
    echo "Usage: $0 -g <genomes_directory> -c <cpus> -e <extension> -d <database directory> -o <output_directory>" #$0 is the name of the script
    echo "Options:"
    echo "  -g   Genomes directory (default: ./)"
    echo "  -c   Number of CPUs (default: 35)"
    echo "  -e   File extension (default: fasta)"
    echo "  -d   Path to the GTDBTK DATA"
    echo "  -o   Output directory (default: out_dir)"
    exit 1
}

while getopts ":g:c:e:d:o:" opt
do 
    case $opt in
        g) 
            genomes="$OPTARG"
            ;;
        c)
            cpus="$OPTARG"
            ;;
        e)  
            extension="$OPTARG"
            ;;
        d)
            db_dir="$OPTARG"
            ;;
        o) 
            out_dir="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac 
done

if [ "$#" -eq 0 ]
then 
    display_help
fi 

if [ ! -d "$out_dir" ]
then 
    mkdir -p "$out_dir"
fi 

# Checking installations
echo "Checking GTDBtk installations ... "

source "$(conda info --base)/etc/profile.d/conda.sh"

######## create the download file from https://github.com/bioconda/bioconda-recipes/blob/master/recipes/gtdbtk/download-db.sh#####

cat << 'EOF' > "${out_dir}"/download-db.sh

set -e

# Configuration
N_FILES_IN_TAR=241860
DB_URL="https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
TARGET_TAR_NAME="gtdbtk_r220_data.tar.gz"

# Script variables (no need to configure)
TARGET_DIR=${1:-$GTDBTK_DATA_PATH}
TARGET_TAR="${TARGET_DIR}/${TARGET_TAR_NAME}"

# Check if this is overriding an existing version
mkdir -p "$TARGET_DIR"
n_folders=$(find "$TARGET_DIR" -maxdepth 1 -type d | wc -l)
if [ "$n_folders" -gt 1 ]; then
  echo "[ERROR] - The GTDB-Tk database directory must be empty, please empty it: $TARGET_DIR"
  exit 1
fi

# Start the download process
# Note: When this URL is updated, ensure that the "--total" flag of TQDM below is also updated
echo "[INFO] - Downloading the GTDB-Tk database to: ${TARGET_DIR}"
wget $DB_URL -O "$TARGET_TAR"

# Uncompress and pipe output to TQDM
echo "[INFO] - Extracting archive..."
tar xvzf "$TARGET_TAR" -C "${TARGET_DIR}" --strip 1 | tqdm --unit=file --total=$N_FILES_IN_TAR --smoothing=0.1 >/dev/null

# Remove the file after successful extraction
rm "$TARGET_TAR"
echo "[INFO] - The GTDB-Tk database has been successfully downloaded and extracted."

# Set the environment variable
if conda env config vars set GTDBTK_DATA_PATH="$TARGET_DIR"; then
  echo "[INFO] - Added GTDBTK_DATA_PATH ($TARGET_DIR) to the GTDB-Tk conda environment."
else
  echo "[INFO] - Conda not found in PATH, please be sure to set the GTDBTK_DATA_PATH envrionment variable"
  echo "export GTDBTK_DATA_PATH=$TARGET_DIR before running GTDB-Tk. "
fi

exit 0
EOF

chmod +x "${out_dir}"/download-db.sh


conda activate bactflow
export GTDBTK_DATA_PATH="${db_dir}"
# Setting up the database
path_f=$(which gtdbtk)
v_f=$(gtdbtk -v)
echo "I found gtdbtk version ${v_f} in ${path_f}"
#finding location of the database
db_loc=$(echo $GTDBTK_DATA_PATH)
db_loc_chk=$(echo $db_loc | grep -o "[a-zA-Z]" | wc -l)

if [ "$db_loc_chk" -gt 0 ]
then
    echo "GTDBtk database is properly setup at $db_loc"
else
    echo "Database is not properly setup, working on it ..."
    if [ ! -d "${out_dir}"/databases/gtdb_db ]
    then
        mkdir -p "${out_dir}"/databases/gtdb_db
        echo "Downloading GTDBtk database, takes a while ..."
        ./"${out_dir}"/download-db.sh "${out_dir}"/databases/gtdb_db
    fi
fi

# Upgrading
python -m pip install gtdbtk --upgrade

# Gene calling
echo "Executing gene calling..."
gtdbtk identify --genome_dir "${genomes}" --out_dir "${out_dir}/identify" --cpus "${cpus}" --extension "${extension}"

# Aligning genome 
echo "Executing aligning..."
gtdbtk align --identify_dir "${out_dir}/identify" --out_dir "${out_dir}/align" --cpus "${cpus}"

# Classification
echo "Executing classification..."
gtdbtk classify --genome_dir "${genomes}" --align_dir "${out_dir}/align" --out_dir "${out_dir}/classify" -x "${extension}" --cpus "${cpus}" --skip_ani_screen
