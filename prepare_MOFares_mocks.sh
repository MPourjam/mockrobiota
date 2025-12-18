#!/bin/bash
set -e  # Exit immediately if a command fails

# ==========================================
# CONFIGURATION
# ==========================================
S3_BASE="https://s3-us-west-2.amazonaws.com/mockrobiota/latest"
GITHUB_BASE="https://raw.githubusercontent.com/caporaso-lab/mockrobiota/master/data"

PRIMER_FWD_STD="GTGCCAGCMGCCGCGGTAA"
PRIMER_REV_STD="GGACTACHVGGGTWTCTAAT"
PRIMER_FWD_MOD="CCGTGCCAGCMGCCGCGGTAA"
PRIMER_REV_MOD="GGACTACHVGGGTWTCTAAT"

# ==========================================
# HELPER: Download & Verify
# ==========================================
download_data() {
    local MOCK_ID=$1
    local HAS_INDEX=$2 

    echo "  > Checking Data for $MOCK_ID..."
    mkdir -p "$MOCK_ID/raw_data"

    # --- Internal Function to Download & Validate ---
    fetch_and_validate() {
        local FILE_PATH=$1
        local URL=$2

        # 1. Delete if file is empty (0 bytes)
        if [ -f "$FILE_PATH" ] && [ ! -s "$FILE_PATH" ]; then
            echo "    ! Found empty file. Deleting: $FILE_PATH"
            rm "$FILE_PATH"
        fi

        # 2. Download if missing
        if [ ! -f "$FILE_PATH" ]; then
             echo "    - Downloading $(basename $FILE_PATH)..."
             wget -q -O "$FILE_PATH" "$URL" || echo "      ! Download failed for $URL"
        fi

        # 3. Verify & Fix Compression (Handle uncompressed downloads)
        if [ -s "$FILE_PATH" ]; then
            if file "$FILE_PATH" | grep -q "gzip compressed data"; then
                : # It is good
            elif file "$FILE_PATH" | grep -q "empty"; then
                echo "    ! Error: File is still empty after download: $FILE_PATH"
                rm "$FILE_PATH" # Remove so it tries again next time
                exit 1
            else
                echo "    ! Warning: File was not gzipped. Compressing..."
                mv "$FILE_PATH" "${FILE_PATH}.tmp"
                gzip -c "${FILE_PATH}.tmp" > "$FILE_PATH"
                rm "${FILE_PATH}.tmp"
            fi
        fi
    }

    # --- Execute Downloads ---
    # Metadata
    if [ ! -s "$MOCK_ID/sample-metadata.tsv" ]; then
        wget -q -O "$MOCK_ID/sample-metadata.tsv" "$GITHUB_BASE/$MOCK_ID/sample-metadata.tsv"
    fi

    # Sequence Files
    fetch_and_validate "$MOCK_ID/raw_data/mock-forward-read.fastq.gz" "$S3_BASE/$MOCK_ID/mock-forward-read.fastq.gz"
    fetch_and_validate "$MOCK_ID/raw_data/mock-reverse-read.fastq.gz" "$S3_BASE/$MOCK_ID/mock-reverse-read.fastq.gz"

    if [ "$HAS_INDEX" == "true" ]; then
        fetch_and_validate "$MOCK_ID/raw_data/mock-index-read.fastq.gz" "$S3_BASE/$MOCK_ID/mock-index-read.fastq.gz"
    fi
}

# ==========================================
# PROCESS: Manifest Method
# ==========================================
process_manifest_mock() {
    local DIR_NAME=$1
    local SAMPLE_ID=$2
    local FWD_PRIMER=$3
    local REV_PRIMER=$4

    echo "----------------------------------------------------------------"
    echo "Processing $DIR_NAME (Manifest)..."

    download_data "$DIR_NAME" "false"
    cd "$DIR_NAME"

    # Create Manifest
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv
    printf "${SAMPLE_ID}\t$PWD/raw_data/mock-forward-read.fastq.gz\t$PWD/raw_data/mock-reverse-read.fastq.gz\n" >> manifest.tsv

    # Import
    if [ ! -f "demux.qza" ]; then
        qiime tools import \
          --type 'SampleData[PairedEndSequencesWithQuality]' \
          --input-path manifest.tsv \
          --output-path demux.qza \
          --input-format PairedEndFastqManifestPhred33V2
    fi

    # Export Untrimmed
    if [ ! -d "untrimmed_fastq" ]; then
        qiime tools export --input-path demux.qza --output-path untrimmed_fastq
    fi

    # Trim
    if [ ! -f "demux-trimmed.qza" ]; then
        qiime cutadapt trim-paired \
          --i-demultiplexed-sequences demux.qza \
          --p-front-f "$FWD_PRIMER" \
          --p-front-r "$REV_PRIMER" \
          --o-trimmed-sequences demux-trimmed.qza
    fi

    # Export Trimmed
    if [ ! -d "trimmed_fastq" ]; then
        qiime tools export --input-path demux-trimmed.qza --output-path trimmed_fastq
    fi

    cd ..
}

# ==========================================
# PROCESS: EMP Method
# ==========================================
process_emp_mock() {
    local DIR_NAME=$1
    echo "----------------------------------------------------------------"
    echo "Processing $DIR_NAME (EMP)..."

    download_data "$DIR_NAME" "true"
    cd "$DIR_NAME"

    mkdir -p emp-paired-end-sequences
    # Force overwrite links to ensure they point to valid files
    ln -sf ../raw_data/mock-forward-read.fastq.gz emp-paired-end-sequences/forward.fastq.gz
    ln -sf ../raw_data/mock-reverse-read.fastq.gz emp-paired-end-sequences/reverse.fastq.gz
    ln -sf ../raw_data/mock-index-read.fastq.gz   emp-paired-end-sequences/barcodes.fastq.gz

    if [ ! -f "emp-seqs.qza" ]; then
        qiime tools import \
          --type EMPPairedEndSequences \
          --input-path emp-paired-end-sequences \
          --output-path emp-seqs.qza
    fi

    if [ ! -f "demux.qza" ]; then
        qiime demux emp-paired \
          --i-seqs emp-seqs.qza \
          --m-barcodes-file sample-metadata.tsv \
          --m-barcodes-column BarcodeSequence \
          --p-rev-comp-mapping-barcodes \
          --o-per-sample-sequences demux.qza \
          --o-error-correction-details demux-details.qza
    fi

    if [ ! -d "untrimmed_fastq" ]; then
        qiime tools export --input-path demux.qza --output-path untrimmed_fastq
    fi

    if [ ! -f "demux-trimmed.qza" ]; then
        qiime cutadapt trim-paired \
          --i-demultiplexed-sequences demux.qza \
          --p-front-f "$PRIMER_FWD_STD" \
          --p-front-r "$PRIMER_REV_STD" \
          --o-trimmed-sequences demux-trimmed.qza
    fi

    if [ ! -d "trimmed_fastq" ]; then
        qiime tools export --input-path demux-trimmed.qza --output-path trimmed_fastq
    fi

    cd ..
}

# ==========================================
# EXECUTION
# ==========================================

# Group 1: Multiplexed
process_emp_mock "mock-4"
process_emp_mock "mock-5"

# Group 2: Standard Primer
process_manifest_mock "mock-12" "Extreme.1" "$PRIMER_FWD_STD" "$PRIMER_REV_STD"
process_manifest_mock "mock-13" "mock-13"   "$PRIMER_FWD_STD" "$PRIMER_REV_STD"
process_manifest_mock "mock-14" "mock-14"   "$PRIMER_FWD_STD" "$PRIMER_REV_STD"
process_manifest_mock "mock-15" "mock-15"   "$PRIMER_FWD_STD" "$PRIMER_REV_STD"
process_manifest_mock "mock-16" "mock-16"   "$PRIMER_FWD_STD" "$PRIMER_REV_STD"

# Group 3: Modified Primer
process_manifest_mock "mock-18" "mock-18"   "$PRIMER_FWD_MOD" "$PRIMER_REV_MOD"
process_manifest_mock "mock-19" "mock-19"   "$PRIMER_FWD_MOD" "$PRIMER_REV_MOD"
process_manifest_mock "mock-20" "mock-20"   "$PRIMER_FWD_MOD" "$PRIMER_REV_MOD"
process_manifest_mock "mock-21" "mock-21"   "$PRIMER_FWD_MOD" "$PRIMER_REV_MOD"
process_manifest_mock "mock-22" "mock-22"   "$PRIMER_FWD_MOD" "$PRIMER_REV_MOD"
process_manifest_mock "mock-23" "mock-23"   "$PRIMER_FWD_MOD" "$PRIMER_REV_MOD"

echo "=========================================="
echo "All Mock Communities Checked & Processed."
