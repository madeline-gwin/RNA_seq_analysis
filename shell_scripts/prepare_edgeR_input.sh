#!/bin/bash

SPECIES="Marigold"
files=(/scratch/magwin/genome/read_mapping/$SPECIES/*ReadsPerGene.out.tab)
output="/scratch/magwin/callus_RNA/"$SPECIES"_edgeR_input.tab"

if [ ${#files[@]} -eq 0 ]; then
    echo "❌ No ReadsPerGene.out.tab files found."
    exit 1
fi

tmpdir=$(mktemp -d)

echo "🔍 Checking files and extracting data..."

# Extract Gene IDs from the first file, skip first 4 lines
tail -n +5 "${files[0]}" | awk '{print $1}' > "$tmpdir/gene_ids.txt"
expected_rows=$(wc -l < "$tmpdir/gene_ids.txt")
echo "✅ ${files[0]} has $expected_rows valid gene rows."

# Start header
header="GeneID"
count_files=()

# Loop through files
for file in "${files[@]}"; do
    sample_name=$(basename "$file" ReadsPerGene.out.tab)
    header+="\t${sample_name}"

    tmp_count_file="$tmpdir/${sample_name}.counts"
    tail -n +5 "$file" | awk '{print $2}' > "$tmp_count_file"

    # Check row count
    actual_rows=$(wc -l < "$tmp_count_file")
    if [ "$actual_rows" -ne "$expected_rows" ]; then
        echo "❌ ERROR: $file has $actual_rows rows, expected $expected_rows."
        echo "🧹 Cleaning up..."
        rm -r "$tmpdir"
        exit 1
    else
        echo "✅ $file has $actual_rows rows."
    fi

    count_files+=("$tmp_count_file")
done

# Write header
echo -e "$header" > "$output"

# Paste gene IDs and counts
paste "$tmpdir/gene_ids.txt" "${count_files[@]}" >> "$output"

# Clean up
rm -r "$tmpdir"
echo "✅ edgeR_input.tab created successfully!"

