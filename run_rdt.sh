#!/bin/bash

cd /mnt/c/Users/Illia/Desktop/Thesis/test/ || { echo "Failed to cd into base directory"; exit 1; }

for rxn_folder in reaction_intermediates4/*; do
    echo "Processing: $rxn_folder"

    if [ -d "$rxn_folder" ]; then
        cd "$rxn_folder" || { echo "Failed to cd into $rxn_folder"; exit 1; }
    else
        echo "Skipping: $rxn_folder is not a directory"
        continue
    fi


    # Check if from_species_with_cmp or to_species_with_cmp exists and is not empty
    if [ ! -s from_species_with_cmp ] || [ ! -s to_species_with_cmp ]; then
        echo "Skipping: Either from_species_with_cmp or to_species_with_cmp is empty or missing in $rxn_folder"
        cd - || { echo "Failed to return to previous directory"; exit 1; }
        continue
    fi

    # first: run rdt
    smiles=$(cat rxn.smiles)
    java -jar /mnt/c/Users/Illia/Desktop/Thesis/test/rdt.jar -Q SMI -q "$smiles" -g -c -j AAM -f TEXT


    # now split RDT output into mol-files
    rm MOL_*
    csplit -f MOL_ ECBLAST_smiles_AAM.rxn '/$MOL/' {*}

    if [ $? -ne 0 ]; then
        echo "Error: Failed to split ECBLAST_smiles_AAM.rxn into MOL files"
        continue
    fi

    if [ ! -e MOL_00 ]; then
        echo "Error: No MOL_ files created. Check ECBLAST_smiles_AAM.rxn content."
        continue
    fi

    # handle the stuff before the first $MOL
    from_to=$(grep -m1 -B1 '$MOL' ECBLAST_smiles_AAM.rxn | head -n 1)
    from_num=$(echo "$from_to" | head -c3 | sed 's/ //g')
    to_num=$(echo "$from_to" | tail -c+4 | sed 's/ //g')
    rm MOL_00

    rm mapping_lines.txt
    counter=1

    # Debug: Check if mapping_lines.txt is empty or missing
    echo "Checking if mapping_lines.txt exists and is populated before processing..."
    if [ -e mapping_lines.txt ]; then
        echo "Existing mapping_lines.txt file found."
    else
        echo "No mapping_lines.txt file found yet."
    fi

# Create mapping_lines.txt file
> mapping_lines.txt  # Clear the file first if it exists

counter=1
last_element=""

# Loop through each molecule in the rdt_index file
sort -u from_species_with_cmp -o from_species_with_cmp
sort -u to_species_with_cmp -o to_species_with_cmp

for fn2 in MOL_??; do
    if [ ! -f "$fn2" ]; then
        echo "Skipping: $fn2 does not exist"
        continue
    fi

    tail -n +2 $fn2 | grep -v '^M  CHG' > ${fn2}.mdl

    awk '(NF==16) { print $4 "\t" $14; }; (NF==15) { print $4 "\t" (0+$13); }' ${fn2}.mdl > ${fn2}.rdt_index

    obabel -i mdl ${fn2}.mdl -o inchi -xa -xT/nochg -O ${fn2}.inchi
    obabel -i mdl ${fn2}.mdl -oinchikey -O ${fn2}.inchikey

    if [ -s ${fn2}.inchikey ]; then

        declare -A element_counts=()
        #species_id_without_cmp=$(grep $(cat ${fn2}.inchikey) species_id_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
        #species_id_without_cmp=$(grep -w "$(cat ${fn2}.inchikey)" species_id_inchikey.txt | cut -f1 | sed 's/_DASH_/-/g')
        species_id_without_cmp=$(grep -F -w "$(tr -d '\r' < ${fn2}.inchikey)" species_id_inchikey.txt | head -1 | cut -f1 | sed 's/_DASH_/-/g')

        echo "DEBUG: InChIKey from ${fn2}.inchikey: $(cat ${fn2}.inchikey)"
        echo "DEBUG: species_id_without_cmp: $species_id_without_cmp"


        inchikey_clean=$(tr -d '\r' < ${fn2}.inchikey)
        species_id="$species_id_without_cmp"

        if [ -z "$species_id" ]; then
            echo "❌ Error: No species match for InChIKey $inchikey_clean"
            continue
        fi
        : '
        if [ $counter -le $from_num ]; then

            # Get the correct line from from_species_with_cmp
            # species_id=$(sed -n "${counter}p" from_species_with_cmp | cut -d ':' -f1 | tr -d '\r')
            mapping_side="from"
            mapping_end="="
        else
            # offset=$((counter - from_num))
            # species_id=$(sed -n "${offset}p" to_species_with_cmp | cut -d ':' -f1 | tr -d '\r')
            mapping_side="to"
            mapping_end=","
        fi
        '

        if grep -Fxq "$species_id" from_species_with_cmp; then
            mapping_side="from"
            mapping_end="="
        elif grep -Fxq "$species_id" to_species_with_cmp; then
            mapping_side="to"
            mapping_end=","
        else
            echo "❌ species_id $species_id not found in from/to species lists"
            continue
        fi

        #echo $species_id > ${fn2}.species_id

        if [[ "$species_id" != "$last_species_id" ]]; then
            last_species_id="$species_id"
            echo "$species_id" > "${fn2}.species_id"
        fi

        atom_counter=0
        last_element=""

        # Handle InChI index to fetch mapping
        if grep -q 'AuxInfo.*/N:' ${fn2}.inchi; then
            inchi_index=$(grep 'AuxInfo' ${fn2}.inchi | sed 's/^.*\/N://; s/\/.*$//; s/,/ /g')
        else
            inchi_index=1
        fi

        #declare -A element_counts

        # Loop through each InChI index entry
        for rdt_line in $inchi_index; do
            element_and_index=$(tail -n +$rdt_line ${fn2}.rdt_index | head -n 1)
            element=$(echo "$element_and_index" | cut -f1)
            mapping_index=$(echo "$element_and_index" | cut -f2)

            # Clean line endings
            species_id=$(echo "$species_id" | tr -d '\r')
            element=$(echo "$element" | tr -d '\r')
            mapping_index=$(echo "$mapping_index" | tr -d '\r')
            mapping_end=$(echo "$mapping_end" | tr -d '\r')
            mapping_side=$(echo "$mapping_side" | tr -d '\r')

            # Unique key per species + element to count properly
            key="${species_id}:${element}"
            count=${element_counts[$key]:-0}
            count=$((count + 1))
            element_counts[$key]=$count

            # Now write with correct element count (e.g., C#1, C#2)
            echo "${mapping_index}	${mapping_side}	${species_id}:${element}#${count}${mapping_end}" >> mapping_lines.txt
        done
    fi

    counter=$((counter + 1))
    done

    # Now ensure that the file is created properly and formatted according to your needs.
    if [ ! -s mapping_lines.txt ]; then
        echo "Error: mapping_lines.txt is empty or missing. Ensure proper input files and processing."
        continue
    fi

        # Create mapping.txt from mapping_lines.txt
        echo "Creating mapping.txt from mapping_lines.txt"
        sort -n mapping_lines.txt | grep -v ':H#' | cut -f3 | tr '\n' ' ' | sed 's/ //g; s/,$//' > mapping.txt
        #sort -n mapping_lines.txt | grep -v ':H#' | cut -f3 | tr '\n' ' ' | sed 's/ //g; s/,$//' > mapping.txt

        if [ ! -s mapping.txt ]; then
            echo "Error: mapping.txt was not created or is empty. Check mapping_lines.txt."
            continue
        fi

        # Check if mapping.txt exists
        if [ -e mapping.txt ]; then
            echo "Successfully created mapping.txt"
        else
            echo "Failed to create mapping.txt."
        fi

        cd - || { echo "Failed to return to previous directory"; exit 1; }
done
