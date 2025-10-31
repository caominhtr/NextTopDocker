#!/bin/bash


cat << 'EOF'
#######################################################################################
#	 _   _           _   _______           ____              _                    #
#	| \ | | _____  _| |_|__   __|___  ___ |  _ \  ___   ___ | | __ ___ ____       #
#	|  \| |/ _ \ \/ / __|  | |  / _ \| _ \| | | |/ _ \ / __\| |/ // _ \| _ \      #
#	| |\  |  __/>  <| |_   | | | (_) ||_)|| |_| | (_) | (__ |   < | __/|(_ |      #
#	|_| \_|\___/_/\_\\__|  |_|  \___/| __/|____/ \___/ \___/|_|\_\\ __/|_|\|      #
#			           	 | |                                          #
#					 |_|                                          #
#       								              #
######################################################################################              
EOF

echo "This work was carried out at the Unité de Biologie Fonctionnelle et Adaptative (BFA), INSERM U1133, CNRS UMR8251, Université Paris Cité, France."
echo "For more information, please visit: https://github.com/caominhtr/NextTopDocker"
echo "If you find it useful, please cite appropriately:                  "
echo ""
echo "----------------------------------------------------------------------------------"
echo "Select task to perform:"
echo "1) Docking"
echo "Docking requires: ID_ligand, ID_protein, ID_ref. Recommended in MOL2 format."
echo "In case of redocking, ID_ref and ID_ligand must be the same. When performing docking of new ligand, ID_ref must be the co-crystal ligand of your protein."
echo ""
echo "2) Prediction only"
echo "This task can only be implemented when redocking results were created. This allows you to change models for pose selection."
echo ""
read -p "Enter choice [1-2]: " choice

merge_file="example/result/merged_score.csv"
output_file="example/result/prediction_full.csv"

case $choice in
    1)
	echo "-------------------------------------------"
	echo "Specify the model(s) you want to use for prediction."
	echo "You can choose multiple models (separated by spaces)."
	echo "Valid options: 10 20 30 40 50 60 70 80 90 100"
	echo ""
	read -p "Enter model numbers: " -a models
	echo ""

	if [ "${#models[@]}" -lt 1 ]; then
		echo "Error: No model numbers provided."
		echo "Example: 10 20 30 40"
		exit 1
	fi

	for m in "${models[@]}"; do
		if ! [[ "$m" =~ ^(10|20|30|40|50|60|70|80|90|100)$ ]]; then
		    echo "Invalid model number: $m"
		    echo "Valid options: 10 20 30 40 50 60 70 80 90 100"
		    exit 1
	fi
	done
	echo "-------------------------------------------"
	echo "Running Docking..."
	temp_file="example/temp.csv"
	> "$temp_file" 

	for ID in example/*; do
	    [ -d "$ID" ] || continue
	    ID_base=$(basename "$ID")
	    lig_path=$(find "$ID" -maxdepth 1 -name "*ligand*" | head -n1)
	    pro_path=$(find "$ID" -maxdepth 1 -name "*protein*" | head -n1)
	    ref_path=$(find "$ID" -maxdepth 1 -name "*ref*" | head -n1)
	    echo "$ID_base,$lig_path,$pro_path,$ref_path" >> "$temp_file"
	done


	#Smina posse sampling 
	echo "Running SMINA docking..."
	stt=1
	len=$(wc -l ./example/temp.csv | awk '{print $1}')

	while IFS=, read -r ID ligand protein ref; do		
		mkdir -p example/${ID}/smina_result
		for size in 10 20 30 40; do
			./smina -r "$protein" -l "$ligand" \
			      --autobox_ligand "$ref" \
			      --size_x "$size" --size_y "$size" --size_z "$size" \
			      --exhaustiveness 8 --num_modes 20 --seed 30602001 \
			      -o "example/${ID}/${ID}_ligand_${size}.mol2" > "example/${ID}/${ID}_ligand_${size}.txt"    2>/dev/null

			mkdir -p "example/${ID}/split${size}"

			awk -v prefix="example/${ID}/split${size}/${ID}_ligand_${size}_" '/^@<TRIPOS>MOLECULE/ {i++; fname = prefix sprintf("%02d", i) ".mol2"} {print > fname}' "example/${ID}/${ID}_ligand_${size}.mol2"

			
			mv example/${ID}/${ID}_ligand_* example/${ID}/smina_result
			mv example/${ID}/split* example/${ID}/smina_result
			
		done
		
		echo "$stt/$len"
		((stt++))
	done < example/temp.csv

	#Gnina 1.3 rescore
	echo "Running GNINA 1.3 rescoring..."
	stt=1

	while IFS=, read -r ID ligand protein ref; do
		mkdir -p example/${ID}/gnina_result
		
		find example/${ID}/smina_result -type f -path "*/split*/*.mol2" | sort >> example/${ID}/${ID}_temp.csv
		
		while read pose; do
			pose_base=$(basename $pose)
			pose_name=${pose_base%.mol2}
			./gnina -l $pose -r $protein --score_only > example/${ID}/gnina_result/${pose_name}_gnina.txt 2>/dev/null
		done < example/${ID}/${ID}_temp.csv
		
		echo "$stt/$len"
		((stt++))
	done < example/temp.csv 
		
	#Extracting Smina scores
	echo "Extracting Smina scores..."
	while IFS=, read -r ID ligand protein ref; do
	    for size in 10 20 30 40; do
		file="example/${ID}/smina_result/${ID}_ligand_${size}.txt"
		[ -f "$file" ] || continue  

		start_line=$(grep -n '^--' "$file" | head -n1 | cut -d: -f1)
		end_line=$(grep -n '^Refine' "$file" | head -n1 | cut -d: -f1)

		[ -z "$start_line" ] || [ -z "$end_line" ] || [ "$start_line" -eq "$end_line" ] && continue

		awk -v prefix="${ID}_ligand_${size}" -v start="$start_line" -v end="$end_line" '
		    NR > start && NR < end {
		        if(NF) {
		            score=$2
		            printf "%s_%02d,%s\n", prefix, ++i, score
		            if(i>=20) exit
		        }
		    }
		' "$file"
	    done
	done < example/temp.csv   >> example/1_smina_score.csv


	#Extracting Gnina scores
	echo "Extracting GNINA 1.3 scores..."
	while IFS=, read -r ID ligand protein ref; do
		while read pose; do
			pose_base=$(basename $pose)
			pose_name=${pose_base%.mol2}
			score=$(grep '^CNNscore' ./example/${ID}/gnina_result/${pose_name}_gnina.txt | cut -d ':' -f2 | xargs)
			echo "$pose_name,$score" >> example/2_gnina_score.csv
		done < example/${ID}/${ID}_temp.csv
		
	done < example/temp.csv

	#Combine
	awk -F, 'NR==FNR {gnina[$1]=$2; next} $1 in gnina {pdb=substr($1,1,4); print pdb","$1","$2","gnina[$1]}' \
	example/2_gnina_score.csv example/1_smina_score.csv > example/3_merged_score.csv

	rm -rf ./example/*/*_temp.csv

	mkdir -p example/result
	mv -f example/*_score.csv example/result
	
	merge_file="example/result/3_merged_score.csv"
	output_file="example/result/4_prediction_full.csv"

	python tools/logreg.py "${models[@]}" "$merge_file" "$output_file"
	status=$?

	if [ $status -ne 0 ]; then
	    echo "Error: Please provide a valid model"
	    exit 1
	fi
	echo "Prediction complete. Results saved to: $output_file"
        ;;
    2)
	echo "-------------------------------------------"
	echo "Specify the model(s) you want to use for prediction."
	echo "You can choose multiple models (separated by spaces)."
	echo "Valid options: 10 20 30 40 50 60 70 80 90 100"
	echo ""

	read -p "Enter model numbers: " -a models
	echo ""

	if [ "${#models[@]}" -lt 1 ]; then
		echo "Error: No model numbers provided."
		echo "Example: 10 20 30 40"
		exit 1
	fi

	for m in "${models[@]}"; do
		if ! [[ "$m" =~ ^(10|20|30|40|50|60|70|80|90|100)$ ]]; then
		    echo "Invalid model number: $m"
		    echo "Valid options: 10 20 30 40 50 60 70 80 90 100"
		    exit 1
		fi
	done

	merge_file="example/result/3_merged_score.csv"
	output_file="example/result/4_prediction_full.csv"

	echo "Running LogReg model(s): ${models[*]}"
	python tools/logreg.py "${models[@]}" "$merge_file" "$output_file"
	status=$?

	if [ $status -ne 0 ]; then
		echo "Error: Prediction failed. Please check model files or inputs."
		exit 1
	fi

	echo ""
	echo "Prediction complete. Results saved to: $output_file"
	;;

    *)
        echo "Invalid option. Exiting."
        exit 1
        ;;
esac

echo "-----------------------------------------"
echo "Done!"
