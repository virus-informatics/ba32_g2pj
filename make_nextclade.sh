!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 18 ##18
#$ -l s_vmem=10G
#$ -l lmem
#$ -q '!mjobs_rerun.q'
#$ #-o /dev/null

nextclade="/Users/okumurakaho/nextclade" #Gomadango2

####Set and forget.
country_info_file="country_info.txt"
py_script="gisaid_recent_filter_ver3.py"
py_script_mut="summarize_mut_info.ver2.py"


####Please edit every time. Probably just the dates.
download_date="2026-07-05" #data download date (actually upload date). Use hyphens.
out_prefix=$(echo ${download_date} | sed -e "s/-/_/g") #path for the output prefix. Date should match Gisaid file name.
working_dir="/output/${out_prefix}"
sequences_fasta="/output/${out_prefix}/sequences_fasta_${out_prefix}.tar.xz"
gisaid_file="/output/${out_prefix}/metadata_tsv_${out_prefix}" #path to the GISAID metadata file
gisaid_metadata="metadata.tsv"
metadata_mut_long="metadata.mut_long.tsv" #path for the output tsv file

####Command
cd ${working_dir}
# << COMMENTOUT
python3 ${py_script} \
       ${sequences_fasta} \
       > ${working_dir}/filtered.fasta

${nextclade} run \
       -d sars-cov-2 -j 6 \
       --cds-selection S \
       --output-tsv=${working_dir}/nextclade.tsv \
       --output-translations=${working_dir}/gene_{cds}.translation.fasta \
       ${working_dir}/filtered.fasta

python3 ${py_script_mut} \
        ${gisaid_file}/${gisaid_metadata} \
        > ${gisaid_file}/${metadata_mut_long}




