# Viral_genomes
NCBI_viral_genomes

anthil psw

7#59XitmPR#

7#59XitmPR#
#From the original sample we took 100,000 without replacement:

cat VT_1A_good_out.fasta |\awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' |\shuf |\head -n 100000 |\awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' |\awk '/^>/{print ">" ++i; next}{print}' > VT_1A_good_out_100000sub.fasta


#Loop sin replacement

for filename in *.fasta; do cat "${filename}" |\awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' |\shuf |\head -n 100000 |\awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' |\awk '/^>/{print ">" ++i; next}{print}' > "${filename}"_100000sub.fasta; done

###Loop bueno


for filename in *.fasta
do  
	cat "${filename}" |\awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t"
,$0); next;} {printf("%s",$0);} END { printf("\n");}' |\shuf |\head -n 10 |\awk 
'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' |\awk '/^>/{print ">" ++i; next}{prin
t}' > "${filename}"_100000sub.fasta
done


 
 #Then we duplicate the file 
 for f in *_100000sub.fasta ; do cp $f duplicates/copy_$f ; done

#Duplicar en mismo folder

for f in *_100000sub.fasta ; do cp $f copy_$f ; done 

#Run find_overlapp

nohup python3 find_overlaps_sept2021.py copy_FM3_2A_good_out_100000sub.fasta copy_FM3_2A_good_out_100000sub.fasta > 2_overlapp_copy_FM3_2A_good_out_100000sub.fasta -e 10 &

#Cut column 1,2

cut -f1,2 2_overlapp_copy_VM3_2A_good_out_100000sub.fasta > new_2_overlapp_copy_VM3_2A_good_out_100000sub.fasta

#Cut column loop

for filename in *.fasta; do 
cut -f1,2 "${filename}" > "col_1and2_${filename}"; done





#Replace \t with a ,
sed -i 's/\t/,/g' new_2_overlapp_copy_*


 #To run contig spectra en loop

for filename in *.fasta; do 
 
 awk -F, '
   {
     if($1 in A && $2 in A) {
       if(A[$1]!=A[$2]) {
         m=A[$2]
         for(i in A)
           if(A[i]==m)
             A[i]=A[$1]
       }
     }
     else if($1 in A)
       A[$2]=A[$1]
     else if($2 in A)
       A[$1]=A[$2]
     else
       A[$1]=A[$2]=++c
   }
   END {
     for(i=1; i<=c; i++) {
       s=x
       for(j in A)
         if(i==A[j])
           s=s (s?FS:x) j
       if(s)
         print "{" s "}"
     }
   }
 ' "${filename}" > "contig_spectra_${filename}"; done


 #continua 

for filename in *.fasta; do 
awk -F\, '{print NF-1}' "${filename}" > "${filename}_con_abundancia"; done


for filename in contig_spectra_*; do 
awk -F\, '{print NF-1}' "${filename}" > "${filename}_con_abundancia"; done

#
for filename in *.fasta; do
pr -mts' ' Contig_spectra_FM3_2C_new_python_col1_counts.txt Contig_spectra_FM3_2C_new_python_col1.txt > final_FM3_2C.txt

for filename in *_con_abundancia; do
pr -mts' ' "${filename}" "${filename}" > "$final_{filename}"; done

for filename in *_con_abundancia; do
pr -mts' ' Contig_spectra_FM3_2C_new_python_col1_counts.txt Contig_spectra_FM3_2C_new_python_col1.txt > final_FM3_2C.txt

Ejemplo usado:
pr -mts' ' 	contig_spectra_MT_3E_good_out_100000sub.fasta_con_abundancia	contig_spectra_MT_3E_good_out_100000sub.fasta	>	Final_MT_3E.fasta	&

#Esto es para contar los elementos desde el all_Samples_contig_spectra

cut -d':' -f1 All_samples_contig_spectra_script_corrected.txt | sort | uniq | xargs -I{} sh -c "echo {} | tr '\n' '\t' ; grep {} All_samples_contig_spectra_script_corrected.txt | cut -d '{' -f2- | tr -d '}' | tr ',' '\n' | sort | uniq | wc -l "


#For loop para comprimir

for filename in *out.fasta_100000sub.fasta; do
gzip -9 "${filename}"; done


for filename in *_good_out.fasta_100000v2_output.fasta; do
gzip -9 "${filename}"; done &

for filename in *_good_out.fasta; do
gzip -9 "${filename}"; done &


#Overlap_spectra: se necesita la primer columna separada por comas del find_overlapp script


Luego se corre este ES EL bueno:
for filename in *.fasta; do awk -F[,] '{c[$1]++; c[$2]++} END {for(i in c) print i, c[i] | "sort"}' "${filename}" > ../Overlap_spectra/"merge_columns${filename}".txt; done 


#Cut columns from the find_overlapp script output

for file in *.fasta;
do
    awk '{ total += $6; count++ } END { print total/count }' "$file" > Name_and_Length/"$(basename "$file")_sorted.txt"
done

for filename in *.fasta;
do
cut -f1,2,5 "${filename}" > Name_and_Length/"merge_columns1_2_5${filename}".txt; done

for filename in *.fasta;
do
awk 'NF == 7 { print $1"_"$2","$5 }
     NF != 7 { print }' "${filename}" > tst_2/"merge_columns1_2_5${filename}".txt; done


Esto es para conseguir las diferentes columnas del find overlapp script output con un loop (En este caso la longitud):

for filename in *.fasta; do awk 'NF == 7 { print $1"_"$2","$5 }
     NF != 7 { print }' "${filename}" > tst_2/"merge_columns1_2_5${filename}".txt; done

grep "." merge_columns1_2_5Entropy10_2_overlapp_copy_* > All_samples_name_length.txt

Esto es para conseguir las diferentes columnas del find overlapp script output con un loop (En este caso la entropia-se corre en la carpeta con los outputs de ese find overlap script):

for filename in *.fasta; do awk 'NF == 7 { print $1"_"$2","$7 }
     NF != 7 { print }' "${filename}" > tst_2/"merge_columns1_2_7${filename}".txt; done



#This is for the rank ocurrences per sample:

cut -d':' -f1 All_samples_overlapp_spectra.txt | sort | uniq | xargs -I{} bash -c "echo {} | tr '\n' '\t'; grep {} All_samples_overlapp_spectra.txt | cut -d' ' -f2 | tr '\n' '\t' ; echo" > test.txt

#This is to transpose the table of rank ocurrences:
transpose the rank of occurrences tables : ~/transpose --limit 10000000x10000000 -t test.txt > tranpose2_redo.txt



#Pasar files a Mark desde edwards
scp  quinto@edwards-data.sdsu.edu:/home3/quinto/frap_maker2/FRAP_Real_Data_results/samples_wo_replacement/Corrected_Script/Col1_2/Contig_spectra_corrected/Final_contig_spectra/All_samples_contig_spectra_script_corrected.txt  mlittle@edwards-data.sdsu.edu:/home1/mlittle/Andres_scp_files/Corrected_Script

##Files para mi desde Mark
scp mlittle@anthill.sdsu.edu:/home1/mlittle/Andres_scp_files/virus_2022/*.fasta quinto@anthill.sdsu.edu:/home3/quinto/FRAP2/FRAP-basic/Marks_help

scp quinto@anthill.sdsu.edu:/home3/quinto/FRAP2/FRAP-basic/* mlittle@anthill.sdsu.edu:/home1/mlittle/Andres_scp_files/ASQ_frap/

2243mark

Para extraer de la tabla:

awk -v s=2 '{print ($1-2)/s, ($2-2)/s}' tabla_probar


awk '{print ($1-2)/s, ($2-2)/s}' tabla_probar


FM3_2A FM3_2C	FM3_2E	FM3_3A
2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	36	4	2	2	58	2	2	2	2	2	2	2	2	42	2	2	2	2	2	6	2	2	2	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	12	2	2	2	2	2	2	2	2	22	2	2	2	2	2	2	2	2	2	2	2	2	2	2	88	2	2	22

FM3_2A	FM3_2C	FM3_2E	FM3_3A
4	5	22	36
42		6	58
22		22	1
8		22	88
			6


2	2	2	2	2


Col1 Col2 Col3 Col4 Col5 Col6
A  1  A  -9  0  -9
B  2   T  -9  -9  -9
C  3   D  -9  1   -9
D  4   R  -9  2   -9



This is to ignore the header
n=2
awk -v c=$n '/_/ && i++ < c {next};1' tabla_probar2 >tabla_probar2_sin_header
From that output you can do the formula
awk '{print ($1-2)/2, ($2-2)/2}'

Remove zeros


49.5 54
51.5 74
98.5 149
221.5 0
5 6
5.5 0
4 7

Sample1 	Sample2
49.5 	54
51.5 	74
98.5 	149
221.5 	0
5 	6
5.5 	0
4 	7
3 	6
10 	2

#To create the Ben's table

cut -d' ' -f2 merge_columnscol_1and2_MM_2C_find_overlap_output.txt | sort -n | uniq -c

for filename in *.txt; do cut -d' ' -f2 merge_columnscol_1and2_MM_2C_find_overlap_output.txt | sort -n | uniq -c "${filename}" > Ben_table/"Ben_table${filename}"; done


Este es el bueno

for filename in *.txt; do cut -d' ' -f2 "${filename}"| sort -n | uniq -c > Ben_table_200/"${filename}"_salida; done



##remove reads less than 200 bp
less_200_script.pl
for filename in *_only_200.fasta; do perl less200bp.pl 200 "${filename}" > "${filename}"_less_200.fasta; done

##cut reads larger than 200bp to 200bp
##loop
for filename in *.fasta; do cut -c1-200  "${filename}" | grep -A1 "^>" | grep -v "^-" >  "${filename}"_only_200.fasta; done






##Script/loop for removing seqs less than 200

for filename in *_good_out.fasta; do perl 200bp_round/less_200_script.pl 200 "${filename}" > 200bp_round/"${filename}"_less_200.fasta; done


for filename in *.fasta; do cut -c1-200  "${filename}" | grep -A1 "^>" | grep -v "^-" >  "${filename}"_only_200.fasta;  done


merge_columnscol_1and2_MM_3A_find_overlap_output_only_200_Entropy.txt.txt.txt_salida:      8 52



scp -P 7010 quinto@edwards-data.sdsu.edu:/home3/quinto/Andres/frap_maker2/FRAP_Real_Data_results/200bp_round/Only_200bp/Subsamples/2_overlap_200bp/Col1_2/Overlap_spectra_only_200/Ben_table_200/All_samples_overlapp_spectra_200_table.txt .

#Para cambiar el nombre de muchos archivos

for i in *; do
  mv "$i" "`echo $i | sed "s/regex/replace_text/"`";
done


for filename in *.fasta; do awk 'NF == 6 { print $1"_"$2","$5 }
     NF != 6 { print }' "${filename}" > Name_Length/"merge_columns1_2_5${filename}".txt; done


#Para correr make_clusters.py en un loop para todas las muestras en un nuevo folder

for filename in *.fasta; do python3 make_clusters_oct.py "${filename}" > Clusters_sample_vs_itself/"${filename}"_clusters.fasta; done

grep "." *.fasta > All_reverse_complement_vs_itself_contig_spectra.txt



#Para contig_spectra

cut -d' ' -f1 Find_overlap_vs_itself_200_only_FM3_2A_only200.fasta_clusters.fasta | sort -n | uniq -c | head

 
grep "." *.fasta > All_sample_clusters.txt

cut -d' ' -f1 All_sample_clusters | sort -n | uniq -c > All_samples_number_of_contig_spectra.txt


Para contig_spectra:
cut -d' ' -f1 All_samples_for_contig_spectra.txt | sort -n | uniq -c > All_samples_number_of_contig_spectra.txt 

################################################################################

ORDEN DESPUES DEL FIND OVERLAP:

CREA CARPETA PARA CLUSTERS
#Para correr make_clusters.py en un loop para todas las muestras en un nuevo folder

for filename in *.fasta; do python3 make_clusters_oct.py "${filename}" > Clusters_sample_vs_itself/"${filename}"_clusters.fasta; done

Para ordernar el numero de clusters de mayor a menor

for filename in *txt; do cut -d' ' -f1 "${filename}" | sort -n | uniq -c > ../Clusters_sorted/"${filename}"_clusters_sorted.txt; done

for filename in *txt; do cut -d' ' -f1 "${filename}" | sort -n | uniq -c > ../Clusters_sorted/"${filename}"_clusters_sorted.txt; done &

#Luego juntalos en un file

grep "." *.fasta > All_samples_vs_reverse_complement_contig_spectra.txt

#Luego cuenta el numeros de clusters por muestra del archivo agrupado

cut -d' ' -f1 All_samples_vs_reverse_complement_contig_spectra.txt | sort -n | uniq -c > All_samples_vs_reverse_complement_number_of_contig_spectra.txt

LISTO!!!

OVERLAP_SPECTRA

En la carpeta donde estan el output del find_overlap crea una carpeta nueva para los col1_2

for filename in *.fasta; do 
cut -f1,2 "${filename}" > caerpeta_nueva/"col_1and2_${filename}"; done

Ejemplo:
for filename in *.fasta; do cut -f1,2 "${filename}" > Col1_2/"col_1and2_${filename}"; done

Se reemplazan los tabs por commas en todos los files

sed -i 's/\t/,/g' col_1and2_Find_overlap_normal_reverse_complement_without_r_*

Se crea nueva carpeta para los overlap_spectra y se corre:

for filename in *.fasta; do awk -F[,] '{c[$1]++; c[$2]++} END {for(i in c) print i, c[i] | "sort"}' "${filename}" > Overlap_spectra/"merge_columns${filename}".txt; done 

Creas una nueva carpeta para la tabla de Ben

for filename in *.txt; do cut -d' ' -f2 "${filename}"| sort -n | uniq -c > Ben_table_sample_vs_reverse_complement/"${filename}"_salida; done

Se agrupan los files en uno solo y ya esta

grep "." *_sample_vs_reverse_complement > All_samples_vs_reverse_comement_overlap_spectra.txt


##Script for randomizing the sequences:

cat combination_sample_and_reverse_numerado_FM3_2A_only200.fasta |awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' |shuf |head -n 200000 |awk '{printf("%s\n%s\n",$1,$2)}' > combination_sample_and_reverse_numerado_FM3_2A_only200_with_random_sequences.fasta

 nano make_pairs.sh 
more make_pairs.sh
 #!/bin/bash
for file1 in *.fasta      ### Outer for loop ###
do
    for file2 in *.fasta ### Inner for loop ###
    do
          echo "$file1 $file2"
    done
done
 #Aqui acaba ese file
chmod 777 make_pairs.sh

 ./make_pairs.sh 
 ./make_pairs.sh > all_combinations.txt
 nano submit_multi.sh 
 more submit_multi.sh
 #!/bin/bash 

#$ -S /bin/bash
#$ -cwd

#$ -pe smp 4

#$ -t 1-4225

SEEDFILE=all_combinations.txt
SEED=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

#$ -o sge/out-$JOB_ID-$TASK_ID
#$ -e sge/err-$JOB_ID-$TASK_ID


python3 find_overlaps_nov_2021_NO_BUG.py -r ${SEED} > "${SEED}.out"
 #Aqui acaba ese file
 qsub -q default submit_multi.sh
 qstat
 

for filename in *._output.txt; do 
grep -E "${filename}" test1.txt > "repeats_${filename}"; done


for filename in repeats*; do 
grep -E '(AAT|TTA)\1{9,}' "${filename}" > more_repeats_"${filename}"test2.txt

for filename in *_output.txt; do 
python3 make_clusters_nov.py "${filename}" > clusters_"${filename}"; done


for filename in *_output.txt; do
awk '/99999/{y=1;next}y' "${filename}" > output_1_sampling2_part_"${filename}"; done &


#Para imprimir todo despues de un patron
sed -n '/>100001/,$p' VT_4D_good_out_200000sub.fasta > test_after.txt &

#Para imprimir todo antes de un patron
sed '/>100001/q' VT_4D_good_out_200000sub.fasta > test_before.txt &

#Para borrar la linea despues del patron
sed '/>100001/d' test_before.txt > no_last_line_test_before.txt &


for filename in *.fasta; do
sed -n '/>100001/,$p' "${filename}" > after_"${filename}" ; done &

for filename in *.fasta; do
sed '/>100001/Q' "${filename}" > before_"${filename}".txt ; done &

for filename in *.fasta; do
sed '/>100001/d' "${filename}" > good_before_"${filename}" ; done &

for filename in *.fasta; do cut -d' ' -f1 "${filename}" > "${filename}"_first_column.txt ; done &

for filename in *.fastq; do seqtk seq -a "${filename}" > salida_fasta_"${filename}"; done &

awk '/^>/{print ">" ++i; next}{print}' 
< Prochlorococcus_metagenome_bac_and_viral_OUTPUT.fa > Prochlorococcus_metagenome_bac_and_viral_OUTPUT_numerado.fa


for filename in find_overlaps_output*; do cut -f1 "${filename}" > col1/col1_"${filename}".txt ; done &
for filename in find_overlaps_output*; do cut -f2 "${filename}" > col2/col2_"${filename}".txt ; done &


for filename in col1_find_overlaps*; do cat "${filename}" | sort | uniq -c > sorted_"${filename}"; done &

for filename in col2_find_overlaps*; do cat "${filename}" | sort | uniq -c > sorted_"${filename}"; done &

cat col1_find_overlaps_outputFM3_2A_vs_FM3_5E_output.txt.txt | sort | uniq -c

find_overlaps_output*

 sed -i 's/ /,/g' sorted_col1_find_overlaps_output* &
 sed -i 's/,,,,,,,//g' *.txt &
 sed -i 's/,,,,,,//g' *.txt &
 sed -i 's/,,,,,//g' *.txt &
 sed -i 's/,,,,//g' *.txt &
 sed -i 's/,,,//g' *.txt &
 sed -i 's/,,//g' *.txt &


cut -d"," -f1 sorted_col1_find_overlaps_outputMT_5A_vs_MM_2E_output.txt.txt | sort | uniq -c

for filename in sorted_col1_*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_"${filename}"; done &

cat "${filename}" | sort | uniq -c > sorted_"${filename}"; done &

for filename in sorted_col1_*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_"${filename}"; done &

 sed -i 's/ /,/g' sorted_col2_find_overlaps_output* &
 sed -i 's/,,,,,,//g' *.txt &
 sed -i 's/,,,,,//g' *.txt &
 sed -i 's/,,,,//g' *.txt &
 sed -i 's/,,,//g' *.txt &
 sed -i 's/,,//g' *.txt &

for filename in sorted_col2_*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_"${filename}"; done &

grep "." overlap_spectra_sorted_col2_find_overlaps_output* > All_samples_col2_overlap_spectra.txt

cat find_overlaps_outputVT_4D_vs_VT_4D_output.txt | cut -f7

#Esto es para cortar la columna de las sequencias en un folder:

for filename in find_overlaps_*; do cat "${filename}" | cut -f7 > overlaps_to_compare_col7/col7_sequence_"${filename}"; done &


###Segundo round para Ben, sample vs itself

for filename in Find_overlpas_vs_itself_1subsamples_*; do cut -f1 "${filename}" > col1/col1_"${filename}".txt ; done &

sed -i 's/ /,/g' sorted_col2_find_overlaps_output* &
 sed -i 's/,,,,,,//g' *.txt &
 sed -i 's/,,,,,//g' *.txt &
 sed -i 's/,,,,//g' *.txt &
 sed -i 's/,,,//g' *.txt &
 sed -i 's/,,//g' *.txt &

 for filename in sorted_col1_*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_"${filename}"; done &

grep "." overlap_spectra_sorted_col1_Find_overlpas_vs_itself_1subsamples_* > All_samples_col1_overlap_spectra.txt

for filename in Find_overlpas_vs_itself_1subsamples_*; do cut -f2 "${filename}" > col2/col2_"${filename}".txt ; done &

sed -i 's/ /,/g' sorted_col2_find_overlaps_output* &
 sed -i 's/,,,,,,//g' *.txt &
 sed -i 's/,,,,,//g' *.txt &
 sed -i 's/,,,,//g' *.txt &
 sed -i 's/,,,//g' *.txt &
 sed -i 's/,,//g' *.txt &


for filename in sorted_col2_*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_"${filename}"; done &

##Segundo round sample vs other
for filename in find_overlaps_output*; do cut -f1 "${filename}" > col1_bueno/col1_"${filename}".txt ; done &

for filename in col1_find_overlaps*; do cat "${filename}" | sort | uniq -c > sorted_"${filename}"; done &
sed -i 's/ /,/g' sorted_col2_find_overlaps_output* &
 sed -i 's/,,,,,,//g' *.txt &
 sed -i 's/,,,,,//g' *.txt &
 sed -i 's/,,,,//g' *.txt &
 sed -i 's/,,,//g' *.txt &
 sed -i 's/,,//g' *.txt &



for filename in sorted_col1_*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_"${filename}"; done &

grep "." overlap_spectra_sorted_col1_Find_overlpas_vs_itself_1subsamples_* > All_samples_col1_overlap_spectra.txt



for filename in find_overlaps_output*; do cut -f2 "${filename}" > col2_bueno/col2_"${filename}".txt ; done &
for filename in col2_find_overlaps*; do cat "${filename}" | sort | uniq -c > sorted_"${filename}"; done &



##################Workflow final:
###cortar columna primera, si fuera la segunda seria cut -f2
for filename in find_overlaps_output*; do cut -f1 "${filename}" > col1/col1_"${filename}".txt ; done &

### numerar y contar la primer columna
for filename in col1_find_overlaps*; do cat "${filename}" | sort | uniq -c > sorted_"${filename}"; done &

###Quitar los espacios por ,
sed -i 's/ /,/g' sorted_col1_find_overlaps_output* &
sed -i 's/,,,,,,,,,,//g' *.txt &
sed -i 's/,,,,,,,,,//g' *.txt &
sed -i 's/,,,,,,,,//g' *.txt &
sed -i 's/,,,,,,,//g' *.txt &
sed -i 's/,,,,,,//g' *.txt &
sed -i 's/,,,,,//g' *.txt &
sed -i 's/,,,,//g' *.txt &
sed -i 's/,,,//g' *.txt &
sed -i 's/,,//g' *.txt &

###sort y contar los unicos
for filename in sorted_col1_find_overlaps_output*; do cut -d"," -f1 "${filename}" | sort -n | uniq -c > overlap_spectra_col1_"${filename}"; done &

### Agrupar los archivos
grep "." overlap_spectra_col1_sorted_col1_find_overlaps_output* > All_samples_col1_overlap_spectra.txt


###Ultimo path de trabajo

cd /home3/quinto/frap_maker2/FRAP_Real_Data_results/samples_wo_replacement/Refresh_ASQ

0) ####For loop para correr todo en una carpeta

 for FN1 in *.fasta; do for FN2 in *.fasta; do python3 ../find_overlaps_july_2022.py $FN1 $FN2 > Output_py_${FN1}_Output_py_${FN2}.txt; done; done &

###Del find overlaps, se cortan la columna de las secuencias en un file

Esto es para cortar la columna de las sequencias del output del find_overlps
cut -f7 Output_FM5_1A_vs_FM3_2A_with_no_r.txt > test_for_id_seqs/Col7_Output_FM5_1A_vs_FM3_2A_with_no_r.txt

1) for filename in *txt; do cut -f7 "${filename}" > ../Seqs_col/"col7_${filename}"; done &

Para agregar ">" y 1 para hacerlo formato fasta

sed 's/^/>1\n/' Col7_Output_FM5_1A_vs_FM3_2A_with_no_r.txt > Col7_Output_FM5_1A_vs_FM3_2A_with_no_r_con header.txt

###########el bueno para sed

2) for filename in *.txt; do sed 's/^/>1\n/' "${filename}" > "${filename}".fasta; done &


awk '/^>/{print ">" ++i; next}{print}' < seqs_test.fa > seqs_test_seqs_numerado.fa

checar si jala: for filename in *.fasta; do awk '/^>/{print ">" ++i; next}{print}' < "${filename}" > "${filename}"_numerado &

#######Este si jalo

3)for filename in *.fasta; do awk '/^>/{print ">" ++i; next}{print}' < "${filename}" > "${filename}"_numerado; done


##Luego se usa cd-hit para reducir las secuencias repetidas

cd-hit -i seqs_test_seqs_numerado.fa -o cd_hit_seqs_test_seqs_numerado_dos.fa

4) for filename in *_numerado; do cd-hit -i "${filename}"  -o "${filename}"_cdhit -M 0 -c 1; done &
- sino corre el cdhit con nohup usa este:
 nohup sh -c 'for filename in *_numerado; do cd-hit -i "${filename}"  -o "${filename}"_cdhit -M 0 -c 1; done' &

5) for filename in *_cdhit; do grep -v ">" "${filename}" > "${filename}"_sin_header.txt; done &

6) cp *_sin_header.txt ../Porites_output_sin_header/

7) python3 compare_seqs_2022_BUENO.py Pocillopora_output_sin_header/ -o all_seqs_Pocillopora_cd_hit_output_table.txt &


################Ultima parte otra vez

Esto es para cortar con loop

for file in *.fasta; do cut -f6 "${file}" > cd_hit_outputs/col6_$file; done

for file in *.fasta; do sed's/^/>1\n/' "${file}" > cd_hit_outputs/col6_con_header_$file; done

Para agregar ">" y 1 para hacerlo formato fasta
sed 's/^/>1\n/' Col7_Output_FM5_1A_vs_FM3_2A_with_no_r.txt > Col7_Output_FM5_1A_vs_FM3_2A_with_no_r_con header.txt


awk '/^>/{print ">" ++i; next}{print}' < seqs_test.fa > seqs_test_seqs_numerado.fa

##Esto es para todas las secs, se corta la columna con las secs, se agrupan en un archivo y se hace el cd-hit

##Una vez cortada la col de las secs se agrupan
grep "." col6_Oct_script_with_reverse_complement_* > All_samples_seqs.txt

###se reemplaza el primer caracter con un >
sed -i 's/col6/>col6/g' All_samples_seqs.txt

###Se reemplazan los : por un salto de linea
sed -i 's/:/\n/g' All_samples_seqs.txt
sed -i 's/col7/>col7/g' All_samples_test_seqs.txt

##Si se hace desde la carpeta con todos los files se hace esto:

grep "." *_numerado > ../All_samples_seqs_numeradas.txt &
sed -i 's/:>/_/g' All_samples_seqs_numeradas.txt
sed -i 's/col7_FM3_2A_only200_vs_FM3_2A_only200.txt.fasta_numerado://g' All_samples_seqs_numeradas.txt
sed -i 's/col7_/>col7_/g' All_samples_seqs_numeradas.txt

sed 's/.*://' All_samples_seqs_numeradas_round4.txt > All_samples_seqs_numeradas_round4v2.txt

##Segundo round Para hacer formato fasta de todas las secuencias separadas

sed -i 's/:>/_/g' All_samples_seqs_numeradas_round4.txt
sed 's/.*://' All_samples_seqs_numeradas_round4.txt > All_samples_seqs_numeradas_round4v2.txt
sed -i 's/col7_/>col7_/g' All_samples_seqs_numeradas_round4v2.txt



###Se corre el cd-hit (deben de ser secuencias fasta con >)

cd-hit -i All_samples_seqs.txt -o output_cd_hit_All_samples_seqs_round2.txt -M 0 &

for filename in *_numerado; do cd-hit -i "${filename}"  -o "${filename}"_cdhit -M 0; done &


for filename in *.txt; do grep -v ">" "${filename}" > "${filename}"_sin_header.txt; done &

##Despues de correr el cd_hit en todos los archivos, quitas el header de los reads:

for filename in *_cdhit; do grep -v ">" "${filename}" > "${filename}"_sin_header.txt; done &

##Los mueves a una nueva carpeta
mkdir ../Cd_hit_output_sin_header
mv *_sin_header.txt ../Cd_hit_output_sin_header/

##Corres el compare seq
nohup python3 compare_seq.py Cd_hit_output_sin_header/ -o all_seqs_cd_hit_output_sin_header_table.txt &



##Quitar header

grep -v ">" lambda_8.fa > sin_header_lambda8.fa

cut -f1 infile | grep -of- infile | sort | uniq -D


Se usa el transpose table

awk -f transpose.awk compare_test.txt > test_transpose


#cp *.fasta ../Refresh_ASQ/Coral_Dataset_200bp_Aug2022/


Para crear todos vs todos:

python3 ../make.py > ../test2
\n sleep 10m \n

###Asi corri los ultimos cd-hit

nohup cd-hit -i All_files_only_seqs_TEST2.txt -o All_files_only_seqs_TEST_OUTPUT_TEST_withc1.txt -c 1 -M 0 -T 10 &
nohup cd-hit -i All_files_only_seqs_TEST2.txt -o All_files_only_seqs_TEST_OUTPUT_TEST.txt -G 1 -M 0 -T 10 &


for filename in *.m4v; do newname=`echo $filename | sed 's/S05/S04/g'`; mv $filename $newname; done

###Para bajar genomas de ncbi usando references

wget --content-disposition "https://download.mozilla.org/?product=firefox-latest-ssl&os=linux64&lang=en-US"


nohup  ls NC_* | xargs -I{} wget -O {}.fasta 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text' &


https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference
https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/refseq_category:reference


###For jellyfish

for filename in *_10_sub_prueba_forloop.fasta; do jellyfish count -m 21 -s 100M -t 10 -C "${filename}" -o "${filename}"_jelly_output; done

for filename in *jelly_output_100mer;  do jellyfish stats "${filename}" -o "${filename}"100merstats; done

###For jelly fish bueno

###how to run jellyfish and count distinct k-mers and output stats into txt files
for filename in *.fasta; do jellyfish count -m 21 -s 100M -t 20 -C "${filename}" -o "${filename}"_jelly_output_21; done
for filename in *_jelly_output_21; do jellyfish stats "${filename}" -o "${filename}"stats; done
for filename in *_jelly_output_21stats; do head -2 "${filename}" | tail -1 > "${filename}"_distinct.txt; done

seqkit stats All_gzips_files/*.gz > All_Euks_lengths.txt &


##Para correr jellyfish comprimido BUENO
for filename in *.fna.gz; do jellyfish count -m 21 -s 100M -t 20 -C <(zcat "${filename}") -o "${filename}"_jelly_output_21; done &

##Luego se corren los demás normal

####For loop para correr todo en una carpeta

 for FN1 in *.fasta; do for FN2 in *.fasta; do python3 ../find_overlaps_july_2022.py $FN1 $FN2 > Output_py_${FN1}_Output_py_${FN2}.txt; done; done &

grep -v ">" lambda_8.fa > sin_header_lambda8.fa

###Para correr FRAP
source ~/anaconda3/etc/profile.d/conda.sh
conda activate my_env
conda activate frap_maker
snakemake --configfile sample.json -p -F


##Para correr el compare

python3 ../compare_seqs_2022_BUENO.py Output_General/ -o Tabla_General.txt &

python3 ncbi-genome-download-runner.py --formats fasta --assembly-levels complete bacteria -o Salida_bacts &

https://www.ncbi.nlm.nih.gov/nuccore/9626243

echo NC_001416 | xargs -I{} wget -O {}.fasta 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'

###Descomprimir varios a la vez:

for file in ls GCF_0*/*fna.gz; do gzip $file -d echo $file | cut -d . -f 1; done

efetch -db nucleotide -id CP007518.2 -mode text -format fasta > CP007518.2.fasta


##Este es el que funciona para toda una lista de Accesion numbers

cat test_file | xargs -I{} wget -O {}.fasta 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'

Este es para las bacterias
cat ../txid2_nucl_acc.txt | xargs -I{} wget -O {}.fasta 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'

####Para obtener los genomas de eucariotes
 zcat assembly_summary_genbank.txt.gz | grep -v "#" | cut -f20 | sed 's/\([^/]*$\)/\1\/\1_genomic.fna.gz/' | head -n10000 | xargs -I{} wget "{}" &


##crear lista de accesion numbers
esearch -db nuccore -query 'txid2 [Orgn]'|efetch -format acc > txid2_nucl_acc.txt &


##Para correr jellyfish comprimido BUENO
for filename in *.fna.gz; do jellyfish count -m 21 -s 100M -t 20 -C <(zcat "${filename}") -o "${filename}"_jelly_output_21; done &


for filename in *.fna; do jellyfish count -m 21 -s 100M -t 20 -C "${filename}" -o "${filename}"_jelly_output_21; done &

for filename in *.fna; do jellyfish count -m 31 -s 100M -t 20 -C "${filename}" -o "${filename}"_jelly_output_31; done &

for filename in *.fna; do jellyfish count -m 51 -s 100M -t 20 -C "${filename}" -o "${filename}"_jelly_output_51; done &

for filename in *.fna; do jellyfish count -m 71 -s 100M -t 20 -C "${filename}" -o "${filename}"_jelly_output_71; done &

for filename in *.fna; do jellyfish count -m 101 -s 100M -t 20 -C "${filename}" -o "${filename}"_jelly_output_101; done &


d = 0 # dimension of matrix
m = {} # the matrix
with open("infile") as fin:
    for (a, b) in [ln.strip().split(",") for ln in fin]:
        r, c = int(a), int(b)
        m[(r, c)] = 1
        d = max(r, c, d)
# printing aligned row & col no. omitted
for r in range(1, d+1):
    for c in range(1, d+1):
        # dict.get(unknown_key) returns None
        print(m.get((r, c)) or m.get((c, r)) or 0, end=" ")
        # or, by converting bool to int
        #print(int((r, c) in m or (c, r) in m), end=" ")
    print()
    

    python3 crear_matrix_FM3_vs_itself_2.py > 
    
##Para sumar los valores de la fila y columna:

awk '{ sumrows=0;
    for (i=1; i<=NF; i++) {
        sumcols[i]+= $i; sumrows+= $i 
    }; print $0, sumrows+0
}
END { for (x in sumcols)
         { printf SEP sumcols[x]+0; SEP=OFS };
     print ""
}' matrix_FM3_vs_itself_with_all_matches.txt > sum_col_and_rows_matrix_FM3_vs_itself_with_all_matches.txt


Para ver ultima columna:

cat sum_col_and_rows_mock_matrix4.txt | awk '{print ($NF)}'

Para ver ultima fila:

tail -1 sum_col_and_rows_mock_matrix4.txt

###Para cambiar de fila a columnas:

tr -s ' '  '\n'< infile > outfile

##Para cortar ultima columna, restar (-1) y dividir entre (/2)

echo -e '1 2 3 4 5 6 7' | awk '{print ($NF-1)/2}'

cat sum_col_and_rows_matrix_FM3_vs_itself_with_all_matches.txt | awk '{print ($NF-1)/2}' > Sumatoria_ultima_col_menos1_entre2_sum_col_and_rows_matrix_FM3_vs_itself_with_all_matches.txt &

Para cortar primer columna y sortearla

cut -c1 sum_col_and_rows_matrix_FM3_vs_itself_with_all_matches.txt | sort -nr | head

sort -nr Sumatoria_ultima_col_menos1_entre2_sum_col_and_rows_matrix_FM3_vs_itself_with_all_matches.txt | head -n 50

##Para unir 2 archivos en un con 2 columnas

paste inputfile1.txt inputfile2.txt > outputfile.txt

###Para checar un valor arriba de cierto valor por columna

awk '{if($4 > 0.5) printf $4"\n"}' File



Para contar el tamaño del genoma con archivos comprimidos
bioawk/./bioawk -c fastx '{ print $name, length($seq) }' < All_gzips_files/GCA_000340645.1_Trep_dent_H-22_V1_genomic.fna.gz > ../conteo.txt

for filename in *.gz; do ../bioawk/./bioawk -c fastx '{ print $name, length($seq) }' < All_gzips_files/"${filename}" > All_Euks_conteo/"${filename}"_count.txt; done &

##Este es el que usa Mark

seqkit stats All_gzips_files/*.gz > All_Euks_lengths.txt &


##Otra forma de bajar genomas:
./genome_updater.sh -o "fung_refseq_cg" -d "refseq" -g "fungi" -l "complete genome" -f "genomic.fna.gz" -t 12

./genome_updater.sh -o "fung_all_genomes" -d "refseq,genbank" -g "fungi" -f "genomic.fna.gz" -t 12

./genome_updater.sh -o "plant_all_genomes" -d "refseq,genbank" -g "plant" -f "genomic.fna.gz" -t 12

./genome_updater.sh -o "mammalian_all_genomes" -d "refseq,genbank" -g "vertebrate_mammalian" -f "genomic.fna.gz" -t 12 &

./genome_updater.sh -o "Archaea_genbank" -d "genbank" -g "archaea" -f "genomic.fna.gz" -t 12 &

###How to move the X first files

mv `ls | head -500` ./subfolder1/


###Correr frap
bash fragplot2.sh /home3/quinto/FRAP2/FRAP-basic/DB_all_bacteria_1000/ Bacteria_DB.fasta /home3/quinto/FRAP2/FRAP-basic/DS/ /home3/quinto/FRAP2/FRAP-basic/results_bact/ 

###Este es el comando que jalo de FRAP

perl jmf4.pl /home/curso/frap/FRAP/FRAP-basic/DB/NitrogenBacteria.fasta /home/curso/frap/FRAP/FRAP-basic/DS/ /home/curso/frap/results/ smalt 5000000

perl jmf4.pl /home/curso/frap/DB_bact/Bacteria_DB.fasta /home/curso/frap/DS/ /home/curso/frap/results_test/ smalt 5000000

rm GCA_903993795.1_10wheat_assembly_jagger_genomic.fna.gz_jelly_output_7*
  180  for filename in *CA_903993795.1_10wheat_assembly_jagger_genomic.fna.gz; do jellyfish count -m 71 -s 100M -t 20 -C <(zcat "${filename}") -o "${filename}"_jelly_output_71; done &




listo y checar-falta stats
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_51
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_510
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_511
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_512
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_51stats
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_710
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_711
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_712
GCA_903994175.1_10wheat_assembly_mace_genomic.fna.gz_jelly_output_713

listo y checar-falta stats
GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna.gz_jelly_output_510
GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna.gz_jelly_output_511
GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna.gz_jelly_output_710
GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna.gz_jelly_output_711
GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna.gz_jelly_output_712
GCA_903994185.1_10wheat_assembly_sy_mattis_genomic.fna.gz_jelly_output_713

listo y checar falta stats
GCA_903994195.1_10wheat_assembly_julius_genomic.fna.gz_jelly_output_510
GCA_903994195.1_10wheat_assembly_julius_genomic.fna.gz_jelly_output_511
GCA_903994195.1_10wheat_assembly_julius_genomic.fna.gz_jelly_output_710
GCA_903994195.1_10wheat_assembly_julius_genomic.fna.gz_jelly_output_711
GCA_903994195.1_10wheat_assembly_julius_genomic.fna.gz_jelly_output_712
GCA_903994195.1_10wheat_assembly_julius_genomic.fna.gz_jelly_output_713

listo y checar falta stats
GCA_903995565.1_10wheat_assembly_landmark1_genomic.fna.gz_jelly_output_510
GCA_903995565.1_10wheat_assembly_landmark1_genomic.fna.gz_jelly_output_511


##From the histo files

awk '{ print $1, $1 * $2 }' NC_063383.fasta.gz_jelly_output_71_histo > output.txt
awk '$0=$1" "$1*$2' NC_063383.fasta.gz_jelly_output_71_histo

##Segundo comando bueno
awk 'NR==FNR{sum+= $2; next} FNR==0{print $1; next} {printf("%s %0f\n",$0,$1/sum)}' row_mult_NC_063383.fasta.gz_jelly_output_71_histo row_mult_NC_063383.fasta.gz_jelly_output_71_histo

awk ' { $0=$1" "$1*$2 } NR==FNR{sum+= $2; next} FNR==0{print $1; next} {printf("%s %0f\n",$0,$1/sum)}' NC_063383.fasta.gz_jelly_output_71_histo NC_063383.fasta.gz_jelly_output_71_histo

##Tercero comando bueno
Para sacar el pi
awk -F" " '{a = log($3)/log(10); printf("%0.4f\n", a)}' output

Este es el del log
awk -F" " '{a = (log($3)/log(10)*$3); printf("%0f\n", a)}' output

##Primmer comando bueno
awk 'NR==FNR{ print $0, $1 * $2}' NC_063383.fasta.gz_jelly_output_71_histo

awk 'NR==FNR{ print $0, $1 * $2} NR==FNR{sum+= $2; next} FNR==0{print $1; next} {printf("%s %0f\n",$0,$1/sum)}'

###Este saca el log y las primeras 3 cols

awk -F" " '{a = (log($3)/log(10)*$3)} NR==FNR{ print $0, a}' output



1 184338
2 6295
6 7
8 9
10 2
11 7
```
My first command is:
```
awk 'NR==FNR{ print $0, $1 * $2}' input.file
```
which will multiply the first for the second column adding the result in a new column:
```
1 184338 184338 0.000005
2 6295 12590 0.000010
6 7 42 0.000030
8 9 72 0.000041
10 2 20 0.000051
11 7 77 0.000056
```
Then the second command I want to add is:
```
awk 'NR==FNR{sum+= $2; next} FNR==0{print $1; next} {printf("%s %0f\n",$0,$1/sum)}'
```
Which will summarize the second column and divide the value of the first column for the total of the sum in col3. This will create a new column with and output like:
```
1 184338 0.000005
2 12590 0.000010
6 42 0.000030
8 72 0.000041
10 20 0.000051
11 77 0.000056
```
The last command is 
```
awk -F" " '{a = (log($3)/log(10)*$3)} NR==FNR{ print $0, a}' output
```


1 184338 184338 0.000005 -2.65051e-05
2 6295 12590 0.000010 -5e-05
6 7 42 0.000030 -0.000135686
8 9 72 0.000041 -0.000179876
10 2 20 0.000051 -0.000218914
11 7 77 0.000056 -0.000238101



########Este es el comando bueno para el awk loop

for filename in *histo; do awk -f test.awk "${filename}" "${filename}" > "${filename}"_salida; done

#########more test.awk 

  NR==FNR{sum+=$2; next}
  {q=$1/sum; printf "%s %s %0f %s\n", $0, $1*$2, q, log(q)/log(10)*q}


Nuevo y corregido con log base 10
awk '
  NR==FNR{sum+=$2; next}
  {q=$1/sum;lq=log(q)/log(10)*q;  printf "%s %s %0f %s %0f\n", $0, $1*$2, q, lq,lq*$2}
' NC_063383.fasta.gz_jelly_output_71_histo NC_063383.fasta.gz_jelly_output_71_histo

##con logaritmo natural:
awk '
  NR==FNR{sum+=$2; next}
  {q=$1/sum;lq=log(q)*q;  printf "%s %s %0f %s %0f\n", $0, $1*$2, q, lq,lq*$2}
' NC_063383.fasta.gz_jelly_output_71_histo NC_063383.fasta.gz_jelly_output_71_histo

##Cambiando el signo de pi ln pi
awk '
  NR==FNR{sum+=$2; next}
  {q=$1/sum;lq=log(q)*q*-1;  printf "%s %s %0f %s %0f\n", $0, $1*$2, q, lq,lq*$2}
' NC_063383.fasta.gz_jelly_output_71_histo NC_063383.fasta.gz_jelly_output_71_histo

awk -F' ' '{sum+=$6;} END{print sum;}' file.txt

###############Pasos para sacar la entropia, necesitas sacar los histo de los archivos:
for filename in *_jelly_output_*; do jellyfish histo "${filename}" -o "${filename}"_histo2; done &

1.-more test_bueno.awk 

NR==FNR{sum+=$2; sum1+=$1*$2; next}
  {q=$1/sum1;lq=log(q)*q*-1;  printf "%s %s %0f %s %0f\n", $0, $1*$2, q, lq,lq*$2}


2)

for filename in *histo; do awk -f TEST_BUENO_usar_este.awk "${filename}" "${filename}" > "${filename}"_salida; done

Este es el script:

more TEST_BUENO_usar_este.awk
NR==FNR{sum+=$2; sum1+=$1*$2; next}
  {q=$1/sum1;lq=log(q)*q*-1;  printf "%s %s %0f %s %0f\n", $0, $1*$2, q, lq,lq*$2}

2.- head salida2:
1 184338 184338 0.000005 6.1843e-05 11.400012
2 12590 25180 0.000010 0.000116654 1.468673
6 42 252 0.000030 0.000316525 0.013294
8 72 576 0.000041 0.000410359 0.029546
10 20 200 0.000051 0.00050163 0.010033
11 77 847 0.000056 0.000546475 0.042079

3.-

for filename in *salida2; do awk -F' ' '{sum+=$6;} END{print sum;}' "${filename}" > "${filename}"_valor; done 

==> NC_063383.fasta.gz_jelly_output_71_histo_salida2_valor <==
12.5234



Nuevo pipeline

1)

for filename in *_21; do jellyfish histo "${filename}" -o "${filename}"histo; done
for filename in *_jelly_output_*; do jellyfish histo "${filename}" -o "${filename}"_histo; done &
2)

for filename in *_101; do jellyfish histo "${filename}" -o "${filename}"_histo; done &
for filename in *_11; do jellyfish histo "${filename}" -o "${filename}"_histo; done &
for filename in *_21; do jellyfish histo "${filename}" -o "${filename}"_histo; done &
for filename in *_31; do jellyfish histo "${filename}" -o "${filename}"_histo; done &
for filename in *_51; do jellyfish histo "${filename}" -o "${filename}"_histo; done &
for filename in *_71; do jellyfish histo "${filename}" -o "${filename}"_histo; done &

3)

for filename in *histo; do awk -f TEST_BUENO_usar_este.awk "${filename}" "${filename}" > "${filename}"_salida; done &

4)

for filename in *salida; do awk -F' ' '{sum+=$6;} END{print sum;}' "${filename}" > "${filename}"_valor; done 

5)

grep "." *valor > sample_valor.txt


####Para bajar sequencias de la compañia:

wget -nH -m --ftp-user=kcarter --ftp-password=M0HH0x7Z8W ftp://igm-storage.ucsd.edu/230106_A01535_0254_BHNWW3DRX2/

##este funciona para uno:
trimmomatic SE -phred33 10_S5_L001_R1_001.fastq.gz 10_S5_L001_R1_001_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


trimmomatic-0.39.jar SE -threads 40 -phred33 /media/drogon/narnia/mrojas/MG_COVID19 /singles_raw/$1 "`basename $1 .fastq`_trim.fastq" ILLUMINACLIP:/media/drogon/narnia/mlittle/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:36

ls -1 *.fastq > SE_files.txt

more run_all_trim.sh 

#!/usr/bin/bash
trimmomatic SE -phred33 $1 "`basename $1 .fastq.gz`.trimmomatic_out.fastq" ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


cat SE_files.txt | xargs -n 1 nohup bash run_all_trim.sh


###Para FRAP
source ~/.bashrc

conda activate frap_maker
Hay que cambiar el nombre de los archivos:
find . -name '*fastq' -exec bash -c ' mv $0 ${0/\_001/.sra_1}' {} \;

snakemake --configfile sample2.json -p -F --cores 16

Replace names:
for file in *.fastq; do mv "$file" "${file/fastq.gz_trim/sra_1}"; done

###Este es el que uso (hay que cambiar el sample.json por tu version)

if [ ! -e sge_err ]; then mkdir -p sge_err sge_out;fi; snakemake -s Snakefile --configfile sample2.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 6 --cores 1000 --latency-wait 60 -F

if [ ! -e sge_err ]; then mkdir -p sge_err sge_out;fi; snakemake -s Snakefile --configfile sample3.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 16 --cores 1000 --latency-wait 60 -F

##Este si jala
if [ ! -e sge_err ]; then mkdir -p sge_err sge_out; fi; snakemake -s Snakefile6 --configfile test_holly2.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 6 --cores 1000 --latency-wait 60 -F

###No lo he probado
if [ ! -e sge_err ]; then mkdir -p sge_err sge_out;fi; snakemake -s Snakefile --configfile sample3.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 6 --cores 1000 --latency-wait 60 -F


if [ ! -e sge_err ]; then mkdir -p sge_err sge_out;fi; snakemake -s Snakefile4 --configfile sample5.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 16 --cores 1000 --latency-wait 60 -F

Para checarlo  se usa:
qstat -u \*



if [ ! -e sge_err ]; then mkdir -p sge_err sge_out;fi; snakemake -s Snakefile4 --configfile sample5_1.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 16 --cores 1000 --latency-wait 60 -F

El 5_1 es con la nueva base de datos.


###Este si jalo:
if [ ! -e sge_err ]; then mkdir -p sge_err sge_out;fi; snakemake -s Snakefile4 --configfile sample6.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 16 --cores 1000 --latency-wait 60 -F

##Para correr el de Tapir:
if [ ! -e sge_err ]; then mkdir -p sge_err sge_out; fi; snakemake -s Snakefile5 --configfile sample1_Tapir_virome.json --cluster 'qsub -cwd -e sge_err -o sge_out -V -q default -pe smp {threads}' --local-cores 6 --cores 1000 --latency-wait 60 -F



###Comando para quitar espacios de bases de datos:

awk -vRS='\n\n' -vORS='' 1 all_virus_ML_80char.fna > all_virus_ML_80char_sin_espacios.fna &


###Formato databases:
Despues de quitar espacios se corre:
awk '
  /^>/ {printf("%s%s\n", !c?"":RS, $0);c=0;next} 
  {printf("%s", $0);c++}
  END {printf ORS}' myInputFile

  taraoceanvirome_v1.fna

  awk '
  /^>/ {printf("%s%s\n", !c?"":RS, $0);c=0;next} 
  {printf("%s", $0);c++}
  END {printf ORS}' taraoceanvirome_v1.fna > taraoceanvirome_v1_w_format.fna

  test_virus_refseq.fasta > test_virus_refseq_w_format.fasta


##Para cambiar caracteres en file:
perl -pe 's/\\n/\n/g' file
>_all_sequences_A


###Controls info:
Genome ID	562.50451
Genome Name	Escherichia coli strain SMEc189
Taxon ID	562
Sequence ID	562.50451.con.0278
Accession	BJKC01000277
Sequence Type	contig
Topology	linear
Description	BJKC01000277.1
GC Content	47.52
Length	9291
Release Date	5/16/2019
Version	1
Insert Date	9/30/2019
Last Modified	9/30/2019

####Este es para crear 200bp de genomas

##se descomprimen genomas:
gzip -d *.gz

##Se quitan los headers:

for FILE in *.fna; do grep -v ">" "$FILE" > "${FILE%.fna}_sin_header.fna"; done

##Se corre el script

for FILE in *.fasta; do cat "$FILE" | tr -d '\n' | dd cbs=200 conv=unblock | nl -n ln | tr '\t' '\n' | paste -d '>\n' /dev/null - - > "${FILE%.fa}_substring_last.fa"; done &


########New FRAP:
awk '                                       
FNR == NR       {S[$1]
                 next
                }
/>/             {HD = substr ($0, 2)
                 L = 0
                 next
                }
                {
                 for (s in S)   {IX = index($0, s)
                                 if (IX) print HD
                                }
                 L += length
                }
' seq3.txt all.txt > rudi_approach3.txt

sed -i 's/ /_/g' rudi_approach3.txt > rudi_approach3_no_spaces.txt

sort rudi_approach3.txt | uniq -c > rudi_approach3_hits.txt

awk '{s+=$1}END{print s}' rudi_approach3_hits.txt

###Este es para lo de arriba en script:

#!/usr/bin/bash 

awk '                                       
FNR == NR       {S[$1]
                 next
                }
/>/             {HD = substr ($0, 2)
                 L = 0
                 next
                }
                {
                 for (s in S)   {IX = index($0, s)
                                 if (IX) print HD
                                }
                 L += length
                }
' seq3.txt all2.txt | sed  's/ /_/g' | sort | uniq -c >  test_gui5.txt

awk '                    
        NR > max { max=NR }
        { tot+=$1; v[NR]=$1; d[NR]=$2 }
        END { for (i=1; i<=max; i++) { print v[i]*100/tot,v[i],d[i] } }

' test_gui5.txt > test_gui6.txt

###Esta es otra version:

awk '
FNR==NR { pat[$1]; next }
function prt_if_pat() {
  if (n++) for (p in pat) if (x=index(buf,p)) { print name }
}
/^>/ { prt_if_pat(); name=$0; buf=""; next }
{ buf=(buf $0) }
END { prt_if_pat() }
' seq4.txt all2.txt | sed  's/ /_/g' | sort | uniq -c





####Para crear 100,000 de Jason seqs
for FILE in *.fasta; do grep -v ">" "$FILE" > "${FILE%.fasta}_sin_header.fasta"; done &
 
for FILE in *_sin_header.fasta; do cat "$FILE" | tr -d '\n' | dd cbs=200 conv=unblock | nl -n ln | tr '\t' '\n' | paste -d '>\n' /dev/null - - > "${FILE%.fasta}_substring_last.fa"; done &
  
 for filename in *_substring_last.fa; do cat "${filename}" |\awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' |\shuf |\head -n 100000 |\awk 'BEGIN{FS="\t"}{printf("%s\n%s\n",$1,$2)}' |\awk '/^>/{print ">" ++i; next}{print}' > "${filename}"_100000sub.fasta; done &

###Para crear random chunks de una base de datos (all5.txt):
#Quitar headers
awk '!/^>/' all5.txt > all5.txt_sinheaders.fa
#Quitar saltos de linea
tr -d '\n' < all5.txt_sinheaders.fa > all5.txt_sinheaders_nisaltos.fa

#Se corre:
for N in {1..10000}; do bash random_chunks.sh >> Exp2_salida_largav6.txt | cat Exp2_salida_largav6.txt  | tr -d '\n' | dd cbs=200 conv=unblock | nl -n ln | tr '\t' '\n' | paste -d '>\n' /dev/null - - > salida_larga_chunks_numerado_10000v6.txt | sed -n '/>10002/q;p'; done  &

more random_chunks:
n=$(stat -c "%s" newfile_lambas_bueno.fa)
r=$(shuf -i1-"$((n-200+1))" -n1)
< newfile_lambas_bueno.fa tail -c+"$r" | head -c200

#Para contar un patron repetido:
grep -EIho 'ggttatatttg' E_coli_with_only_first11bp.fa | wc -l

#Para contar las diferentes pares de bases en un genoma
 sed 's/\(.\)/\1\n/g' Escherichia_coli_strain_SMEc189_ONLY_ATCG_sin_header.fasta | sort | uniq -c

 #Contar todo lo que no sea AGTC
  grep -v ">" Escherichia_coli_strain_SMEc189_ONLY_ATCG.fasta | grep -E -i -o -v "G|C|T|A|N" | wc -l

  
  #######Otra version del nuevo FRAP (es mas rapida)

  #!/bin/bash

awk '
  FNR==NR {
     re = re ((re) ? "|" : "") $0
     next
 }
 FNR>1 && $2~re {print $1}
' 45_S45_L002_R1_good_out_no_headers_subsample2000.fasta RS='>' FS='\n' all4.txt


###Para controles, cambiar caracteres en alguna posicion:

grep -v ">" Escherichia_coli_strain_SMEc189.fa > Escherichia_coli_strain_SMEc189_sin_header.fa

more Escherichia_coli_strain_SMEc189_sin_header.fa
   awk '{
  for (i = 0; ++i <= NF;)
    ++c % n || $i = v 
  }1' OFS= FS= n=10 v=n Escherichia_coli_strain_SMEc189_sin_header.fa > E_coli_original_+_n_every_10bp_sin_header.fa

  ##Es el mismo del de arriba para cambiar characteres en alguna posicion:
  
  awk '{
  for (i = 0; ++i <= NF;)
    ++c % n || $i = v 
  }1' OFS= FS= n=3 v=A infile



##Este for loop cambia el formato de los genomas a fasta:
for filename in *.fasta; do awk '                      
  /^>/ {printf("%s%s\n", !c?"":RS, $0);c=0;next} 
  {printf("%s", $0);c++}
  END {printf ORS}' "${filename}" > "Fasta_formato_${filename}".txt; done

for filename in *.fna; do awk '                      
  /^>/ {printf("%s%s\n", !c?"":RS, $0);c=0;next} 
  {printf("%s", $0);c++}
  END {printf ORS}' "${filename}" > "Fasta_formato_${filename}".txt; done

###Este es el bueno para contar todos los caracteres de un genoma:

grep -v ">" NC_000852.fasta | sed 's/\(.\)/\1\n/g' | sort | uniq -c



###Este es el bueno para contar las bases que no son ACGT
grep -v ">" file.fasta  | tr -d 'ATCG' | tr -d '\n'| wc -c

###Usando un loop
for FILE in *.txt; do grep -v ">" "$FILE" | tr -d 'atcgATCG' | tr -d '\n'| wc -c > "${FILE%.fasta}_counts.txt"; done &

###Este es para sacar el porcentage de las bases raras:

IGNORECASE=1 gawk '!/^>/ { c=gsub(/[^atcgATCG]/, "&"); printf("%d out of %d = %.2f%\n", c, length, (c/length)*100)}' Fasta_formato_NC_000852.fasta.txt

##Haciendo loop
for filename in Fasta_formato*; do IGNORECASE=1 gawk '!/^>/ { c=gsub(/[^atcgATCG]/, "&"); printf("%d out of %d = %.2f%\n", c, length, (c/length)*100)}' "${filename}" > "Porcentaje_${filename}".txt; done

##Scrutinier

for filename in Fasta_formato*; do awk -F'[^acgtACGT]' '!/^>/{n=NF-1; printf("%d out of %d = %.2f%\n", n, length, 100*n/length)}' "${filename}" > "Scrutiner_${filename}".txt; done


grep "." *.fasta.txt.txt > All_percentage_non_DNA.txt
grep "." Porcentaje_Fasta_formato_NC_0* > All_virus_porcentaje.txt &
grep "." Scrutiner_Fasta_formato* > All_virus_Scrutiner.txt &
more Fasta_formato_NC_038319.fasta.txt_counts
grep "." *.fasta.txt_counts > All_virus_non_dna.txt &



##Control para comparar:
./sequence-stats/src/sequence-stats -b NC_000866.fasta
./sequence-stats/src/sequence-stats -b NC_000874.fasta


#Quitar filas con zeros:

#opcion 1
awk '$2' Archivos_y_valor_Porcentaje.txt > Archivos_y_valor_Porcentaje_sin_zeros.txt

# opcion 2
cut -f2  Archivos_y_valor_Porcentaje.txt | grep -v -e "^0$" > Archivos_y_valor_Porcentaje_sin_zeros.txt


Para contar posiciones de los characteres que no son acgt
grep -v ">" E_coli_with_y_every20bp_dos_veces.fa |  tr -d 'ATCGatcg' | tr -d '\n'| grep -o . |  grep -n ''


#####Script para encontrar las bases ambiguas en genomas
for filename in *.fasta; do IGNORECASE=1 gawk '
    !/^>/ { len=length
       printf("Seq%d\n", ++seq)
       while (match($0, /([^atcg])/, a)) {
         printf("%s position %d out of %d \n", a[1], a[1, "start"], len)
         sub(a[1], "a")
       }
}' "${filename}" > "Positions_${filename}".txt; done

for filename in *.txt; do IGNORECASE=1 gawk '
    !/^>/ { len=length
       printf("Seq%d\n", ++seq)
       while (match($0, /([^atcg])/, a)) {
         printf("%s position %d out of %d \n", a[1], a[1, "start"], len)
         sub(a[1], "a")
       }
}' "${filename}" > "Positions_${filename}".txt; done &


###Last path
/home3/quinto/Para_genomas_All/controls/Genome_controls_w_and_wo_C


##Para quitar sequencias mas grandes de 100bp

cat input.fa| \
awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])<=100) {printf("%s\n%s\n",L[0],L[1]);}}'

for filename in *.fna.txt; do awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])<=100) {printf("%s\n%s\n",L[0],L[1]);}}' "${filename}" > "less_than_100_${filename}".txt; done &

###Para borrar archivos vacios
find . -type f -empty -print -delete


###Comandos para bajar todos los genomas

./genome_updater.sh -o "Viral_genbank_refseq_complete_genomes" -d "refseq,genbank" -g "viral" -l "complete genome" -f "genomic.fna.gz" -t 12 &

./genome_updater.sh -o "Viral_genbank_refseq_complete_genomes" -d "refseq,genbank" -g "viral" -l "complete genome" -f "genomic.fna.gz" -t 12 &

./genome_updater.sh -o "Bacteria_genbank_refseq_complete_genomes" -d "refseq,genbank" -g "bacteria" -l "complete genome" -f "genomic.fna.gz" -t 16 &

