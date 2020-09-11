
zenity --info --title="Consenso all'analisi" --window-icon=/dati/server/Software/appGEN.png --text=" Vuoi lanciare il programma appGEN_AMPLICON? 
programma scritto da Davide Gentilini 
         gentilini.davide@gmail.com" --ok-label="Certo che si" 
descrizione=$(zenity --forms --title="Nome Esperimento" --text="Assegna un Codice Esperimento" \
   --add-entry="NOME" \
   --add-entry="COGNOME" \
   --add-entry="Username" \
   --add-list="Insert your choice" --list-values 'Pannello1|Pannello2' \
   --add-calendar="data esperimento")





#################################################################################################
PATHDATA=$(zenity --file-selection --directory --title="****INDICA PATH-DATA***"  --text="INDICA PATH-DATA")
echo $PATHDATA >> ISTRUZIONI.txt

PATHPIPE=$(zenity --file-selection --directory --title="****INDICA PATH-appGEN***"  --text="INDICA PATH-appGEN")
echo $PATHPIPE >> ISTRUZIONI.txt

TARGET=$(zenity --file-selection --title="*** SELEZIONA FILE .BED TARGET ***" )
echo $TARGET >> ISTRUZIONI.txt

GENREF=$(zenity --file-selection --title="*** SELEZIONA FILE GENOMA DI RIFERIMENTO --- genome.fa --- ***" )
echo $GENREF >> ISTRUZIONI.txt

PATHVCF_DATABASE=$(zenity --file-selection --directory --title="*** Select Directory for VCF-History-DATABASE ***" )
echo $PATHVCF_DATABASE >> ISTRUZIONI.txt


THRESHOLD=$(zenity --forms --title="INDICA THRESHOLD" --text="INDICA THRESHOLD" \
   --add-entry="THRESHOLD")

ANNOT=$(zenity --file-selection --title="Select a VCF FILE with ANNOTATIONs")
echo $ANNOT >> ISTRUZIONI.txt


###################################################################################################
######################################################################################################

cd $PATHDATA/
echo ' Rename Files '
for i in `seq 10`; do rename 's/-/_/' *; done
echo ' only underscores  '


echo 'Alignement step starting'
for file in *R1_001.fastq.gz; 
do
$PATHPIPE/bwa-0.7.10/bwa mem -R "@RG\tID:IDa\tSM:SM\tPL:Illumina" -t 4 $GENREF $PATHDATA/"$file" $PATHDATA/"$(basename "$file" _R1_001.fastq.gz)_R2_001.fastq.gz" > $PATHDATA/"$(basename "$file" _R1_001.fastq.gz).PRE.sam" ;done

#######################################################################################################
for file in *.PRE.sam; 
do
java -jar -Xms4g -Xmx4g $PATHPIPE/picard-tools-1.119/SortSam.jar INPUT=$PATHDATA/"$file" OUTPUT=$PATHDATA/"$(basename "$file" .PRE.sam).SORTED.bam" SORT_ORDER=coordinate ;done

rm *.PRE.sam


for file in *.SORTED.bam; 
do
java -jar -Xms4g -Xmx4g $PATHPIPE/picard-tools-1.119/CollectAlignmentSummaryMetrics.jar R=$GENREF I=$PATHDATA/"$file" O=$PATHDATA/"$(basename "$file" .SORTED.bam).METRIX.txt" ;done


for file in *.SORTED.bam; 
do
java -jar -Xms4g -Xmx4g $PATHPIPE/picard-tools-1.119/BuildBamIndex.jar INPUT=$PATHDATA/"$file" ;done


for file in *.SORTED.bam; 
do
java -d64 -Xms4g -Xmx4g -jar $PATHPIPE/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $GENREF -I $PATHDATA/"$file" -L:capture,BED $TARGET -o $PATHDATA/"$(basename "$file" .SORTED.bam).targets.list" ;done

echo 'INDEL-Realigning step starting'
for file in *.SORTED.bam; 
do
java -d64 -Xms4g -Xmx4g -jar $PATHPIPE/GenomeAnalysisTK.jar -T IndelRealigner -I $PATHDATA/"$file" -targetIntervals $PATHDATA/"$(basename "$file" .SORTED.bam).targets.list" -R $GENREF -known $PATHPIPE/MANIFEST-BED-DATABASE/1000G_phase1.indels.hg19.vcf -known $PATHPIPE/MANIFEST-BED-DATABASE/dbsnp_135.hg19.vcf --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4 -o $PATHDATA/"$(basename "$file" .SORTED.bam).aligned.bam" ;done
echo 'INDEL-Realigned step is done!!!'


echo 'BaseRecalibrator step starting'

for file in *.aligned.bam; 
do
java -d64 -Xms4g -Xmx4g -jar $PATHPIPE/GenomeAnalysisTK.jar -T BaseRecalibrator --maximum_cycle_value 600 -I $PATHDATA/"$file" -R $GENREF -knownSites $PATHPIPE/MANIFEST-BED-DATABASE/dbsnp_135.hg19.vcf -cov QualityScoreCovariate -cov CycleCovariate -o $PATHDATA/"$(basename "$file" .aligned.bam).grp" ;done



for file in *.aligned.bam; 
do 
java -d64 -Xms4g -Xmx4g -jar $PATHPIPE/GenomeAnalysisTK.jar -T PrintReads -I $PATHDATA/"$file" -R $GENREF -BQSR $PATHDATA/"$(basename "$file" .aligned.bam).grp" -o $PATHDATA/"$(basename "$file" .aligned.bam).bam" ;done

rm -f, yes | rm *.SORTED*
rm -f, yes | rm *.aligned*
rm -f, yes | rm *.REALIGNED*
rm -f, yes | rm *.DEDUP*



for file in *.bam; 
do
java -jar -Xms4g -Xmx4g $PATHPIPE/GenomeAnalysisTK.jar -T UnifiedGenotyper -I $PATHDATA/"$file" -R $GENREF --dbsnp $PATHPIPE/MANIFEST-BED-DATABASE/dbsnp_135.hg19.vcf -glm BOTH -o $PATHDATA/"$(basename "$file" .bam).row.vcf" -L:capture,BED $TARGET -stand_call_conf 30.0 -stand_emit_conf 10.0;done 


echo 'Variant Call step is done!!!'


for file in *.row.vcf; 
do
java -jar $PATHPIPE/GenomeAnalysisTK.jar -T VariantFiltration -R $GENREF -V $PATHDATA/"$file" --filterExpression "QD < 2.0" --filterName "LowQualityByDepth_min" --filterExpression "MQRankSum < 12.5" --filterName "LowMQRankSum" --filterExpression "MQRankSum < 12.5" --filterName "LowMQRankSum" --filterExpression "DP < 20" --filterName "COPERTURA_MIN_20" --filterExpression "GQ < 30" --filterName "LowGQX" -o $PATHDATA/"$(basename "$file" .row.vcf).vcf" ;done 

rm -f, yes | rm *.row.vcf*

#####################################################################################################
# creazione vcf multiplo

rm -f, yes | rm MULTIPLE.vcf
rm -f, yes | rm ANNOTATION_INDEX.vcf
for file in *.vcf; do                                                                                                                                                                                                                                                                      grep -v -h '^#' $file > "$(basename "$file" .vcf).ID"; 
done

ls *.ID|xargs -I% sed -i 's/$/\t%/' %

cat *.ID > temp0

awk -F'\t' '{print $1"-"$2"-"$4"-"$5}' temp0 > temp1

awk -F '\t' '{print $1}' temp1 | sort | uniq -c | sort -nr > temp1x

awk -F' ' '{print $2, $1}' temp1x > temp1xx


sort -k1 temp1 > temp1.1
sort -k1 temp1xx > temp1xx.1


join -a 1 -a 1 temp1.1 temp1xx.1 > temp1.2a
awk -F"\t" '!_[$1,$2,$4,$5]++' temp1.2a > temp1.2
paste temp1 temp0 > temp3

sort -k1 temp3 > temp4
sort -k1 temp1.2 > temp5

join -a 1 -a 1 temp4 temp5 | awk -F' ' '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$1}' > temp1000
awk -v OFS="\t" '$1=$1' temp1000 >  MULTIPLE.DGV


awk -F"\t" '!_[$1,$2,$4,$5]++' temp0 > ANNOTATION_INDEX.DGV

rm -f, yes | rm *.ID
rm -f, yes | rm temp*


##########################################################################################################################
#####################################################################################################
#####################################################################################################
echo 'Creation of Historical Database of Genomic variants'

cp $PATHDATA/*.vcf $PATHVCF_DATABASE/

rm $PATHVCF_DATABASE/MULTIPLE.vcf
rm $PATHVCF_DATABASE/COLLAPSED_MULTIPLE.vcf

cd $PATHVCF_DATABASE/


rm -f, yes | rm MULTIPLE.vcf
rm -f, yes | rm ANNOTATION_INDEX.vcf
for file in *.vcf; do                                                                                                                                                                                                                                                                      grep -v -h '^#' $file > "$(basename "$file" .vcf).ID"; 
done


ls *.ID|xargs -I% sed -i 's/$/\t%/' %

cat *.ID > temp0

awk -F'\t' '{print $1"-"$2"-"$4"-"$5}' temp0 > temp1

awk -F '\t' '{print $1}' temp1 | sort | uniq -c | sort -nr > temp1x

awk -F' ' '{print $2, $1}' temp1x > temp1xx




sort -k1 temp1 > temp1.1
sort -k1 temp1xx > temp1xx.1



join -a 1 -a 1 temp1.1 temp1xx.1 > temp1.2a
awk -F"\t" '!_[$1,$2,$4,$5]++' temp1.2a > temp1.2
paste temp1 temp0 > temp3

sort -k1 temp3 > temp4
sort -k1 temp1.2 > temp5

join -a 1 -a 1 temp4 temp5 | awk -F' ' '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$1}' > temp1000
awk -v OFS="\t" '$1=$1' temp1000 >  MULTIPLE.DGV


awk -F"\t" '!_[$1,$2,$4,$5]++' temp0 > ANNOTATION_INDEX.DGV

rm -f, yes | rm *.ID
rm -f, yes | rm temp*


############################################################################################################################



cd $PATHDATA


############################################################################################################################
mkdir RECUPERI

echo 'COVERAGE Analysis step starting'


zenity progress bar dialog
(
echo "20" ; sleep 1
echo "# CALCOLO COPERTURE E PROFONDITA'" ; sleep 1

cd $PATHDATA
depth=$(ls | grep .deep.txt)
if [ -z "${depth}" ];
then
for file in *.bam; 
do
samtools depth $PATHDATA/"$file" > $PATHDATA/"$(basename "$file" .bam).deep.txt" ;done
	
else
echo "$depth found."
fi

echo "40" ; sleep 1
echo "#under threshold--calculation " ; sleep 1

for file in *.deep.txt; 
do

awk -F"\t" '$3 < '$THRESHOLD' { print $1"\t"$2"\t"$2 }' $PATHDATA/"$file" > $PATHDATA/"$(basename "$file" .deep.txt).underT.bed" ;done

echo "50" ; sleep 1
echo "#CALCOLO COORDINATE REGIONI LETTE " ; sleep 1

TOTBED=$(ls | grep .TOT.bed)
if [ -z "${TOTBED}" ];
then
for file in *.bam; 
do
bamToBed -i $PATHDATA/"$file" > $PATHDATA/"$(basename "$file" .bam).TOT.bed" ;done
	
else
echo "$TOTBED found."
fi

echo "80" ; sleep 1
echo "# region in target calculation" ; sleep 1
for file in *.TOT.bed; 
do
subtractBed -a $TARGET -b $PATHDATA/"$file" > $PATHDATA/"$(basename "$file" .TOT.bed).zero.bed" ;done

for file in *.underT.bed; 
do
intersectBed -a $TARGET -b $PATHDATA/"$file" > $PATHDATA/"$(basename "$file" .underT.bed).underTHRESHOLD.bed" ;done


for file in *.underTHRESHOLD.bed; 
do
cat $PATHDATA/"$file" $PATHDATA/"$(basename "$file" .underTHRESHOLD.bed).zero.bed"> $PATHDATA/RECUPERI/"$(basename "$file" .underTHRESHOLD.bed).RECOVER.bed" ;done
echo "100" ; sleep 1
) |
zenity --progress --auto-close "PREMI OK PER PROSEGUIRE CON ANALISI" 
--text "CARICAMENTO DEI DATI ESEGUITO PROCEDO CON ANALISI" --percentage=0
rm -f, yes | rm *.deep*
rm -f, yes | rm *.underTHRESHOLD.bed
rm -f, yes | rm *.TOT.bed
rm -f, yes | rm *.zero.bed
rm -f, yes | rm *.underT.bed
rm -f, yes | rm *.grp
rm -f, yes | rm *.targets.list
rm -f, yes | rm *.METRIX.txt


mkdir SNPEFF

for file in *.DGV; 
do
java -jar $PATHPIPE/snpEff/snpEff.jar -v hg19 $PATHDATA/"$file" > $PATHDATA/SNPEFF/"$(basename "$file" .vcf).SNPEFF.DGV";done 



for file in *.bam; 
do
java -Xms2g -Xmx4g -cp  $PATHPIPE/NGSRICH/NGSrich_0.7.8/bin NGSrich evaluate -r $PATHDATA/"$file" -u hg19 -a  $PATHPIPE/MANIFEST-BED-DATABASE/refGene.txt -t $TARGET -o $PATHDATA/outputCOVERAGE -p 20 -h 200 --no-details; done 


mkdir BAM
mkdir VCF


mv -f, yes | mv *.bam $PATHDATA/BAM/
mv -f, yes | mv *.bai $PATHDATA/BAM/
mv -f, yes | mv *.vcf $PATHDATA/VCF/
mv -f, yes | mv *.DGV $PATHDATA/VCF/
cd $PATHDATA/VCF/



                                                                                                                                                                                                                                                                grep -v -h '^#' $ANNOT > temp0x; 



awk -F'\t' '{print $1"-"$2"-"$4"-"$5}' temp0x > temp0

paste temp0 temp0x > temp

awk '!seen[$1]++ >= 1' temp > TEMPO1



grep -v -h '^#' MULTIPLE.DGV > temp1x; 



awk -F'\t' '{print $1"-"$2"-"$4"-"$5}' temp1x > temp1
paste temp1 temp1x > temp2


sort -k1 temp2 > temp5
sort -k1 TEMPO1 > temp6

join -a 1 -a 1 temp5 temp6 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.13,2.14,2.15,2.16 > temp10

sed "s/ /\t/g" temp10 > temp11


sed "s/ /\t/g" temp11 > "$(basename "MULTIPLE" .vcf).COMMENTED.vcf"

rm -f, yes | rm temp*
rm -f, yes | rm TEMPO1
rm -f, yes | rm risultato

for i in `seq 1`; do rename 's/.DGV/.vcf/' *; done


zenity --info --title="Consenso all'analisi" --window-icon=/dati/server/Software/appGEN.png --text=" SEI SODDISFATTO DELL'ANALISI? 
Ti auguro una buona giornata 
         gentilini.davide@gmail.com" --ok-label="CHIUDIMI" 

		 
		












