#######
######
######    First I am going to move into the Level_3 folder with all of the raw count files
######    Then I am going to take the FILE_SAMPLE_MAP.txt file and rearrange the columns
######
######    Added 4-4-2016
######    I have had problems with passing the variables into sed, and awk
######    I am going to instead just write this for the Sking Cancer dataset and use sed to change the future scripts
######
echo -e "\n\tThis is a script written to parse the RNAseqV2 data from TCGA and put the data into a spreadsheet format that can be used in Excel."
echo -e "\n\tThis script was written in April of 2016. Please be advised that the structure of data output from TCGA may have changed and this script may need to be modified."
cd RNASeqV2/*RNASeqV2/Level_3
cp ../../../FILE_SAMPLE_MAP.txt .
awk '{print $2 "\t" $1}' FILE_SAMPLE_MAP.txt > make.headers.sh
sed -i '/filename/d' make.headers.sh
sed -i 's/^/echo -e "GENE\\t/g' make.headers.sh
sed -i 's/\t/" > /g' make.headers.sh
sed -i 's/$/.HEADER/g' make.headers.sh
sh make.headers.sh
######
######
######
######   Now we have the header files, so we need to cat the count files with the headers
echo -e "\n\tSample name is being added to the top line of each file using information from FILE_SAMPLE_MAP.txt provided."
grep "junction_quantification" FILE_SAMPLE_MAP.txt > name.junction_quantification.samples.txt
grep "rsem.genes.results" FILE_SAMPLE_MAP.txt > name.rsem.genes.results.samples.txt
grep "rsem.isoforms.results" FILE_SAMPLE_MAP.txt > name.rsem.isoforms.results.samples.txt
grep "rsem.genes.normalized_results" FILE_SAMPLE_MAP.txt > name.rsem.genes.normalized_results.samples.txt
grep "rsem.isoforms.normalized_results" FILE_SAMPLE_MAP.txt > name.rsem.isoforms.normalized_results.samples.txt
grep "bt.exon_quantification.txt" FILE_SAMPLE_MAP.txt > name.bt.exon_quantification.samples.txt
awk '{print "cat " $1 ".HEADER " $1 " > " $2 ".junction_quantification.txt"}' name.junction_quantification.samples.txt > name.junction_quantification.sh
sh name.junction_quantification.sh
awk '{print "cat " $1 ".HEADER " $1 " > " $2 ".rsem.genes.results"}' name.rsem.genes.results.samples.txt > name.rsem.genes.results.sh
sh name.rsem.genes.results.sh
awk '{print "cat " $1 ".HEADER " $1 " > " $2 ".rsem.isoforms.results"}' name.rsem.isoforms.results.samples.txt > name.rsem.isoforms.results.sh
sh name.rsem.isoforms.results.sh
awk '{print "cat " $1 ".HEADER " $1 " > " $2 ".rsem.genes.normalized_results"}' name.rsem.genes.normalized_results.samples.txt > name.rsem.genes.normalized_results.sh
sh name.rsem.genes.normalized_results.sh
awk '{print "cat " $1 ".HEADER " $1 " > " $2 ".rsem.isoforms.normalized_results"}' name.rsem.isoforms.normalized_results.samples.txt > name.rsem.isoforms.normalized_results.sh
sh name.rsem.isoforms.normalized_results.sh
awk '{print "cat " $1 ".HEADER " $1 " > " $2 ".bt.exon_quantification.txt"}' name.bt.exon_quantification.samples.txt > name.bt.exon_quantification.samples.sh
sh name.bt.exon_quantification.samples.sh
echo -e "\n\tFile Labeling complete. Organizing data"
#######
######
######    Now I am going to make folders and move the labeled samples just for simplicity
######    I am also going to delete the header files because they are annoying.
rm *.HEADER
mkdir Labeled_Files
mv TCGA* Labeled_Files
mkdir Unlabeled_Files
cd Unlabeled_Files
mkdir Junction_Quantification
mkdir RSEM_Genes_Results
mkdir RSEM_Isoforms_Results
mkdir RSEM_Genes_NORMALIZED_Results
mkdir Isoforms_Normalized_Results
mkdir Exon_Quantification
cd ../
mv *junction_quantification.txt Unlabeled_Files/Junction_Quantification/
mv *rsem.genes.results Unlabeled_Files/RSEM_Genes_Results/
mv *rsem.isoforms.results Unlabeled_Files/RSEM_Isoforms_Results/
mv *rsem.genes.normalized_results Unlabeled_Files/RSEM_Genes_NORMALIZED_Results/
mv *rsem.isoforms.normalized_results Unlabeled_Files/Isoforms_Normalized_Results/
mv *bt.exon_quantification.txt Unlabeled_Files/Exon_Quantification/
mkdir tmp
mv name* tmp/
mv make.headers.sh tmp/
cd Labeled_Files
mkdir Junction_Quantification
mkdir RSEM_Genes_Results
mkdir RSEM_Isoforms_Results
mkdir RSEM_Genes_NORMALIZED_Results
mkdir Isoforms_Normalized_Results
mkdir Exon_Quantification
mv *junction_quantification.txt Junction_Quantification/
mv *rsem.genes.results RSEM_Genes_Results/
mv *rsem.isoforms.results RSEM_Isoforms_Results/
mv *.rsem.genes.normalized_results RSEM_Genes_NORMALIZED_Results/
mv *.rsem.isoforms.normalized_results Isoforms_Normalized_Results/
mv *.bt.exon_quantification.txt Exon_Quantification/
####### For Housekeeping, I am goin to jump up a level and make a tmp file to hold all of the other scripts that I made that will not be needed later, but do not junk things up.
##### This step will probably be moved to the end when I am finished.
## cd ..
## mkdir tmp
## mv name* tmp/
## mv make.headers.sh tmp/
######
######     New steps of analysis added 4-1-2016
######
######     Now I have made all of the labeled files, but I need to combine them.
######
######
######     At this point I have moved into the Labeled_Files directory
###### First make the file of the gene names
cd RSEM_Genes_NORMALIZED_Results/
ls *.rsem.genes.normalized_results > tmp1.rsem.genes.normalized_results.txt
head -n1 tmp1.rsem.genes.normalized_results.txt > grab.genes.sh
sed -i 's/^/cut -f1 /g' grab.genes.sh
sed -i 's/$/ > genes.txt/g' grab.genes.sh
sh grab.genes.sh
##### Now cut the counts from every single file in the directory
paste tmp1.rsem.genes.normalized_results.txt tmp1.rsem.genes.normalized_results.txt > tmp2.rsem.genes.normalized_results.sh
sed -i 's/^/cut -f2 /g' tmp2.rsem.genes.normalized_results.sh
sed -i 's/\t/ > /g' tmp2.rsem.genes.normalized_results.sh
sed -i 's/$/.counts.txt/g' tmp2.rsem.genes.normalized_results.sh
sh tmp2.rsem.genes.normalized_results.sh
######  Now I need to paste all of the temp count files together into one giant file.
ls *counts.txt > tmp3.rsem.genes.normalized_results.txt
tr "\n" "\t" < tmp3.rsem.genes.normalized_results.txt > paste.rsem.genes.normalized_results.sh
sed -i 's/^/paste genes.txt /g' paste.rsem.genes.normalized_results.sh
sed -i 's/$/ > BRCA.Combined.rsem.genes.normalized_results/g' paste.rsem.genes.normalized_results.sh
sh paste.rsem.genes.normalized_results.sh
mkdir ../../../../../Results
cp BRCA.Combined.rsem.genes.normalized_results ../../../../../Results/BRCA.Combined.rsem.genes.normalized_results
rm *counts.txt
echo -e "\n\trsem.genes.normalized_results file compiled."
######
######
###### That was the most important file, but I am going to do the same thing for the other files too.
######
cd ../Exon_Quantification
ls *bt.exon_quantification.txt > tmp1.bt.exon_quantification.txt
head -n1 tmp1.bt.exon_quantification.txt > grab.genes.sh
sed -i 's/^/cut -f1 /g' grab.genes.sh
sed -i 's/$/ > genes.txt/g' grab.genes.sh
sh grab.genes.sh
##### Now cut the counts from every single file in the directory
paste tmp1.bt.exon_quantification.txt tmp1.bt.exon_quantification.txt > tmp2.bt.exon_quantification.sh
sed -i 's/^/cut -f2 /g' tmp2.bt.exon_quantification.sh
sed -i 's/\t/ > /g' tmp2.bt.exon_quantification.sh
sed -i 's/$/.counts.txt/g' tmp2.bt.exon_quantification.sh
sh tmp2.bt.exon_quantification.sh
######  Now I need to paste all of the temp count files together into one giant file.
ls *counts.txt > tmp3.bt.exon_quantification.txt
tr "\n" "\t" < tmp3.bt.exon_quantification.txt > paste.bt.exon_quantification.sh
sed -i 's/^/paste genes.txt /g' paste.bt.exon_quantification.sh
sed -i 's/$/ > BRCA.Combined.bt.exon_quantification.txt/g' paste.bt.exon_quantification.sh
sh paste.bt.exon_quantification.sh
cp BRCA.Combined.bt.exon_quantification.txt ../../../../../Results/BRCA.Combined.bt.exon_quantification.txt
rm *counts.txt
######
######   Now for the Isoforms_Normalized_Results files
######
echo -e "\n\tbt.exon_quantification file compiled."
cd ../Isoforms_Normalized_Results
ls *rsem.isoforms.normalized_results > tmp1.rsem.isoforms.normalized_results.txt
head -n1 tmp1.rsem.isoforms.normalized_results.txt > grab.genes.sh
sed -i 's/^/cut -f1 /g' grab.genes.sh
sed -i 's/$/ > genes.txt/g' grab.genes.sh
sh grab.genes.sh
##### Now cut the counts from every single file in the directory
paste tmp1.rsem.isoforms.normalized_results.txt tmp1.rsem.isoforms.normalized_results.txt > tmp2.rsem.isoforms.normalized_results.sh
sed -i 's/^/cut -f2 /g' tmp2.rsem.isoforms.normalized_results.sh
sed -i 's/\t/ > /g' tmp2.rsem.isoforms.normalized_results.sh
sed -i 's/$/.counts.txt/g' tmp2.rsem.isoforms.normalized_results.sh
sh tmp2.rsem.isoforms.normalized_results.sh
######  Now I need to paste all of the temp count files together into one giant file.
ls *counts.txt > tmp3.rsem.isoforms.normalized_results.txt
tr "\n" "\t" < tmp3.rsem.isoforms.normalized_results.txt > paste.rsem.isoforms.normalized_results.sh
sed -i 's/^/paste genes.txt /g' paste.rsem.isoforms.normalized_results.sh
sed -i 's/$/ > BRCA.Combined.rsem.isoforms.normalized_results.txt/g' paste.rsem.isoforms.normalized_results.sh
sh paste.rsem.isoforms.normalized_results.sh
cp BRCA.Combined.rsem.isoforms.normalized_results.txt ../../../../../Results/BRCA.Combined.rsem.isoforms.normalized_results.txt
rm *counts.txt
######
######
######   Now for the Junction_Quantification
######
echo -e "\n\trsem.isoforms.normalized_results file compiled."
cd ../Junction_Quantification
ls *junction_quantification.txt > tmp1.junction_quantification.txt
head -n1 tmp1.junction_quantification.txt > grab.genes.sh
sed -i 's/^/cut -f1 /g' grab.genes.sh
sed -i 's/$/ > genes.txt/g' grab.genes.sh
sh grab.genes.sh
##### Now cut the counts from every single file in the directory
paste tmp1.junction_quantification.txt tmp1.junction_quantification.txt > tmp2.junction_quantification.sh
sed -i 's/^/cut -f2 /g' tmp2.junction_quantification.sh
sed -i 's/\t/ > /g' tmp2.junction_quantification.sh
sed -i 's/$/.counts.txt/g' tmp2.junction_quantification.sh
sh tmp2.junction_quantification.sh
######  Now I need to paste all of the temp count files together into one giant file.
ls *counts.txt > tmp3.junction_quantification.txt
tr "\n" "\t" < tmp3.junction_quantification.txt > paste.junction_quantification.sh
sed -i 's/^/paste genes.txt /g' paste.junction_quantification.sh
sed -i 's/$/ > BRCA.Combined.junction_quantification.txt/g' paste.junction_quantification.sh
sh paste.junction_quantification.sh
cp BRCA.Combined.junction_quantification.txt ../../../../../Results/BRCA.Combined.junction_quantification.txt
rm *.counts.txt
######
######
######   Now for the RSEM_Genes_Results UN-normalized
######
echo -e "\n\tjunction_quantification file compiled."
cd ../RSEM_Genes_Results
ls *rsem.genes.results > tmp1.rsem.genes.results.txt
head -n1 tmp1.rsem.genes.results.txt > grab.genes.sh
sed -i 's/^/cut -f1 /g' grab.genes.sh
sed -i 's/$/ > genes.txt/g' grab.genes.sh
sh grab.genes.sh
##### Now cut the counts from every single file in the directory
paste tmp1.rsem.genes.results.txt tmp1.rsem.genes.results.txt > tmp2.rsem.genes.results.sh
sed -i 's/^/cut -f2 /g' tmp2.rsem.genes.results.sh
sed -i 's/\t/ > /g' tmp2.rsem.genes.results.sh
sed -i 's/$/.counts.txt/g' tmp2.rsem.genes.results.sh
sh tmp2.rsem.genes.results.sh
######  Now I need to paste all of the temp count files together into one giant file.
ls *counts.txt > tmp3.rsem.genes.results.txt
tr "\n" "\t" < tmp3.rsem.genes.results.txt > paste.rsem.genes.results.sh
sed -i 's/^/paste genes.txt /g' paste.rsem.genes.results.sh
sed -i 's/$/ > BRCA.Combined.rsem.genes.results.txt/g' paste.rsem.genes.results.sh
sh paste.rsem.genes.results.sh
cp BRCA.Combined.rsem.genes.results.txt ../../../../../Results/BRCA.Combined.rsem.genes.results.txt
rm *.counts.txt
#####
#####
#####    I think this if the final round
#####
echo -e "\n\trsem.genes.results file compiled."
cd ../RSEM_Isoforms_Results
ls *rsem.isoforms.results > tmp1.rsem.isoforms.results.txt
head -n1 tmp1.rsem.isoforms.results.txt > grab.genes.sh
sed -i 's/^/cut -f1 /g' grab.genes.sh
sed -i 's/$/ > genes.txt/g' grab.genes.sh
sh grab.genes.sh
##### Now cut the counts from every single file in the directory
paste tmp1.rsem.isoforms.results.txt tmp1.rsem.isoforms.results.txt > tmp2.rsem.isoforms.results.sh
sed -i 's/^/cut -f2 /g' tmp2.rsem.isoforms.results.sh
sed -i 's/\t/ > /g' tmp2.rsem.isoforms.results.sh
sed -i 's/$/.counts.txt/g' tmp2.rsem.isoforms.results.sh
sh tmp2.rsem.isoforms.results.sh
######  Now I need to paste all of the temp count files together into one giant file.
ls *counts.txt > tmp3.rsem.isoforms.results.txt
tr "\n" "\t" < tmp3.rsem.isoforms.results.txt > paste.rsem.isoforms.results.sh
sed -i 's/^/paste genes.txt /g' paste.rsem.isoforms.results.sh
sed -i 's/$/ > BRCA.Combined.rsem.isoforms.results.txt/g' paste.rsem.isoforms.results.sh
sh paste.rsem.isoforms.results.sh
cp BRCA.Combined.rsem.isoforms.results.txt ../../../../../Results/BRCA.Combined.rsem.isoforms.results.txt
rm *.counts.txt
cd ../../../../../Results
echo -e "\n\trsem.isoforms.results file compiled."
echo -e "\n\tMoving onto statistical analysis of Hypoxia genes."
######
######    Now I am going to try to make the scatter plots of HIF vs CA9, HIF vs PGK1, and CA9 vs PGK1
######    First I will grab the hypoxia genes from the normalized file, BRCA.Combined.rsem.genes.normalized_results
head -n1 BRCA.Combined.rsem.genes.normalized_results >> BRCA.hypoxia.genes.txt
awk '$1 == "HIF1A|3091"' BRCA.Combined.rsem.genes.normalized_results >> BRCA.hypoxia.genes.txt
awk '$1 == "CA9|768"' BRCA.Combined.rsem.genes.normalized_results >> BRCA.hypoxia.genes.txt
awk '$1 == "PGK1|5230"' BRCA.Combined.rsem.genes.normalized_results >> BRCA.hypoxia.genes.txt
#####   Now I need to turn it sideways using gawk
awk '
{
        for (i=1; i<=NF; i++) {
                a[NR,i] = $i
        }
}
NF>p { p = NF}
END {
        for(j=1; j<=p; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){
                        str=str"\t"a[i,j];
                }
                print str
        }
}' BRCA.hypoxia.genes.txt > BRCA.hypoxia.genes.Flipped.txt
sed -i 's/HIF1A|3091/HIF1A/g' BRCA.hypoxia.genes.Flipped.txt
sed -i 's/CA9|768/CA9/g' BRCA.hypoxia.genes.Flipped.txt
sed -i 's/PGK1|5230/PGK1/g' BRCA.hypoxia.genes.Flipped.txt
######
######    Now I am going to use echo commands to build the R submission file.
######
######
echo "library(ggplot2)" >> BRCA.hypoxia.analysis.R
echo "hypoxia_genes <- read.table(\"BRCA.hypoxia.genes.Flipped.txt\", header = TRUE, row.names = 1)" >> BRCA.hypoxia.analysis.R
echo "png(\"HIF1AvCA9_BRCA.png\", width=8, height=8, units = \"in\", res=600)" >> BRCA.hypoxia.analysis.R
echo "eqn_hifvca9 <- function(hypoxia_genes){" >> BRCA.hypoxia.analysis.R
echo "  m <- lm(CA9~HIF1A, data=hypoxia_genes);" >> BRCA.hypoxia.analysis.R
echo "  eq <- substitute(italic(Y) == a ~ + ~ b %*% italic(X)*\",\" ~~italic(r)^2~\"=\"~r2," >> BRCA.hypoxia.analysis.R
echo "                   list(a = format(coef(m)[1], digits = 4)," >> BRCA.hypoxia.analysis.R
echo "                        b = format(coef(m)[2], digits = 4)," >> BRCA.hypoxia.analysis.R
echo "                        r2 = format(summary(m)\$r.squared, digits = 4)))" >> BRCA.hypoxia.analysis.R
echo "  as.character(as.expression(eq));" >> BRCA.hypoxia.analysis.R
echo "}" >> BRCA.hypoxia.analysis.R
#### Calculate the regression again, but this time save it as function so it can be graphed.
echo "reg_hifvca9 <- lm(CA9~HIF1A, data=hypoxia_genes)" >> BRCA.hypoxia.analysis.R
#### Create a ggplot vector that we will use for graphing
echo "hifvca9 <- ggplot(hypoxia_genes, aes(HIF1A, CA9))" >> BRCA.hypoxia.analysis.R
#### Now graph the points
echo "hifvca9 + geom_point(colour = \"blue\") +" >> BRCA.hypoxia.analysis.R
echo "  theme_bw() +" >> BRCA.hypoxia.analysis.R
echo "  geom_abline(intercept = reg_hifvca9\$coefficients[1], slope = reg_hifvca9\$coefficients[2]) +" >> BRCA.hypoxia.analysis.R
echo "  geom_text(x = 9000, y = 10000, label = eqn_hifvca9(hypoxia_genes), parse = TRUE, size = 5) +" >> BRCA.hypoxia.analysis.R
echo "  geom_text(x = 9000, y = 12000, label = \"CA9 as a function of HIF1A in BRCA\", size = 5)" >> BRCA.hypoxia.analysis.R
#### color = blue is changing the color of the dots
#### theme_bw() changes the background to white instead of gray
#### geom_text writes something on the graph
##### size = 5 is the font size
#### X = and Y = are the locations on the graph. This will change for each graph.
#### Now you need to dev.off() to turn off the device taking into into the PNG. It will save it as a PNG.
echo "dev.off()" >> BRCA.hypoxia.analysis.R
#####
#####
#### Below is a failed attempt at combining everything into one line.
#echo "library(ggplot2)" >> BRCA.hypoxia.analysis.R
#echo "hypoxia_genes <- read.table(\"BRCA.hypoxia.genes.Flipped.txt\", header = TRUE, row.names = 1)" >> BRCA.hypoxia.analysis.R
#echo "hifvca9 <- ggplot(hypoxia_genes, aes(HIF1A, CA9))" >> BRCA.hypoxia.analysis.R
#echo "png(\"HIF1AvCA9_BRCA.png\", width=8, height=8, units = \"in\", res=600)" >> BRCA.hypoxia.analysis.R
#echo "eqn_hifvca9 <- function(hypoxia_genes){ m <- lm(CA9~HIF1A, data=hypoxia_genes); eq <- substitute(italic(Y) == a ~ + ~ b %*% italic(X)*\",\" ~~italic(r)^2~\"=\"~r2, list(a = format(coef(m)[1], digits = 4), b = format(coef(m)[2], digits = 4), r2 = format(summary(m)\$r.squared, digits = 4))) as.character(as.expression(eq))}" >> BRCA.hypoxia.analysis.R
#echo "reg_hifvca9 <- lm(CA9~HIF1A, data=hypoxia_genes)" >> BRCA.hypoxia.analysis.R
#echo "hifvca9 + geom_point(colour = \"blue\") + theme_bw() + geom_abline(intercept = reg_hifvca9\$coefficients[1], slope = reg_hifvca9\$coefficients[2]) + geom_text(x = 9000, y = 5000, label = eqn_hifvca9(hypoxia_genes), parse = TRUE, size = 5) + geom_text(x = 9000, y = 5500, label = \"CA9 as a function of HIF1A in BRCA\", size = 5)" >> BRCA.hypoxia.analysis.R
#echo "dev.off()" >> BRCA.hypoxia.analysis.R
######
######
######
######   I am going to try it with keeping the spaces and see what happens.
######
######
#echo "library(ggplot2)" >> BRCA.hypoxia.analysis.R
echo "hypoxia_genes <- read.table(\"BRCA.hypoxia.genes.Flipped.txt\", header = TRUE, row.names = 1)" >> BRCA.hypoxia.analysis.R
echo "png(\"HIF1AvPGK1_BRCA.png\", width=8, height=8, units = \"in\", res=600)" >> BRCA.hypoxia.analysis.R
echo "eqn_hifvpgk1 <- function(hypoxia_genes){" >> BRCA.hypoxia.analysis.R
echo "  m <- lm(PGK1~HIF1A, data=hypoxia_genes);" >> BRCA.hypoxia.analysis.R
echo "  eq <- substitute(italic(Y) == a ~ + ~ b %*% italic(X)*\",\" ~~italic(r)^2~\"=\"~r2," >> BRCA.hypoxia.analysis.R
echo "                   list(a = format(coef(m)[1], digits = 4)," >> BRCA.hypoxia.analysis.R
echo "                        b = format(coef(m)[2], digits = 4)," >> BRCA.hypoxia.analysis.R
echo "                        r2 = format(summary(m)\$r.squared, digits = 4)))" >> BRCA.hypoxia.analysis.R
echo "  as.character(as.expression(eq));" >> BRCA.hypoxia.analysis.R
echo "}" >> BRCA.hypoxia.analysis.R
#### Calculate the regression again, but this time save it as function so it can be graphed.
echo "reg_hifvpgk1 <- lm(PGK1~HIF1A, data=hypoxia_genes)" >> BRCA.hypoxia.analysis.R
#### Create a ggplot vector that we will use for graphing
echo "hifvpgk1 <- ggplot(hypoxia_genes, aes(HIF1A, PGK1))" >> BRCA.hypoxia.analysis.R
#### Now graph the points
echo "hifvpgk1 + geom_point(colour = \"blue\") +" >> BRCA.hypoxia.analysis.R
echo "  theme_bw() +" >> BRCA.hypoxia.analysis.R
echo "  geom_abline(intercept = reg_hifvpgk1\$coefficients[1], slope = reg_hifvpgk1\$coefficients[2]) +" >> BRCA.hypoxia.analysis.R
echo "  geom_text(x = 9000, y = 42000, label = eqn_hifvpgk1(hypoxia_genes), parse = TRUE, size = 5) +" >> BRCA.hypoxia.analysis.R
echo "  geom_text(x = 9000, y = 46000, label = \"PGK1 as a function of HIF1A in BRCA\", size = 5)" >> BRCA.hypoxia.analysis.R
#### color = blue is changing the color of the dots
#### theme_bw() changes the background to white instead of gray
#### geom_text writes something on the graph
##### size = 5 is the font size
#### X = and Y = are the locations on the graph. This will change for each graph.
#### Now you need to dev.off() to turn off the device taking into into the PNG. It will save it as a PNG.
echo "dev.off()" >> BRCA.hypoxia.analysis.R
echo "png(\"CA9vPGK1_BRCA.png\", width=8, height=8, units = \"in\", res=600)" >> BRCA.hypoxia.analysis.R
echo "eqn_ca9vpgk1 <- function(hypoxia_genes){" >> BRCA.hypoxia.analysis.R
echo "  m <- lm(CA9~PGK1, data=hypoxia_genes);" >> BRCA.hypoxia.analysis.R
echo "  eq <- substitute(italic(Y) == a ~ + ~ b %*% italic(X)*\",\" ~~italic(r)^2~\"=\"~r2," >> BRCA.hypoxia.analysis.R
echo "                   list(a = format(coef(m)[1], digits = 4)," >> BRCA.hypoxia.analysis.R
echo "                        b = format(coef(m)[2], digits = 4)," >> BRCA.hypoxia.analysis.R
echo "                        r2 = format(summary(m)\$r.squared, digits = 4)))" >> BRCA.hypoxia.analysis.R
echo "  as.character(as.expression(eq));" >> BRCA.hypoxia.analysis.R
echo "}" >> BRCA.hypoxia.analysis.R
#### Calculate the regression again, but this time save it as function so it can be graphed.
echo "reg_ca9vpgk1 <- lm(CA9~PGK1, data=hypoxia_genes)" >> BRCA.hypoxia.analysis.R
#### Create a ggplot vector that we will use for graphing
echo "ca9vpgk1 <- ggplot(hypoxia_genes, aes(PGK1, CA9))" >> BRCA.hypoxia.analysis.R
#### Now graph the points
echo "ca9vpgk1 + geom_point(colour = \"blue\") +" >> BRCA.hypoxia.analysis.R
echo "  theme_bw() +" >> BRCA.hypoxia.analysis.R
echo "  geom_abline(intercept = reg_ca9vpgk1\$coefficients[1], slope = reg_ca9vpgk1\$coefficients[2]) +" >> BRCA.hypoxia.analysis.R
echo "  geom_text(x = 27000, y = 4770, label = eqn_ca9vpgk1(hypoxia_genes), parse = TRUE, size = 5) +" >> BRCA.hypoxia.analysis.R
echo "  geom_text(x = 27000, y = 5500, label = \"CA9 as a function of PGK1 in BRCA\", size = 5)" >> BRCA.hypoxia.analysis.R
#### color = blue is changing the color of the dots
#### theme_bw() changes the background to white instead of gray
#### geom_text writes something on the graph
##### size = 5 is the font size
#### X = and Y = are the locations on the graph. This will change for each graph.
#### Now you need to dev.off() to turn off the device taking into into the PNG. It will save it as a PNG.
echo "dev.off()" >> BRCA.hypoxia.analysis.R
R CMD BATCH BRCA.hypoxia.analysis.R
