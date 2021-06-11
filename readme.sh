

Mehdi,
 
Here is the excel document. Check the last 2 sheets: Low Freq neg EUR and the one that says 0.08%.
 
Thanks for the help
 
Gabriel


in analysis/02_merge/
merge.688.nomissing.withRSID.vcf.gz

need to remove 5 duplicate samples:


51-07536	27
51-07540	31
51-08273	301
51-16513	648
51-16520	655

zcat merge.688.nomissing.withRSID.vcf.gz | head -99 | grep ^# | bgzip -c > merge.683.nomissing.withRSID.vcf.gz
zcat merge.688.nomissing.withRSID.vcf.gz | sed '1,99d' |  cut -f -26,28-30,32-300,302-647,649-654,656- | bgzip -c  >> merge.683.nomissing.withRSID.vcf.gz

#################################################
#
# Haplotype Analysis
#
#################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --snps /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR_snp_ids.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.vcf.gz
bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.vcf.gz

tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

shapeit4 --input /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap --region chr1 --output /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased.vcf.gz

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased

awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased.ped > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased.set.missing.zero.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased.map | awk '{print $2,$4}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased.set.missing.zero.map

java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_Low_freq_neg_EUR.recode.phased

#################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --snps /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08_snpids.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.vcf.gz
bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

shapeit4 --input /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap --region chr1 --output /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.vcf.gz

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased

awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.ped > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.recode.phased.set.missing.zero.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.map | awk '{print $2,$4}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.recode.phased.set.missing.zero.map

java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.ped --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683-0.08.recode.phased.ped

#################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/249_snp_id_positions.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

shapeit4 --input /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap --region chr1 --output /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased.vcf.gz

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recode12 --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased
plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeAD --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeAD.phased
plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeA --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs//PKLR683_249_snps.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.ped

awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.ped

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.map | awk '{print $2,$4}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.map

cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.map

java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased

#################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/179_snp_id_biallelic_positions.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

shapeit4 --input /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap --region chr1 --output /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased.vcf.gz

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased

awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased.ped > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased.set.missing.zero.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased.map | awk '{print $2,$4}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased.set.missing.zero.map

java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_179_snps_biallelic.recode.phased

#################################################

#################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz \
	--remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt \
	--positions /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/69_AA_variants.txt \
	--recode \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--out /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps

bgzip -c /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.vcf > /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

shapeit4 --input /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.vcf.gz \
	--map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap \
	--region chr1 \
	--log /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.phased.shapeit4.log \
	--output /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.phased.vcf.gz

plink --vcf /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.phased.vcf.gz \
	--make-bed \
	--biallelic-only strict \
	--out /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.phased

plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.recode.phased \
	--recode \
	--out /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.plink.recode.phased

# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recode12 --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeAD --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeAD.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeA --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs//PKLR683_249_snps.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.ped

awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.plink.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.plink.recode.phased.set.missing.zero.ped

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.map | awk '{print $2,$4}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.map

cat /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.plink.recode.phased.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /mnt/lab-thein/Fayaz/04-Gabriel/05-PKLR683_69_SNPs/PKLR683_69_snps.plink.recode.phased.set.missing.zero.map

java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased

#################################################

#################################################################################
#
# Find unique haplotypes out of 214 phased variants
#
#################################################################################

cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased.raw | sed '1,1d' | cut -d" " -f7- > tmp

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' tmp > tmp1

cat tmp1 | sed 's/ //g' > tmp2 
cat tmp2 | sort -k1 | uniq -c > tmp3
cat tmp3 | sort -rk1 > tmp4
cat tmp4 | wc -l #188 unique haplotypes


java -Xms64g -Xmx256g -jar /home/seifuddinft/splitstree5/jars/splitstree5.jar

#################################################################################
#
# Haplotype Analysis using haplo.stats
#
#################################################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/69_AA_variants.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.vcf.gz

shapeit4 --input /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.vcf.gz \
	--map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap \
	--region chr1 \
	--log /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.phased.shapeit4.log \
	--output /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.phased.vcf.gz

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_69_snps.recode.phased.vcf.gz \
	--make-bed \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.recode.phased

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.recode.phased \
	--recode \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.plink.recode.phased

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.plink.recode.phased.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.plink.recode.phased.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.plink.recode.phased.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/06-PKLR_69_SNPs_haplostats/PKLR683_variants.plink.recode.phased.haplo.stats.map

#################################################################################
#
# Haplotype Analysis using haplo.stats, 1000 genomes, phase3, 69 variants
#
#################################################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/69_AA_variants.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz

perl /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/vcf_to_ped_convert.pl \
	-vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz \
	-sample_panel_file /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel \
	-region 1:155262400-155262414 \
	-population GBR \
	-output_ped ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.ped \
	-output_info ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.map \
	-output_dir /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3 \
	-min_maf 0 \
	-base_format 'letters'

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz \
	--make-bed \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="EUR"){print $1,$1}}' > integrated_call_samples_v3.20130502.EUR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.EUR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.EUR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.EUR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.EUR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.EUR.haplo.stats.input.map



echo -e "

REQUIRED ARGUMENTS
            -vcf                Path to a locally or remotely accessible vcf file.
                                The vcf file must be compressed by bgzip and indexed by tabix if it is a remote file.
                                The vcf format is a tab format for presenting variation sites and 
                                genotypes data and is described at http://vcftools.sourceforge.net/specs.html.
                                This tool takes both vcf4.0 and vcf4.1 format files.
            -sample_panel_file  Path to a locally or remotely accessible sample panel file, listing all individuals (first column)
                                and their population (second column)
            -region             Chromosomal region in the format of chr:start-end (e.g. 1:1000000-100500) or chr (e.g. 1)
            -population         A population name, which must appear in the second column of the sample panel file.
                                Can be specified more than once for multiple populations.

OPTIONAL ARGUMENTS
            -tabix              Path to the tabix executable; default is to search the path for 'tabix'
                                tabix is not required if the vcf file is uncompressed and locally accessible
            -output_ped         Name of the output ped file (linkage pedigree file);
                                default is region.ped (e.g. 1_100000-100500.ped)
            -output_info        Name of the output info file (marker information file);
                                default is region.info (e.g. 1_1000000-100500.info)
            -output_dir         Name of a directory to place the output_ped and output_info files
            -min_maf            Only include variations with a minor allele_frequency greater than or equal to this value
            -max_maf            Only include variations with a minor allele_frequency less than or equal to this value
            -base_format        Either 'letter' or 'number'. Genotypes in the ped file can be coded either ACGT or 1-4
                                where 1=A, 2=C, 3=G, T=4.  Default is 'number'.
            -help               Print out help menu

OUTPUT FILES
            The file formats of the linkage pedigree and marker information files are described at
" > /dev/null

#################################################################################
#
# Haplotype Analysis using haplo.stats, 1000 genomes, phase3, 66 SNPs
#
#################################################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/66_AA_variants.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz

perl /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/vcf_to_ped_convert.pl \
	-vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz \
	-sample_panel_file /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel \
	-region 1:155262400-155262414 \
	-population GBR \
	-output_ped ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.ped \
	-output_info ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.map \
	-output_dir /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3 \
	-min_maf 0 \
	-base_format 'letters'

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz \
	--make-bed \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode

###AFR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="AFR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.haplo.stats.input.map

###EUR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="EUR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.haplo.stats.input.map



echo -e "

REQUIRED ARGUMENTS
            -vcf                Path to a locally or remotely accessible vcf file.
                                The vcf file must be compressed by bgzip and indexed by tabix if it is a remote file.
                                The vcf format is a tab format for presenting variation sites and 
                                genotypes data and is described at http://vcftools.sourceforge.net/specs.html.
                                This tool takes both vcf4.0 and vcf4.1 format files.
            -sample_panel_file  Path to a locally or remotely accessible sample panel file, listing all individuals (first column)
                                and their population (second column)
            -region             Chromosomal region in the format of chr:start-end (e.g. 1:1000000-100500) or chr (e.g. 1)
            -population         A population name, which must appear in the second column of the sample panel file.
                                Can be specified more than once for multiple populations.

OPTIONAL ARGUMENTS
            -tabix              Path to the tabix executable; default is to search the path for 'tabix'
                                tabix is not required if the vcf file is uncompressed and locally accessible
            -output_ped         Name of the output ped file (linkage pedigree file);
                                default is region.ped (e.g. 1_100000-100500.ped)
            -output_info        Name of the output info file (marker information file);
                                default is region.info (e.g. 1_1000000-100500.info)
            -output_dir         Name of a directory to place the output_ped and output_info files
            -min_maf            Only include variations with a minor allele_frequency greater than or equal to this value
            -max_maf            Only include variations with a minor allele_frequency less than or equal to this value
            -base_format        Either 'letter' or 'number'. Genotypes in the ped file can be coded either ACGT or 1-4
                                where 1=A, 2=C, 3=G, T=4.  Default is 'number'.
            -help               Print out help menu

OUTPUT FILES
            The file formats of the linkage pedigree and marker information files are described at
" > /dev/null

#################################################################################
#
# Haplotype Analysis using haplo.stats, 255 variants, 683 sickle-cell subjects
#
#################################################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/255_PKLR_variants_positions.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

# shapeit4 phases SNPs and small INDELs only

shapeit4 --input /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/01-genetic_maps_b37_by_chromosome/chr1.b37.gmap --region chr1 --log /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.phased.shapeit4.log --sequencing --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --pbwt-modulo 0.000000005 --pbwt-mdr 0.0 --pbwt-mac 0 --ibd2-maf 0.0 --ibd2-mdr 0.0 --ibd2-count 0 --output /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.phased.vcf.gz

echo -e "
		Basic options:
		  --help                                Produce help message
		  --seed arg (=15052011)                Seed of the random number generator
		  -T [ --thread ] arg (=1)              Number of thread used

		Input files:
		  -I [ --input ] arg                    Genotypes to be phased in VCF/BCF 
				                        format
		  -H [ --reference ] arg                Reference panel of haplotypes in 
				                        VCF/BCF format
		  -S [ --scaffold ] arg                 Scaffold of haplotypes in VCF/BCF 
				                        format
		  -M [ --map ] arg                      Genetic map
		  -R [ --region ] arg                   Target region
		  --use-PS arg                          Informs phasing using PS field from 
				                        read based phasing
		  --sequencing                          Default parameter setting for 
				                        sequencing data (e.g. high variant 
				                        density)

		MCMC parameters:
		  --mcmc-iterations arg (=5b,1p,1b,1p,1b,1p,5m)
				                        Iteration scheme of the MCMC
		  --mcmc-prune arg (=0.999)             Pruning threshold in genotype graphs

		PBWT parameters:
		  --pbwt-modulo arg (=0.02)             Storage frequency of PBWT indexes in cM
				                        (i.e. storage every 0.02 cM by default)
		  --pbwt-depth arg (=4)                 Depth of PBWT indexes to condition on
		  --pbwt-mac arg (=2)                   Minimal Minor Allele Count at which 
				                        PBWT is evaluated
		  --pbwt-mdr arg (=0.5)                 Maximal Missing Data Rate at which PBWT
				                        is evaluated

		IBD2 parameters:
		  --ibd2-length arg (=3)                Minimal size of IBD2 tracks for 
				                        building copying constraints
		  --ibd2-maf arg (=0.01)                Minimal Minor Allele Frequency for 
				                        variants to be considered in the IBD2 
				                        mapping
		  --ibd2-mdr arg (=0.5)                 Maximal Missing data rate for variants 
				                        to be considered in the IBD2 mapping
		  --ibd2-count arg (=150)               Minimal number of filtered variants in 
				                        IBD2 tracks
		  --ibd2-output arg                     Output all IBD2 constraints in the 
				                        specified file (useful for debugging!)

		HMM parameters:
		  -W [ --window ] arg (=2.5)            Minimal size of the phasing window in 
				                        cM
		  --effective-size arg (=15000)         Effective size of the population

		Output files:
		  -O [ --output ] arg                   Phased haplotypes in VCF/BCF format
		  --log arg                             Log file
" > /dev/null

# plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
# plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recode12 --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeAD --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeAD.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeA --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs//PKLR683_249_snps.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.ped

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.ped

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.map | awk '{print $2,$4}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.map

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.map

# java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

# plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_255_variants.recode.phased.vcf.gz \
	--make-bed \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased \
	--recode \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased.genotypes

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased.genotypes.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased.genotypes.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased.genotypes.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/11-PKLR_255_variants_haplostats/PKLR683_220_variants.recode.phased.genotypes.haplo.stats.input.map

#################################################

#################################################################################
#
# Haplotype Analysis using haplo.stats, 1000 genomes, phase3, 220 variants
#
#################################################################################

vcftools --gzvcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/12-PKLR_220_variants_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/12-PKLR_220_variants_1000genomes_phase3/PKLR_220_variants_positions.txt --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/12-PKLR_220_variants_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_220_variants_b37

bgzip -c /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz

perl /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/vcf_to_ped_convert.pl \
	-vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz \
	-sample_panel_file /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel \
	-region 1:155262400-155262414 \
	-population GBR \
	-output_ped ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.ped \
	-output_info ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.map \
	-output_dir /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3 \
	-min_maf 0 \
	-base_format 'letters'

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz \
	--make-bed \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode

###AFR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="AFR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.haplo.stats.input.map

###EUR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="EUR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.haplo.stats.input.map



echo -e "

REQUIRED ARGUMENTS
            -vcf                Path to a locally or remotely accessible vcf file.
                                The vcf file must be compressed by bgzip and indexed by tabix if it is a remote file.
                                The vcf format is a tab format for presenting variation sites and 
                                genotypes data and is described at http://vcftools.sourceforge.net/specs.html.
                                This tool takes both vcf4.0 and vcf4.1 format files.
            -sample_panel_file  Path to a locally or remotely accessible sample panel file, listing all individuals (first column)
                                and their population (second column)
            -region             Chromosomal region in the format of chr:start-end (e.g. 1:1000000-100500) or chr (e.g. 1)
            -population         A population name, which must appear in the second column of the sample panel file.
                                Can be specified more than once for multiple populations.

OPTIONAL ARGUMENTS
            -tabix              Path to the tabix executable; default is to search the path for 'tabix'
                                tabix is not required if the vcf file is uncompressed and locally accessible
            -output_ped         Name of the output ped file (linkage pedigree file);
                                default is region.ped (e.g. 1_100000-100500.ped)
            -output_info        Name of the output info file (marker information file);
                                default is region.info (e.g. 1_1000000-100500.info)
            -output_dir         Name of a directory to place the output_ped and output_info files
            -min_maf            Only include variations with a minor allele_frequency greater than or equal to this value
            -max_maf            Only include variations with a minor allele_frequency less than or equal to this value
            -base_format        Either 'letter' or 'number'. Genotypes in the ped file can be coded either ACGT or 1-4
                                where 1=A, 2=C, 3=G, T=4.  Default is 'number'.
            -help               Print out help menu

OUTPUT FILES
            The file formats of the linkage pedigree and marker information files are described at
" > /dev/null

#################################################################################
#
# Haplotype Analysis using haplo.stats, 99 variants, 683 sickle-cell subjects
# DARC/ACKR1 gene
#	
#################################################################################

workingdir="/home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes"

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --positions $workingdir/DARC_ACKR1_99_variants_by_position_hg19.txt --recode --out $workingdir/DARC_683_hg19

bgzip -c $workingdir/DARC_683_hg19.recode.vcf > $workingdir/DARC_683_hg19.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf $workingdir/DARC_683_hg19.recode.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

# shapeit4 phases SNPs and small INDELs only

shapeit4 --input $workingdir/DARC_683_hg19.recode.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/01-genetic_maps_b37_by_chromosome/chr1.b37.gmap --region chr1 --log $workingdir/DARC_683_hg19.recode.phased.shapeit4.log --sequencing --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --pbwt-modulo 0.000000005 --pbwt-mdr 0.0 --pbwt-mac 0 --ibd2-maf 0.0 --ibd2-mdr 0.0 --ibd2-count 0 --output $workingdir/DARC_683_hg19.recode.phased.vcf.gz

echo -e "
		Basic options:
		  --help                                Produce help message
		  --seed arg (=15052011)                Seed of the random number generator
		  -T [ --thread ] arg (=1)              Number of thread used

		Input files:
		  -I [ --input ] arg                    Genotypes to be phased in VCF/BCF 
				                        format
		  -H [ --reference ] arg                Reference panel of haplotypes in 
				                        VCF/BCF format
		  -S [ --scaffold ] arg                 Scaffold of haplotypes in VCF/BCF 
				                        format
		  -M [ --map ] arg                      Genetic map
		  -R [ --region ] arg                   Target region
		  --use-PS arg                          Informs phasing using PS field from 
				                        read based phasing
		  --sequencing                          Default parameter setting for 
				                        sequencing data (e.g. high variant 
				                        density)

		MCMC parameters:
		  --mcmc-iterations arg (=5b,1p,1b,1p,1b,1p,5m)
				                        Iteration scheme of the MCMC
		  --mcmc-prune arg (=0.999)             Pruning threshold in genotype graphs

		PBWT parameters:
		  --pbwt-modulo arg (=0.02)             Storage frequency of PBWT indexes in cM
				                        (i.e. storage every 0.02 cM by default)
		  --pbwt-depth arg (=4)                 Depth of PBWT indexes to condition on
		  --pbwt-mac arg (=2)                   Minimal Minor Allele Count at which 
				                        PBWT is evaluated
		  --pbwt-mdr arg (=0.5)                 Maximal Missing Data Rate at which PBWT
				                        is evaluated

		IBD2 parameters:
		  --ibd2-length arg (=3)                Minimal size of IBD2 tracks for 
				                        building copying constraints
		  --ibd2-maf arg (=0.01)                Minimal Minor Allele Frequency for 
				                        variants to be considered in the IBD2 
				                        mapping
		  --ibd2-mdr arg (=0.5)                 Maximal Missing data rate for variants 
				                        to be considered in the IBD2 mapping
		  --ibd2-count arg (=150)               Minimal number of filtered variants in 
				                        IBD2 tracks
		  --ibd2-output arg                     Output all IBD2 constraints in the 
				                        specified file (useful for debugging!)

		HMM parameters:
		  -W [ --window ] arg (=2.5)            Minimal size of the phasing window in 
				                        cM
		  --effective-size arg (=15000)         Effective size of the population

		Output files:
		  -O [ --output ] arg                   Phased haplotypes in VCF/BCF format
		  --log arg                             Log file
" > /dev/null

# plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
# plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recode12 --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeAD --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeAD.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeA --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs//PKLR683_249_snps.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.ped

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.ped

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.map | awk '{print $2,$4}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.map

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.map

# java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

# plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased

plink --vcf $workingdir/DARC_683_hg19.recode.phased.vcf.gz \
	--make-bed \
	--out $workingdir/DARC_683_hg19.recode.phased

plink --bfile $workingdir/DARC_683_hg19.recode.phased \
	--recode \
	--out $workingdir/DARC_683_hg19.recode.phased.genotypes

cat $workingdir/DARC_683_hg19.recode.phased.genotypes.ped | cut -d" " -f7- > $workingdir/DARC_683_hg19.recode.phased.genotypes.haplo.stats.input.ped

cat $workingdir/DARC_683_hg19.recode.phased.genotypes.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > $workingdir/DARC_683_hg19.recode.phased.genotypes.haplo.stats.input.map

#################################################################################
#
# Haplotype Analysis using haplo.stats, 81 variants, 683 sickle-cell subjects
# DARC/ACKR1 gene
# redo by excluding variants with MAF=0 (initially called using 688 subjects)
# redo by excluding some multiallelic indels	
# redo by splitting multiallelic
#
#################################################################################

workingdir="/home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants"

vcftools --gzvcf /home/seifuddinft/lab-thein/Analysis/02_merge/merge.688.nomissing.withRSID.vcf.gz --remove /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/exclude_5_individuals_wgs_vcf.txt --positions $workingdir/DARC_ACKR1_81_variants_by_position_hg19.txt --recode --out $workingdir/DARC_683_hg19

bgzip -c $workingdir/DARC_683_hg19.recode.vcf > $workingdir/DARC_683_hg19.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf $workingdir/DARC_683_hg19.recode.vcf.gz

/home/seifuddinft/bin/bcftools norm $workingdir/DARC_683_hg19.recode.vcf.gz --multiallelics - both --output $workingdir/DARC_683_hg19.recode.multisplit.vcf.gz --output-type z

tabix -p vcf $workingdir/DARC_683_hg19.recode.multisplit.vcf.gz

# tar -zxvf /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/genetic_maps.b38.tar.gz
# gunzip /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap.gz
# cat /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.gmap | sed '1,1d' |awk 'BEGIN{print "pos\tchr\tcM"};{print $1"\t""chr"$2"\t"$3}' > /mnt/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/chr1.b38.chr.gmap

# shapeit4 phases SNPs and small INDELs only

shapeit4 --input $workingdir/DARC_683_hg19.recode.multisplit.vcf.gz --map /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/shapeit4/maps/01-genetic_maps_b37_by_chromosome/chr1.b37.gmap --region chr1 --log $workingdir/DARC_683_hg19.recode.phased.shapeit4.log --sequencing --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --pbwt-modulo 0.000000005 --pbwt-mdr 0.0 --pbwt-mac 0 --ibd2-maf 0.0 --ibd2-mdr 0.0 --ibd2-count 0 --output $workingdir/DARC_683_hg19.recode.multisplit.phased.vcf.gz

tabix -p vcf $workingdir/DARC_683_hg19.recode.multisplit.phased.vcf.gz

echo -e "
		Basic options:
		  --help                                Produce help message
		  --seed arg (=15052011)                Seed of the random number generator
		  -T [ --thread ] arg (=1)              Number of thread used

		Input files:
		  -I [ --input ] arg                    Genotypes to be phased in VCF/BCF 
				                        format
		  -H [ --reference ] arg                Reference panel of haplotypes in 
				                        VCF/BCF format
		  -S [ --scaffold ] arg                 Scaffold of haplotypes in VCF/BCF 
				                        format
		  -M [ --map ] arg                      Genetic map
		  -R [ --region ] arg                   Target region
		  --use-PS arg                          Informs phasing using PS field from 
				                        read based phasing
		  --sequencing                          Default parameter setting for 
				                        sequencing data (e.g. high variant 
				                        density)

		MCMC parameters:
		  --mcmc-iterations arg (=5b,1p,1b,1p,1b,1p,5m)
				                        Iteration scheme of the MCMC
		  --mcmc-prune arg (=0.999)             Pruning threshold in genotype graphs

		PBWT parameters:
		  --pbwt-modulo arg (=0.02)             Storage frequency of PBWT indexes in cM
				                        (i.e. storage every 0.02 cM by default)
		  --pbwt-depth arg (=4)                 Depth of PBWT indexes to condition on
		  --pbwt-mac arg (=2)                   Minimal Minor Allele Count at which 
				                        PBWT is evaluated
		  --pbwt-mdr arg (=0.5)                 Maximal Missing Data Rate at which PBWT
				                        is evaluated

		IBD2 parameters:
		  --ibd2-length arg (=3)                Minimal size of IBD2 tracks for 
				                        building copying constraints
		  --ibd2-maf arg (=0.01)                Minimal Minor Allele Frequency for 
				                        variants to be considered in the IBD2 
				                        mapping
		  --ibd2-mdr arg (=0.5)                 Maximal Missing data rate for variants 
				                        to be considered in the IBD2 mapping
		  --ibd2-count arg (=150)               Minimal number of filtered variants in 
				                        IBD2 tracks
		  --ibd2-output arg                     Output all IBD2 constraints in the 
				                        specified file (useful for debugging!)

		HMM parameters:
		  -W [ --window ] arg (=2.5)            Minimal size of the phasing window in 
				                        cM
		  --effective-size arg (=15000)         Effective size of the population

		Output files:
		  -O [ --output ] arg                   Phased haplotypes in VCF/BCF format
		  --log arg                             Log file
" > /dev/null

# plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased.vcf.gz --make-bed --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
# plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased --recode --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.recode.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recode12 --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeAD --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeAD.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeA --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs//PKLR683_249_snps.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.ped

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.ped

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.map | awk '{print $2,$4}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.map

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased.set.missing.zero.map

# java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

# plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased

plink --vcf $workingdir/DARC_683_hg19.recode.multisplit.phased.vcf.gz \
	--make-bed \
	--out $workingdir/DARC_683_hg19.recode.multisplit.phased

plink --bfile $workingdir/DARC_683_hg19.recode.multisplit.phased \
	--recode \
	--out $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes

cat $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes.ped | cut -d" " -f7- > $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes.haplo.stats.input.ped

cat $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes.haplo.stats.input.map

###Haploview analysis

plink --bfile $workingdir/DARC_683_hg19.recode.multisplit.phased \
	--recode 12 \
	--out $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes.recode12

plink2 --bfile $workingdir/DARC_683_hg19.recode.multisplit.phased \
	--recode HV \
	--biallelic-only strict \
	--out $workingdir/DARC_683_hg19.recode.multisplit.phased.genotypes.recode12

plink --vcf $workingdir/DARC_683_hg19.recode.multisplit.phased.vcf.gz \
	--make-bed \
	--biallelic-only strict \
	--snps-only \
	--out $workingdir/DARC_683_hg19.recode.multisplit.biallelic.strict.only.snps.phased

plink2 --bfile $workingdir/DARC_683_hg19.recode.multisplit.biallelic.strict.only.snps.phased \
	--recode HV \
	--out $workingdir/DARC_683_hg19.recode.multisplit.biallelic.strict.only.snps.phased.haploview.input

# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recode12 --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode12.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeAD --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeAD.phased
# plink --bfile /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased --recodeA --out /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recodeA.phased

# awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs//PKLR683_249_snps.recode.phased.ped > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.ped

awk -F' ' '{ if ($6==-9) $6=0;}1' OFS=' ' /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/DARC_683_hg19.recode.multisplit.phased.genotypes.recode12.ped > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/DARC_683_hg19.recode.multisplit.phased.genotypes.recode12.set.missing.pheno.to.zero.haploview.input.ped

# cat /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.map | awk '{print $2,$4}' > /mnt/lab-thein/Fayaz/04-Gabriel/03-PKLR683_249_SNPs/PKLR683_249_snps.recode.phased.set.missing.zero.map

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/DARC_683_hg19.recode.multisplit.phased.genotypes.recode12.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/DARC_683_hg19.recode.multisplit.phased.genotypes.recode12.set.missing.pheno.to.zero.haploview.input.map

java -Xms64g -Xmx256g -jar /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/Haploview.jar

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased --blocks no-pheno-req --out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.phased


#################################################################################
#
# Haplotype Analysis using haplo.stats, 1000 genomes, phase3, 81 variants
# DARC/ACKR1 gene
# 
#
#################################################################################

workingdir="/home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_haplotypes_1000genomesAFR_phase3"

vcftools --gzvcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/12-PKLR_220_variants_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --positions /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_haplotypes_1000genomesAFR_phase3/DARC_ACKR1_81_variants_by_position_hg19.txt --recode --out $workingdir/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_DARC_variants_b37

bgzip -c $workingdir/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_DARC_variants_b37.recode.vcf > $workingdir/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_DARC_variants_b37.recode.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf $workingdir/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_DARC_variants_b37.recode.vcf.gz

perl /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/vcf_to_ped_convert.pl \
	-vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.vcf.gz \
	-sample_panel_file /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel \
	-region 1:155262400-155262414 \
	-population GBR \
	-output_ped ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.ped \
	-output_info ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_69_AA_variants_b37.recode.map \
	-output_dir /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3 \
	-min_maf 0 \
	-base_format 'letters'

plink --vcf /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.vcf.gz \
	--make-bed \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode

###AFR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="AFR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.AFR.haplo.stats.input.map

###EUR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="EUR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel

plink --bfile /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.EUR.panel \
	--out /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.ped | cut -d" " -f7- > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.haplo.stats.input.ped

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_66_AA_variants_b37.recode.EUR.haplo.stats.input.map


#################################################################################
#
# Merging 683WGS, 1kgenomesphase3ALL, 81 variants
# DARC/ACKR1 gene
# 
#
#################################################################################

workingdir="/home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_merged_VCF_683WGS_1kgenomesphase3ALL"

/home/seifuddinft/bin/bcftools merge \
								--missing-to-ref \
								--merge both \
								--output $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.vcf \
								--output-type v  \
								/home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/DARC_683_hg19.recode.multisplit.phased.vcf.gz /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_haplotypes_1000genomesAFR_phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_DARC_variants_b37.recode.vcf.gz

bgzip -c $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.vcf > $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.vcf.gz
# bcftools index /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/PKLR683_249_snps.vcf.gz
tabix -p vcf $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.vcf.gz

/home/seifuddinft/bin/bcftools norm \
								$workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.vcf.gz \
								--multiallelics - both \
								--output $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti.vcf.gz \
								--output-type z

tabix -p vcf $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti.vcf.gz

echo -e "

About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.
         Note that only records from different files can be merged, never from the same file. For
         \"vertical\" merge take a look at \"bcftools norm\" instead.
Usage:   bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]

Options:
        --force-samples                resolve duplicate sample names
        --print-header                 print only the merged header and exit
        --use-header <file>            use the provided header
    -0  --missing-to-ref               assume genotypes at missing sites are 0/0
    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. \"PASS,.\")
    -F, --filter-logic <x|+>           remove filters if some input is PASS (\"x\"), or apply all filters (\"+\") [+]
    -g, --gvcf <-|ref.fa>              merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max
    -i, --info-rules <tag:method,..>   rules for merging INFO fields (method is one of sum,avg,min,max,join) or \"-\" to turn off the default [DP:sum,DP4:sum]
    -l, --file-list <file>             read file names from the file
    -m, --merge <string>               allow multiallelic records for <snps|indels|both|all|none|id>, see man page for details [both]
        --no-version                   do not append version and command line to the header
    -o, --output <file>                write output to a file [standard output]
    -O, --output-type <b|u|z|v>        'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
    -r, --regions <region>             restrict to comma-separated list of regions
    -R, --regions-file <file>          restrict to regions listed in a file
        --threads <int>                use multithreading with <int> worker threads [0]

" > /dev/null

plink --vcf $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti.vcf.gz \
	--make-bed \
	--out $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti


###AFR
cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/07-PKLR_69_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel | awk '{if($3=="AFR"){print $1,$1}}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel

plink --bfile $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/08-PKLR_66_SNPs_haplostats_1000genomes_phase3/integrated_call_samples_v3.20130502.AFR.panel \
	--out $workingdir/merged_VCF_683WGS_1kgenomesphase3AFRonly_DARC_ACKR1.splitmulti

cat $workingdir/merged_VCF_683WGS_1kgenomesphase3AFRonly_DARC_ACKR1.splitmulti.ped | cut -d" " -f7- > $workingdir/merged_VCF_683WGS_1kgenomesphase3AFRonly_DARC_ACKR1.splitmulti.haplo.stats.input.ped

cat $workingdir/merged_VCF_683WGS_1kgenomesphase3AFRonly_DARC_ACKR1.splitmulti.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > $workingdir/merged_VCF_683WGS_1kgenomesphase3AFRonly_DARC_ACKR1.splitmulti.haplo.stats.input.map

###683WGS

cat /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/DARC_683_hg19.recode.multisplit.biallelic.strict.only.snps.phased.fam | awk '{print $1,$1}' > /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/683WGSsubjects.txt

plink --bfile $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti \
	--recode \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/683WGSsubjects.txt \
	--out $workingdir/merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti

cat $workingdir/merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.ped | cut -d" " -f7- > $workingdir/merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.haplo.stats.input.ped

cat $workingdir/merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.map | awk '{if($2=="."){print"chr1_"$4,$4}else{print $2,$4}}' > $workingdir/merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.haplo.stats.input.map


###Sliptstree clustering analysis

workingdir="/home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_merged_VCF_683WGS_1kgenomesphase3ALL"

plink --bfile $workingdir/merged_VCF_683WGS_1kgenomesphase3ALL_DARC_ACKR1.splitmulti \
	--recode A \
	--keep /home/seifuddinft/lab-thein/Fayaz/04-Gabriel/10-DARC_ACKR1_Duffy_Antigen_haplotypes/00-redo-haplotypes-using-filtered-DARC-variants/683WGSsubjects.txt \
	--out $workingdir/merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.recodeA

# cat merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.recodeA.raw | sed '1,1d' | cut -d" " -f7- | sed 's/ //g' > tmp
# cat merged_VCF_683WGS_1kgenomesphase3WGS683only_DARC_ACKR1.splitmulti.recodeA.raw | sed '1,1d' | cut -d" " -f1 > tmp1
# paste tmp1 tmp > tmp2

