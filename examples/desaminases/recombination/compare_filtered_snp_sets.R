# Test to check if there is difference in mutation number between different days

data <- read.csv("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_with_RUN7/per_sample_vcf_number_add_info.stat", sep = "\t", header = TRUE)

#PmCDA1
PmCDA1 <- subset(data, deaminase=='PmCDA1' & genotype=='wt')

kruskal.test(number_of_variants~day, PmCDA1)
#Kruskal-Wallis chi-squared = 5.1671, df = 2, p-value = 0.07551
wilcox.test(number_of_variants~day, PmCDA1, day == '1' | day == '3')
#W = 30, p-value = 0.6165
wilcox.test(number_of_variants~day, PmCDA1, day == '3' | day == '6')
#W = 85, p-value = 0.02782
wilcox.test(number_of_variants~day, PmCDA1, day == '1' | day == '6')
#W = 39, p-value = 0.181

anova_results <- aov(number_of_variants~day, data=PmCDA1) 
summary(anova_results)
summary(lm(anova_results))

#PmCDA1_no_suspicious samples
PmCDA1_no_suspicious_samples <- subset(data, deaminase=='PmCDA1' & genotype=='wt' & number_of_variants > 108 )
kruskal.test(number_of_variants~day, PmCDA1_no_suspicious_samples)
#Kruskal-Wallis chi-squared = 1.485, df = 2, p-value = 0.4759
wilcox.test(number_of_variants~day, PmCDA1_no_suspicious_samples, day == '1' | day == '3')
#W = 30, p-value = 0.6165
wilcox.test(number_of_variants~day, PmCDA1_no_suspicious_samples, day == '3' | day == '6')
#W = 49, p-value = 0.2496
wilcox.test(number_of_variants~day, PmCDA1_no_suspicious_samples, day == '1' | day == '6')
#W = 21, p-value = 0.6991

#PmCDA1_sub1
PmCDA1_sub1 <- subset(data, deaminase=='PmCDA1' & genotype=='sub1')

wilcox.test(number_of_variants~day, PmCDA1_sub1)
#W = 8, p-value = 0.2


#AID
AID <- subset(data, deaminase=='AID' & genotype=='wt')

kruskal.test(number_of_variants~day, AID)
#Kruskal-Wallis chi-squared = 1.8094, df = 2, p-value = 0.4047
wilcox.test(number_of_variants~day, AID, day == '1' | day == '3')
#W = 4, p-value = 0.8571
wilcox.test(number_of_variants~day, AID, day == '3' | day == '6')
#W = 23, p-value = 0.1699  WARNING
wilcox.test(number_of_variants~day, AID, day == '1' | day == '6')
#W = 6, p-value = 1  WARNING

#A1
A1 <- subset(data, deaminase=='A1' & genotype=='wt')

kruskal.test(number_of_variants~day, A1)
#Kruskal-Wallis chi-squared = 1.7064, df = 2, p-value = 0.426
wilcox.test(number_of_variants~day, A1, day == '1' | day == '3')
#W = 19, p-value = 0.2635
wilcox.test(number_of_variants~day, A1, day == '3' | day == '6')
#W = 43.5, p-value = 0.4344 WARNING
wilcox.test(number_of_variants~day, A1, day == '1' | day == '6')
#W = 17, p-value = 0.6166 WARNING


