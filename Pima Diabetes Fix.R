# Library
library(ggplot2) # Membuat diagram batang
library(RColorBrewer) # Memberi warna pada diagram batang
library(dplyr) # Membuat data frame di uji normal multivariat
library(biotools) # Menghitung nilai chi di uji kesamaan matriks kovarian
library(MASS) # Membuat aturan klasifikasi
library(rrcov) # Memuat MCD estimator
library(mvoutlier) # Mendeteksi outlier

# Import Dataset
Diabetes <- read.csv("C:/Users/khans/Documents/diabetes.csv")

# Data Pre-processing
Diabetes$Outcome[Diabetes$Outcome == "1"] <- "2"
Diabetes$Outcome[Diabetes$Outcome == "0"] <- "1" # Mengganti label sesuai dengan ruang sampel

Diabetes <- Diabetes %>%
  filter(
    !if_any(-Pregnancies, ~ . == 0 | is.na(.)) 
  ) # Menghapus missing value pada data

rownames(Diabetes) <- NULL

# Exploratory Data Analysis
colors <- brewer.pal(2, "Set2")

ggplot(Diabetes, aes(x = Outcome, fill = Outcome)) +
  geom_bar(position = "dodge", width = 0.5, color = "black") +
  scale_fill_manual(
    values = brewer.pal(2, "Set2"),
    labels = c("1" = "Negatif Diabetes", "2" = "Positif Diabetes")
  ) +
  labs(
    x = "Kondisi Diabetes", 
    y = "Jumlah Penduduk", 
    fill = "Kondisi",
    title = "Distribusi Kondisi Diabetes Berdasarkan Hasil Pemeriksaan"
  ) +
  geom_text(
    stat = "count", 
    aes(label = ..count..), 
    position = position_dodge(width = 0.5), 
    vjust = -0.5, 
    color = "black", size = 3
  ) +
  theme_minimal() # Diagram batang variabel respon

descriptive_stats <- Diabetes %>%
  group_by(Outcome) %>%                     
  summarise(across(everything(),            
                   list(mean = ~mean(.),    
                        sd = ~sd(.)),       
                   .names = "{col}_{fn}"))

View(descriptive_stats) # Statistika deskriptif variabel penjelas

# Comparing Two Mean Vectors
grup1 <- Diabetes[Diabetes$Outcome == "1", c("Pregnancies", "Glucose","BloodPressure","SkinThickness","Insulin","BMI","DiabetesPedigreeFunction","Age")]
grup2 <- Diabetes[Diabetes$Outcome == "2", c("Pregnancies", "Glucose","BloodPressure","SkinThickness","Insulin","BMI","DiabetesPedigreeFunction","Age")]

mean_grup1 <- colMeans(grup1)
mean_grup2 <- colMeans(grup2)

cov_grup1 <- cov(grup1)
cov_grup2 <- cov(grup2)

n1 <- nrow(grup1)
n2 <- nrow(grup2)
spooled <- ((n1 - 1) * cov_grup1 + (n2 - 1) * cov_grup2) / (n1 + n2 - 2) # Calculate the pooled covariance matrix

mean_diff <- as.matrix(mean_grup1 - mean_grup2)

t2 <- (n1 * n2) / (n1 + n2) * t(mean_diff) %*% solve(spooled) %*% mean_diff
t2 # T^2 hotelling

F_critical <- qf(0.95, 8, 383)

T2_critical <- (8 * (392 - 2)) / (392 - 8 - 1) * F_critical
T2_critical #T table

# Assessing Assumptions
md_Diabetes<-mahalanobis(Diabetes[1:8],colMeans(Diabetes[1:8]),cov(Diabetes[1:8])) # namanya aja mahalanobis, padahal hasilnya jarak kuadrat

chikuadrat_Diabetes <- qchisq((nrow(Diabetes[1:8]) - seq_len(nrow(Diabetes[1:8])) + 0.5) / nrow(Diabetes[1:8]), df = ncol(Diabetes[1:8]))

qqplot(chikuadrat_Diabetes,md_Diabetes,main="Q-Q Plot",ylab="Jarak Kuadrat",xlab="Chi-square")
abline(a=0,b=1) # Plot normal multivariat

normalitas_Diabetes <- data.frame(md_Diabetes)
normalitas_Diabetes$chikuadrat_Diabetes <- chikuadrat_Diabetes
normalitas_Diabetes$results<-ifelse(normalitas_Diabetes$md_Diabetes <= normalitas_Diabetes$chikuadrat_Diabetes, 'True', 
                                    ifelse(normalitas_Diabetes$chikuadrat_Diabetes > normalitas_Diabetes$md_Diabetes, 'No', 'None'))

sum_md_less_than_chi_Diabetes<-length(which(normalitas_Diabetes$results=="True"))
n_Diabetes<-nrow(normalitas_Diabetes)
(sum_md_less_than_chi_Diabetes/n_Diabetes)*100 # Normal multivariat

boxM <- boxM(as.matrix(Diabetes[1:8]), Diabetes$Outcome)
boxM

v<-0.5*ncol(Diabetes[1:8])*(ncol(Diabetes[1:8])+1)*(2-1)
chisq<-qchisq(c(0.05),df=v,lower.tail=FALSE)
chisq # Kesamaan matriks kovarian

# Detecting Outliers
outlier <- dd.plot(Diabetes[, 1:8], quan = 0.5114796)
outlier_detect <- outlier$outliers

Diabetes$outlier <- outlier_detect

Outliers <- subset(Diabetes, outlier == 'TRUE')
No_Outliers <- subset(Diabetes, outlier == 'FALSE')
str(Outliers)
str(No_Outliers)

# Splitting dataset
n <- round(nrow(No_Outliers)*0.705)
set.seed(12345)
sample <- sample(seq_len(nrow(No_Outliers)), size=n)
sample

pre_train_diabetes <- No_Outliers[sample, ]
train_diabetes <- rbind(Outliers, pre_train_diabetes)
test_diabetes <- No_Outliers[-sample, ]

group1_bv <- train_diabetes[train_diabetes$Outcome == "1", c("Pregnancies", "Glucose","BloodPressure","SkinThickness","Insulin","BMI","DiabetesPedigreeFunction","Age", "outlier")]
group2_bv <- train_diabetes[train_diabetes$Outcome == "2", c("Pregnancies", "Glucose","BloodPressure","SkinThickness","Insulin","BMI","DiabetesPedigreeFunction","Age", "outlier")]

table(group1_bv$outlier) # banyak outlier di kelompok negatif diabetes
table(group2_bv$outlier) # banyak outlier di kelompok positif diabetes

# Analisis Diskriminan Kuadratik
group1 <- group1_bv[1:8]
group2 <- group2_bv[1:8]

mean_group1 <- colMeans(group1)
mean_group2 <- colMeans(group2) # Vektor rata-rata sampel

cov_group1 <- cov(group1)
cov_group2 <- cov(group2) # Matriks kovarian sampel

det_Sigma1 <- det(cov_group1)
det_Sigma2 <- det(cov_group2) # Determinan matriks kovarian

inv_Sigma1 <- solve(cov_group1) 
inv_Sigma2 <- solve(cov_group2) # Inverse matriks kovarian

term1 <- 0.5 * log(det_Sigma1 / det_Sigma2)
term2 <- 0.5 * (t(mean_group1) %*% inv_Sigma1 %*% mean_group1 - t(mean_group2) %*% inv_Sigma2 %*% mean_group2)
k <- term1 + term2

compute_score <- function(row) {
  lhs <- -0.5 * t(row) %*% (inv_Sigma1 - inv_Sigma2) %*% row +
    (t(mean_group1) %*% inv_Sigma1 - t(mean_group2) %*% inv_Sigma2) %*% row - 
    k
  return(lhs)
} # Aturan klasifikasi analisis diskriminan kuadratik

data_matrix <- as.matrix(test_diabetes[1:8])
scores <- apply(data_matrix, 1, compute_score)

predicted_classes <- ifelse(scores >= 0, "1", "2")

testing_adk <- data.frame(
  score = as.numeric(scores),
  decision_boundary = 0,
  predicted_class = predicted_classes,
  actual_class = test_diabetes$Outcome
)

head(testing_adk)

# APER
confusion_matrix <- table(Predicted = testing_adk$predicted_class, Actual = testing_adk$actual_class)
print(confusion_matrix) # Tabel klasifikasi

n11 <- confusion_matrix["1", "1"]
n12 <- confusion_matrix["1", "2"]
n21 <- confusion_matrix["2", "1"]
n22 <- confusion_matrix["2", "2"]

n1 <- n11 + n12 # Total in group 1
n2 <- n21 + n22 # Total in group 2

APER <- ((n12 + n21) / (n1 + n2))* 100
print(paste("Apparent Error Rate (APER):", APER,"%")) # Nilai APER

CCR <- (n11 + n22) / (n1 + n2) * 100
print(paste("Apparent Correct Classification Rate (CCR):", CCR,"%")) # Nilai ACCR

# Analisis diskriminan kuadratik robust
mcd1 <- CovMcd(group1, alpha = 0.5227273, nsamp = 500) # alpha = h1
mcd2 <- CovMcd(group2, alpha = 0.5387931, nsamp = 500)

mean_mcd_group1 <- mcd1@center
mean_mcd_group2 <- mcd2@center # Vektor rata-rata sampel robust

cov_mcd_group1 <- mcd1@cov
cov_mcd_group2 <- mcd2@cov # Matriks kovarian sampel robust

det_mcd_Sigma1 <- det(cov_mcd_group1)
det_mcd_Sigma2 <- det(cov_mcd_group2) # Determinan matriks kovarian sampel robust

inv_mcd_Sigma1 <- solve(cov_mcd_group1)
inv_mcd_Sigma2 <- solve(cov_mcd_group2) # Invers matriks kovarian sampel robust

term1_mcd <- 0.5 * log(det_mcd_Sigma1 / det_mcd_Sigma2)
term2_mcd <- 0.5 * (t(mean_mcd_group1) %*% inv_mcd_Sigma1 %*% mean_mcd_group1 - t(mean_mcd_group2) %*% inv_mcd_Sigma2 %*% mean_mcd_group2)
k_mcd <- term1_mcd + term2_mcd

compute_score_mcd <- function(row) {
  lhs <- -0.5 * t(row) %*% (inv_mcd_Sigma1 - inv_mcd_Sigma2) %*% row +
    (t(mean_mcd_group1) %*% inv_mcd_Sigma1 - t(mean_mcd_group2) %*% inv_mcd_Sigma2) %*% row -
    k_mcd
  return(lhs)
} # Aturan klasifikasi diskriminan kuadratik robust

scores_mcd <- apply(data_matrix, 1, compute_score_mcd)

predicted_classes_mcd <- ifelse(scores_mcd >= 0, "1", "2")

result_mcd <- data.frame(
  score = as.numeric(scores_mcd),
  decision_boundary = 0,
  predicted_class = predicted_classes_mcd,
  actual_class = test_diabetes$Outcome
)

View(result_mcd)

# APER
confusion_matrix_mcd <- table(Predicted = result_mcd$predicted_class, Actual = result_mcd$actual_class)
print(confusion_matrix_mcd) # Tabel klasifikasi mcd

n11_mcd <- confusion_matrix_mcd["1", "1"]
n12_mcd <- confusion_matrix_mcd["1", "2"]
n21_mcd <- confusion_matrix_mcd["2", "1"]
n22_mcd <- confusion_matrix_mcd["2", "2"]

n1_mcd <- n11_mcd + n12_mcd
n2_mcd <- n21_mcd + n22_mcd

APER_mcd <- (n12_mcd + n21_mcd) / (n1_mcd + n2_mcd)
print(paste("Apparent Error Rate (APER):", APER_mcd)) # Nilai APER MCD

CCR_mcd <- (n11_mcd + n22_mcd) / (n1_mcd + n2_mcd)
print(paste("Apparent Correct Classification Rate (CCR):", CCR_mcd)) # Nilai ACCR MCD