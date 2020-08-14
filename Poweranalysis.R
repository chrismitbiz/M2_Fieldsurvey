#power analysis


metadata<- readr::read_tsv("~/Documents/Study/LaTrobe/Research/phD/Milestone2_Field_Survey/M2_Bacteria/M2_Qiime_output_Bacteria2020run/metadata.tsv") 

metadata2 <- metadata[c(2:nrow(metadata)),c(1:length(metadata))] %>%
  tibble::rownames_to_column("spl") %>%
  dplyr::mutate_all(type.convert) %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  tibble::column_to_rownames("spl")

data <- metadata2 %>% group_by(Soil, Paddock) %>% summarise(Dlossmean = mean(Dloss), CNmean = mean(CN))
shapiro.test(data$Dlossmean)
rcompanion::groupwiseMean(Dlossmean ~ 1, data = data, conf = 0.95, digits = 3,
              boot = TRUE, R = 1000)

rcompanion::groupwiseMean(D1517 ~ 1, data = metadata2, conf = 0.95, digits = 3,
                          boot = TRUE, R = 1000)

m <- lm(Dlossmean ~ CNmean, data = data)
summary(m)

# Calculate the number of samples needed to model Dlossmean (12 samples) with CN as predictor
# https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html
# The effect size is the R2 of the model. the coefficient of determination, aka the “proportion of variance explained”
pwr.f2.test(u = 1, f2 = 0.5/(1 - 0.5), sig.level = 0.001, power = 0.8)
#n = v + u + 1
# n = 21 + 1 + 1 = 22 paddocks needed at an effec size of 0.5


# samples (n) needed for a rubust mean
#http://r-video-tutorial.blogspot.com/2017/07/power-analysis-and-sample-size.html
# n = ((2 x STD) / SE )^2
data <- metadata2 %>% group_by(Paddock) %>% summarise(Dlossmean = mean(Dloss), CNmean = mean(CN))
sum_stat <- data %>% rstatix::get_summary_stats(Dlossmean, show = c("mean", "sd", "se"))
n = ( (2 * sum_stat$sd) / sum_stat$se )^2
n
# n = 48 samples needed for the mean value to be robust mean of dieldrin concentrations


# Sample number required if we want to test the effect of different paddocks
pwr.anova.test(k=12, f=0.25, sig.level=0.05, power=0.8)  
# at an medium effect size of 0.25 we need 28 samples per paddock to have an 80% power. 


# Sample number required if we want to correlate dieldrin where the correlation is predicted to be 0.5
cohen.ES(test = "r", size = "medium")
pwr.r.test(r = 0.5, sig.level = 0.05, power = 0.8, alternative = "two.sided")
# We need 28 samples in that case




