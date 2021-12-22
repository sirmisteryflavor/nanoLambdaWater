
# Load required libraries
library(caret)
library(tidyverse)
library(ggplot2)
library(ggcorrplot)

# Read in the two datasets and merge them.
spectra <- read.csv("projects/nanolambda/Averaged_spectra.csv")
labResults <- read.csv("projects/nanolambda/Water_quality_data.csv")
waterQuality <- merge(spectra, labResults, by="Name")


################################################################################
#                                                                              #
#                               Pre-processing                                 #
#                                                                              #
################################################################################


# Rename the spectra column names to contain only the spectra value
newColNames <- c()

for(col in colnames(waterQuality)){
  if(length(grep("X[0-9]+", col)) > 0){
    newCol <- gsub("X", "", col)
    newCol <- tolower(newCol)
    newColNames <- append(newColNames, newCol)
  } else {
    newColNames <- append(newColNames, tolower(col))
  }
}
colnames(waterQuality) <- newColNames

# Also identify and label non-numeric column names for later. 
nonSpectraColNames <- c()
for(col in newColNames){
  if (length(grep("[a-zA-Z]", col)) > 0) {
    nonSpectraColNames <- c(col, nonSpectraColNames)
  }
}

# Convert to long format so that we can plot it easily. 
waterQualityLonger <- 
  pivot_longer(waterQuality, 
               !nonSpectraColNames, 
               names_to = "wavelength", 
               values_to = "absorption") %>%
  mutate(absorption = as.numeric(absorption))

# Convert to tibble
spectraOnly <- waterQuality[, grepl("[0-9]", names(waterQuality))]


################################################################################
#                                                                              #
#                                 Exploratory                                  #
#                                                                              #
################################################################################

# Plot the individual wavelengths from each sample.
waterQualityLonger %>%
  ggplot(aes(x = as.numeric(wavelength), y = absorption, group = name)) + 
  geom_line(aes(alpha = .1)) + 
  xlab("Wavelength") +
  ylab("Absorption") +
  theme_bw() + 
  theme(legend.position = "none") +
  ggtitle("Absorption Measurements From India Water Samples")

# Visualize independent variable correlations.
waterQualityLonger %>%
  select(nonSpectraColNames) %>%
  pivot_longer(!name, 
               names_to = "property", 
               values_to = "measurement") %>%
  ggplot(aes(x = measurement)) + 
  geom_density() + 
  facet_wrap(~ property, scales = "free") + 
  xlab("Measurement") +
  ylab("Density") +
  theme_bw() + 
  ggtitle("Property Measurements From India Water Samples")

# plot the correlation matrix as a plot
ggcorrplot(cor(spectraOnly))
  
# Run a PCA on the spectra data
pca <- prcomp(spectraOnly, scale = TRUE)

# Plot principle component 1 against principle component 2
plot(pca$x[, 1], pca$x[, 2])

# See how much variation each component accounts for
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var) * 100, 1)
barplot(pca.var.per[1:10], main = "Scree Plot", xlab = "Principle Component", ylab = "Percent Variation")

pca.data <- data.frame(pc1 = pca$x[, 1], pc2 = pca$x[, 2], sampleName = waterQuality$name)

ggplot(data = pca.data, aes(x = pc1, y = pc2, label = sampleName)) +
    geom_text() + 
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
    theme_bw() + 
    ggtitle("PCA of Soil Sample")

# Retrieve the loading scores for pc1
loading_scores_pc1 <- pca$rotation[, 1]
loading_scores_pc2 <- pca$rotation[, 2]

# See the wavelength that has the highest influence
spectra_scores_pc1 <- abs(loading_scores_pc1)
spectra_score_ranked_pc1 <- sort(spectra_scores_pc1, decreasing = TRUE)
spectra_scores_pc2 <- abs(loading_scores_pc2)
spectra_score_ranked_pc2 <- sort(spectra_scores_pc2, decreasing = TRUE)

# Top 10 wavelength measurements
top_10_spectras_pc1 <- names(spectra_score_ranked_pc1[1:10])
top_10_spectras_pc2 <- names(spectra_score_ranked_pc2[1:10])

# Combine the two sets and use it as the selected independent variable for
# fitting into a logistic regression.
selected_variables <- union(top_10_spectras_pc1, top_10_spectras_pc2)

################################################################################
#                                                                              #
#                                  Modeling                                    #
#                                                                              #
################################################################################

# Choose a dependent variable we want to model (nitrate) and discretize to 
# turn it into a classification problem.
waterQuality <- 
  waterQuality %>%
  select(all_of(selected_variables), nitrate) %>%
  mutate(safeNitrate = as.factor(ifelse(nitrate <= 10, 1, 0))) %>%
  select(-nitrate)

set.seed(321321)

# Split into Training and Testing
obsIndex <- seq(nrow(waterQuality))
trainingSize <- round(length(obsIndex) * (.75))
trainIndex <- sample(obsIndex, trainingSize)
waterQualityTrain <- waterQuality[trainIndex, ]
waterQualityTest <- waterQuality[-trainIndex, ]

# Fit a simple logistic regression and identify the most predictive spectra
glm.fits <- glm(safeNitrate ~ ., data = waterQualityTrain, family = binomial)
summary(glm.fits)

# Predict probabilities from model
glm.probs <- predict(glm.fits, waterQualityTest, type="response")
glm.pred = rep(1, length(glm.probs))
glm.pred[glm.probs < 0.5] = 0
table(glm.pred, waterQualityTest$safeNitrate)

# Compute the error rate
mean(glm.pred != waterQualityTest$safeNitrate)

################################################################################
#                                                                              #
#                              Cross Validation                                #
#                                                                              #
################################################################################

# Using the caret package
train_control <- trainControl(method = "cv", number = 5, repeats = 3)

set.seed(123)
model_1 <- train(safeNitrate ~ ., 
               data = waterQualityTrain, 
               preProcess=c("scale", "pca"),
               trControl = train_control,
               method = "glm",
               family = binomial())

# Print model summary and results.
summary(model_1)
model_1$results

# Compare to the test accuracy of our model.
model_1_predictions <- predict(model_1, newdata = waterQualityTest)
confusionMatrix(data = model_1_predictions, waterQualityTest$safeNitrate)