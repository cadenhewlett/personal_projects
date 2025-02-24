library(ARTool)
library(emmeans)
library(WRS2)
set.seed(1928)
# load data
goggles$index <- 1:nrow(goggles)
indices <- sample(goggles$index, size = nrow(goggles) * 0.75)
train <- goggles[indices, ]
test  <- goggles[-indices, ]

#----- Aligned Rank Transform -----#
model_art <- art(attractiveness ~ gender * alcohol, data = train)
art_nova <- anova(model_art)
# the linear model
art_lm_model <- artlm(model_art, term = "gender:alcohol")

#----- Trimmed Means -----#
# ANOVA for trimmed means, M-estimators, and medians.
t2w_results <- t2way(attractiveness ~ gender * alcohol,
                     data = train,
                     tr = 0.2)
# post-hoc
coef_info <- mcp2a(attractiveness ~ gender * alcohol, data = train)
#----- Classical ANOVA -----#
lm_train <- train
lm_train$gender <- factor(lm_train$gender)
lm_train$alcohol <- factor(lm_train$alcohol)
lm_train$gender_alcohol <- interaction(lm_train$gender, lm_train$alcohol)
lm_train <- na.omit(lm_train)
lm_results <-  lm(attractiveness  ~ gender * alcohol, data = lm_train)
lm_anova <- anova(lm_results)


# --- Comparison of p-Values --- #
names <- rownames(lm_anova)[-4]
# extract
lm_p <- ( data.frame(lm_anova)[1:3, 5] )
tm_p <- c(t2w_results$A.p.value, t2w_results$B.p.value, t2w_results$AB.p.value)
ar_p <- data.frame( art_nova )$`Pr..F.`

all_results <- data.frame(rbind(lm_p, tm_p, ar_p))
colnames(all_results) <- names

# --- Prediction --- #
rbind(
  lm_results$coefficients,
  art_lm_model$coefficients
)

