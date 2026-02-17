## ============================================================
## LabPhon 2026: PB52 /ɪ/ vs /ɛ/ comparison
## ============================================================

set.seed(2026)

# Core packages
library(phonJSD)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(phonTools)
library(MASS)   # LDA (optional)

# If running from source, you can use:
# devtools::load_all()

## ============================================================
## 1. Load Petersen & Barney (1952) data
## ============================================================

data(pb52, package = "phonTools")

eh_ih <- pb52 %>%
  filter(vowel %in% c("I", "E")) %>%
  droplevels()

features     <- c("f1", "f2")
category_col <- "vowel"

## ============================================================
## 2. Global metrics (pooled)
## ============================================================

global_jsd <- estimate_jsd(
  data         = eh_ih,
  features     = features,
  category_col = category_col,
  do_boot      = TRUE,
  n_boot       = 500,
  min_tokens   = 20
)

global_pillai <- estimate_pillai(
  data         = eh_ih,
  features     = features,
  category_col = category_col
)

global_bhatt <- estimate_bhatt(
  data         = eh_ih,
  features     = features,
  category_col = category_col
)

global_overlap <- estimate_overlap(
  data         = eh_ih,
  features     = features,
  category_col = category_col
)

global_jsd
global_pillai
global_bhatt
global_overlap

## ============================================================
## 3. Grouped metrics (sex, type)
## ============================================================

compare_by_group <- function(df,
                             group_col,
                             features,
                             category_col,
                             n_boot     = 200,
                             min_tokens = 20) {
  jsd <- estimate_jsd(
    data         = df,
    features     = features,
    category_col = category_col,
    group_col    = group_col,
    do_boot      = TRUE,
    n_boot       = n_boot,
    min_tokens   = min_tokens
  )

  pillai <- estimate_pillai(
    data         = df,
    features     = features,
    category_col = category_col,
    group_col    = group_col,
    min_tokens   = min_tokens
  )

  bhatt <- estimate_bhatt(
    data         = df,
    features     = features,
    category_col = category_col,
    group_col    = group_col,
    min_tokens   = min_tokens
  )

  overlap <- estimate_overlap(
    data         = df,
    features     = features,
    category_col = category_col,
    group_col    = group_col,
    min_tokens   = min_tokens
  )

  jsd %>%
    left_join(pillai,  by = c("group", "n_tokens")) %>%
    left_join(bhatt,   by = c("group", "n_tokens")) %>%
    left_join(overlap, by = c("group", "n_tokens"))
}

by_sex <- compare_by_group(
  df           = eh_ih,
  group_col    = "sex",
  features     = features,
  category_col = category_col
)

by_type <- compare_by_group(
  df           = eh_ih,
  group_col    = "type",
  features     = features,
  category_col = category_col
)

by_sex
by_type

## ============================================================
## 4. Visualization (F1/F2)
## ============================================================

plot_df <- eh_ih %>%
  mutate(
    vowel_ipa = case_when(
      vowel == "E" ~ "/ɛ/",
      vowel == "I" ~ "/ɪ/",
      TRUE ~ as.character(vowel)
    )
  )

ggplot(plot_df, aes(x = f2, y = f1, color = vowel_ipa)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(level = 0.68) +
  scale_x_reverse() +
  scale_y_reverse() +
  labs(
    x = "F2 (Hz)",
    y = "F1 (Hz)",
    color = "Vowel",
    title = "PB52 /ɛ/ vs /ɪ/ distributions"
  ) +
  theme_minimal()

## ============================================================
## 5. Optional: LDA accuracy by group (perceptual proxy)
## ============================================================

lda_accuracy <- function(df,
                         features,
                         category_col = "vowel",
                         test_prop    = 0.3) {
  df <- df %>% mutate(id_row = dplyr::row_number())
  n  <- nrow(df)
  n_test <- max(2, floor(test_prop * n))

  test_ids  <- sample(df$id_row, n_test)
  train_ids <- setdiff(df$id_row, test_ids)

  train <- df %>% filter(id_row %in% train_ids)
  test  <- df %>% filter(id_row %in% test_ids)

  form <- as.formula(
    paste(category_col, "~", paste(features, collapse = " + "))
  )

  lda_fit <- tryCatch(
    MASS::lda(form, data = train),
    error = function(e) NULL
  )
  if (is.null(lda_fit)) return(NA_real_)

  pred <- tryCatch(
    predict(lda_fit, newdata = test),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NA_real_)

  mean(pred$class == test[[category_col]])
}

lda_by_sex <- eh_ih %>%
  group_by(sex) %>%
  group_modify(~ tibble(lda_acc = lda_accuracy(.x, features, category_col))) %>%
  ungroup()

lda_by_type <- eh_ih %>%
  group_by(type) %>%
  group_modify(~ tibble(lda_acc = lda_accuracy(.x, features, category_col))) %>%
  ungroup()

lda_by_sex
lda_by_type
