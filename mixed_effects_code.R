### Beginning ####

rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lme4, marginaleffects, optimx, broom.mixed)

anes <- readRDS("./data/anes/anes_recoded.rds"); str(anes)

### Random intercept ####

ri.fit <- lmer(ft_dem ~ 1 + age + female + race +
                   (1 | state),
               data = anes)
summary(ri.fit, 
        correlation = FALSE) # correlation of FE, of little interest

pacman::p_load(lmerTest)

# we need to re-fit the model after loading lmerTest
ri.fit <- lmer(ft_dem ~ 1 + age + female + race +
                   (1 | state),
               data = anes)
summary(ri.fit, 
        correlation = FALSE)

# extracting random effects

ranef(ri.fit, 
      drop = TRUE,     # to have them as a vector
      condVar = FALSE) # I don't want to see their conditional variance

# plotting random effects

broom.mixed::tidy(ri.fit, effects = "ran_vals") %>%
    ggplot(aes(x = reorder(level, -estimate), y = estimate)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_point(color = "cyan4") +
    geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                  width = 0, color = "cyan4") +
    labs(x = "State", y = "RE estimate") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_classic() +
    coord_flip()

### Crossed random effects ####

crossed.re.fit <- lmer(ft_dem ~ 1 + age + female + race +
                           (1 | state) + (1 | cohort),
                       data = anes)
summary(crossed.re.fit, 
        correlation = FALSE)

# plotting random effects

bind_rows(
    predictions(crossed.re.fit, 
                by = "state") %>% 
        dplyr::select(state, estimate, conf.low, conf.high) %>% 
        mutate(re = "State") %>% 
        rename(value = state),
    predictions(crossed.re.fit, 
                by = "cohort") %>% 
        dplyr::select(cohort, estimate, conf.low, conf.high) %>% 
        mutate(re = "Cohort") %>% 
        rename(value = cohort)
) %>% 
    ggplot(aes(x = reorder(value, -estimate), y = estimate)) +
    geom_point(aes(col = re)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high, col = re),
                  width = 0) +
    labs(x = "Random effect", y = "Predicted Democratic FT score") +
    scale_color_manual("", values = c("cyan4", "royalblue4")) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_classic() +
    theme(legend.position = "none") +
    facet_wrap(. ~ re,
               scales = "free_y") +
    coord_flip()

### Random slope ####

rs.fit <- lmer(ft_dem ~ 1 + age + female + race +
                   (1 + age | state),
               data = anes)

# scaling age

anes = anes %>%
    mutate(age.s = (age - mean(age)) / sd(age))

# re-estimating model

rs.fit <- lmer(ft_dem ~ 1 + age.s + female + race +
                   (1 + age.s | state),
               data = anes)
summary(rs.fit,
        correlation = FALSE)

# trying different optimizers

optims <- allFit(rs.fit,
                 verbose = FALSE) # don't report progress

optims_OK <- optims[sapply(optims, is, "merMod")]
lapply(optims_OK, function(x) x@optinfo$conv$lme4$messages)
# taken from: https://joshua-nugent.github.io/allFit/

# re-fitting model with optimizer of choice

rs.fit.bob <- lmer(ft_dem ~ 1 + age.s + female + race +
                       (1 + age.s | state),
                   control = lmerControl(optimizer = "bobyqa"),
                   data = anes)
summary(rs.fit.bob, 
        correlation = FALSE)

# spaghetti plot

data.preds = expand.grid(age.s = unique(anes$age.s),
                         state = unique(anes$state),
                         female = 1,
                         race = "hispanic")

data.preds$preds = predict(rs.fit.bob,
                           newdata = data.preds)

ggplot(data.preds, aes(x = age.s, y = preds)) +
    geom_line(aes(color = state)) +
    labs(x = "Age (z-score)", y = "Predicted Democratic FT score") +
    theme_classic() +
    theme(legend.position = "none")

# plotting joint distribution of random effects

data.frame(ranef(rs.fit.bob)) %>% 
    dplyr::select(-c(grpvar, condsd)) %>% 
    pivot_wider(names_from = term, values_from = condval) %>% 
    rename(state = grp,
           intercept = `(Intercept)`) %>% 
    ggplot(aes(x = age.s, y = intercept)) +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_point(color = "cyan4") +
    labs(x = "RE (age)", y = "RE (intercept)") +
    theme_classic()

### Uncorrelated random effects ####

rs.fit.bob.uncorr <- lmer(ft_dem ~ 1 + age.s + female + race +
                              (1 | state) + 
                              (0 + age.s | state), # calling zero to explicitly remove the intercept
                          control = lmerControl(optimizer = "bobyqa"),
                          data = anes)
summary(rs.fit.bob.uncorr, 
        correlation = FALSE)

# plotting joint distribution of random effects

data.preds$preds.uncorr = predict(rs.fit.bob.uncorr,
                                  newdata = data.preds)

ggplot(data.preds, aes(x = age.s, y = preds.uncorr)) +
    geom_line(aes(color = state)) +
    scale_color_manual(values = viridis(5)) +
    labs(x = "Age (z-score)", y = "Predicted Democratic FT score") +
    theme_classic() +
    theme(legend.position = "none")

### Likelihood ratio test ####

# random intercept vs random slope

# re-fit models with ML
ri.fit.ml <- update(ri.fit, REML = FALSE)
rs.fit.bob <- update(rs.fit.bob, REML = FALSE)

anova(ri.fit, rs.fit.bob) # put the simpler model first

# uncorrelated random slope vs correlated random slope
rs.fit.bob <- update(rs.fit.bob, REML = FALSE)
rs.fit.bob.uncorr <- update(rs.fit.bob.uncorr, REML = FALSE)

anova(rs.fit.bob.uncorr, rs.fit.bob) # put the simpler model first

