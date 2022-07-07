require(tidyverse)
theme_set(theme_light())

#### Sample size ####

# expected conversion rate in arm 1
exp1 = 0.2

# expected conversion rate in arm 2
exp2 = 0.03

# size of each arm
max_ss_tot <- 300
ss_step <- 10
n_seq <- seq(25, max_ss_tot/2, ss_step)

# number of iterations
niter <- 1000

# alpha error
alpha <- 0.05

#char prev
char_prev   <- 0.33
char_effect <- 2
char_tested <- 50

# sa
effect_sa  <- 3
subgroup_n <- 4

# dropouts
dropouts <- 0.2


bug_fix_low_threshold <- 0


set.seed(1234)





power_db <- map_df(n_seq, function(n) {
  
  res <- data.frame(iteration = 1:niter, p_pa = NA, sa_p = NA, mc_p = NA, sa_p_eter = NA)
  
  # iterate
  for (i in 1:niter) {
    
    # conversion simulation
    conv1  <- rbinom(n, 1, exp1)
    conv2  <- rbinom(n, 1, exp2)
    sim_db <- data.frame(id = 1:n*2, 
                         group = c(rep("exp1", n), rep("exp2", n)), 
                         conv = c(conv1, conv2))
    converters    <- sim_db$conv == 1
    res[i,]$p_pa  <- chisq.test(sim_db$group, sim_db$conv)$p.value
    
    # multiple comparisons
    sim_db$char1               <- NA
    sim_db[converters,]$char1  <- rbinom(sum(converters), 1, (char_prev*char_effect) )
    sim_db[!converters,]$char1 <- rbinom(sum(!converters), 1, (char_prev/char_effect) )
    with_char                  <- sim_db$char1 == 1
    mdl <- glm(conv ~ char1 + group, data = sim_db, family = "binomial")
    res[i,]$mc_p <- summary(mdl)$coefficients[,'Pr(>|z|)']['char1']
    # res[i,]$mc_p               <- ifelse(n < bug_fix_low_threshold, NA, chisq.test(sim_db$char1, sim_db$conv)$p.value)
    
    # subgroup analysis
    sim_db$con_type              <- NA
    effect_single_sa             <- (100/subgroup_n)*effect_sa
    effect_others_sa             <- (100-effect_single_sa)/(subgroup_n-1)
    
    sim_db[converters & with_char,]$con_type <- sample(1:subgroup_n,
      sum(converters & with_char), 
      replace = T, 
      prob = c(effect_single_sa, rep(effect_others_sa, (subgroup_n-1))))
    
    sim_db[converters & !with_char,]$con_type <- sample(1:subgroup_n,
      sum(converters & !with_char), 
      replace = T)
    
    # res[i,]$sa_p_eter <- ifelse(n < bug_fix_low_threshold, NA, 
    #                        chisq.test(table(sim_db[converters,]$char1, 
    #                                         sim_db[converters,]$con_type))$p.value)
    # 
    # res[i,]$sa_p <- ifelse(n < bug_fix_low_threshold, NA, 
    #                        chisq.test(table(sim_db$char1, 
    #                                         case_when(sim_db$con_type == 1 ~ 1, 
    #                                                   sim_db$con_type == 2 ~ 2,
    #                                                   T ~ NA_real_)))$p.value)
  }
  
  power_pa        <- mean(res$p_pa < alpha)
  power_sa_p_eter <- mean(res$sa_p_eter < alpha)
  power_sa        <- mean(res$sa_p < alpha/(subgroup_n-1))
  power_mc        <- mean(res$mc_p < (alpha/char_tested))
  
  smry <- data.frame(n               = 2*n*(1+dropouts), 
                     power_pa        = power_pa, 
                     power_sa_p_eter = power_sa_p_eter, 
                     power_sa        = power_sa, 
                     power_mc        = power_mc)
  
})


power_db_vert <- power_db %>% 
  gather(type, power, -n) %>% 
  filter(!is.na(power)) %>% 
  mutate(type = case_when(type == "power_pa" ~ "\n\nBetween group conversion comparison\n\n",
                          type == "power_mc" ~ "Factors association to conversion\nindependently from group\n(Bonferroni adjusted)"))

  


ggplot(power_db_vert, aes(n, power, col = type)) +
  geom_point(alpha = 0.5) +
  geom_line() +
  geom_vline(xintercept = 198, lty = 3) +
  geom_hline(yintercept = c(0.80, 0.85, 0.90, 0.95), lty = 2, alpha = c(.4, .7, .9, 1)) +
  scale_x_continuous(name = "\nTotal n. of patients\n",
                     breaks = seq(0, max_ss_tot, 100)) +
  scale_y_continuous(name = "Power\n",
                     breaks = seq(0, 1, .1), 
                     labels = function(x) { scales::percent(x, accuracy = 1) },
                     limits = c(0, 1)) +
  ggsci::scale_color_aaas(name = "Analysis")

# ggsave("power_curve.png", width = 10, height = 5)
  
  








