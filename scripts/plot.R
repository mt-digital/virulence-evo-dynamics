require(ggplot2)
require(latex2exp)
require(dplyr)
require(readr)
require(purrr)
require(reshape2)
require(stringr)

mytheme <- theme(axis.line = element_line(), legend.key=element_rect(fill = NA),
                text = element_text(size=22), 
                panel.background = element_rect(fill = "white"))


plot_series <- function(csv_file, write_path="tmp/output.pdf") {

  df <- read.csv(csv_file)

  # dflim <- df[,!names(df) %in% c("susceptible")]
  # melted <- melt(dflim, id.vars = "step")
  melted <- melt(df, id.vars = "step")

  ggplot(melted, aes(x=step, y=value, col=variable)) + geom_line() + mytheme

  ggsave(write_path, width=6.75, height=5)
}


total_infections <- function(
    csv_files = c("50_trials_3_mrates_1.csv", "50_trials_3_mrates_2.csv",
                  "50_trials_3_mrates_vinit=1.0:0.1:2.0_1.csv", 
                  "50_trials_3_mrates_vinit=1.0:0.1:2.0_2.csv"),
    write_file = "figures/total_infections.pdf"
    ){
    
  df <- read_csv_series(csv_files) 
  df$mutation_rate <- factor(df$mutation_rate) 

  p <- ggplot(df, aes(x=virulence_init, y=total_infected, color=mutation_rate)) + 
    stat_summary(geom = "line", size=1.25, fun = mean, aes(group=mutation_rate)) +
    stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3, aes(fill=mutation_rate, group=mutation_rate)) +
    xlab(TeX("Initial virulence, $v_0$")) + ylab("Total infected\n(mean with 95% CI)") +
    scale_color_discrete(TeX("Mutation rate, $\\mu$")) + scale_fill_discrete(TeX("Mutation rate, $\\mu$")) +
    mytheme
  
  ggsave(write_file, p, width=8.25, height = 4.5)
  
  return (p)
}

rescue_factor <- function(
    csv_files = c("50_trials_3_mrates_1.csv", "50_trials_3_mrates_2.csv",
                  "50_trials_3_mrates_vinit=1.0:0.1:2.0_1.csv", 
                  "50_trials_3_mrates_vinit=1.0:0.1:2.0_2.csv"),
    write_file = "figures/evolutionary_rescue.pdf"
) {
  df <- read_csv_series(csv_files) 
  df$virulence_init <- factor(df$virulence_init) 
  df$mutation_rate <- factor(df$mutation_rate) 
  
  aggdf <- group_by(df, virulence_init, mutation_rate) %>%
    summarize(mean_total_infected = mean(total_infected))
  
  # Separate different virulence inits for calculating rescue factor, log of diff
  # between no evolution and that particular evo rate setting. The following is
  # not a smart way to do it, but good enough for a first pass, optimize later.
  aggdf_mu0 <- filter(aggdf, mutation_rate == 0.0)
  
  # Calculate log of difference for virulence 0.4. Automate later if neccessary.
  aggdf_mu4 <- filter(aggdf, mutation_rate == 0.4)
  head(aggdf_mu4)
  
  # Calculate log of difference for virulence 0.8. Automate later if neccessary.
  aggdf_mu8 <- filter(aggdf, mutation_rate == 0.8)
  head(aggdf_mu8)
  
  # print(aggdf_mu8$mean_total_infected - aggdf_mu0$mean_total_infected)
  
  aggdf_mu4$log_diff_total <- (aggdf_mu4$mean_total_infected - aggdf_mu0$mean_total_infected) / aggdf_mu0$mean_total_infected
  print(aggdf_mu4$log_diff_total)
  aggdf_mu8$log_diff_total <- (aggdf_mu8$mean_total_infected - aggdf_mu0$mean_total_infected) / aggdf_mu0$mean_total_infected
  
  logdiffdf <- rbind(aggdf_mu4, aggdf_mu8)
  
  p <- ggplot(logdiffdf, aes(x=virulence_init, y=log_diff_total, color=mutation_rate)) +
    geom_line(aes(group=mutation_rate), size=1.25) +
    xlab(TeX("Initial virulence, $v_0$")) + ylab(TeX("Rescue factor, $\\rho$")) + 
    scale_color_manual(TeX("Mutation rate, $\\mu$"), breaks=c(0.4, 0.8), 
                       values=c("#00b938", "#609cff")) +
    mytheme
  
  ggsave(write_file, p, width=8.25, height=4.5)
  
  return (p)
}

read_csv_series <- function(csv_files) {
  return (
    csv_files %>%
      map_df(~read_csv(., show_col_types = FALSE))
  )
}
