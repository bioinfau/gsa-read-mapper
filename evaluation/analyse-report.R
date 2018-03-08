
check_package <- function(x) {
    if (!require(x, character.only = TRUE))
    {
    	install.packages(x, dep=TRUE)
    	if(!require(x,character.only = TRUE))
    		stop("Package not found")
    }
 }

check_package("ggplot2")
check_package("dplyr")

args <- commandArgs(TRUE)
if (length(args) < 1) {
	stop("The evaluation report needs to be provided as the first argument.")
}
time_results <- read.table(args[1], header = TRUE)

ordered <- time_results %>%
    group_by(mapper) %>%
    summarise(mean_time = mean(time)) %>%
    arrange(mean_time)
ordered_mappers <- factor(time_results$mapper, levels = ordered$mapper)
fastest <- ordered %>% top_n(n = -1, wt = mean_time) %>% as.list

normalised <- time_results %>%
    mutate(time = (time / fastest$mean),
           mapper = factor(mapper, levels = ordered$mapper))
ordered_normalised <- normalised %>%
    group_by(mapper) %>%
    summarise(mean_time = mean(time))

ggplot(normalised, aes(x = mapper, y = time)) +
    geom_hline(
        aes(yintercept = mean_time),
        color = "lightgray", linetype="dashed",
        data = ordered_normalised
    ) +
    geom_boxplot(fill = "#E7B800", alpha = 0.7) +
    geom_jitter(width = 0.1, colour = "#fc6721", size = 2) +
    scale_y_continuous(
        labels = paste(fastest$mapper, "x", round(ordered_normalised$mean_time)),
        breaks = ordered_normalised$mean_time
    ) +
    theme_classic() +
    xlab("Read-mapper") +
    ylab("Normalised running time")
ggsave(paste0(args[1], ".png"))
