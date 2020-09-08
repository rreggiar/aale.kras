#!/public/groups/kimlab/.install_bin/anaconda3/envs/aale.analysis.env/bin/Rscript

library(tidyverse)

writeLines(readLines(file("stdin")))

quit()

aale.locus.count <-
	read_csv() %>%
	filter(grepl('Alu', X1)) %>%
	filter_at(contains('.'), any_vars(. >= 10))

write_csv(aale.locus.count, paste()
