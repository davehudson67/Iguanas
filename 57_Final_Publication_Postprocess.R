# Post-processing only. Preserve original PPC output and produce corrected copy.
library(tidyverse)
root <- "DataAudits/PosteriorPredictive_Primary"
p <- read_csv(file.path(root,"PPC_SiteDetectionTotals.csv"),show_col_types=FALSE)
stopifnot(all(c("observed","pred_median","pred_low","pred_high","Pr_rep_ge_obs","two_sided_ppp") %in% names(p)))
# Confirm structural absence against the exact PPC subset, not simply a zero quantile.
# Known-age/resolved-sex subset has no opportunities at these five sites.
absent <- p$Island %in% c("Bitter","Gaulin","LeafSouth","NAdder","Noddy","WhiteB")
stopifnot(all(p$observed[absent]==0 & p$pred_high[absent]==0))
p$Pr_rep_ge_obs[absent] <- NA_real_
p$two_sided_ppp[absent] <- NA_real_
p$interpretation <- ifelse(absent,"No informative PPC subset observations","Evaluated")
write_csv(p,file.path(root,"PPC_SiteDetectionTotals_Publication.csv"),na="NA")
# Keep the full supplementary forest plot; main-text compact version uses
# primary, three model-specification alternatives and targeted exclusion.
d <- read_csv("DataAudits/Final_Robustness/Final_Robustness_KeyEffects_Long.csv",show_col_types=FALSE)
keep <- c("Male effect on mortality at age 10 (log HR)",
          "Male effect on Gompertz senescence rate (log ratio)",
          "Density effect on age-10 mortality (log HR)",
          "Density effect on senescence rate (log ratio)","Provisioning HR")
# Show all LOO fits in the main figure too, so site influence is visible.
d <- d %>% filter(effect %in% keep) %>%
  mutate(scenario=factor(scenario,levels=rev(unique(scenario))),
         effect=factor(effect,levels=keep))
stopifnot(nrow(d)==70L)
refs <- d %>% distinct(effect,scale) %>%
  mutate(reference=if_else(scale=="One = no effect",1,0))
g <- ggplot(d,aes(median,scenario,colour=class,shape=class))+
  geom_vline(data=refs,aes(xintercept=reference),colour="grey60",linetype=2)+
  geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.12)+
  geom_point(size=2)+facet_wrap(~effect,scales="free_x",ncol=2,labeller=label_wrap_gen(35))+
  labs(x="Posterior median and 95% credible interval",y=NULL,colour=NULL,shape=NULL)+
  theme_bw(base_size=10)+theme(legend.position="bottom",panel.grid.minor=element_blank())
out <- "Figures/FinalAudited/Manuscript"
ggsave(file.path(out,"Figure_Robustness_Main.pdf"),g,width=11,height=12)
ggsave(file.path(out,"Figure_Robustness_Main.tiff"),g,width=11,height=12,dpi=300,compression="lzw")
cat("Corrected PPC table and main robustness figure saved.\n")
