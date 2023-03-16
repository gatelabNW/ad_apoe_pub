# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
# 
# Date: 10-06-2022
# Written by: Abhi Ramakrishnan
# Summary: Custom violin plotting function using ggplot2
# Use: Change custom parameters as needed below, source this script and add custom 
# parameters to base ggplot function
# E.g.: ggplot(x,y,fill) + geom_violin() + scale_fill_manual() + custom_violin
# 
#-------------------------------------------------------------------------------

# Make a list of parameters that can be added to ggplot
custom_violin_upright <-list(stat_summary(fun=median, geom="point", size=2, color="black"),
                             scale_x_discrete(limits=rev),
                             guides(colour = guide_legend(override.aes = list(size = 3, linetype=0, shape= 16))),
                             theme(plot.title = element_text(hjust = 0.5)))
                             
custom_violin <-list(coord_flip(),
                             stat_summary(fun=median, geom="point", size=2, color="black"),
                             scale_x_discrete(limits=rev),
                             guides(colour = guide_legend(override.aes = list(size = 3, linetype=0, shape= 16))),
                             theme(plot.title = element_text(hjust = 0.5)))            
#-------------------------------------------------------------------------------  
# Extras
#-------------------------------------------------------------------------------  

custom_pal <- colorRampPalette(colors=c("darkorchid1","bisque", "red"))
my_palette <- custom_pal(50)          

diagnosis_order <- c("Healthy Control",
                     "Alzheimers Disease")
group_order <- c("Healthy Control_E3/E3",
                 "Healthy Control_E3/E4",
                 "Healthy Control_E4/E4",
                 "Alzheimers Disease_E3/E3",
                 "Alzheimers Disease_E3/E4",
                 "Alzheimers Disease_E4/E4")
