---
title: "qPCR_analysis_tutorial_in_R_Tyron"
output: 
  html_document:
    theme: cerulean
    toc: true
    toc_float: true
date: "2023-09-22"
---
### qPCR analysis in R
This R Markdown is for tutorial on how to perform qPCR analysis.

The steps are described below:

* Set working directory and load packages and raw data(excel_file)
* Data wrangling and reshape dataframe.
* Select a gene you want to analyze and create a new dataframe (e.g. df_GeneX).
* Calculation of delta Ct and expression values for qPCR analysis
* Plot a bargraph with ggplot2

### Set working directory and load packages and raw data(excel_file)
```{r}
#First copy everything from your result sheet and paste it into a new excel file!
#Here it is called test.xlsx.

#Set up your working directory
getwd()
setwd('/Users/S157501/Desktop/R_script')

#Load your packages
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggpubr)
df<-read_excel('test.xlsx',skip = 42) #Remove rows above columns and import excel raw data.
df<-df[-(37:41),] #Remove bottom rows that don't belong to dataframe 
colnames(df)<-sub(' ','_',colnames(df))#Remove space in the dataframe 
df$Sample_Name<-str_replace_all(df$Sample_Name,' ','_')#Remove space in the Sample_Name column 
df #check your dataframe

```
### Data wrangling and reshape dataframe.

```{r}
# Select the columns you want to use for analysis

df<-df%>%
  select(Sample_Name,Target_Name,CT)

# Add new columns with replicate #s and group certain columns
new_df <-df %>%
  group_by(Sample_Name, Target_Name) %>%
  mutate(Replicate = row_number())
new_df #check the new dataframe. Note here you add a new column called Relicate!

```
### Select a gene you want to analyze and create a new dataframe (e.g. df_GeneX).

```{r}
# Reshape the dataframe

final_df_raw <- new_df %>%
  pivot_wider(
    names_from = c(Sample_Name, Target_Name),
    values_from = CT
  )
dim(final_df_raw) #Check the structure of the dataframe. 3 rows, and 13 columns!
final_df_raw #Check the reshaped dataframe.Note Replicate column moves to the first one, and the rest of columns are grouped into  Sample_Name+Target_Name!

```

```{r}
df_IRF1<-final_df_raw%>% #Create a new dataframe for IRF1 gene grown in glucose.#Here the treatment is IFNg and untreated.
  select(GLU_UTX_CRFK_IRF1,GLU_IFNG_CRFK_IRF1,
         GLU_UTX_CRFK_ACTIN,GLU_IFNG_CRFK_ACTIN)
colnames(df_IRF1) #check column names of your new IRF1 dataframe!
```

### Calculation of delta Ct and expression values for qPCR analysis

```{r}
###dCt calculation
dCt_IRF1_glu_utx<-df_IRF1$GLU_UTX_CRFK_IRF1-df_IRF1$GLU_UTX_CRFK_ACTIN
dCt_IRF1_glu_ifng<-df_IRF1$GLU_IFNG_CRFK_IRF1-df_IRF1$GLU_IFNG_CRFK_ACTIN
df_IRF1<-cbind(df_IRF1,dCt_IRF1_glu_utx)#add to the IRF1 dataframe
df_IRF1<-cbind(df_IRF1,dCt_IRF1_glu_ifng)#add to the IRF1 dataframe

colnames(df_IRF1)

###ddCt calculation
dCt_IRF1_glu_utx_avg<-mean(dCt_IRF1_glu_utx)
ddCt_IRF1_glu_utx<-dCt_IRF1_glu_utx-dCt_IRF1_glu_utx_avg
ddCt_IRF1_glu_ifng<-dCt_IRF1_glu_ifng-dCt_IRF1_glu_utx_avg
df_IRF1<-cbind(df_IRF1,ddCt_IRF1_glu_utx)#add to the IRF1 dataframe
df_IRF1<-cbind(df_IRF1,ddCt_IRF1_glu_ifng)#add to the IRF1 dataframe

colnames(df_IRF1)

##expression and standard deviation calculation (final step)
expression_IRF1_glu_utx<-2^-ddCt_IRF1_glu_utx
expression_IRF1_glu_ifng<-2^-ddCt_IRF1_glu_ifng
expression_IRF1_glu_utx_sd<-sd(expression_IRF1_glu_utx)#standard deviation calculation
expression_IRF1_glu_ifng_sd<-sd(expression_IRF1_glu_ifng)#standard deviation calculation

df_IRF1<-cbind(df_IRF1,expression_IRF1_glu_utx)#add to the IRF1 dataframe
df_IRF1<-cbind(df_IRF1,expression_IRF1_glu_ifng)#add to the IRF1 dataframe
df_IRF1<-cbind(df_IRF1,expression_IRF1_glu_utx_sd)#add to the IRF1 dataframe
df_IRF1<-cbind(df_IRF1,expression_IRF1_glu_ifng_sd)#add to the IRF1 dataframe

df_IRF1 #check the new finalized dataframe!

#Save the file with all your raw and calculated data. I save it as a CSV or excel file
library(writexl)
write_xlsx(df_IRF1,"qPCR_IRF1_raw_and_processed_data.xlsx")# save the file as excel file!

```


### Plot a bargraph with ggplot2

```{r}
# Prepare a new dataframe for graph

qpcr_df_irf1_glu<-data.frame(expression_IRF1_glu_utx,#change the column names and only selected expression values!
                               expression_IRF1_glu_ifng)
colnames(qpcr_df_irf1_glu)<-c('Untreated','IFNg')
colnames(qpcr_df_irf1_glu)#check colunm names!


qpcr_df_irf1_glu <- data.frame(
  Treatment = c(rep('Untreated', times = 3), rep('IFNg', times = 3)),
  Expression = c(expression_IRF1_glu_utx, expression_IRF1_glu_ifng)
)# Reshape the dataframe in a way where the two colunmns names changed into Treatment and Expression.

qpcr_df_irf1_glu #check the dataframe!

# Calculate mean and SD for each Treatment---- This step is essential for plotting!!!
summary_df <- qpcr_df_irf1_glu %>%
  group_by(Treatment) %>%
  summarise(mean_Expression = mean(Expression), sd_Expression = sd(Expression))

# Create the plot
qpcr_irf1_glu_plot <- ggplot(summary_df, aes(x = Treatment, y = mean_Expression, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_Expression - sd_Expression, ymax = mean_Expression + sd_Expression), width = 0.2, position = position_dodge(0.9)) +
  labs(title="Expression of IRF1 in glucose", x = "Treatment", y = "Expression") +
  scale_fill_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))

qpcr_irf1_glu_plot#view your diagram

ggsave("qpcr_irf1_glu_plot.png")# save the diagram
```

