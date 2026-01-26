##Code for plotting StairwayPlot output

# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _MBS/MBS_stairway/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p1<- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("Magdalena (shallow)")+
  xlim(c(0, 60000))+
  ylab("log (Ne)")

p1
# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _BA_done/BA_stairway/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p2 <- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("Asuncion (shallow)")+
  xlim(c(0, 60000))+
  ylab("log (Ne)")


p2

# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _SQD/SQD_stairway/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p3 <- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("San Quintin (deep)")+
  xlim(c(0, 60000))

p3

# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _SQS/SQS_stairway/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p4<- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("San Quintin (shallow)")+
  xlim(c(0, 60000))

p4

# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _EDreal/stairway_plot_es/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p5<- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("Ensenada (deep)")+
  xlim(c(0, 60000))
  

p5


# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _ES/ED_stairway/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p6<- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("Ensenada (shallow)")+
  xlim(c(0, 60000))
  
  
p6

# Read the file (adjust the file name/path and separator as needed)
sp_data <- read.table("./SFS_analysis3/stairway_plot_v2.1.2 _DS_done/ES_stairway/sample_population Stairway Plot.final.summary", header = TRUE, sep = "\t")

# Convert column names to valid R identifiers
names(sp_data) <- make.names(names(sp_data))
head(sp_data)  # inspect the data

# Create the ggplot object
p7<- ggplot(sp_data, aes(x = year, y = log10(Ne_median))) +
  geom_line(color = "blue") +
  ylim(2.5, 5.5)+
  theme_classic()+
  geom_ribbon(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), alpha = 0.2)+
  ggtitle("Danver's (shallow)")+
  xlim(c(0, 60000))+
  ylab("log (Ne)")

p7

cowplot::plot_grid(p1,p2,p3,p4, p5, p6, p7, p_present, align = "h")

##

# Define a helper function that reads a file, extracts the present estimate,
# and adds a column for the population name.
get_present_est <- function(file, pop_name) {
  sp_data <- read.table(file, header = TRUE, sep = "\t")
  names(sp_data) <- make.names(names(sp_data))
  # Ensure the 'year' column is numeric.
  sp_data$year <- as.numeric(as.character(sp_data$year))
  # Extract the row corresponding to the present (year == 0). If not found, take the smallest year.
  present <- sp_data[sp_data$year == 0, ]
  if(nrow(present) == 0) {
    present <- sp_data[which.min(sp_data$year), ]
  }
  present$pop <- pop_name
  return(present)
}

# Read each file and extract the present estimate with the desired population names.
pop_data <- bind_rows(
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _MBS/MBS_stairway/sample_population Stairway Plot.final.summary", "Magdalena"),
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _BA_done/BA_stairway/sample_population Stairway Plot.final.summary", "Asuncion"),
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _SQD/SQD_stairway/sample_population Stairway Plot.final.summary", "San Quintin deep"),
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _SQS/SQS_stairway/sample_population Stairway Plot.final.summary", "San Quintin shallow"),
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _EDreal/stairway_plot_es/sample_population Stairway Plot.final.summary", "Ensenada deep"),
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _ES/ED_stairway/sample_population Stairway Plot.final.summary", "Ensenada shallow"),
  get_present_est("./SFS_analysis3/stairway_plot_v2.1.2 _DS_done/ES_stairway/sample_population Stairway Plot.final.summary", "Danvers")
)

# Set the factor levels for the population names in the desired order.
pop_data$pop <- factor(pop_data$pop, levels = c("Magdalena", "Asuncion", "San Quintin deep", "San Quintin shallow", "Ensenada deep", "Ensenada shallow", "Danvers"))

# Create a plot that compares the present Ne estimates across populations.
p_present <- ggplot(pop_data, aes(x = pop, y = log10(Ne_median))) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = log10(Ne_2.5.), ymax = log10(Ne_97.5.)), width = 0.2, color = "blue") +
  theme_classic() +
  xlab("Population") +
  ylab("Effective Population Size (Ne)") +
  ggtitle("Effective Population Size (Ne) at Present") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(p_present)


cowplot::plot_grid(p1,p2,p3,p4, p5, p6, p7, p_present)

