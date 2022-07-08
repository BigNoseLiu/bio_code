args=commandArgs(T)
xlim_c <- args[1]
ylim_c <- args[2]
output_file_name <- args[3]
input_file_name <- args[4]
png(file = output_file_name)
data <- read.table(file = input_file_name)
data2 <- as.matrix(data)
density_data <- density(rep(data2[,1],data2[,2]))
density_data
plot(density_data,xlim=c(0,  as.numeric(xlim_c)))
