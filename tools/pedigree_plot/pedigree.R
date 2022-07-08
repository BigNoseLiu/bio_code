##-------确认所依赖的包已正确安装-----------------------
args=commandArgs(T)
check_bioconductor_pkg <- function(pkgs) {
    if (require(pkgs, character.only = TRUE)) {
        print(paste(pkgs, "is loaded correctly"))
    }
    else {
        print(paste("Trying to install ", pkgs))
        update.packages(ask = FALSE)
        source("http://bioconductor.org/biocLite.R")
        biocLite(pkgs)
        ## Re-check
        if (require(pkgs, character.only = TRUE)) {
            print(paste(pkgs, "installed and loaded!"))
        }
        else {
            stop(paste("Error: could not install", pkgs))
        }
    }
}

check_pkg <- function(pkgs) {
    if (require(pkgs , character.only = TRUE)) {
        print(paste(pkgs, "is loaded correctly"))
    }
    else {
        print(paste("Trying to install ", pkgs))
        update.packages(ask = FALSE)
        install.packages(pkgs)
        ## Re-check
        if (require(pkgs, character.only = TRUE)) {
            print(paste(pkgs, "installed and loaded!"))
        }
        else {
            stop(paste("Error: could not install", pkgs))
        }
    }
}
##---------------------------------------------

check_pkg("kinship2")
check_pkg("reshape2")

##---------------------------------------------
current_time <- Sys.time()
file_name_tag <- format(current_time, format = "%Y%m%d%H%M%S")
input_file_name <- args[1]
output_file_name <- args[2]
#output_file_name <- "plot_pedigree"
#input_file_name <- "kinship.csv"

##-----------------------------------------------------------------------------
ped_data_frame <- read.csv(file = input_file_name)

relationship_data_frame <-
    ped_data_frame[c("id", "relation.id2", "relation")]
relationship_data_frame <-
    relationship_data_frame[which(relationship_data_frame$relation > -1), ]
relationship_matrix <- as.matrix(relationship_data_frame)

ped_names <- names(ped_data_frame)
affected_indices <- grep("affected", ped_names)
affected_data_frame <- ped_data_frame[affected_indices]
names(affected_data_frame) <- sub("affected.", "", names(affected_data_frame))
affected_matrix <- as.matrix(affected_data_frame)

label_indices <- grep("label", ped_names)
label_data_frame <- ped_data_frame[label_indices]

labels <- apply(label_data_frame,1,paste,collapse="\n")

ped <- pedigree(
    id = ped_data_frame$id,
    dadid = ped_data_frame$father,
    momid = ped_data_frame$mother,
    sex = ped_data_frame$sex,
    status = ped_data_frame$death_or_not,
    affected = affected_matrix,
    relation = relationship_matrix
)
ped_data_frame$sample.avail

ped_data_frame$sample.avail <- ifelse(is.na(ped_data_frame$sample.avail), 0, ped_data_frame$sample.avail)
ped_data_frame$sample.avail

pdf(file = paste(output_file_name, ".pdf", sep = ""))
plot(ped, id = labels,
     col=ifelse(ped_data_frame$sample.avail, 2, 1))
title(main = args[3])
pedigree.legend(ped, location="topright", radius=0.2)
dev.off()

png(file = paste(output_file_name, ".png", sep = ""))
plot(ped, 
     id = labels,
     col=ifelse(ped_data_frame$sample.avail, 2, 1))
title(main = args[3])
pedigree.legend(ped, location="topright", radius=.3)
dev.off()
