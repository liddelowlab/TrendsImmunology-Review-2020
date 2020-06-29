#demultiplexing pipeline written by Taitea Dyskra
#purpose: to prepare barcode, feature, and matrix files for each donor in Mathys et al. (2019) aggregated dataset

library(data.table)
library(Matrix)

setwd("/gpfs/data/liddelowlab/data/Mathys_2019")

barcode.names <- read.csv(file="notfiltered_column_metadata.txt", header=TRUE, sep="\t")

barcode <- barcode.names[grepl(".1", barcode.names$TAG), ]
barcode <- barcode[-(grep(".41", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS13/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".2", barcode.names$TAG), ]
barcode <- barcode[-(grep(".42", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS45/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".3", barcode.names$TAG), ]
barcode <- barcode[-(grep(".43", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS14/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".4", barcode.names$TAG), ]
barcode <- barcode[-(grep(".48", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS38/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".5", barcode.names$TAG), ]
barcode <- barcode[-(grep(".45", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS1/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".6", barcode.names$TAG), ]
barcode <- barcode[-(grep(".46", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS34/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".7", barcode.names$TAG), ]
barcode <- barcode[-(grep(".47", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS15/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".8", barcode.names$TAG), ]
barcode <- barcode[-(grep(".48", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS37/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".9", barcode.names$TAG), ]
barcode <- barcode[-(grep(".39", barcode$TAG)),]
barcode <- barcode$TAG
write.table(barcode, file='./ROS2/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".10", barcode.names$TAG),1]
write.table(barcode, file='./ROS32/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".11", barcode.names$TAG),1]
write.table(barcode, file='./ROS3/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".12", barcode.names$TAG),1]
write.table(barcode, file='./ROS29/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".13", barcode.names$TAG),1]
write.table(barcode, file='./ROS16/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".14", barcode.names$TAG),1]
write.table(barcode, file='./ROS41/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".15", barcode.names$TAG),1]
write.table(barcode, file='./ROS17/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".16", barcode.names$TAG),1]
write.table(barcode, file='./ROS43/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".17", barcode.names$TAG),1]
write.table(barcode, file='./ROS4/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".18", barcode.names$TAG),1]
write.table(barcode, file='./ROS33/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".19", barcode.names$TAG),1]
write.table(barcode, file='./ROS18/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".20", barcode.names$TAG),1]
write.table(barcode, file='./ROS39/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".21", barcode.names$TAG),1]
write.table(barcode, file='./ROS5/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".22", barcode.names$TAG),1]
write.table(barcode, file='./ROS25/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".23", barcode.names$TAG),1]
write.table(barcode, file='./ROS6/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".24", barcode.names$TAG),1]
write.table(barcode, file='./ROS31/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".25", barcode.names$TAG),1]
write.table(barcode, file='./ROS11/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".26", barcode.names$TAG),1]
write.table(barcode, file='./ROS26/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".27", barcode.names$TAG),1]
write.table(barcode, file='./ROS12/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".28", barcode.names$TAG),1]
write.table(barcode, file='./ROS35/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".29", barcode.names$TAG),1]
write.table(barcode, file='./ROS19/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".30", barcode.names$TAG),1]
write.table(barcode, file='./ROS44/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".31", barcode.names$TAG),1]
write.table(barcode, file='./ROS7/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".32", barcode.names$TAG),1]
write.table(barcode, file='./ROS36/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".33", barcode.names$TAG),1]
write.table(barcode, file='./ROS8/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".34", barcode.names$TAG),1]
write.table(barcode, file='./ROS30/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".35", barcode.names$TAG),1]
write.table(barcode, file='./ROS20/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".36", barcode.names$TAG),1]
write.table(barcode, file='./ROS48/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".37", barcode.names$TAG),1]
write.table(barcode, file='./ROS9/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".38", barcode.names$TAG),1]
write.table(barcode, file='./ROS28/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".39", barcode.names$TAG),1]
write.table(barcode, file='./ROS21/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".40", barcode.names$TAG),1]
write.table(barcode, file='./ROS42/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".41", barcode.names$TAG),1]
write.table(barcode, file='./ROS22/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".42", barcode.names$TAG),1]
write.table(barcode, file='./ROS46/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".43", barcode.names$TAG),1]
write.table(barcode, file='./ROS10/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".44", barcode.names$TAG),1]
write.table(barcode, file='./ROS27/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".45", barcode.names$TAG),1]
write.table(barcode, file='./ROS23/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".46", barcode.names$TAG),1]
write.table(barcode, file='./ROS47/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".47", barcode.names$TAG),1]
write.table(barcode, file='./ROS24/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcode <- barcode.names[grepl(".48", barcode.names$TAG),1]
write.table(barcode, file='./ROS40/barcodes.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)


genes <- read.csv(file="notfiltered_gene_row_names.txt", header=FALSE, sep="\t")
genes$V3 <- c("Gene Expression")
write.table(genes, file='./ROS1/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS2/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS3/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS4/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS5/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS6/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS7/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS8/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS9/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS10/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS11/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS12/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS13/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS14/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS15/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS16/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS17/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS18/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS19/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS20/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS21/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS22/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS23/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS24/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS25/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS26/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS27/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS28/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS29/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS30/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS31/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS32/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS33/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS34/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS35/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS36/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS37/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS38/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS39/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS40/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS41/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS42/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS43/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS44/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS45/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS46/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS47/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genes, file='./ROS48/features.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)


mat <- readMM(file = "notfiltered_count_matrix.mtx")

barcode <- read.csv(file='./ROS13/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,1:512]
writeMM(obj = submat, file="./ROS13/mat.mtx")

barcode <- read.csv(file='./ROS45/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,513:898]
writeMM(obj = submat, file="./ROS45/mat.mtx")

barcode <- read.csv(file='./ROS14/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,899:2131]
writeMM(obj = submat, file="./ROS14/mat.mtx")

barcode <- read.csv(file='./ROS38/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,2132:4340]
writeMM(obj = submat, file="./ROS38/mat.mtx")

barcode <- read.csv(file='./ROS1/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,4341:7433]
writeMM(obj = submat, file="./ROS1/mat.mtx")

barcode <- read.csv(file='./ROS34/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,7434:9163]
writeMM(obj = submat, file="./ROS34/mat.mtx")

barcode <- read.csv(file='./ROS15/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,9164:11380]
writeMM(obj = submat, file="./ROS15/mat.mtx")

barcode <- read.csv(file='./ROS37/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,11381:14845]
writeMM(obj = submat, file="./ROS37/mat.mtx")

barcode <- read.csv(file='./ROS2/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,14846:16383]
writeMM(obj = submat, file="./ROS2/mat.mtx")

barcode <- read.csv(file='./ROS32/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,16384:19607]
writeMM(obj = submat, file="./ROS32/mat.mtx")

barcode <- read.csv(file='./ROS3/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,19608:21052]
writeMM(obj = submat, file="./ROS3/mat.mtx")

barcode <- read.csv(file='./ROS29/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,21053:21976]
writeMM(obj = submat, file="./ROS29/mat.mtx")

barcode <- read.csv(file='./ROS16/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,21977:23809]
writeMM(obj = submat, file="./ROS16/mat.mtx")

barcode <- read.csv(file='./ROS41/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,23810:25351]
writeMM(obj = submat, file="./ROS41/mat.mtx")

barcode <- read.csv(file='./ROS17/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,25352:28842]
writeMM(obj = submat, file="./ROS17/mat.mtx")

barcode <- read.csv(file='./ROS43/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,28843:30472]
writeMM(obj = submat, file="./ROS43/mat.mtx")

barcode <- read.csv(file='./ROS4/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,30473:31800]
writeMM(obj = submat, file="./ROS4/mat.mtx")

barcode <- read.csv(file='./ROS33/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,31801:33629]
writeMM(obj = submat, file="./ROS33/mat.mtx")

barcode <- read.csv(file='./ROS18/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,33630:35242]
writeMM(obj = submat, file="./ROS18/mat.mtx")

barcode <- read.csv(file='./ROS39/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,35243:36493]
writeMM(obj = submat, file="./ROS39/mat.mtx")

barcode <- read.csv(file='./ROS5/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,36494:38330]
writeMM(obj = submat, file="./ROS5/mat.mtx")

barcode <- read.csv(file='./ROS25/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,38331:40035]
writeMM(obj = submat, file="./ROS25/mat.mtx")

barcode <- read.csv(file='./ROS6/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,40036:40821]
writeMM(obj = submat, file="./ROS6/mat.mtx")

barcode <- read.csv(file='./ROS31/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,40822:42300]
writeMM(obj = submat, file="./ROS31/mat.mtx")

barcode <- read.csv(file='./ROS11/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,42301:44022]
writeMM(obj = submat, file="./ROS11/mat.mtx")

barcode <- read.csv(file='./ROS26/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,44023:46429]
writeMM(obj = submat, file="./ROS26/mat.mtx")

barcode <- read.csv(file='./ROS12/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,46430:48672]
writeMM(obj = submat, file="./ROS12/mat.mtx")

barcode <- read.csv(file='./ROS35/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,48673:50681]
writeMM(obj = submat, file="./ROS35/mat.mtx")

barcode <- read.csv(file='./ROS19/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,50682:52428]
writeMM(obj = submat, file="./ROS19/mat.mtx")

barcode <- read.csv(file='./ROS44/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,52429:53167]
writeMM(obj = submat, file="./ROS44/mat.mtx")

barcode <- read.csv(file='./ROS7/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,53168:54130]
writeMM(obj = submat, file="./ROS7/mat.mtx")

barcode <- read.csv(file='./ROS36/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,54131:56593]
writeMM(obj = submat, file="./ROS36/mat.mtx")

barcode <- read.csv(file='./ROS8/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,56594:57821]
writeMM(obj = submat, file="./ROS8/mat.mtx")

barcode <- read.csv(file='./ROS30/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,57822:58888]
writeMM(obj = submat, file="./ROS30/mat.mtx")

barcode <- read.csv(file='./ROS20/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,58889:60465]
writeMM(obj = submat, file="./ROS20/mat.mtx")

barcode <- read.csv(file='./ROS48/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,60466:61771]
writeMM(obj = submat, file="./ROS48/mat.mtx")

barcode <- read.csv(file='./ROS9/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,61772:62659]
writeMM(obj = submat, file="./ROS9/mat.mtx")

barcode <- read.csv(file='./ROS28/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,62660:64768]
writeMM(obj = submat, file="./ROS28/mat.mtx")

barcode <- read.csv(file='./ROS21/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,64769:65626]
writeMM(obj = submat, file="./ROS21/mat.mtx")

barcode <- read.csv(file='./ROS42/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,65627:67620]
writeMM(obj = submat, file="./ROS42/mat.mtx")

barcode <- read.csv(file='./ROS22/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,67621:69550]
writeMM(obj = submat, file="./ROS22/mat.mtx")

barcode <- read.csv(file='./ROS46/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,69551:71015]
writeMM(obj = submat, file="./ROS46/mat.mtx")

barcode <- read.csv(file='./ROS10/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,71016:72949]
writeMM(obj = submat, file="./ROS10/mat.mtx")

barcode <- read.csv(file='./ROS27/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,72950:74450]
writeMM(obj = submat, file="./ROS27/mat.mtx")

barcode <- read.csv(file='./ROS23/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,74451:76260]
writeMM(obj = submat, file="./ROS23/mat.mtx")

barcode <- read.csv(file='./ROS47/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,76261:77940]
writeMM(obj = submat, file="./ROS47/mat.mtx")

barcode <- read.csv(file='./ROS24/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,77941:80273]
writeMM(obj = submat, file="./ROS24/mat.mtx")

barcode <- read.csv(file='./ROS40/barcodes.tsv', sep='\t', header=FALSE)
submat <- mat[ ,80274:80660]
writeMM(obj = submat, file="./ROS40/mat.mtx")
