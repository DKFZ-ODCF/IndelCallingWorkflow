library(tidyverse)
library(optparse)
library(grid)
library(gridExtra)
library(Canopy)
library(jsonlite)

### Reading in arguments

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="All base coverage file"),
  make_option(c("-oP", "--oPlot"), type="character", default=NULL, help="Output png file path"),
  make_option(c("-oF", "--oFile"), type="character", default=NULL, help="Output table file path"),
  make_option(c("-p", "--pid"), type="character", default=NULL, help="Name of the pid"),
  make_option(c("-c", "--chrLength"), type="character", default=NULL, help="Chromosomes length file"),
  make_option(c("-s", "--cFunction"), type="character", default=NULL, help="Updated canopy function"),
  make_option(c("-t", "--SeqType"), type="character", default = NULL, help="WES or WGS"),
  make_option(c("-r", "--rightBorder"), type="character", default = NULL, help="Maximum control AF"),
  make_option(c("-b", "--bottomBorder"), type="character", default = NULL, help="Minimum tumor AF")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$file)) {
  print_help(opt_parser)
  stop("Input file name missing\n", call.=F)
} else if (is.null(opt$oPlot)) {
    print_help(opt_parser)
    stop("ouput png file name missing\n", call.=F)
} else if (is.null(opt$oFile)) {
    print_help(opt_parser)
    stop("ouput table file name missing\n", call.=F)
} else if(is.null(opt$pid)) {
    print_help(opt_parser)
    stop("Name of the PID missing\n", call.=F)
} else if(is.null(opt$cFunction)) {
    print_help(opt_parser)
    stop("Canopy updated function not provided\n", call.=F)
} else if(is.null(opt$SeqType)) {
    print_help(opt_parser)
    stop("Sequence type not provided\n", call.=F)
}

## 
source(opt$cFunction)
##### Data Analysis
# read the Rare.txt file and chromosome left file
dat<-read.delim(opt$file, header=T, sep="\t")
chr.length <- read.table(opt$chrLength, header=T)

# Initial cluster centroid
clusterCentroid <- function (seqType, maxControl=0.45, minTumor=0.01){
  if(seqType == "WGS") {
    mu.init <- cbind(c(0.5, 0.95, 0.5,  0.5,  0.5,  0.5,  0.02, 0.02, 0.02, 0.02, 0.10), 
                     c(0.5, 0.95, 0.25, 0.75, 0.95, 0.05, 0.30, 0.5,  0.95, 0.10, 0.10))
    numberCluster <- 11
  } else if(seqType == "WES") {
    mu.init <- cbind(c(0.5, 0.02, 0.02, 0.1, 0.02), 
                     c(0.5, 0.30, 0.5,  0.1, 0.10))
    numberCluster <- 5
    maxControl <- 0.35 # since few points to form a stable cluster 
  }
  return(list("mu.init" = mu.init, "numberCluster"= numberCluster, 
              "maxControl" = maxControl, "minTumor" = minTumor))
}

centroid <- clusterCentroid(opt$SeqType, maxControl = as.numeric(opt$rightBorder),
                            minTumor = as.numeric(opt$bottomBorder))

# Running Canopy
R <-as.matrix(dat[,c(7,9)])
X <-as.matrix(dat[,c(8,10)])

# Canopy run and assigning centers
canopy.clust <- tryCatch(
  {
  canopy.cluster(R, X, num_cluster = centroid$numberCluster, 
                             num_run = 1, Mu.init = centroid$mu.init)
  }, error = function(e){
    if(opt$SeqType == "WGS") {
    centroid <- clusterCentroid("WES")
    print("Error catched!")
    canopy.cluster(R, X, num_cluster = centroid$numberCluster, 
                   num_run = 1, Mu.init = centroid$mu.init)
    } else if(opt$SeqType == "WES") {
      print("Error in centriod assignment in WES data")
    }
  }
)

dat$canopyCluster<-canopy.clust$sna_cluster

## Select the TiN cluster
somaticClass <- dat %>%  
  mutate(squareRescue = Control_AF < centroid$maxControl & 
           Tumor_AF > centroid$minTumor) %>% 
  mutate(diagonalRescue = Control_AF < Tumor_AF) %>% 
  group_by(canopyCluster) %>%
  dplyr::summarise(prop1 = mean(squareRescue == T), prop2 = mean(diagonalRescue==T)) %>% 
  filter(prop1 > 0.85 & prop2 > 0.85) %>% 
  select(canopyCluster) %>% collect() %>% .[[1]]

### Removing unusual clusters
for(cl in somaticClass){
  cluster.size <- nrow(dat[dat$canopyCluster == cl,])
  control.max.line <- dat[dat$canopyCluster == cl,][which.max(dat[dat$canopyCluster == cl,]$Control_AF),]
  
  count <- nrow(dat[dat$Control_AF <=  control.max.line$Control_AF & 
                      dat$Tumor_AF > control.max.line$Tumor_AF & 
                      !(dat$canopyCluster %in% somaticClass),])
  
  if(count > cluster.size/4){
    somaticClass <- somaticClass[!somaticClass == cl]
  }
}

## Trying to rescue the homozygous cluster left alone
if(length(somaticClass) != 0) {
  dat %>% group_by(canopyCluster) %>% 
      summarise(max_C = max(Control_AF), med_T = quantile(Tumor_AF, 0.80)) %>% 
      filter(canopyCluster %in% somaticClass ) %>% filter(med_T == max(med_T)) -> homo.threshold
  
  # Couting the total number of TiN variants if we rescue the homozygous variants we well
  dat %>% filter((Tumor_AF > homo.threshold$med_T & 
                   Control_AF < homo.threshold$max_C) |
                   canopyCluster %in% somaticClass) %>% nrow() -> TiN.homo.updated
  # Only the variants from Tin Clusters
  dat %>% filter(canopyCluster %in% somaticClass) %>% nrow() -> TiN.not.homo.updated
  
  # The Rescued TiN should only add 25% more variants
  if(TiN.homo.updated <= (TiN.not.homo.updated + (TiN.not.homo.updated/4))) {
    dat %>% 
      mutate(TiN_Class = ifelse((Tumor_AF > homo.threshold$med_T & 
                                 Control_AF < homo.threshold$max_C) |
                                  Control_AF == 0 | 
                                canopyCluster %in% somaticClass, "Somatic_Rescue", "Germline")) -> dat
  } else {
    dat %>% mutate(TiN_Class = ifelse(Control_AF == 0 | 
                                        canopyCluster %in% somaticClass, "Somatic_Rescue", "Germline")) -> dat
  } 
} else {
  dat %>% mutate(TiN_Class = ifelse(Control_AF == 0, "Somatic_Rescue", "Germline")) -> dat
}

## Exome removing common variants from somatic rescue
dat %>% filter(grepl("Somatic_Rescue", TiN_Class)) %>% 
  mutate(TiN_Class = ifelse(grepl("Common", Rareness) & 
                              grepl("Somatic_Rescue", TiN_Class), "Germline", "Somatic_Rescue")) -> somRes

dat %>% filter(grepl("Germline", TiN_Class)) %>% 
  rbind(somRes) -> dat

#dat %>% group_by(Rareness, TiN_Class) %>% summarise(count=n())

# Plot 1 with canopy cluster
poly.df <- data.frame(x=c(0, 0, centroid$maxControl, centroid$maxControl), 
                      y=c(centroid$minTumor, 1, 1, centroid$maxControl))

p1 <- ggplot() + geom_point(aes(Control_AF, Tumor_AF, color=factor(canopyCluster)), alpha=0.5, data=dat) + 
  theme_bw() + theme(text = element_text(size=15), legend.position="bottom") + 
  xlab("Control VAF") + ylab("Tumor VAF") + 
  xlim(0,1) + ylim(0,1) +
  guides(color=guide_legend("Canopy clusters")) + 
  ggtitle(paste0("Clusters from Canopy")) + 
  geom_polygon(data=poly.df, aes(x, y), alpha=0.2, fill="gold")

# Plot 2 with TiN cluster
p2 <- ggplot() + geom_polygon(data=poly.df, aes(x, y), alpha=0.2, fill="#d8161688") + 
  geom_point(aes(Control_AF, Tumor_AF, color=TiN_Class), alpha=0.3, data=dat) + 
  theme_bw() + theme(text = element_text(size=15), legend.position="bottom") + 
  xlab("Control VAF") + ylab("Tumor VAF") + 
  xlim(0,1) + ylim(0,1) +
  guides(color=guide_legend("TiN clusters")) + 
  ggtitle(paste0("TiN clusters")) + 
  geom_polygon(data=poly.df, aes(x, y), alpha=0.2, fill="gold")
  
## function to plot linear chromosome
plotGenome_ggplot <- function(data, Y, chr.length, colorCol) {
  chr.length$cumShiftLength <- cumsum(as.numeric(chr.length$shiftLength))
  chr.length$cumLength <- cumsum(as.numeric(chr.length$Length))
  chrTotalLength <- sum(as.numeric(chr.length$Length))
  
  data %>% mutate(CHR=as.factor(CHR)) %>% left_join(chr.length, by="CHR") %>% 
    mutate(cumPOS = POS + cumShiftLength) %>% 
    ggplot() + geom_point(aes_string(x="cumPOS", y=Y, color=colorCol), shape=124) + theme_bw() + 
    scale_x_continuous(breaks = chr.length$cumLength - (chr.length$Length/2),
                       labels = chr.length$CHR,
                       minor_breaks = chr.length$cumShiftLength,
                       expand = c(0, 0)) + 
    theme(axis.title.x = element_blank(), text = element_text(size=15),
          panel.grid.major.x = element_line(color="lightgrey", linetype = 0),
          panel.grid.minor.x = element_line(color="grey")) 
}

## Plot whole genome
p3<-plotGenome_ggplot(dat, 'Control_AF', chr.length, 'TiN_Class')
p4<-plotGenome_ggplot(dat, 'Tumor_AF', chr.length, 'TiN_Class')

#### multi plotting
# Blank region 
#blank <- grid.rect(gp=gpar(col="white"))

# Rescue info table
#rescueInfo<-as.data.frame(table(dat$TiN_Class))
#colnames(rescueInfo)<-c("Reclassification", "Counts")

rescueInfo <- dat %>% 
  group_by(TiN_Class) %>% 
  summarise(Count =n(), Median_Control_VAF = formatC(median(Control_AF), digits=5, format="f"),
                        Median_Tumor_VAF = formatC(median(Tumor_AF), digits=5, format="f"))

rescueInfo.toFile <- rescueInfo
if("Somatic_Rescue" %in% rescueInfo.toFile$TiN_Class) {
  rescueInfo.toFile$Pid<-opt$pid
} else {
  rescueInfo.toFile<-rbind(rescueInfo.toFile, c("Somatic_Rescue", 0, 0, 0))
  rescueInfo.toFile$Pid<-opt$pid
}
write_json(rescueInfo.toFile, path=paste0(opt$oFile, "_summary.json"))

TableTheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 1, hjust=1, x=0.95)),
  colhead = list(fg_params=list(cex = 1, hjust=1, x=0.95)),
  rowhead = list(fg_params=list(cex = 1, hjust=1, x=0.95)))

TableAnn <-tableGrob(rescueInfo, rows = c(), theme=TableTheme)

PlotLayout <-rbind(c(1,2,3),
                   c(1,2,3),
                   c(4,4,4),
                   c(5,5,5))

### Writing as png file


png(file = opt$oPlot, width=1500, height=800)
grid.arrange(p1, p2, TableAnn, p3, p4, 
             layout_matrix = PlotLayout,
             top=textGrob(paste0('TiN Detection Analysis - TiNDA : ', opt$pid), gp=gpar(cex=2)))

dev.off()

# Saving the rescue table file 
write.table(dat, file=opt$oFile, sep="\t", row.names = F, quote = F)

## 
#reg.finalizer(environment(), cleanup, onexit = FALSE)
