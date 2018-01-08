library(ggplot2)
library(data.table)
library(sqldf)

# [clear all variables]
rm(list=ls())
cutoff <- 5

# [set working dir]
setwd("/media/panisa/BKWORK/BK_Thesis_260815/IS_work/")


# [loading data to table]
dat <- read.table("output_B_15/info_bwa_8000847_7132693", sep="\t", header=TRUE)


# [set working dir]
setwd("/media/panisa/BKWORK/BK_Thesis_260815/IS_work/output_B_15")



# [define variables]
ref_id <- dat$ref_id
read_name <- dat$read_name
read_start <-dat$read_start
read_end <- dat$read_end
read_length <- dat$read_length
ref_start <- dat$ref_start
ref_end <- dat$ref_end
ref_length <- dat$ref_length
cigar <- dat$cigar
clip_start <-dat$clip_start
clip_end <- dat$clip_end
strand <- dat$strand
unmapped <- dat$unmapped
seq <- dat$seq
genome <- 'bwa_8000847_7132693'



# [create table group by position of ref_start]
DT <- data.table(dat)
GN <- data.table(genome)
# DT_start <- DT[, sum(clip_start), by = list(ref_id,ref_start)]
# DT_end <- DT[, sum(clip_end), by = list(ref_id,ref_end)]

DT_start <- sqldf('select ref_id, ref_start, strand, sum(clip_start) as num_clip_start, (ref_start+15) as ref_start_upper
	from DT 
	where clip_start = 1
	group by ref_id, ref_start 
	having num_clip_start >= 1
	order by ref_id')

DT_end <- sqldf('select ref_id, ref_end, strand, sum(clip_end) as num_clip_end, (ref_end-15) as ref_end_lower 
	from DT 
	where clip_end = 1
	group by ref_id, ref_end
	having num_clip_end >= 1
	order by ref_id')

DT_SE <- sqldf('select * 
	from DT_start DS 
	natural join (select * from DT_end) DE 
	where (DS.ref_id = DE.ref_id) and (DE.ref_end between DS.ref_start and DS.ref_start_upper) and DS.num_clip_start >= 5 and DE.num_clip_end >= 5')

DT_ES <- sqldf('select *
	from DT_end E
	natural join (select * from DT_start) S
	where (E.ref_id = S.ref_id) and (S.ref_start between E.ref_end_lower and E.ref_end) and E.num_clip_end >= 5 and S.num_clip_start >= 5')

DT_detail_list <- sqldf('select ref_id, ref_start, ref_end, strand, GN.genome as genome, num_clip_start, num_clip_end
	from DT_SE
	natural join (select * from GN) GN')


# [write files]
write.table(DT_start, "count_start_pos.xls", sep="\t", col.names = NA, row.names = TRUE)
write.table(DT_end, "count_end_pos.xls", sep="\t", col.names = NA, row.names = TRUE)
write.table(DT_SE, "count_start_end_5.xls", sep="\t", col.names = NA, row.names = TRUE)
write.table(DT_ES, "count_end_start_5.xls", sep="\t", col.names = NA, row.names = TRUE)
write.csv(DT_detail_list, "position/detail_list.csv")



-----

# [plot histograms]
DT_start_hist <- ggplot(DT_start, aes(x=ref_start, fill=ref_id)) + geom_histogram(binwidth=1000, alpha=1, position="dodge") + theme(legend.position="bottom")
DT_end_hist <- ggplot(DT_end, aes(x=ref_end, fill=ref_id)) + geom_histogram(binwidth=100, alpha=1, position="dodge") + theme(legend.position="bottom")


# [save file]
png(filename="hist_start_ref.png", width=18, height=12, units="in", res=300)
DT_start_hist + ylab("Freq. of clip_start")
dev.off()

png(filename="hist_end_ref.png", width=18, height=12, units="in", res=300)
DT_end_hist + ylab("Freq. of clip_end")
dev.off()


