#This function finds overlapping marks in two different
#chip tables. It assumes exact matches rather than true overlap

	find.chip.overlap <- function(chip.table1, chip.table2){
		chr <- unique(chip.table1[,1])
		
		common.mark.table <- NULL
		for(ch in 1:length(chr)){
			#first find all the peaks in the tables on the
			#same chromosome
			chr.locale1 <- which(chip.table1[,1] == chr[ch])
			subtable1 <- chip.table1[chr.locale1,]
			chr.locale2 <- which(chip.table2[,1] == chr[ch])
			subtable2 <- chip.table2[chr.locale2,]
			
			common.start <- intersect(subtable1[,"start"], subtable2[,"start"])
			
			common.locale1 <- match(common.start, subtable1[,"start"])
			common.marks1 <- subtable1[common.locale1,]
			
			common.locale2 <- match(common.start, subtable2[,"start"])
			common.marks2 <- subtable2[common.locale2,]
			
			result <- cbind(subtable1[common.locale1,], subtable2[common.locale2,4:9])
			
			common.mark.table <- rbind(common.mark.table, result)
			}
		return(common.mark.table)
		}

