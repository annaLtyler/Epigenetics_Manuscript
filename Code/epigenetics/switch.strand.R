#When a gene is identifies as being on the 
#negative strand, this function switches
#the start and end feature names so other
#scripts can look in the correct place for
#gene features

switch.strand <- function(start.feature, end.feature){
	feature.names <- c(start.feature, end.feature)
	start.locale <- grep("start", feature.names)
	end.locale <- grep("end", feature.names)

	if(length(start.locale) > 0){
		feature.names[start.locale] <- gsub("start", "end", feature.names[start.locale])
		}
	if(length(end.locale) > 0){
		feature.names[end.locale] <- gsub("end", "start", feature.names[end.locale])
		}
	return(feature.names)	
	}
