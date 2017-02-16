library(ggplot2)
library(tidyr)
read.outfile.admb<-function(path.string)
{
	tt.chars <- scan(path.string, sep = "", what = "",quiet=TRUE)
	tt.nums <- as.numeric(tt.chars)
	list.length <- length(tt.nums[is.na(tt.nums)])
	out <- vector("list", list.length)
	object.counter <- 1
	element.counter <- 1
	while(object.counter <= list.length) {
		names(out[object.counter]) <- tt.chars[element.counter]
		print(tt.chars[element.counter])
		element.counter <- element.counter + 1
		tt.rows <- tt.nums[element.counter]
		element.counter <- element.counter + 1
		tt.cols <- tt.nums[element.counter]
		element.counter <- element.counter + 1
		if(tt.rows > 1) {
			out[[object.counter]] <- matrix(data = tt.nums[element.counter:(element.counter + tt.rows * tt.cols - 1)], nrow = 
				tt.rows, ncol = tt.cols, byrow = T)
			element.counter <- element.counter + (tt.rows * tt.cols)
		}
		if(tt.rows == 1) {
			out[[object.counter]] <- c(tt.nums[element.counter:(element.counter + tt.cols - 1)])
			element.counter <- element.counter + tt.cols
		}
		object.counter <- object.counter + 1
	}
	out.names  <- tt.chars[is.na(as.numeric(tt.chars))]
	names(out) <- out.names
	return(out)
}
#' Write and run model configu
#' can do retrospectives
#' TODO: Insert more information here.
#'
#' @param Nretro how many peels to go back
#' @param yrs_sel_change vector of years that selectivity is changed
#' @param Lambda vector penalty weights
#' @return list of results for each retrospective
#' @export
Do_Run<-function(Nretro=10,Lambda=c(200,1,1,1,1,.1,.1,12.5,1,1),yrs_sel_change,rn="Retro_",nselages=8,run=TRUE,keeporig=TRUE)
{
res <-list()
if (keeporig) file.copy("ssas.dat","ssas_orig.dat",overwrite = TRUE)
for (i in 0:Nretro){
  if (run){
    cat(file="ssas.dat",
    "#Lambdas \n", 
    Lambda,
    "#Retro number \n", 
    i,"\n 
    # Datafile \n 
    ssas_1.dat \n 
    #Routput file \n ",
    paste0(rn,i,".rep \n"),
    "#Number of ages to estimate fishery selectivity \n",
    nselages, " \n
    #Number of yrs fishery sel changes \n",
    length(yrs_sel_change), " \n
    #Yrs of change \n",
    yrs_sel_change," \n" )
    system("./ssas -nox -iprint 100")
  }
  if (Nretro>0)
    res[[i+1]] <- read.outfile.admb(paste0(rn,i,".rep"))
  else
    res <- read.outfile.admb(paste0(rn,i,".rep"))
}
if (keeporig) file.copy("ssas_orig.dat","ssas.dat",overwrite = TRUE)
return(res)
}
mytheme <- theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank() )
mytheme <- mytheme + theme(text=element_text(size=18)) + theme(axis.title.x=element_text(size=24) ,axis.title.y=element_text(size=24))
mytheme <- mytheme + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_line(colour="grey60", linetype="dashed"), panel.grid.major.y = element_blank() )
mytheme <- mytheme + theme( panel.background = element_rect(fill="white"), panel.border = element_rect(colour="black", fill=NA, size=1))

