#!/usr/bin/Rscript


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17


####################################################################################
# Barplot creation

create_barplot<-function(Nb,filename,db)
{
  png(filename)
  
  bplt<-barplot(Nb,
          ylim=c(0, max(Nb)+10000),
          main=paste("Editing Sites (",db,")"),
          xlab=("Differenciation stages"),
          ylab=("Editing Sites Number"),
          col="black"
  )
  
  text(x= bplt, y = Nb, label= Nb, pos= 3, cex = 0.8, col = "red")
  dev.off()
}

#####################################################################################
#main

#WARNING install.packages("stringr", dependencies=TRUE)
library(stringr)

#Argument recovery
args = commandArgs(trailingOnly=TRUE)
number=args[1]
NAME=args[2]
filename=args[3]
db=args[4]

Nb<-str_split(number,",")

c(for (n in Nb){
	number<-as.numeric(n)
})


name<-str_split(NAME,",")

c(for (n in name){
	list_name<-n
})
  names(number)=list_name

#Barplot creation
create_barplot(number,filename,db)

