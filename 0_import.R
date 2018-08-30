# file import -------------------------------------------------------------

  rm(list = ls())
  setwd("//homer.uit.no/mse041/Desktop/ms/ham")


# positif ve negatif olarak mass dosyalari - listele , oku, rename

  positives<-list.files(pattern = "pos") 
  negatives<-list.files(pattern = "neg")
  pos_data<-lapply(positives, read.csv)
  names(pos_data)<-stringr::str_replace(positives,pattern = ".csv",replacement = "")
  neg_data<-lapply(negatives, read.csv)
  names(neg_data)<-stringr::str_replace(negatives,pattern = ".csv",replacement = "")

# esi kolonu ekleme
  
  pos_data<-lapply(pos_data, function(x) cbind(x,esi=1))
  neg_data<-lapply(neg_data, function(x) cbind(x,esi=-1))

# sample kolonu ekleme
  
  pos_st<-names(pos_data)
  pos_st<-as.numeric(stringr::str_replace(pos_st,pattern = "pos",replacement = ""))
  for (i in 1:length(pos_data)){
    pos_data[[i]]$sample<-pos_st[i]
    
  }
  

  neg_st<-names(neg_data)
  neg_st<-as.numeric(stringr::str_replace(neg_st,pattern = "neg", replacement = ""))
  for (i in 1:length(neg_data)){
    neg_data[[i]]$sample<-neg_st[i]
  }
  
# Na ve Cl colon ismi - nacl olarak degistirme

  for (i in 1:length(pos_data)){
    colnames(pos_data[[i]])[colnames(pos_data[[i]])=="Na"]<-"nacl"
  }
  
  for (i in 1:length(neg_data)){
    colnames(neg_data[[i]])[colnames(neg_data[[i]])=="Cl"]<-"nacl"
  }
  
# Index colonu olusturma

  for (i in 1:length(pos_data)){
   pos_data[[i]]$index<-paste(pos_data[[i]]$C,pos_data[[i]]$H,pos_data[[i]]$N,
                             pos_data[[i]]$O,pos_data[[i]]$S,pos_data[[i]]$nacl,sep=" ")
  }
  
  for (i in 1:length(neg_data)){
    neg_data[[i]]$index<-paste(neg_data[[i]]$C,neg_data[[i]]$H,neg_data[[i]]$N,
                               neg_data[[i]]$O,neg_data[[i]]$S,neg_data[[i]]$nacl,sep=" ")
  }
  
  
# Aromaticity index ve category
  
  for (i in 1:length(pos_data)){
    pos_data[[i]]$AI <- (1+pos_data[[i]]$C-pos_data[[i]]$O-pos_data[[i]]$S-(0.5*(pos_data[[i]]$H)))/
      (pos_data[[i]]$C-pos_data[[i]]$O-pos_data[[i]]$S-pos_data[[i]]$N)
  }
  
  for (i in 1:length(neg_data)){
    neg_data[[i]]$AI <- (1+neg_data[[i]]$C-neg_data[[i]]$O-neg_data[[i]]$S-(0.5*(neg_data[[i]]$H)))/
      (neg_data[[i]]$C-neg_data[[i]]$O-neg_data[[i]]$S-neg_data[[i]]$N)
  }

  
  for (i in 1:length(pos_data)){
    pos_data[[i]]$AIC <- ifelse(pos_data[[i]]$AI > 0.5 & pos_data[[i]]$AI < 0.67,"A",
                                ifelse(pos_data[[i]]$AI<1 & pos_data[[i]]$AI >= 0.67,"CA",
                                       "NA"))
  }
  
  for (i in 1:length(neg_data)){
    neg_data[[i]]$AIC <- ifelse(neg_data[[i]]$AI > 0.5 & neg_data[[i]]$AI < 0.67,"A",
                                ifelse(neg_data[[i]]$AI<1 & neg_data[[i]]$AI >= 0.67,"CA",
                                       "NA"))
  }
  
  
# CHO kolonu olusturma
  
  for (i in 1:length(pos_data)){
    pos_data[[i]]$CHO <- ifelse(pos_data[[i]]$N==0 & pos_data[[i]]$S ==0,"CHO",
                         ifelse(pos_data[[i]]$N==1 & pos_data[[i]]$S == 0,"CHON",
                         ifelse(pos_data[[i]]$N==2 & pos_data[[i]]$S == 0,"CHON2",
                         ifelse(pos_data[[i]]$N==1 & pos_data[[i]]$S == 1,"CHONS",
                                                     "CHON2S"))))
  }
  
  
  for (i in 1:length(neg_data)){
    neg_data[[i]]$CHO <- ifelse(neg_data[[i]]$N==0 & neg_data[[i]]$S ==0,"CHO",
                         ifelse(neg_data[[i]]$N==1 & neg_data[[i]]$S == 0,"CHON",
                         ifelse(neg_data[[i]]$N==2 & neg_data[[i]]$S == 0,"CHON2",
                         ifelse(neg_data[[i]]$N==1 & neg_data[[i]]$S == 1,"CHONS",
                                                     "CHON2S"))))
  }
  


# station kolonu ekleme
  
  for (i in 1:length(pos_data)) {
    pos_data[[i]]$station <-ifelse(pos_data[[i]]$sample %in% c(1, 2, 3, 4), 769,
        ifelse(pos_data[[i]]$sample %in% c(5, 6, 7, 8), 770 ,
        ifelse(pos_data[[i]]$sample %in% c(9, 10, 11, 12),780 ,
        ifelse(pos_data[[i]]$sample %in% c(13, 14, 15, 16),782 ,
        ifelse(pos_data[[i]]$sample %in% c(17, 18, 19, 20), 806 ,
        ifelse(pos_data[[i]]$sample %in% c(21, 22, 23, 24), 821,
        ifelse(pos_data[[i]]$sample %in% c(25, 26, 27, 28), 827 ,
        ifelse(pos_data[[i]]$sample %in% c(29, 30, 31, 32), 830,
        ifelse(pos_data[[i]]$sample %in% c(33, 34, 35, 36, 37, 38),834 ,
        ifelse(pos_data[[i]]$sample %in% c(39, 40, 41, 42, 43),842,
        ifelse(pos_data[[i]]$sample %in% c(44, 45, 46, 47, 48),863,
        ifelse(pos_data[[i]]$sample %in% c(49, 50, 51, 52, 53),867,
        ifelse(pos_data[[i]]$sample %in% c(54, 55, 56, 57),907 ,
        ifelse(pos_data[[i]]$sample %in% c(58, 59, 60, 61), 909,
        ifelse(pos_data[[i]]$sample %in% c(62, 63, 64, 65),912,
        ifelse(pos_data[[i]]$sample %in% c(66, 67, 68, 69), 936,
        ifelse(pos_data[[i]]$sample %in% c(70, 71, 72, 73, 74, 75), 944,
        ifelse(pos_data[[i]]$sample %in% c(76, 77, 78, 79, 80), 949 ,
                             "NA"))))))))))))))))))
                                      
  }
  
  for (i in 1:length(neg_data)) {
    neg_data[[i]]$station <-ifelse(neg_data[[i]]$sample %in% c(1, 2, 3, 4), 769,
       ifelse(neg_data[[i]]$sample %in% c(5, 6, 7, 8), 770 ,
       ifelse(neg_data[[i]]$sample %in% c(9, 10, 11, 12),780 ,
       ifelse(neg_data[[i]]$sample %in% c(13, 14, 15, 16),782 ,
       ifelse(neg_data[[i]]$sample %in% c(17, 18, 19, 20), 806 ,
       ifelse(neg_data[[i]]$sample %in% c(21, 22, 23, 24), 821,
       ifelse(neg_data[[i]]$sample %in% c(25, 26, 27, 28), 827 ,
       ifelse(neg_data[[i]]$sample %in% c(29, 30, 31, 32), 830,
       ifelse(neg_data[[i]]$sample %in% c(33, 34, 35, 36, 37, 38),834 ,
       ifelse(neg_data[[i]]$sample %in% c(39, 40, 41, 42, 43),842,
       ifelse(neg_data[[i]]$sample %in% c(44, 45, 46, 47, 48),863,
       ifelse(neg_data[[i]]$sample %in% c(49, 50, 51, 52, 53),867,
       ifelse(neg_data[[i]]$sample %in% c(54, 55, 56, 57),907 ,
       ifelse(neg_data[[i]]$sample %in% c(58, 59, 60, 61), 909,
       ifelse(neg_data[[i]]$sample %in% c(62, 63, 64, 65),912,
       ifelse(neg_data[[i]]$sample %in% c(66, 67, 68, 69), 936,
       ifelse(neg_data[[i]]$sample %in% c(70, 71, 72, 73, 74, 75), 944,
       ifelse(neg_data[[i]]$sample %in% c(76, 77, 78, 79, 80), 949 ,
                             "NA"))))))))))))))))))
    
  }
  
  
  
  
  
  
  
  
# positif ve negatif birlestirme
 
 
  neg_data2<- neg_data[-c(1,50)]
  all_data<-mapply(rbind,neg_data2,pos_data, SIMPLIFY = FALSE) # positif negatif birlestirme
  
  com_data<-lapply(all_data,function(x) x[!duplicated(x$index),]) # duplicate olanalarin atilmasi
  
  positives2<-stringr::str_replace(positives,pattern ="pos", replace="com")
  coms<-stringr::str_replace(positives2,pattern =".csv", replace="")
  names(com_data)<-coms
  
  
  # setwd("//homer.uit.no/mse041/Desktop/ms/uni")
  # lapply(1:length(uni_data), function(i) write.csv(uni_data[[i]], 
  #                                                 file = paste0(names(uni_data[i]), ".csv"),
  #                                                 row.names = TRUE))
  
  
# remove all except lists  
  
  rm(list=setdiff(ls(), c("neg_data","pos_data","com_data")))
