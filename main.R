library(arules)
library(arulesSequences)
library(ggplot2)
library(dplyr)

download.file('http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/diab_trans.data', destfile = 'data/diab_trans.data')
diab.df <- read.csv("data/diab_trans.data", header=TRUE, stringsAsFactors = FALSE)

head(diab.df)

diab.df$code_id = as.numeric(substring(diab.df$code, 4,5))
description.df <- read.csv('data/description.txt', header=TRUE)

input.df <- inner_join(diab.df,description.df)
head(select(input.df, code_id, description_PL))


to_chunk.df = input.df %>% 
  group_by(code_id) %>%
  distinct(value) %>%
  summarise(count=n()) %>%
  filter(count>1)
to_chunk.df$measure = TRUE 

base.df <- left_join(input.df, to_chunk.df)
dawka = c(33, 34, 35)
base.df <- base.df %>% mutate(type = ifelse(measure==T, ifelse(code_id %in% dawka, 'dawka', 'pomiar'), 'wydarzenie'))
to_chunk_temp.df = inner_join(description.df, to_chunk.df)
as.vector(to_chunk_temp.df$description_PL)



to_hist_pomiar.df <- base.df %>% na.omit() %>% filter(measure==TRUE, type=='pomiar')
ggplot(to_hist_pomiar.df, aes(x=value, fill=description_PL))+geom_histogram()
to_hist_dawka.df <- base.df %>% na.omit() %>% filter(measure==TRUE, type=='dawka')
ggplot(to_hist_dawka.df, aes(x=value, fill=description_PL))+geom_histogram(binwidth = 1)+scale_x_continuous(limits = c(2,40))

to_hist.df <- base.df %>% na.omit() %>% filter(measure==TRUE, code_id==62)

q <- quantile(to_hist.df$value, c(0.20, 0.80))

ggplot(to_hist.df, aes(x=value))+geom_histogram()+geom_vline(xintercept = c(q[1], q[2]))

to_divide.df <- base.df %>% na.omit() %>% filter(measure==TRUE)
for_loop <- to_divide.df %>% distinct(code_id)
qa = c()
qb = c()
code = c()

for(i in for_loop$code_id){
  to_divide_tmp.df <- base.df %>% na.omit() %>% filter(measure==TRUE, code_id==i)
  qa = c(qa,  quantile(to_divide_tmp.df$value, c(0.20)))
  qb = c(qb,  quantile(to_divide_tmp.df$value, c(0.80)))
  code = c(code, i)
  
}
df <- data.frame(code,qa,qb)

cluster <- function(code_id, value){

      qa=df[df$code==code_id,"qa"]
      qb=df[df$code==code_id,"qb"]
      if(length(qa)==0 || length(qb)==0 )
        return(NA)

      if(is.na(value))
        return(NA)
      if(is.null(value))
        return(NA)
      
      if(value<qa){
        return('niski')
        }
      else if(value>qb){
        return('wysoki')
        }
      else if(value>=qa & value<=qb){
        return('normalny')
        }
      else{
        return(NA)}
}

base.df$value.level <- mapply(cluster, base.df$code_id, base.df$value)
df <-base.df %>%
  mutate(description = if_else(!is.na(value.level),paste(description_PL,"poziom", value.level), as.character(description_PL)))  %>%
  select(patient_id, time_sek, description)
new_id <- df %>%arrange(description) %>% distinct(description)
new_id$id <- seq.int(nrow(new_id))
df = inner_join(df, new_id)
length(unique(df$patient_id))
View(df)

#zapisanie data.frame do pliku aby wygodnie wczyta³o siê g jako sewkencjê
write.table(df[,c(1,2,3)], file = "data/clean_data.csv", row.names=FALSE,col.names=FALSE, sep=",")

#wczytanie z pliku w postaci sekwencyjnej 
dseq <- read_baskets(con = "data/clean_data.csv", sep =",", info = c("sequenceID","eventID"))
summary(dseq)
#test, czy wczytanie jest poprawne
frameS =   as(dseq,"data.frame")
View(frameS)

#liczba ró¿nych zdarzeñ
nitems(dseq)

#czêstoœæ zdarzeñ 
freqItem = itemFrequency(dseq)
freqItem = sort(freqItem, decreasing = TRUE )
freqItem
head(freqItem,30)
#wsparcie objawów hipoglikemii to ~0,11 dlatego zaczniemy rozwa¿ania od wsparcia = 0.1

#ustawienie parametrów
#parametry: 
#support=0.05: wsparcie ustawiamy doœæ nisko, ¿eby "z³apaæ" tak¿e regu³y zwi¹zane z wydarzeniami, które wydarzaj¹ siê rzadko
#maxsize = 1 wszystkie "koszyki" to jednorazowe zdarzenia, dlatego mo¿na zastosowaæ maxsize=1
#maxgap = 28800 w tym wypadku interesuj¹ nas zdarzenia, które wystêpuj¹ w zakresie 8 godzin
#maxlen = 4 zdarzenia w obrêbie max 8 godzin daj¹ ³¹czny zakres objawów w ramach doby, co wydaje siê rozs¹dn¹ granic¹ dla hipoglikemii 
#maxwin - w algorytmie cspade pakietu arulesSequences ma status "disabled", dlatego nie mo¿emy go u¿yæ.
parametry = new ("SPparameter",support = 0.05, maxsize = 1, maxgap = 28800, maxlen=4)
wzorce= cspade(dseq, parametry, control = list(verbose = TRUE, tidLists = FALSE, summary= TRUE))
#odkrycie regu³
reguly = ruleInduction(wzorce,confidence = 0.8)

#podsumowanie, przyk³ad
length(reguly)
inspect(head(reguly,30))

#wyszukanie regu³ interesuj¹cych - zdarzeñ 
hipoglik=subset(reguly, !(lhs(reguly) %in% c("\"Objawy hipoglikemii\"")) & rhs(reguly) %in% c("\"Objawy hipoglikemii\""))
hipoglik = sort(hipoglik, by = "lift", decreasing=TRUE)
inspect(hipoglik)

#uzyskane regu³y dla parametrów: support = 0.05, maxsize = 1, maxgap = 28800, maxlen = 4 )
'
 1 <{"Dawka insuliny regularnej poziom wysoki"},                                                                                                       
     {"Pomiar zawartoœci glukozy we krwi przed obiadem poziom niski"},                                                                                  
{"Pomiar zawartoœci glukozy we krwi przed kolacj¹ poziom niski"}>     => <{"Objawy hipoglikemii"}>                                              0.06060606  1.0000000 1.783784 
2 <{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Dawka insuliny NPH poziom normalny"},                                                                                                            
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.09090909  1.0000000 1.783784 
3 <{"Pomiar zawartoœci glukozy we krwi po kolacji poziom normalny"},                                                                                  
{"Dawka insuliny NPH poziom normalny"},                                                                                                            
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.07575758  1.0000000 1.783784 
4 <{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Dawka insuliny regularnej poziom niski"},                                                                                                        
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.07575758  1.0000000 1.783784 
5 <{"Pomiar zawartoœci glukozy we krwi po kolacji poziom normalny"},                                                                                  
{"Dawka insuliny regularnej poziom niski"},                                                                                                        
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.09090909  1.0000000 1.783784 
6 <{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Dawka insuliny regularnej poziom normalny"},                                                                                                     
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.07575758  1.0000000 1.783784 
7 <{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Pomiar zawartoœci glukozy we krwi po kolacji poziom normalny"},                                                                                  
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.06060606  1.0000000 1.783784 
8 <{"Dawka insuliny regularnej poziom normalny"},                                                                                                     
{"Pomiar zawartoœci glukozy we krwi po kolacji poziom normalny"},                                                                                  
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.07575758  1.0000000 1.783784 
9 <{"Pomiar zawartoœci glukozy we krwi przed kolacj¹ poziom normalny"},                                                                               
{"Pomiar zawartoœci glukozy we krwi po kolacji poziom normalny"},                                                                                  
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.06060606  1.0000000 1.783784 
10 <{"Pomiar zawartoœci glukozy we krwi po kolacji poziom normalny"},                                                                                  
{"Nieokreœlony pomiar zawartoœci glukozy we krwi poziom niski"}>      => <{"Objawy hipoglikemii"}>                                              0.09090909  0.8571429 1.528958 
11 <{"Pomiar zawartoœci glukozy we krwi przed kolacj¹ poziom normalny"},                                                                               
{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Pomiar zawartoœci glukozy we krwi przed przek¹sk¹ poziom niski"}>   => <{"Objawy hipoglikemii"}>                                              0.06060606  0.8000000 1.427027 
12 <{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Pomiar zawartoœci glukozy we krwi przed obiadem poziom niski"}>     => <{"Objawy hipoglikemii"}>                                              0.06060606  0.8000000 1.427027 
13 <{"Dawka insuliny regularnej poziom normalny"},                                                                                                     
{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Pomiar zawartoœci glukozy we krwi przed obiadem poziom niski"}>     => <{"Objawy hipoglikemii"}>                                              0.06060606  0.8000000 1.427027 
14 <{"Pomiar zawartoœci glukozy we krwi przed œniadaniem poziom wysoki"},                                                                              
{"Dawka insuliny NPH poziom niski"},                                                                                                               
{"Pomiar zawartoœci glukozy we krwi przed obiadem poziom niski"}>     => <{"Objawy hipoglikemii"}>                                              0.06060606  0.8000000 1.427027 
15 <{"Spo¿ycie posi³ku wiêkszego ni¿ zazwyczaj"},                                                                                                      
{"Pomiar zawartoœci glukozy we krwi przed obiadem poziom niski"},                                                                                  
{"Pomiar zawartoœci glukozy we krwi przed kolacj¹ poziom niski"}>     => <{"Objawy hipoglikemii"}>                                              0.06060606  0.8000000 1.427027 
''