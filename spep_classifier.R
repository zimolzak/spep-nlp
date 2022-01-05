library(tm)
library(Matrix)
library(e1071)
library(ROCR)
library(stringr)
library(ggplot2)

source("spep_functions.R")
options(stringsAsFactors = FALSE)




#### Read data and preprocess

Xa = read.delim("spep_6sta3n_fix_anno3.txt")
cat("\n")
cat('Distribution of annotations:\n')
table(Xa$monoclonal_01, exclude=NULL)
cat('\n')

Xa$labpanelcomment = gsub("_", " ", Xa$labpanelcomment)

Xa$labpanelcomment = gsub(".", "", Xa$labpanelcomment, fixed = TRUE)
Xa$labpanelcomment = gsub(",", "", Xa$labpanelcomment, fixed=TRUE)
Xa$labpanelcomment = gsub(":", "", Xa$labpanelcomment, fixed=TRUE)
Xa$labpanelcomment = gsub(";", "", Xa$labpanelcomment, fixed=TRUE)
Xa$labpanelcomment = gsub("(", "", Xa$labpanelcomment, fixed=TRUE)
Xa$labpanelcomment = gsub(")", "", Xa$labpanelcomment, fixed=TRUE)
Xa$labpanelcomment = gsub("PEP-", "", Xa$labpanelcomment, fixed=TRUE)
Xa$labpanelcomment = gsub("immunofix-", "", Xa$labpanelcomment, 
                          ignore.case=TRUE)

skey = data.frame(sta3n=c(523, 534, 583, 621, 662, 692), 
                  name=c("Boston", "Charleston", "Indianapolis", 
                         "Mt Home TN", "San Francisco", "White City OR"))




#### Simple modeling (train on a station, test on other SPEPs from same station).

message("Simple modeling")
par(mfrow=c(1,1))
all_good_stations = as.numeric(names(table(Xa$sta3n))[names(table(Xa$sta3n)) != '692'])
# Full list of stations: 523 534 583 621 662 692
# 692 really no data

# TECHNICALLY shouldn't need to exclude NA if we have annotated our whole CSV but you never know
Xaa = Xa[! is.na(Xa$monoclonal_01), ] 

for(this_station in all_good_stations) {
	  Whole_sdtm = char2sparse(Xaa$labpanelcomment)
        X_sdtm = Whole_sdtm[Xaa$sta3n == this_station, ]  ## needs the right num of columns
        X_annotations = Xaa$monoclonal_01[Xaa$sta3n == this_station]
	  message(paste("  ", this_station))
        M = xy2plot(X_sdtm, X_annotations, 
                paste("Station", this_station, 
                      skey$name[skey$sta3n == this_station] ))
	  line_of_code = paste("model_", this_station, " = M", sep='') ## inelegant, side effects!
	  eval(parse(text = line_of_code))
        readline("waiting\n")
}

# All stations.
shuffled = Xaa[sample(nrow(Xaa)),] # Kind of important
X_sdtm = char2sparse(shuffled$labpanelcomment)
X_annotations = shuffled$monoclonal_01
message("  all")

ytrain_temp = X_annotations[1:(length(X_annotations)/2)]
ytest_temp = X_annotations[(length(X_annotations)/2 + 1):(length(X_annotations))]
cat('Distribution of anno TRAINING SET:\n')
table(ytrain_temp)
cat('Distribution of anno TEST SET:\n')
table(ytest_temp)
cat('\n')

model_all = xy2plot(X_sdtm, X_annotations, "All stations", do_roc=TRUE)
# AUC on global (all stations) is 0.96
cat("\n")




#### Which single term has best separation? (variable importance)

message("Variable importance")
term_performance = data.frame(term=NULL, p=NULL, odds_ratio=NULL, a=NULL, b=NULL, c=NULL, d=NULL)
for (j in 1:dim(X_sdtm)[2]){
	Ta = table(data.frame(xj=(X_sdtm[,j] > 0), y=X_annotations))
	##       x       y
	a = Ta['TRUE',  '1']
	b = Ta['TRUE',  '0']
	c = Ta['FALSE', '1']
	d = Ta['FALSE', '0']
	f = fisher.test(Ta)
	tm = colnames(X_sdtm)[j]
	p = f$p.value
	o = as.numeric(f$estimate)
	term_performance = rbind(term_performance, data.frame(term=tm, p=p, odds_ratio=o, a=a, b=b, c=c, d=d))
}

bonf_threshold = 0.05 / dim(term_performance)[1]

for (thresh in c(bonf_threshold, 1e-10)) {
	cat("P <", thresh, "\n")
	sig = term_performance[term_performance$p < thresh,]
	sig$normal = sig$odds_ratio < 1
	print(sig[order(sig$odds_ratio), ])
	cat("\n")
}

# dx plot, histo of 400 ish term p-vals
qplot(term_performance$p) + geom_vline(xintercept = 0.05) + geom_vline(xintercept = bonf_threshold) + scale_x_log10()

# PUBLICATION FIGURE (histogram of terms odds ratios)
ggplot(aes(x=odds_ratio, fill=(p < bonf_threshold)), data=term_performance) + 
	geom_histogram(show.legend = FALSE) +
	scale_x_log10() +	
	scale_fill_manual(values=c('#838383', 'black')) +
#                                FALSE   TRUE
	geom_vline(xintercept = 1) +
	xlab('Odds ratio') + ylab('Count') +
	annotate('text', label = 'Suggestive of\nmonoclonal', x=10, y=60) + 
	annotate('text', label = "Suggestive of\nnot monoclonal", x=0.1, y=60) +
	theme(axis.title = element_text(size=12))

# axis.text.x
# geom_text given to annotate() understands size=.. as an aesthetic.
# 80 mm * 80 mm, 300 dpi. ==> 945 * 945

ggsave("feat_hist_fig3_8080.tiff", width = 80, height = 80, units="mm") # 2100 x 2096 default

ggsave("feat_hist_fig3_scale05.tiff", scale=0.5)

# jamia seems to be "medium format" page size.
# therefore single column 80mm wide, dbl col 160mm



#### Test models based on complete graph on 5 stations

message("complete graph")
cat('Matrix of optimum ROCC points:\n')
AUC_matrix = data.frame(AUC=NULL, s1=NULL, s2=NULL)
for (i in all_good_stations) {
        message(paste("  ", i))
        for (j in all_good_stations) {
                auc = station_compare(Xaa, i, j)
                AUC_matrix = rbind(AUC_matrix, data.frame(AUC=auc, s1=i, s2=j))
        }
}




#### Make GraphViz directed graph from ROC cutoff (& histogram)

# dx plot, histo of 25 AUCs
qplot(data=AUC_matrix, x=AUC, xlab="Area under ROC curve", 
      main="Performances of 5 classifiers, each tested on all 5 stations")

aucm_selected =  AUC_matrix[AUC_matrix$AUC > 0.8,] ### magic number

dot_lines = 'digraph {'
for (i in 1:length(aucm_selected[,'s1'])){
        src_n = aucm_selected[i,'s1']
        dest_n = aucm_selected[i,'s2']
        src = skey$name[skey$sta3n == src_n]
        dest = skey$name[skey$sta3n == dest_n]
        src = gsub(" ", "_", src)
        dest = gsub(" ", "_", dest)
        if(src != dest){
                dot_lines = c(dot_lines, (paste(src, '->', dest)))
        }
}
dot_lines = c(dot_lines, '}')

fileConn = file("roc_graph.dot")
writeLines(dot_lines, fileConn)
close(fileConn)




#### Heat map of the 5x5 classifier applications

join2 = sta3n_to_city(AUC_matrix, skey)
cat("\nMatrix of AUCs:\n")
print(join2)
cat("\n")

# PUBLICATION FIGURE (heat map AUC)
p = ggplot(join2, aes(x=station1, y=station2, label=round(AUC, 2)))
p + geom_tile(aes(fill=AUC)) + labs(x="Station for training", 
                                    y ="Station for testing") +
					geom_label(fontface='bold') +
					scale_fill_gradient(low='black', high='white') +
 					 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




##### Rules based

message("Rules based")

predictions = c()
for (i in 1:dim(Xaa)[1]) {
  predictions = c(predictions, 
                  rules_predict(Xaa$labpanelcomment[i], 
                                Xaa$sta3n[i])
                  )
}

Comparison = data.frame(y = Xaa$monoclonal_01, yhat = predictions)
cat("Rules-based confusion matrix:\n")
table(Comparison)
cat("\n")
cat("Rules-based accuracy (all stations, proper application):\n")
dim(Comparison[(Comparison$y==1 & Comparison$yhat==1) | 
               (Comparison$y==0 & Comparison$yhat==0), 
               ])[1] / dim(Comparison)[1]
cat("\n")

# Rules, "big model"
# Copy pasted, DRY don't repeat yourself, yeah, yeah....
predictions = c()
for (i in 1:dim(Xaa)[1]) {
  # note different function rules_big_predict() which takes just 1 arg.
  predictions = c(predictions, rules_big_predict(Xaa$labpanelcomment[i]))
}
Comparison = data.frame(y = Xaa$monoclonal_01, yhat = predictions)
cat("Rules-based BIG MODEL confusion matrix:\n")
table(Comparison)
cat("\n")
cat("Rules-based BIG MODEL accuracy (all stations, NO HINTS):\n")
dim(Comparison[(Comparison$y==1 & Comparison$yhat==1) | 
               (Comparison$y==0 & Comparison$yhat==0), 
               ])[1] / dim(Comparison)[1]
cat("\n")
FalsePos = Comparison[(Comparison$y==0 & Comparison$yhat==1), ]
FalsePosTxt = Xaa$labpanelcomment[(Comparison$y==0 & Comparison$yhat==1)]
cat("Just confirm that this looks like sub-cell of the confus matrix:\n")
table(FalsePos)
cat("Here are the offending reports:\n")
print(FalsePosTxt)
cat("\n")


#### Apply rules to *everybody*. (Histogram, then heat map)

Acc_matrix = data.frame(Accuracy=NULL, s1=NULL, s2=NULL)
for (s1 in all_good_stations){
        for (s2 in all_good_stations){
                predictions = c()
                reports_to_test = Xaa$labpanelcomment[Xaa$sta3n == s2]
                for (r in reports_to_test) {
                  predictions = c(predictions, rules_predict(r, s1))
                }
                Comparison = data.frame(y = Xaa$monoclonal_01[Xaa$sta3n == s2],
                                        yhat = predictions)
                A = dim(Comparison[(Comparison$y==1 & Comparison$yhat==1) | 
                               (Comparison$y==0 & Comparison$yhat==0), 
                               ])[1] / dim(Comparison)[1]
                cat("Rules-based accuracy", s1, s2, ":", A, "\n")
                Acc_matrix = rbind(Acc_matrix, data.frame(Accuracy=A, s1=s1, s2=s2))
        }
}

# dx plot, histo of 25 rules accuracies.
qplot(data=Acc_matrix, x=Accuracy, 
      main="Performances of 5 rule systems, each tested on all 5 stations")

Acc_matrix_names = sta3n_to_city(Acc_matrix, skey)

# PUBLICATION FIGURE (heat map of 25 rules Accs)
p = ggplot(Acc_matrix_names, aes(x=station1, y=station2, label=round(Accuracy, 2)))
p + geom_tile(aes(fill=Accuracy)) + labs(x="Station for training", 
                                         y ="Station for testing")+ 
						geom_label(fontface='bold') + 
						scale_fill_gradient(low='black', high='white')+
 					 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())






#### Length analysis

Xa = merge(Xa, skey)   ## come back and fixme

qplot(1:dim(Xa)[1], str_length(Xa$labpanelcomment), color=Xa$name, log="y",
      xlab="SPEP number", ylab="Report length (characters)") + labs(color='Station')

A = Xa[! is.na(Xa$monoclonal_01), ]

qplot(1:dim(A)[1], str_length(A$labpanelcomment), log="y",
      color=as.factor(A$monoclonal_01), # shape=as.factor(A$sta3n), 
      xlab="SPEP number", ylab="Report length (characters)") + labs(color='Monoclonal')
