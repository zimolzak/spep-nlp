char2sparse = function(charvector) {
        MyCorpus <- VCorpus(VectorSource(charvector),
                readerControl = list(language = "en"))
        MyCorpus <- tm_map(MyCorpus, content_transformer(tolower))
        DTM <- DocumentTermMatrix(MyCorpus, 
                                  control = list(bounds = list(global = c(0, Inf)))) 
        sparse_DTM <- sparseMatrix(i = DTM$i, j = DTM$j, x = DTM$v,
                                   dims = dim(DTM),
                                   dimnames = list(rownames(DTM), colnames(DTM)))
        return(sparse_DTM)
}

xy2plot = function(sparse_DTM, annotations, title, do_roc = FALSE) {
	  stopifnot(dim(sparse_DTM)[1] == length(annotations))

        # Don't need to reshuffle because I already shuffled to decide which to annotate.
        xtrain = sparse_DTM[1:(dim(sparse_DTM)[1]/2), ]
        ytrain = annotations[1:(length(annotations)/2)]
        xtest = sparse_DTM[(dim(sparse_DTM)[1]/2 + 1) : (dim(sparse_DTM)[1]), ]
        ytest = annotations[(length(annotations)/2 + 1):(length(annotations))]

        # svm.model <- svm(xtrain, ytrain, cost = 10, gamma = 1)
	  ST = tune.svm(xtrain, ytrain, gamma=10^(-4:2), cost=10^(0:4))
	  svm.model = ST$best.model
	  g = as.numeric(ST$best.parameters['gamma'])
	  c = as.numeric(ST$best.parameters['cost'])
	  p = ggplot(ST$performances, aes(x=gamma, y=cost))
	  p2 = p + geom_tile(aes(fill=log(error))) + geom_vline(xintercept=g) + geom_hline(yintercept=c) + scale_x_log10() + scale_y_log10() + labs(title = title)
	  plot(p2)

        preds = predict(svm.model, xtest)
        rocr_pred = prediction(preds, ytest)
        perf <- performance(rocr_pred, "tpr", "fpr")
        auc_object = performance(rocr_pred, "auc")
        A = as.numeric(slot(auc_object, 'y.values'))

        mymain = paste(title, ": AUC", as.character(round(A, 3)))
	  if(do_roc){
	        plot(perf, colorize=T, lwd= 3, main=mymain)
	  }
        cat(paste(title, ": AUC =", A, "\n"))

	  perf2 =  performance(rocr_pred, 'acc', 'fnr')
	  perf3 = performance(rocr_pred, 'tnr')
	  roc_df = data.frame(	fpr = slot(perf, 'x.values')[[1]], 
					tpr = slot(perf, 'y.values')[[1]],
#					fnr = slot(perf2, 'x.values')[[1]], 
					acc = slot(perf2, 'y.values')[[1]], 
					sp = slot(perf3, 'y.values')[[1]], # tnr
					a = slot(rocr_pred, 'tp')[[1]],
					b = slot(rocr_pred, 'fp')[[1]],
					c = slot(rocr_pred, 'fn')[[1]],
					d = slot(rocr_pred, 'tn')[[1]],
					cutoff = slot(perf, 'alpha.values')[[1]]
				)
	  roc_df$dist = sqrt(roc_df$fpr ^ 2 + (1 - roc_df$tpr) ^ 2)
	  print(roc_df[which.min(roc_df$dist),])
	  return(svm.model)
}

station_compare = function(Xa, s1, s2) {
        ### Note! Assumes names of columns of Xa!
	  ### also going to assume that Xa is fully annotated, no more !is.na() subsetting

	  # Set up test data for station 2
        Whole_sdtm = char2sparse(Xa$labpanelcomment)
        xtest = Whole_sdtm[Xa$sta3n == s2, ]
        ytest = Xa$monoclonal_01[Xa$sta3n == s2]
	  stopifnot(dim(xtest)[1] == length(ytest))

	  # Retrieve our model for station 1
	  line_of_code = paste("svm.model = model_", s1, sep='')
	  eval(parse(text = line_of_code))

	  # Apply model to test data.
        preds = predict(svm.model, xtest)
        rocr_pred = prediction(preds, ytest)


        perf <- performance(rocr_pred, "tpr", "fpr")
	  perf2 =  performance(rocr_pred, 'acc', 'fnr')
	  perf3 = performance(rocr_pred, 'tnr')
	  roc_df = data.frame(	fpr = slot(perf, 'x.values')[[1]], 
					tpr = slot(perf, 'y.values')[[1]],
#					fnr = slot(perf2, 'x.values')[[1]], 
					acc = slot(perf2, 'y.values')[[1]], 
					sp = slot(perf3, 'y.values')[[1]], # tnr
					a = slot(rocr_pred, 'tp')[[1]],
					b = slot(rocr_pred, 'fp')[[1]],
					c = slot(rocr_pred, 'fn')[[1]],
					d = slot(rocr_pred, 'tn')[[1]],
					cutoff = slot(perf, 'alpha.values')[[1]]
				)
	  roc_df$dist = sqrt(roc_df$fpr ^ 2 + (1 - roc_df$tpr) ^ 2)
	  cat('#### ', s1, '->', s2, '\n')
	  print(roc_df[which.min(roc_df$dist),])

        auc_object = performance(rocr_pred, "auc")
        return(as.numeric(slot(auc_object, 'y.values')))
}

sta3n_to_city = function (matrix, key) {
        ## Generate a data frame w/ numeric stations decoded to strings (station1 station2).
        ## matrix is a data.frame with numeric cols s1 and s2.
        ## key is a data.frame that should have numeric col 'sta3n' and string col 'name'.

        join1 =  merge(x = matrix, y = key, by.x = 's1', by.y = 'sta3n')
        join1$station1 = join1$name
        join1$name = NULL

        join2 =  merge(x = join1, y = key, by.x = 's2', by.y = 'sta3n')
        join2$station2 = join2$name
        join2$name = NULL

        return(join2)
}

rules_predict = function (x, sta3n) {
  if (sta3n == 523 & (
          grepl("identified in gamma region", x, ignore.case=TRUE) |
          grepl("monoclonal", x, ignore.case=TRUE) | 
          grepl("SER PARAPROTEIN Ig", x, ignore.case=TRUE)
          )
      ) {
    return (1) 
  }

  if (sta3n == 534 & (
          grepl("abnormal monoclonal", x, ignore.case=TRUE) |
          grepl("M Spike is too small", x, ignore.case=TRUE) |
          grepl("may be due to a monoclonal gammopathy", x, ignore.case=TRUE) |
          grepl("present in the gamma globulin", x, ignore.case=TRUE) |
          grepl("an abnormal monoclonal band", x, ignore.case=TRUE) |
          grepl("monoclonal gammopathy", x, ignore.case=TRUE) |
          grepl("evaluation reveals a restricted", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  if (sta3n == 583 & (
          grepl("abnormal monoclonal band", x, ignore.case=TRUE) |
          grepl("M-spike", x, ignore.case=TRUE) |
          grepl("M spike", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  if (sta3n == 621 & (
          grepl("a monoclonal band", x, ignore.case=TRUE) |
          grepl("faint monoclonal band", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  if (sta3n == 662 & (
          grepl("monoclonal pattern", x, ignore.case=TRUE) |
          grepl("monclonal pattern", x, ignore.case=TRUE) |  ## catching a typo!
          grepl("faint band detected", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  return(0)
}




## ok this is copy pasted

rules_big_predict = function (x) {
  sta3n = 0 # this is a hack so I don't have to delete a bunch of 'if' statements. I'm in a hurry.
  if (sta3n == 0 & (
          grepl("identified in gamma region", x, ignore.case=TRUE) |
          grepl("monoclonal", x, ignore.case=TRUE) | 
          grepl("SER PARAPROTEIN Ig", x, ignore.case=TRUE)
          )
      ) {
    return (1) 
  }

  if (sta3n == 0 & (
          grepl("abnormal monoclonal", x, ignore.case=TRUE) |
          grepl("M Spike is too small", x, ignore.case=TRUE) |
          grepl("may be due to a monoclonal gammopathy", x, ignore.case=TRUE) |
          grepl("present in the gamma globulin", x, ignore.case=TRUE) |
          grepl("an abnormal monoclonal band", x, ignore.case=TRUE) |
          grepl("monoclonal gammopathy", x, ignore.case=TRUE) |
          grepl("evaluation reveals a restricted", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  if (sta3n == 0 & (
          grepl("abnormal monoclonal band", x, ignore.case=TRUE) |
          grepl("M-spike", x, ignore.case=TRUE) |
          grepl("M spike", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  if (sta3n == 0 & (
          grepl("a monoclonal band", x, ignore.case=TRUE) |
          grepl("faint monoclonal band", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  if (sta3n == 0 & (
          grepl("monoclonal pattern", x, ignore.case=TRUE) |
          grepl("monclonal pattern", x, ignore.case=TRUE) |  ## catching a typo!
          grepl("faint band detected", x, ignore.case=TRUE) 
          )
      ) {
    return (1) 
  }
  
  return(0)
}