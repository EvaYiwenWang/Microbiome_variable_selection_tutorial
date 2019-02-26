#######################################################################################################
#######################################################################################################
#
# Read data "Crohn.rda"
#
#######################################################################################################
#######################################################################################################


load("Crohn_data.rda")


#######################################################################################################
#
# Comparison of the 12 variables selected by the 3 methods: selbal, clr+logistic_lasso, coda_logistic_lasso
#
#######################################################################################################


#######################################################################################################
#######################################################################################################
#
# CODA LOGISTIC LASSO 
#
#######################################################################################################
#######################################################################################################

source("functions_coda_penalized_regression.R")

# Main function:
#   coda_logistic_lasso(y,x,lambda)

# y is the binary outcome, can be numerical (values 0 and 1), factor (2 levels) or categorical (2 categories)
#
# x is the matrix of microbiome abundances, either absolute abundances (counts) or relative abundances (proportions)
#   the rows of x are individuals/samples, the columns are taxa

# Imputation of zeros:
# The user can provide x without zeros (after implementation of an external zero imputation method) and no additional transformation will be applied 
# If x contains zeros, 

# Log-transformation: 
# x should not be the matrix of log(counts) or log(proportions). The method itself performs the log-transformation of the abundances.

# lambda is the penalization parameter: the larger the value of lambda the fewer number of variables will be selected.
#
# function rangLambda(y,x,numVar, lambdaIni=1) provides a rang of lambda values corresponding to a given number of variables to be selected (numVar). The default initial lambda is lambdaIni=1.  


y<-y_Crohn 

x<-x_Crohn 

rangLambda(y,x,numVar=12, lambdaIni =0.15)
  
results_codalasso<-coda_logistic_lasso(y,x,lambda=0.19)
results_codalasso
selected_codalasso<-results_codalasso[[4]][-1]

columns_selected_codalasso<-results_codalasso[[3]][-1]

#coef<-results[[5]][-1]

write.csv(data.frame(columns_selected_codalasso,selected_codalasso),"results_codalasso_Crohn12.csv")



#######################################################################################################
#######################################################################################################
#
# CLR + LOGISTIC LASSO
#
#######################################################################################################
#######################################################################################################

y<-y_Crohn_numeric  # glmnet() requires y to be numeric

x<-as.matrix(x_Crohn) 

# CLR transformation Z=log(x)
z<-log(x)
clrx<- apply(z,2,function (x) x-rowMeans(z))
#rowMeans(clrx)  

library(glmnet)


clrlasso <- glmnet(clrx, y, standardize=FALSE , alpha=1,family="binomial", lambda=seq(0.015,2,0.01)) 
plot(clrlasso)
print(clrlasso)
clrlasso_coef<-coef(clrlasso,s=0.10)
sum(abs(clrlasso_coef)>0)
selected_clrlasso<-which(as.numeric(abs(coef(clrlasso, s=0.1)))>0)[-1];  # indices of selected variables (abs(coef)>0)
selected_clrlasso<-selected_clrlasso-1
taxa_id<-colnames(x_Crohn)[selected_clrlasso]

write.csv(data.frame(selected_clrlasso,taxa_id),"results_clrlasso_Crohn12.csv")

#######################################################################################################
#######################################################################################################
#
# selbal: Selection of balances
#
#######################################################################################################
#######################################################################################################

y<-y_Crohn  # for binary outcomes (logistic regression) selbal requires y to be factor, 
#                             if y is numeric selbal implements linear regression

x<-x_Crohn 
colnames(x)<-(1:ncol(x))

selbal_Crohn<-selbal(x = x, y = y, logt=T, maxV=12)

selected_selbal<-as.numeric(c(selbal_Crohn[[6]][,1],selbal_Crohn[[6]][,2]))
id.na<-which(is.na(selected_selbal))
selected_selbal<-selected_selbal[-id.na]
selected_selbal<-as.character(selected_selbal)

columns_selected_selbal<-which(colnames(x)%in% selected_selbal)
taxa_id<-colnames(x_Crohn)[columns_selected_selbal]

write.csv(data.frame(columns_selected_selbal, taxa_id),"results_selbal_Crohn12.csv")



#######################################################################################################
#######################################################################################################
#
# Concordance of variables selected by the three methods
#
#######################################################################################################
#######################################################################################################

require("VennDiagram")
require("gplots")



d_selbal<-read.csv("results_selbal_Crohn12.csv", header = T)
d_clrlasso<-read.csv("results_clrlasso_Crohn12.csv", header = T)
d_codalasso<-read.csv("results_codalasso_Crohn12.csv", header = T)

#(d_selbal[,2], d_clrlasso[,2], d_codalasso[,2])
  
taxa_selected<-list(d_clrlasso[,2], d_codalasso[,2], d_selbal[,2])
taxa.id_selected<-list(d_clrlasso[,3], d_codalasso[,3], d_selbal[,3])


venn.plot <- venn.diagram(taxa_selected , NULL, fill=c("magenta", "blue", "lightblue"), 
                          alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("clr+lasso", "coda-lasso", "selbal"),
                          main="Concordance of selected taxa for Crohn data")
grid.draw(venn.plot)

taxa.id <- venn(taxa.id_selected, show.plot=FALSE)
taxa <- venn(taxa_selected, show.plot=FALSE)

inters_taxa.id <- attr(taxa.id,"intersections")
inters_taxa <- attr(taxa,"intersections")

lapply(inters_taxa, head) 

lapply(inters_taxa.id, head) 

