# Disscetor function

# Author: Stephanie Jurburg


#### Dissector function ####
#' Compare distances between groups of samples
#'
#' Creates a 6 column data frame from a distance matrix. The first two columns correspond to the two samples between which the comparison is being made. The third is the distance value, the fourth and fifth columns are the variable assigned to each of the two samples, and these can be different factors. The sixth column contains "total" for all rows, and is meant to be used as a boxplot for reference. The recommended method for trimming and plotting this dataframe d is as follows:
#' 
#' @param physeq a phyloseq object
#' @param distm a distance metric. Inherits from \code{\link{vegdist}}. Options are "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
#' @param var1 grouping variable, must match one of \code{sample_variable(physeq)}. 
#' @param var2 second grouping variable, must match one of \code{sample_variable(physeq)}, and usually var1=var2. var1!=var2 for very specific cases, for example when comparing the distances between mother-offspring couples and the variables coding for the relationship are named differently for mothers and offspring. 
#'
#' @return 6-column dataframe. The first two columns correspond to the two samples between which the comparison is being made. The third is the distance value, the fourth and fifth columns are the variable assigned to each of the two samples, and these can be different factors. The sixth column contains "total" for all rows, and is meant to be used as a separate boxplot for reference (a null hypothesis). 
#' @note the result must be further trimmed according to user preferences. Further processing according to example below will produce a boxplot and an ANOVA.
#' @export
#' @import vegan
#' @import ggplot2
#' @import phyloseq
#' @author Stephanie Jurburg (Dr. Carrot) \email{s.d.jurburg@gmail.com}
#'
#' @examples
#' d=dissector(physeq, distm, var1, var2)
#' e=d[d$var2 == d$var2, ]
#' #Plot as follows
#' ggplot()+geom_boxplot(data=e, aes(x=var2, y=value))+geom_boxplot(data=d, aes(x=Total, y=value))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
#' anova(lm(value~var2, e))

fun_dissector <- function(physeq, distm, var1, var2) {
  # Args:
  # - physeq: phyloseq object
  # - distm: a dist class object
  # - var1: a category within which to compare samples. For example, if var1= "Time", the distance between samples from the samples of the same time point will be measured.
  if (physeq@otu_table@taxa_are_rows=="FALSE"){
    dist.mat=vegdist((otu_table(physeq)), method="bray")
  } else{
    dist.mat=vegdist(t(otu_table(physeq)), method="bray")
  }
  A <- attr(dist.mat, "Size")
  B <- if (is.null(attr(dist.mat, "Labels"))) sequence(A) else attr(dist.mat, "Labels")
  if (isTRUE(attr(dist.mat, "Diag"))) attr(dist.mat, "Diag") <- FALSE
  if (isTRUE(attr(dist.mat, "Upper"))) attr(dist.mat, "Upper") <- FALSE
  d=data.frame(
    var1.names = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    var2.names = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(dist.mat))
  v.1=as.data.frame(sample_data(physeq)[[var1]], row.names(sample_data(physeq)))
  colnames(v.1)=var1
  v.2=as.data.frame(sample_data(physeq)[[var2]], row.names(sample_data(physeq)))
  d[,var1]=v.1[,var1][match(d$var1.names, row.names(v.1))]
  if (var1==var2){
    colnames(v.2)=paste(var2,".")
    d[,paste(var2,".")]=v.2[,paste(var2,".")][match(d$var2.names, row.names(v.2))]
  } else{
    colnames(v.2)=var2
    d[,var2]=v.2[,var2][match(d$var2.names, row.names(v.2))]
  }
  d$Total="total"
  return(d)
}