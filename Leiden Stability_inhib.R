# Leiden Cluster Stability

# For excitatory and inhibitory cells, cells were clustered with k.param = 11, 
# resolution = to 2.6. For excitatory, very distinct cell types like cholinergic,
# serotoninergic, and dopaminergic cells dominated the variance. Hence, we subset the larger bulk of
# non-super discrete excitatory clusters to resolve this population further. For this subset, k=11 and
# resolution = 2.8.
# To search for stable clusters, we performed a bootstrapped clustering procedure 30 times with a
# random subset of 80% of the cells. For each resulting clustering, we calculated the Pairwise Jaccard
# coefficient between the full and bootstrapped clustering. Clusters with a mean Jaccard coefficient
# of less than .40 were classified as unstable. Cells deemed unstable were then classified into a remaining
# stable cluster using a radial support vector machine classifier that was trained on 
# 200 cells from each stable cluster with 5-fold repeated cross validation using the caret package.
# For inhib and ex.sub the classifier accuracy was approximately 90% when tested on the 
# test set of untrained cells.
# For excitatory cluster 1, we noticed a large amount of diversity in this large cluster and again,
# subdivided this cluster by reclustering with 20pcs, k=15, and resolution =0.4. We repeated the bootstrap
# procedure mentioned above 50 times to define 7 stable clusters, each with discrete markers and spatial
# profiles.


library(scclusteval)

setwd('C:/Users/Eric Vaughn/Dropbox/Sequencing/PAG/MERFISH/December 2020/')
load('C:/Users/Eric Vaughn/Dropbox/Sequencing/PAG/MERFISH/June 2020/inhibitory.RObj')


load('Leiden_Bootstrap_Stability_list_seed11.complete.R')
boot.cl = lapply(boot.cl.11, function(x) x['cl.2.6'])

stab = list()
for(i in 1:30){
  x.int = intersect(Cells(inhib),names(boot.cl[[i]]$cl.2.6))
  x= PairWiseJaccardSets(inhib$pca.snn_res.2.6[x.int],sort(boot.cl[[i]]$cl.2.6)[x.int])
  tmp = list(mat = x, max = sort(apply(x,MARGIN = 1, max), decreasing = T), cells = x.int)
  stab[[i]] = tmp
}


# retrieve the average max jaccard coefficient for each cell (max of each cluster it belongs to) across
# all repeated measures.

x = lapply(stab, function(x) apply(x$`mat`, MARGIN = 1, function(x) sort(x,decreasing = T)[1]))
x = do.call(rbind,x)
y = apply(x, 2, mean) %>% sort()
barplot(sort(y))

# how many are below .5 average Jaccard?

length(which(y<.4)) # 7
bad = names(which(y<.4))
table(inhib$pca.snn_res.2.6)[bad]

i.av = data.frame(t(inhib@assays$MERFISH@scale.data), cl = inhib$pca.snn_res.2.6) %>% group_by(cl) %>% summarise_all(mean) %>% 
  tibble::column_to_rownames('cl') %>% t()

sapply(bad, function(x) topz(x, i.av, n=10), USE.NAMES = T, simplify = F)



df = x %>% reshape2::melt()
df$Var2 = factor(df$Var2)
ggplot(df, aes(x=Var2, y=value, color = Var2)) + geom_violin() + geom_jitter() + 
  xlab('Inhib Cluster') + ylab('Jaccard Index') + NoLegend()

ggsave('Inhib Stability.png')


# Classify Unstable Cells -------------------------------------------------


# make parallel
library(doParallel)
library(parallel)
library(caret)
cl = makeCluster(16)
clusterEvalQ(cl, .libPaths("C:/Users/Eric Vaughn/Documents/R/win-library/3.6"))
registerDoParallel(cl, cores = 16)

Idents(inhib) = 'pca.snn_res.2.6'

set.seed(1)
mat = inhib@assays$MERFISH@scale.data
table(inhib$pca.snn_res.2.6)
cells = WhichCells(inhib, downsample = 200, idents = bad, invert = T)
train.ind = sample(cells, length(cells)*.8)
train = t(mat[, train.ind])
test = t(mat[, setdiff(cells, train.ind)])
test.class = inhib@active.ident[rownames(test)]
train.class = inhib@active.ident[train.ind]
train.class = data.frame(train.class)
colnames(train.class) = 'cl'
train.class$cl = factor(train.class$cl)
train.class$cl = droplevels(train.class$cl)

# train control
TrainCtrl1 <- trainControl(method = "repeatedcv", number = 5, verboseIter = T)

# svm radial --------------------------------------------------------------
inhib.svm <- train(x = train, y = train.class$cl, method="svmRadial", trControl=TrainCtrl1,
                    tuneLength = 10, verbose=T)
# check accuracy of classifier
pred.class = predict(inhib.svm , test)
svm.conf = confusionMatrix(pred.class, test.class)
svm.conf # 0.89 accuracy of model

save(inhib.svm, file = 'inhib.svm.RObj')

bad.cells = WhichCells(inhib, idents = bad)
pred.unstable = predict(inhib.svm, t(mat[,bad.cells]))
table(pred.unstable)
save(bad.cells, file = 'inhib.bad.cells.R')

inhib$stable = inhib$pca.snn_res.2.6
inhib$stable[bad.cells] = pred.unstable

DimPlot(inhib, group.by = 'stable', cells = sample(Cells(inhib),30000), label = T) + NoLegend()
DimPlot(inhib, group.by = 'stable', cells = sample(Cells(inhib),30000), cells.highlight = bad.cells) + NoLegend()
DimPlot(inhib, cells = sample(Cells(inhib),30000), label = T, group.by = 'pca.snn_res.2.6', 
        cells.highlight = WhichCells(inhib, idents = bad)) + NoLegend()


save(inhib, file = 'inhib.RObj')

