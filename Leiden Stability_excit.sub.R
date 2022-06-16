# Leiden Cluster Stability

library(scclusteval)
library(future.apply)


plan(multiprocess)
results = future.apply::future_sapply(c(2.2,2.4,2.6,2.8,3,3.2,3.4), simplify = F,
                                      USE.NAMES = T, FUN = function(res) {

res = 2.8
boot.cl = lapply(boot.cl.ex, function(x) x[paste0('cl.',res)])

stab = list()
for(i in 1:30){
  x.int = intersect(Cells(ex.sub),names(boot.cl[[i]][[paste0('cl.',res)]]))
  cells2use = ex.sub[[paste0('pca.snn_res.',res)]][x.int,]
  names(cells2use) = x.int
  x= PairWiseJaccardSets(cells2use, sort(boot.cl[[i]][[paste0('cl.',res)]])[x.int])
  tmp = list(mat = x, max = sort(apply(x,MARGIN = 1, max), decreasing = T), cells = x.int)
  stab[[i]] = tmp
}


x = lapply(stab, function(x) apply(x$`mat`, MARGIN = 1, function(x) sort(x,decreasing = T)[1]))
x = do.call(rbind,x)
y = apply(x, 2, mean)
barplot(sort(y))

df = x %>% reshape2::melt()
df$Var2 = factor(df$Var2)
ggplot(df, aes(x=Var2, y=value, color = Var2)) + geom_violin() + geom_jitter() + NoLegend()

# how many are below .4 average Jaccard?

length(which(y<.4)) # 7
bad = names(which(y<.4))
table(ex[[paste0('pca.snn_res.',res)]])[bad]

tmp = list(x = x, y = sort(y), bad = bad)

results[[as.character(res)]] = tmp
}
)


AvHeat(as.data.frame(results$`2.8`$x), fill.val = 'expression', gene.order = gtools::mixedsort(colnames(results$`2.8`$x)))

Idents(ex.sub) = 'pca.snn_res.2.8'

bad = names(which(results$`2.8`$y<.4))
results$`2.8`$y[bad]
results$`2.8`$y




table(ex.sub$pca.snn_res.2.8)
table(data.frame(rep = ex.sub$rep_name, cl = ex.sub$pca.snn_res.2.8))

esub.av = data.frame(t(ex.sub@assays$MERFISH@scale.data), cl = ex.sub$pca.snn_res.2.8) %>% group_by(cl) %>% summarise_all(mean) %>% 
  tibble::column_to_rownames('cl') %>% t()


# dropping c( 69, 68, 65, 59, 54, 46, 55)


VlnPlot(ex.sub, pt.size = 0, features = c('nFeature_MERFISH','DAPI','mRNA_counts','blank.sum'), ncol = 1)



# make parallel
library(doParallel)
library(parallel)
library(caret)
cl = makeCluster(16)
clusterEvalQ(cl, .libPaths("C:/Users/Eric Vaughn/Documents/R/win-library/3.6"))
registerDoParallel(cl, cores = 16)

Idents(ex.sub) = 'pca.snn_res.2.8'

bad = c( 69, 68, 65, 54, 46, 55)

set.seed(1)
mat = ex.sub@assays$MERFISH@scale.data
table(ex.sub$pca.snn_res.2.8)
cells = WhichCells(ex.sub, downsample = 200, idents = bad, invert = T)
train.ind = sample(cells, length(cells)*.8)
train = t(mat[, train.ind])
test = t(mat[, setdiff(cells, train.ind)])
test.class = ex.sub@active.ident[rownames(test)]
train.class = ex.sub@active.ident[train.ind]
train.class = data.frame(train.class)
colnames(train.class) = 'cl'
train.class$cl = factor(train.class$cl)
train.class$cl = droplevels(train.class$cl)

# train control
TrainCtrl1 <- trainControl(method = "repeatedcv", number = 5, verboseIter = T)

# svm radial --------------------------------------------------------------
ex.sub.svm <- train(x = train, y = train.class$cl, method="svmRadial", trControl=TrainCtrl1,
                   tuneLength = 10, verbose=T)
# check accuracy of classifier
pred.class = predict(ex.sub.svm , test)
svm.conf = confusionMatrix(pred.class, test.class)
svm.conf


# looks good, predict the remaining unstable cells
bad.cells = WhichCells(ex.sub, idents = bad)
pred.unstable = predict(ex.sub.svm, t(mat[,bad.cells]))
table(pred.unstable)

ex.sub$stable = ex.sub$pca.snn_res.2.8
ex.sub$stable[bad.cells] = pred.unstable

DimPlot(ex.sub, group.by = 'stable', cells = sample(Cells(ex.sub),30000), label = T) + NoLegend()
DimPlot(ex.sub, group.by = 'stable', cells = sample(Cells(ex.sub),30000), cells.highlight = bad.cells) + NoLegend()
DimPlot(ex.sub, cells = sample(Cells(ex.sub),30000), label = T, group.by = 'pca.snn_res.3', 
        cells.highlight = WhichCells(ex.sub, idents = bad)) + NoLegend()


save(ex.sub, file = 'ex.sub.RObj')


# check cluster 1 for more diversity

Idents(ex.sub) = 'stable'
c1 = subset(ex.sub, idents = 1)
c1 = RunPCA(c1, npcs = 20, features = n.genes)
c1 = FindNeighbors(c1, dims = 1:20, k.param = 15)
plan(sequential)
c1 = FindClusters(c1,algorithm = 4, n.start = 10, n.iter = 10, resolution = .4)
plan(multiprocess)
c1 = RunUMAP(c1, dims = 1:20, n.neighbors = 15, min.dist=0.1)
DimPlot(c1, label = T) + NoLegend()

c1.markers = data.frame(t(c1@assays$MERFISH@scale.data), cl = c1$MERFISH_snn_res.0.4) %>% 
  group_by(cl) %>% summarise_all(mean) %>% tibble::column_to_rownames('cl') %>% t()

# stability for c1
set.seed(11)
cells = lapply(1:50, function(x) sample(Cells(c1), ncol(c1)*0.8))
boot.cl.c1 = list()
for (i in 1:50){
  message(paste0('iteration ',i))
  plan(sequential)
  b = subset(c1, cells = cells[[i]]) # sample 1 mil cells
  b = RunPCA(b, npcs = 20, verbose = F, features = n.genes)
  b = FindNeighbors(b, k.param = 15, verbose = T, dims = 1:20, 
                    reduction = 'pca', graph.name = 'pca.snn', 
                    nn.method = 'annoy', annoy.metric = 'cosine')
  plan(multiprocess, workers = 5)
  b = FindClusters(object = b, modularity.fxn = 1, resolution = c(.1,.2,.3,.4,.5), 
                   algorithm = 4, n.start = 20, graph.name = 'pca.snn', method = 'igraph',
                   group.singletons = T)
  boot.cl.c1[[i]] = list(cl.0.1 = b$pca.snn_res.0.1,
                         cl.0.2 = b$pca.snn_res.0.2,
                         cl.0.3 = b$pca.snn_res.0.3,
                         cl.0.4 = b$pca.snn_res.0.4,
                         cl.0.5 = b$pca.snn_res.0.5,
                         cells = Cells(b))
}


stab.c1 = list()
res = 0.4
for(i in 1:50){
  x.int = intersect(Cells(c1),names(boot.cl.c1[[i]][[paste0('cl.',0.5)]]))
  cells2use = c1[[paste0('MERFISH_snn_res.',res)]][x.int,]
  names(cells2use) = x.int
  x= PairWiseJaccardSets(cells2use, sort(boot.cl.c1[[i]][[paste0('cl.',0.5)]])[x.int])
  tmp = list(mat = x, max = sort(apply(x,MARGIN = 1, max), decreasing = T), cells = x.int)
  stab.c1[[i]] = tmp
}


x = lapply(stab.c1, function(x) apply(x$`mat`, MARGIN = 1, function(x) sort(x,decreasing = T)[1]))
x = do.call(rbind,x)
y = apply(x, 2, mean)
barplot(sort(y))
y

df = x %>% reshape2::melt()
df$Var2 = factor(df$Var2)
ggplot(df, aes(x=Var2, y=value, color = Var2)) + geom_violin() + geom_jitter() + 
  xlab('Ex.sub Cluster') + ylab('Jaccard Index') + NoLegend()


# c1 bad cell classifier --------------------------------------------------

# classify remaining
cl = makeCluster(16)
clusterEvalQ(cl, .libPaths("C:/Users/Eric Vaughn/Documents/R/win-library/3.6"))
registerDoParallel(cl, cores = 16)

Idents(c1) = 'MERFISH_snn_res.0.4'

bad = c(8,9)

set.seed(1)
mat = c1@assays$MERFISH@scale.data
table(c1$MERFISH_snn_res.0.4)
cells = WhichCells(c1, downsample = 200, idents = bad, invert = T)
train.ind = sample(cells, length(cells)*.8)
train = t(mat[, train.ind])
test = t(mat[, setdiff(cells, train.ind)])
test.class = c1@active.ident[rownames(test)]
train.class = c1@active.ident[train.ind]
train.class = data.frame(train.class)
colnames(train.class) = 'cl'
train.class$cl = factor(train.class$cl)
train.class$cl = droplevels(train.class$cl)

# train control
TrainCtrl1 <- trainControl(method = "repeatedcv", number = 5, verboseIter = T)

c1.svm <- train(x = train, y = train.class$cl, method="svmRadial", trControl=TrainCtrl1,
                    tuneLength = 10, verbose=T)
# check accuracy of classifier
pred.class = predict(c1.svm , test)
svm.conf = confusionMatrix(pred.class, test.class)
svm.conf # accuracy of .90


bad.cells = WhichCells(c1, idents = bad)
pred.unstable = predict(c1.svm, t(mat[,bad.cells]))
table(pred.unstable)

c1$stable = c1$MERFISH_snn_res.0.4
c1$stable[bad.cells] = pred.unstable

DimPlot(c1, group.by = 'stable', cells = sample(Cells(c1),3000), label = T) + NoLegend()



# Merge c1 with ex.sub ----------------------------------------------------

ex.sub$stable2 = as.numeric(ex.sub$stable) + 6
ex.sub$stable2[Cells(c1)] = c1$stable
ex.sub$stable = factor(ex.sub$stable2)

