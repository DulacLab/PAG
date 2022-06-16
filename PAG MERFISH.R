library(Seurat)
library(dplyr)
library(purrr)
library(ggplot2)
library(plotly)
library(future)
library(future.apply)
library(gridExtra)



# Load data ---------------------------------------------------------------

db = readxl::read_xlsx('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/Animal_database.xlsx')[1:50,]
db = db[!is.na(db$`Analysis done`),] # remove NAs
db = db[!is.na(db$`Date sliced`),] # remove NAs
db = db[is.na(db$`QC metrics`),] # remove dropped datasets; anything that's not NA is 'Dropped'

pag = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/Location_annotated_data.csv', header = T)
pag$rep_name = plyr::mapvalues(pag$Sample, 
                               from = unique(pag$Sample),
                               to = c('Agg_M','Inf_M','Nai_M','Vir_F','Nai_F','Mat_M',
                                      'Mot_F','Fea','Fat_M','Nai_M','Mat_F')
)
pag$rep_name = as.character(pag$rep_name)
pag$rep_name[pag$rep_name=='Fea'] = paste0('Fea_', pag$Sex[pag$rep_name=='Fea'])
pag$rep_numb = plyr::mapvalues(pag$Animal_id, 
                               from = as.numeric(db$Analysis_animal_id),
                               to = as.numeric(db$`Animal Number`))
pag$rep_name = paste0(pag$rep_name,pag$rep_numb)
table(pag$rep_name)


# individually add --------------------------------------------------------


agg.m3.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190821_PAG3_2A.csv', header = T)
agg.m3.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190823_PAG3_2P.csv', header = T)
agg.m3 = rbind(agg.m3.a, agg.m3.p)
agg.m3$rep_name = 'Agg_M3'
agg.m2.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190828_PAG3_1A.csv', header = T)
agg.m2.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190826_PAG3_1P.csv', header = T)
agg.m2 = rbind(agg.m2.a,agg.m2.p)
agg.m2$rep_name = 'Agg_M2'
agg.m4 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191024_PAG3_2.csv', header = T)
agg.m4$rep_name = 'Agg_M4'
agg = list(agg.m3, agg.m2, agg.m4)
conserved.cols = intersect(colnames(agg.m2), colnames(agg.m3))
agg = lapply(agg, function(x) x[,conserved.cols]) # subset to conserved rows (m1 doesn't have mRNA_counts)
agg = do.call("rbind", agg) # bind into one matrix

# naive males
nai.m1.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/190903_PAG3_1A.csv', header = T)
nai.m1.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/190905_PAG3_1P.csv', header = T)
nai.m1 = rbind(nai.m1.a, nai.m1.p)
nai.m1$rep_name = 'Nai_M1'
nai.m3.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191009_PAG3_3A.csv', header = T)
nai.m3.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191011_PAG3_3P.csv', header = T)
nai.m3 = rbind(nai.m3.a, nai.m3.p)
nai.m3$rep_name = 'Nai_M3'
nai.m2.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/191014_PAG3_2A.csv', header = T)
nai.m2.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/191017_PAG3_1P.csv', header = T)
nai.m2 = rbind(nai.m2.a, nai.m2.p)
nai.m2$rep_name = 'Nai_M2'
nai.m5 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/191127_PAG3_1.csv', header = T)
nai.m5$rep_name = 'Nai_M5'
nai.m6 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200129_PAG3_2.csv', header = T)
nai.m6$rep_name = 'Nai_M6'
nai.m7 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200216_PAG3_2.csv', header = T)
nai.m7$rep_name = 'Nai_M7'
# naive females
nai.f1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191107_PAG3_3.csv', header = T)
nai.f1$rep_name = 'Nai_F1'
nai.f2 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191113_PAG3_1.csv', header = T)
nai.f2$rep_name = 'Nai_F2'
nai.f3 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191125_PAG3_2.csv', header = T)
nai.f3$rep_name = 'Nai_F3'
nai.f4 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191203_PAG3_2.csv', header = T)
nai.f4$rep_name = 'Nai_F4'
nai = list(nai.m1, nai.m3, nai.m2, nai.m5, nai.m6,nai.m7, nai.f1, nai.f2, nai.f3, nai.f4)
nai = lapply(nai, function(x) x[,conserved.cols]) # subset to conserved rows 
nai = do.call("rbind", nai)

inf.m1.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190911_PAG3_2A.csv', header = T)
inf.m1.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190909_PAG3_1P.csv', header = T)
inf.m1 = rbind(inf.m1.a, inf.m1.p)
inf.m1$rep_name = 'Inf_M1'
inf.m2.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/190916_PAG3_2A.csv', header = T)
inf.m2.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/190913_PAG3_1P.csv', header = T)
inf.m2 = rbind(inf.m2.a, inf.m2.p)
inf.m2$rep_name = 'Inf_M2'
inf.m3.a = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190920_PAG3_3A.csv', header = T)
inf.m3.p = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/190923_PAG3_2P.csv', header = T)
inf.m3 = rbind(inf.m3.a, inf.m3.p)
inf.m3$rep_name = 'Inf_M3'
inf = list(inf.m3, inf.m2, inf.m1)
inf = lapply(inf, function(x) x[,conserved.cols]) # subset to conserved rows (m1 doesn't have mRNA_counts)
inf = do.call("rbind", inf)

vir.f4 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191104_PAG3_1.csv', header = T)
vir.f4$rep_name = 'Vir_F4'
vir.f5 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/191111_PAG3_1.csv', header = T)
vir.f5$rep_name = 'Vir_F5'
vir.f1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200214_PAG3_2.csv', header = T)
vir.f1$rep_name = 'Vir_F1'
vir = list(vir.f1, vir.f4, vir.f5)
vir = lapply(vir, function(x) x[,conserved.cols]) # subset to conserved rows (m1 doesn't have mRNA_counts)
vir = do.call("rbind", vir)

# parenting
#mother
mot1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191216_PAG3_1.csv', header = T)
mot1$rep_name = 'Mot_F1'
mot1 = mot1[,conserved.cols]
mot2 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200120_PAG3_1.csv', header = T)
mot2$rep_name = 'Mot_F2'
mot2 = mot2[,conserved.cols]
mot3 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200228_PAG3_1.csv', header = T)
mot3$rep_name = 'Mot_F3'
mot3 = mot3[,conserved.cols]
#father
fat1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/Outdated_data_files/191213_PAG3_1.csv', header = T)
fat1$rep_name = 'Fat_M1'
fat1 = fat1[,conserved.cols]
fat2 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200122_PAG3_2.csv', header = T)
fat2$rep_name = 'Fat_M2'
fat2 = fat2[,conserved.cols]
fat3 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200226_PAG3_1.csv', header = T)
fat3$rep_name = 'Fat_M3'
fat3 = fat3[,conserved.cols]
mot = rbind(mot1,mot2,mot3)
fat = rbind(fat1,fat2,fat3)


# mating
#male
mat1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191121_PAG3_3.csv', header = T)
mat1$rep_name = 'Mat_M1'
mat2 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191119_PAG3_2.csv', header = T)
mat2$rep_name = 'Mat_M2'

mat.f1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200316_PAG3_1.csv', header = T)
mat.f1$rep_name = 'Mat_F1'
mat.f2 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200318_PAG3_2.csv', header = T)
mat.f2$rep_name = 'Mat_F2'
mat = list(mat2, mat1, mat.f1, mat.f2)
mat = lapply(mat, function(x) x[,conserved.cols]) # subset to conserved rows (m1 doesn't have mRNA_counts)
mat = do.call("rbind", mat)

# fear
fea.m1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/191218_PAG3_1.csv', header = T)
fea.m1$rep_name = 'Fea_M1'
fea.m2 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200206_PAG3_2.csv', header = T)
fea.m2$rep_name = 'Fea_M2'
fea.f1 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200204_PAG3_1.csv', header = T)
fea.f1$rep_name = 'Fea_F1'
fea.f3 = read.csv('C:/Users/Eric Vaughn/Dropbox/PAG_MERFISH_data/raw_data/200310_PAG3_2.csv', header = T)
fea.f3$rep_name = 'Fea_F3'

fea = list(fea.m1, fea.m2, fea.f1, fea.f3)
fea = lapply(fea, function(x) x[,conserved.cols]) # subset to conserved rows (m1 doesn't have mRNA_counts)
fea = do.call("rbind", fea)


cise = read.csv('A:/Dropbox/Sequencing/PAG/MERFISH/Exercise/210106_PAG3_sample6.csv', header = T)
cise.meta = cise[,1:10]
rownames(cise.meta) = cise.meta$X
cise = cise[,-(1:10)] %>% t()
colnames(cise) = cise.meta$X
cise = CreateSeuratObject(cise, project = 'Exercise', assay = 'MERFISH', 
                          meta.data = cise.meta)



pag = list(agg, nai, inf, vir, mot, fat, mat, fea)

pag = do.call(rbind, pag)




# Filter dataset ----------------------------------------------------------

multi.genes = colnames(pag)[16:273]
ieg.genes = colnames(pag)[301:326]
seq.genes = c('Snap25', 'Hepacam', 'Ucn','Cartpt')
blanks = grep('blank',colnames(pag), value = T)
genes = c(multi.genes, seq.genes, ieg.genes)


# remove cells that have 0 counts for multiplexed genes
good = Matrix::rowSums(pag[,multi.genes])>0

# save meta data
pag.meta = as.data.frame(pag[good,c(colnames(pag)[1:15],'rep_name', 'polyT', 'DAPI')])
blank.sum = rowSums(pag[good,blanks])*10 # multiplied by 10 for convenience of getting out of decimals
pag.meta$blank.sum = blank.sum
rownames(pag.meta) = 1:nrow(pag.meta)

pag = pag[good, genes]
rownames(pag) = 1:nrow(pag)
pag = t(pag) # reduce to just expression data and transpose
all(!is.na(pag)) # should be true if there are no NAs


# Create Seurat object for all cells --------------------------------------
mer = CreateSeuratObject(pag, project = 'pag',assay = 'MERFISH', meta.data = pag.meta)

Idents(mer) = 'rep_name'

# remove low quality cells
mer = subset(mer, subset = mRNA_counts< 1800 & mRNA_counts >20)

# QC
# VlnPlot(mer, 'median_total_density', pt.size = 0)
# VlnPlot(mer, 'mRNA_counts', pt.size = 0)
# VlnPlot(mer, 'Hepacam', pt.size = 0)
# VlnPlot(mer, 'Snap25', pt.size = 0)


# this takes a very long time to run for all cells, so use a reduced pca space to filter glia for now
glia.genes = c('Vtn','Cldn5','Slc1a3','Agt','Arhgap45','Lum','Mog','Pdgfra','Hepacam','Hdc','X1500015O10Rik','Mrc1')
neuron.genes = c('Slc17a6','Gad2','Snap25','Slc32a1','Slc17a7','Slc17a8','Cacna2d1','Cacna2d3','Ntng1', 'Tph2','Th','Chat','Tac1','Penk','Sst')
mer$glia = PercentageFeatureSet(mer, features = setdiff(glia.genes,'Hepacam'))
# remove most glia to ease load on classifier
mer = ScaleData(mer, features = c(multi.genes, seq.genes), do.center = T, scale.max = 10, vars.to.regress = 'rep_name')

mer = subset(mer, glia <10)

# classify non neurons

mer$isNeuron = predict(class.svm, t(mer@assays$MERFISH@scale.data))



# Add other identities ----------------------------------------------------

mer$beh = plyr::mapvalues(mer$rep_name, 
                          from = unique(c(grep('Agg',mer$rep_name, value =T),
                                          grep('Vir',mer$rep_name, value =T),
                                          grep('Nai',mer$rep_name, value =T),
                                          grep('Inf',mer$rep_name, value =T),
                                          grep('Mot',mer$rep_name, value =T),
                                          grep('Fat',mer$rep_name, value =T),
                                          grep('Mat',mer$rep_name, value =T),
                                          grep('Fea',mer$rep_name, value =T))),
                          to = c(rep('Agg', length(unique(grep('Agg',mer$rep_name, value =T)))),
                                 rep('Vir', length(unique(grep('Vir',mer$rep_name, value =T)))),
                                 rep('Nai', length(unique(grep('Nai',mer$rep_name, value =T)))),
                                 rep('Inf', length(unique(grep('Inf',mer$rep_name, value =T)))),
                                 rep('Mot', length(unique(grep('Mot',mer$rep_name, value =T)))),
                                 rep('Fat', length(unique(grep('Fat',mer$rep_name, value =T)))),
                                 rep('Mat', length(unique(grep('Mat',mer$rep_name, value =T)))),
                                 rep('Fea', length(unique(grep('Fea',mer$rep_name, value =T)))))
)

# create set to distinguish each imaging session
half.sets = data.frame(rep = mer$rep_name, slice = mer$Brain_pos) %>%
  group_by(rep) %>% summarise(set = n_distinct(slice)) %>% filter(set<20) %>% select(rep) %>% unlist()
Idents(mer) = 'rep_name'
mer$Brain_pos[WhichCells(mer, idents = half.sets)] = mer$Brain_pos[WhichCells(mer, idents = half.sets)]*2
mer$set = if_else(mer$Brain_pos>12, true = 'P',false = 'A')
mer$set = if_else(mer$rep_name %in% half.sets, true = 'all',false = mer$set)
mer$set = factor(paste(mer$set, mer$rep_name,sep = '_'))



# Nonneural --------------------------------------------------------------------


load('cldn5.cells.R')
load('low.quality.nn.R')

Idents(mer) = 'isNeuron'
nn = subset(mer, idents = 'nonneuronal')
nn = subset(nn, subset = Ucn < 2, slot = 'scale.data')# Ucn and Cartpt cells made it passed the classifier and need to be moved to neuronal
nn = subset(nn, subset = Cartpt<2.5, slot = 'scale.data')
nn = RunPCA(nn, npcs = 20, features =  c(multi.genes,seq.genes))
plan(multiprocess) # use all 16 cores
nn = RunUMAP(object = nn, reduction = "pca", dims = 1:20, 
             n.components = 2L, n.neighbors = 15, min.dist = .05,
             learning.rate = 1, angular.rp.forest = T, metric = 'cosine',
             uwot.sgd = T, n.epochs = 1000, reduction.name = 'umap')
DimPlot(nn, group.by = 'rep_name', cells = sample(Cells(nn),30000)) + NoLegend()
plan(sequential)
nn = FindNeighbors(nn, k.param = 15, verbose = T, dims = 1:20, 
                   reduction = 'pca', graph.name = 'pca.snn',
                   nn.method = 'annoy', annoy.metric = 'cosine')
plan(multiprocess)
nn = FindClusters(object = nn, modularity.fxn = 1, resolution = c(1.2), 
                  algorithm = 3, n.start = 10, graph.name = 'pca.snn')
DimPlot(nn, cells = sample(Cells(nn),30000), label = T, group.by = 'pca.snn_res.1.2') + NoLegend()

save(nn, file = 'nonneuronal.RObj')

nn.av = data.frame(t(nn@assays$MERFISH@scale.data), cl = nn$pca.snn_res.2)
nn.av = nn.av %>% group_by(cl) %>% summarise_all(.funs = mean) %>% tibble::column_to_rownames('cl') %>% t()

# rep breakdown
prop = table(data.frame(nn$rep_name, nn$pca.snn_res.1.2))

# cerebellum cells
Idents(nn) = 'pca.snn_res.2'
to.inhib = WhichCells(nn, idents = 35)
save(to.inhib, file = 'to.inhib.R')

poor = WhichCells(nn, idents = c(9,8,26,24,29,21,19)) # Either low z-score or doublet mixes of Endo/Astro/Oligo mixed



# repeat with poor quality cells removed# Meis1 and Foxo1 appear to be artificially creating new clusters 
# in a way that seems artifactual

nn = subset(nn, cells = c(to.inhib, poor), invert = T)
nn = RunPCA(nn, npcs = 15, features =  setdiff(c(multi.genes,seq.genes),c('Foxo1','Meis1')))
plan(multiprocess) # use all 16 cores
nn = RunUMAP(object = nn, reduction = "pca", dims = 1:15, 
             n.components = 2L, n.neighbors = 15, min.dist = .05,
             learning.rate = 1, angular.rp.forest = T, metric = 'cosine',
             uwot.sgd = T, n.epochs = 1000, reduction.name = 'umap')
DimPlot(nn, group.by = 'rep_name', cells = sample(Cells(nn),30000)) + NoLegend()
plan(sequential)
nn = FindNeighbors(nn, k.param = 15, verbose = T, dims = 1:15, 
                   reduction = 'pca', graph.name = 'pca.snn',
                   nn.method = 'annoy', annoy.metric = 'cosine')
plan(sequential)
nn = FindClusters(object = nn, modularity.fxn = 1, resolution = c(1.2), 
                  algorithm = 1, n.start = 10, graph.name = 'pca.snn')
DimPlot(nn, cells = sample(Cells(nn),30000), label = T, group.by = 'pca.snn_res.1.2') + NoLegend()

nn.av = data.frame(t(nn@assays$MERFISH@scale.data), cl = nn$pca.snn_res.1.2)
nn.av = nn.av %>% group_by(cl) %>% summarise_all(.funs = mean) %>% tibble::column_to_rownames('cl') %>% t()


# Cldn5 = Endothelial
# Vtn = Pericyte, Vasc Smooth Muscle, or FB like
# Pdgfra = Oligo Precursor Cell (OPC) or FB-like
# Mecom = Endothelial or ependymal
# Lef1 = high in Endothelial, low in pericyte
# Sema3c = Endothelial
# Arhgef17 = Pericyte, SMC, or FB-like
# Pde3a = Pericyte or SMC (mural)

types = list(
  Endo = 3,
  Mural = 8,
  Microglia = 12,
  Astrocyte = c(13,0,5),
  Ependymal = 7,
  Radial_Glia_Like = 15,
  OPC = c(6,21),
  Immature_Oligo = 14,
  Fibroblast_like = 17,
  Mature_Oligo = c(2,11,4,1,9,10)
)


nn = subset(nn, idents = c(16,18,19,20,22,23,24,25,26), invert = T) # 16 are poorly defined neurons; throwing out

type.names = unlist(lapply(names(types), function(x) rep(x, length(types[[x]]))))

nn$types = plyr::mapvalues(nn$pca.snn_res.1.2, from = unlist(types), to = type.names)

nn$types = droplevels(nn$types)

DimPlot(nn, cells = sample(Cells(nn),30000), label = T, group.by = 'types') + NoLegend()

nn$stable = nn$types # rename to stable so merging with neurons below is simpler

save(nn, file = 'nonneuronal.RObj')


# Neuron Analysis ---------------------------------------------------------


Idents(mer) = 'isNeuron'
mer.n = subset(mer, idents = 'neuron')
mer.n = merge(mer.n, subset(mer, cells = to.inhib))

# add back Ucn/Cartpt cells that were lost in nonneuronal set
Idents(nn) = 'pca.snn_res.1'
ucn = subset(nn, idents = 20)

mer.n = merge(mer.n, ucn)

mer.n = subset(mer.n, cells = cldn5, invert = T) # remove these cells that made it through

VlnPlot(mer.n, 'mRNA_counts', group.by = 'rep_name', pt.size = 0)

ncol(mer.n)
mer.n = subset(mer.n, subset = nFeature_MERFISH > 50)
mer.n = subset(mer.n, subset = mRNA_counts > 100)
mer.n = subset(mer.n, subset = mRNA_counts < 1500)
mer.n = subset(mer.n, subset = DAPI > 2.38) # 2sd below mean DAPI signal
ncol(mer.n)

plan(multiprocess) # use all 16 cores
mer.n = ScaleData(mer.n, features = c(n.genes,'Cartpt','Ucn'), scale.max = 100,
                  vars.to.regress = c('rep_name','mRNA_counts'), model.use = 'linear',
                  use.umi = F)
mer.n = RunPCA(mer.n, npcs = 30, 
               features =  c(n.genes,'Cartpt','Ucn'))
mer.n = RunUMAP(object = mer.n, reduction = "pca", dims = 1:30, 
                n.components = 2L, n.neighbors = 15, min.dist = .05,
                learning.rate = 1, angular.rp.forest = T, metric = 'cosine',
                uwot.sgd = T, n.epochs = 500, reduction.name = 'umap')
DimPlot(mer.n, group.by = 'rep_name', cells = sample(Cells(mer.n),30000)) + NoLegend()
plan(sequential)
mer.n = FindNeighbors(mer.n, k.param = 15, verbose = T, dims = 1:2, 
                      reduction = 'umap', graph.name = 'umap.snn')
mer.n = FindClusters(object = mer.n, modularity.fxn = 1, resolution = c(.1), 
                     algorithm = 1, n.start = 10, graph.name = 'umap.snn')
DimPlot(mer.n, cells = sample(Cells(mer.n),30000), label = T) + NoLegend()

g = c('Th','Slc6a3','Chat','Dbh','Tph2','Slc32a1','Gad2','Slc17a6','Slc17a7','Slc17a8')
av = data.frame(t(mer.n@assays$MERFISH@scale.data[g,]), cl = mer.n$umap.snn_res.0.1)
av = av %>% group_by(cl) %>% summarise_at(.vars = g, .funs = mean)

i.cl = av %>% filter(Slc32a1>.8) %>% select(cl) %>% unlist()
DimPlot(mer.n, cells = sample(Cells(mer.n),30000), 
        cells.highlight = WhichCells(mer.n, idents = i.cl)) + NoLegend()


# inhibitory --------------------------------------------------------------

# low quality cells were defined after one round of analysis and tossed stated below
load(file = 'low.quality.cells.R')

inhib = subset(mer.n, idents = i.cl)

ex.cells = WhichCells(inhib, expression = Slc32a1 < .2) # mean Slc32a1 is 0.67; Gad2 is 1.72
ex.cells = intersect(ex.cells, WhichCells(inhib, expression = Slc17a6 >1.12)) # .78 is 2 sd above mean for vglut2
ex.cells = intersect(ex.cells, WhichCells(inhib, expression = Cacna2d1 >0.58)) # .58 is 2sd above mean for Cacna2d1 selective to glut

# confirm these should be moved to excitatory
FeatureScatter(inhib, cells = ex.cells, feature1 = 'Slc17a6','Slc32a1')
DimPlot(inhib, cells.highlight = ex.cells)

inhib = subset(inhib, cells = c(ex.cells,low.quality), invert = T)
inhib = subset(inhib, cells = intersect(Cells(inhib),Cells(mer.n))) # update to mrna_count regressed 7/19/20
inhib@assays$MERFISH@scale.data = mer.n@assays$MERFISH@scale.data[,Cells(inhib)]

plan(sequential)

# removing genes from PCA that shouldn't be in GABAergic neurons (infering from knowledge and nucSeq specificity)
inhib = RunPCA(inhib, npcs = 50, features =  setdiff(n.genes,c('Chat','Tph2','Slc17a7','Npw')))
plan(multiprocess, workers = 16) # use all 16 cores
genes = setdiff(n.genes, c('Chat','Tph2','Tfap2d','Slc18a2','Slc17a7','Gad2','Slc32a1','Pou4f1'))
inhib = RunUMAP(object = inhib, dims = 1:50, 
                n.components = 2L, n.neighbors = 11, min.dist = .05,
                learning.rate = 1, angular.rp.forest = T, metric = 'cosine',
                uwot.sgd = T, n.epochs = 1000, seed.use = 1)
DimPlot(inhib, group.by = 'rep_name', cells = sample(Cells(inhib),30000)) + NoLegend()
plan(sequential)
inhib = FindNeighbors(inhib, k.param = 11, verbose = T, dims = 1:50, 
                      reduction = 'pca', graph.name = 'umap.snn',
                      nn.method = 'annoy', annoy.metric = 'cosine')
plan(multiprocess)
inhib = FindClusters(object = inhib, modularity.fxn = 1, resolution = seq(2,3.2,.1), 
                     algorithm = 4, n.start = 40, graph.name = 'umap.snn', method = 'igraph',
                     group.singletons = F)
DimPlot(inhib, group.by = 'pca.snn_res.2.6', cells = sample(Cells(inhib),30000), label = T) + NoLegend()

i.av = data.frame(t(inhib@assays$MERFISH@scale.data), cl = inhib$stable)
i.av = i.av %>% group_by(cl) %>% summarise_at(.vars = setdiff(n.genes,c('Cartpt','Ucn')), .funs = mean) %>% 
  tibble::column_to_rownames(var = 'cl') %>% select(setdiff(n.genes,c('Cartpt','Ucn'))) %>% t()


# cluster 31 appears to be a hybrid cluster
FeaturePlot(inhib, 'Slc17a6', min.cutoff = .3, order = T, cells = WhichCells(inhib, idents = 58))
FeaturePlot(inhib, 'Slc32a1', min.cutoff = .3, order = T, cells = WhichCells(inhib, idents = 58))

# cluster 25 has vglut and vgat, but little else to identify (highest z-score is Spon1= .71)
# has very low mRNA_counts and unstable proportions in replicates and its location. Throwing these cells out
low.quality = WhichCells(inhib, idents = 25)

inhib$beh = gsub('[[:digit:]]', '', inhib$rep_name) # label behaviors with sex

x = sort(unique(inhib$stable))
inhib$stable = plyr::mapvalues(inhib$stable, from = x,
                               to = 1:length(x))

inhib$stable = paste0('i.', inhib$stable)

save(inhib, file = 'inhibitory.RObj')
save(low.quality, file = 'low.quality.cells.R')


# Excitatory and others ---------------------------------------------------

load('cldn5.cells.R') # some cells that got through somehow?
ex = subset(mer.n, cells = c(Cells(inhib), low.quality), invert = T) # includes other cell types such as Chat, Dbh, Th, etc.

ex = RunPCA(ex, npcs = 80, features =  c(n.genes, "Cartpt", "Ucn"))
plan(multiprocess) # use all 16 cores
ex = RunUMAP(object = ex, reduction = "pca", dims = 1:60, 
             n.components = 2L, n.neighbors = 11, min.dist = .05,
             learning.rate = 1, angular.rp.forest = T, metric = 'cosine',
             uwot.sgd = T, n.epochs = 1000, reduction.name = 'umap.30pc')
DimPlot(ex, group.by = 'pca.snn_res.3', cells = sample(Cells(ex),30000), reduction = 'umap') + NoLegend()
plan(sequential)
ex = FindNeighbors(ex, k.param = 11, verbose = T, dims = 1:60, 
                   reduction = 'pca', graph.name = 'pca.snn', 
                   nn.method = 'annoy', annoy.metric = 'cosine')
plan(multiprocess)
ex = FindClusters(object = ex, modularity.fxn = 1, resolution = c(2,2.2,2.4,2.6,2.8,3,3.2,3.4), 
                  algorithm = 4, n.start = 40, graph.name = 'pca.snn', method = 'igraph',
                  group.singletons = F)
DimPlot(ex, cells = sample(Cells(ex),30000), label = T, group.by = 'pca.snn_res.3') + NoLegend()

e.av = data.frame(t(ex@assays$MERFISH@scale.data), cl = ex$pca.snn_res.3)

e.av = e.av %>% group_by(cl) %>% summarise_at(.vars = n.genes, .funs = mean) %>% 
  tibble::column_to_rownames('cl') %>% t()

ex$beh = gsub('[[:digit:]]', '', ex$rep_name) # label behaviors with sex

save(ex, file = 'excitatory.RObj')


Idents(ex) = 'pca.snn_res.3'
region = list(
  Cb = WhichCells(ex, idents =73),
  Troc.nerve = WhichCells(ex, idents = 60),
  Trochlear = WhichCells(ex, idents =),
  Thal = WhichCells(ex, idents =70),
  InfCol = WhichCells(ex, idents =c(5)),
  Oculomotor = WhichCells(ex, idents =25),
  D.Raphe = WhichCells(ex, idents = c(45,1)),
  NorAdr = WhichCells(ex, idents = 67),
  Trigem.nerve = WhichCells(ex, idents = 72),
  Dopa = WhichCells(ex, idents = c(58,57,61))
)

# artifactual cells disproportionately found in one/few replicates
# not due to spatial boundaries
poor.ex = WhichCells(ex, idents = c(62,69)) 

# subset to bulk of PAG excitatory without other sources of high diversity
ex.sub = subset(ex, idents = c(62,69,1,45,54,58,25,57,37,67,42,60,73,5,70,11,22,46), invert = T)


ex.sub = RunPCA(ex.sub, npcs = 60, features =  c(n.genes, "Cartpt"))
plan(multiprocess) # use all 16 cores
ex.sub = RunUMAP(object = ex.sub, reduction = "pca", dims = 1:60, 
                 n.components = 2L, n.neighbors = 11, min.dist = .05,
                 learning.rate = 1, angular.rp.forest = T, metric = 'cosine',
                 uwot.sgd = T, n.epochs = 1000, reduction.name = 'umap')
DimPlot(ex.sub, cells = sample(Cells(ex.sub),30000), label = T, group.by = 'pca.snn_res.2.8') + NoLegend()
plan(sequential)
ex.sub = FindNeighbors(ex.sub, k.param = 11, verbose = T, dims = 1:60, 
                       reduction = 'pca', graph.name = 'pca.snn', 
                       nn.method = 'annoy', annoy.metric = 'cosine')


plan(multiprocess)
ex.sub = FindClusters(object = ex.sub, modularity.fxn = 1, resolution = c(2,2.2,2.4,2.6,2.8,3,3.2,3.4), 
                      algorithm = 4, n.start = 40, graph.name = 'pca.snn', method = 'igraph',
                      group.singletons = T)
DimPlot(ex.sub, cells = sample(Cells(ex.sub),30000), label = T, group.by = 'pca.snn_res.3') + NoLegend()

esub.av = data.frame(t(ex.sub@assays$MERFISH@scale.data), cl = ex.sub$pca.snn_res.2.8)

esub.av = esub.av %>% group_by(cl) %>% summarise_at(.vars = n.genes, .funs = mean) %>% 
  tibble::column_to_rownames('cl') %>% t()

save(ex.sub, file = 'ex.sub.RObj')


# Now add ex.sub back to ex

# throw out any missing cluster numbers and rename from 1:n_clusters
ex.sub$stable2 = plyr::mapvalues(ex.sub$stable2, from = sort(unique(ex.sub$stable2)),
                                 to = 1:length(unique(ex.sub$stable2)))

# remove ex.62, it's artifactual to Agg_M3 (in just about 2 slices everything expressed)
Idents(ex) = 'pca.snn_res.3'
ex = subset(ex, idents = 62, invert = T)

# define super discrete clusters
x = c(69,1,45,54,58,25,57,37,67,42,60,73,5,70,11,22,46)

DimPlot(ex, cells = sample(Cells(ex),30000), label = T, group.by = 'pca.snn_res.3',
        cells.highlight = WhichCells(ex, idents = x)) + NoLegend()


# name these from 70 (there's 69 stable clusters in ex.sub)
ex$new = plyr::mapvalues(ex$pca.snn_res.3, 
                         from = x,
                         to = 70:(70+length(x)-1))

ex$new = as.vector(ex$new)

ex$new[Cells(ex.sub)] = as.numeric(as.character(ex.sub$stable2))

ex$new = factor(ex$new, levels = 1:length(unique(ex$new)))

ex$stable = ex$new

DimPlot(ex, group.by = 'stable', cells = sample(Cells(ex),30000), reduction = 'umap', label = T) + NoLegend()

e.av = data.frame(t(ex@assays$MERFISH@scale.data), cl = ex$stable)

e.av = e.av %>% group_by(cl) %>% summarise_at(.vars = n.genes, .funs = mean) %>% 
  tibble::column_to_rownames('cl') %>% t()

ex$stable = paste0('e.', ex$stable)

# remove cells that exist in nonneuronal set as well
dup = intersect(Cells(ex), Cells(nn))
ex = subset(ex, cells = dup, invert = T)

save(ex, file = 'ex.RObj')



# Merge nonneuronal, inhib, and excit together ----------------------------

# excit stable validated at res 3 k =11
# ex.sub stable validated at resolution 2.8 k=11
# inhib stable validated at inhib resolution 2.6 k =11
# nn types merged but originally clustered res 1.2 k =15

# i.56 and e.82 are cerebellum
# e.84 is probably very posterior end of thalamus 


pag = merge(x = nn, y = list(Excitatory =ex, Inhibitory = inhib))

pag$stable = factor(pag$stable)

pag$types[as.character(ex$X)] = 'Excitatory' 

pag$types[as.character(inhib$X)] = 'Inhibitory' 

pag = RenameCells(pag, new.names = pag$X)
  

# IEG Calculation ------------------------------------------------------------

  pag$beh = plyr::mapvalues(pag$beh, from = unique(pag$beh), to = c('Aggression_M','Infanticide_M',
                                                                    'Naive', 'Virgin Parenting_F', 'Mating','Parenting_F','Fear','Parenting_M','Aggression_M',
                                                                    'Infanticide_M','Naive_M','Virgin Parenting_F','Naive_F','Mating_M','Parenting_F',
                                                                    'Fear_M','Parenting_M','Fear_F','Mating_F'))
  g = c('Fos','Fosb','Fosl2','Egr1','Nr4a3','Nr4a1','Junb')
  pag$ieg.score = Matrix::colSums(pag@assays$MERFISH@scale.data[g,])
  
  # remove clusters that aren't at least 1% present in each replicate
  # these tend to be on the edges of the cut region
  
  fringe = table(data.frame(rep = pag$rep_name, cl = pag$stable)) %>% as.data.frame() %>% group_by(cl) %>%
    mutate(frac = Freq/sum(Freq)) %>% filter(frac<.01) %>% tally %>% filter(n>3) %>% dplyr::select(cl) %>% unlist %>% unique %>% as.vector
  
  
  ieg.av = data.frame(t(pag@assays$MERFISH@scale.data[g,]),rep = pag$rep_name, 
                      cl = pag$stable) %>% group_by(cl,rep) %>% summarise_all(.funs = mean)
  
  i.df = data.frame(score = pag$ieg.score, 
                    cl = pag$stable, rep = pag$rep_name,
                    beh = pag$beh, sex = pag$Sex)
