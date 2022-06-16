# Adjust IEG
# October 5 2021
# Eric Vaughn

# The lowest amount of Fos found in a cell across each replicate is 1: the volume-normalized
# conversion of 1 Fos count should be equivalent across replicates.
# I'll rescale all non-zero counts to a common scale (median min/max per gene across replicates) for each replicate 
# to reduce these batch effects.

load('pag.RObj')

g = c('Fos','Fosb','Fosl2','Egr1','Nr4a3','Nr4a1','Junb')

ieg = data.frame(rep= pag$rep_name, name = pag$X, t(as.matrix(pag@assays$MERFISH@data[g,]))) 
ieg.melt = ieg %>% reshape2::melt()
ieg.stats = ieg.melt %>% group_by(rep, variable) %>% summarise(lo = min(value[value>0]), hi = max(value))
ieg.summary = ieg.stats %>% group_by(variable) %>% summarise(med.lo = median(lo), med.hi = median(hi))

cells = Cells(pag)

# rescale the non-zero distribution to median min and max that exists across all replicates
adjustIEG = 
  sapply(unique(ieg$rep), function(r){
    sapply(g, function(x){
      cells = which(pag@assays$MERFISH@data[x, as.character(ieg$name[ieg$rep==r])]>0) %>% names()
      low.th = ieg.summary$med.lo[ieg.summary$variable==x]
      up.th = ieg.summary$med.hi[ieg.summary$variable==x]
      scales::rescale(
        pag@assays$MERFISH@data[x,cells],
        to = c(low.th, up.th)
      )
    },
    USE.NAMES = T,
    simplify = F)
    
  },
  USE.NAMES = T,
  simplify = F
  )

# fill each nonzero IEG vector with zeros where they exist so I can compose the appropriate matrix
adjustIEG2 =
  sapply(1:31, function(r){
    sapply(g, function(g){
      cells = as.character(ieg$name[ieg$rep== unique(ieg$rep)[r]])
      out= setdiff(cells, names(adjustIEG[[r]][[g]]))
      a = c(rep(0,length(out)), adjustIEG[[r]][[g]])
      names(a) = c(out, names(adjustIEG[[r]][[g]]))
      a = a[cells]
    },
    USE.NAMES = T,
    simplify = F)
  },
  USE.NAMES = T,
  simplify = F
  )

adjustIEG3 = sapply(1:31, function(r){
  do.call(rbind, adjustIEG2[[r]])
},
  simplify = F,
  USE.NAMES = T
)

# replace in object with sparse matrix
pag@assays$MERFISH@data[g,] = Matrix::Matrix(do.call(cbind, adjustIEG3), sparse = T)[,Cells(pag)]

pag = ScaleData(pag, vars.to.regress = c('rep_name'), model.use = 'linear', features = g)
pag$ieg.score = Matrix::colSums(pag@assays$MERFISH@scale.data[g,])

save(pag, file = 'pag.batch corrected IEG 10.5.21.RObj')


# Statistical tests -------------------------------------------------------


i.df = data.frame(score = pag$ieg.score, 
                  cl = pag$stable, rep = pag$rep_name,
                  beh = pag$beh)

pag.cl = grep(pattern = '\\.', x = as.vector(unique(pag$stable)), value = T)

beh.ieg = i.df %>% filter(cl %in% pag.cl)  %>% group_by(cl,beh) %>% 
  summarise(av.ieg = mean(score), n= n())

# filter out exercise clusters that are very low and are likely to be poor representations of IEG signal
beh.ieg = beh.ieg %>% filter(n>10)

fc = beh.ieg %>% group_by() %>% mutate(sc = scales::rescale(av.ieg, c(0,5))) %>%
  group_by(cl) %>% 
  mutate(cl.av.sc = mean(sc), cl.av = mean(av.ieg), nai.av = mean(sc[beh %in% c('Naive_F','Naive_M')])) %>% 
  mutate(fc = (sc/cl.av.sc)-1, diff= sc-cl.av.sc, wt_fc = (sc/cl.av.sc)-1+diff, 
         nai.diff = sc-nai.av, nai.fc = (sc/nai.av)-1, nai.wt_fc = (sc/nai.av)-1 + nai.diff)


cells = data.frame(cl=pag$stable[pag$stable %in% pag.cl], 
                   beh = pag$beh[pag$stable %in% pag.cl], 
                   cellname = Cells(pag)[pag$stable %in% pag.cl])
beh.t = list()
beh.n = table(data.frame(beh = pag$beh, cl = pag$stable))
scores = pag$ieg.score
for(b in unique(fc$beh)){
  message(b)
  beh.t[[b]]= sapply(pag.cl, FUN = function(x){
    message(x)
    if(beh.n[b,x]>0){
      t.test(
        scores[as.character(dplyr::filter(cells, beh==b & cl == x) %>% dplyr::select(cellname) %>% unlist())], 
        scores[as.character(dplyr::filter(cells, beh!=b & cl == x) %>% dplyr::select(cellname) %>% unlist())]) 
    }
    else{
      message(x,' too low')
      return(1)
    }
  },
  simplify = F, USE.NAMES = T)
  beh.t[[b]] = lapply(beh.t[[b]], broom::tidy) # tidy up output
  beh.t[[b]] = do.call(rbind,beh.t[[b]]) # compress into a df   
}

beh.t2 = mapply(cbind, beh.t, "Behavior"=names(beh.t), SIMPLIFY=F)
beh.t2 = sapply(beh.t2, function(x) cbind(x, Cluster=rownames(x)), simplify = F, USE.NAMES = T)

cl.pval = purrr::reduce(beh.t2, full_join)
cl.pval$fdr = p.adjust(cl.pval$p.value, method = 'fdr')

fc = merge(fc, cl.pval[,c(5,11:13)], by.x = c('cl','beh'), by.y = c('Cluster','Behavior'))

# make a new column for significant clusters that are both greater than ieg.score 0.5 from the cluster average and
# have a fdr corrected pvalvue of less than .05

fc = fc %>% group_by(cl) %>% mutate(sig = ifelse(av.ieg-cl.av > .3 & fdr < 0.051 & av.ieg>1, yes = T, no = F))

# Catherine wants parenting first
fc$beh = factor(fc$beh, levels = c('Virgin Parenting_F','Parenting_F','Parenting_M','Infanticide_M',
                                   'Aggression_M','Mating_M','Mating_F','Fear_F','Fear_M',
                                   'Exercise_F','Exercise_M','Naive_M','Naive_F'), ordered = T)



# Rerun ttest but against naive conditions --------------------------------

nai.t = list()
nai.n = table(data.frame(beh = pag$beh, cl = pag$stable))
scores = pag$ieg.score


plan('multisession')
# test all behaviors except naive
bs = unique(fc$beh)[-c(12,13)] %>% as.character()
nai.t = future_sapply(bs, function(b){
  l = sapply(pag.cl, FUN = function(x){
    if(nai.n[b,x]>0){
      t.test(
        scores[as.character(dplyr::filter(cells, beh==b & cl == x) %>% dplyr::select(cellname) %>% unlist())], 
        scores[as.character(dplyr::filter(cells, beh %in% c('Naive_F','Naive_M
                                                            ') & cl == x) %>% dplyr::select(cellname) %>% unlist())]) 
    }
    else{
      return(1)
    }
  },
  simplify = F, USE.NAMES = T)
  l = lapply(l, broom::tidy) # tidy up output
  l = do.call(rbind,l) # compress into a df
  return(l)
  },
  USE.NAMES = T, simplify = F)

nai.t2 = mapply(cbind, nai.t, "Behavior"=names(nai.t), SIMPLIFY=F)
nai.t2 = sapply(nai.t2, function(x) cbind(x, Cluster=rownames(x)), simplify = F, USE.NAMES = T)

cl.pval = purrr::reduce(nai.t2, full_join)
cl.pval$n.fdr = p.adjust(cl.pval$p.value, method = 'fdr')
cl.pval$n.p.value = cl.pval$p.value
#filler spaces for naive animals
tmp = cl.pval[cl.pval$Behavior %in% c('Fear_M','Fear_F'),]
tmp$n.fdr = 1
tmp$n.p.value = 1
tmp$Behavior = plyr::mapvalues(x =tmp$Behavior, from = c('Fear_M','Fear_F'), c('Naive_M','Naive_F'))
cl.pval = rbind(cl.pval,tmp)

fc = merge(fc, cl.pval[,11:14], by.x = c('cl','beh'), by.y = c('Cluster','Behavior'))
fc = fc %>% group_by(cl, beh) %>% slice(1)

# make a new column for significant clusters that are both greater than ieg.score 0.5 from the cluster average and
# have a fdr corrected pvalvue of less than .05

fc = fc %>% group_by(cl) %>% mutate(n.sig = ifelse(nai.diff > .3 & n.fdr < 0.051 & av.ieg>1, yes = T, no = F)) %>% group_by()

# put parenting first
fc$beh = factor(fc$beh, levels = c('Virgin Parenting_F','Parenting_F','Parenting_M','Infanticide_M',
                                   'Aggression_M','Mating_M','Mating_F','Fear_F','Fear_M',
                                   'Exercise_F','Exercise_M','Naive_M','Naive_F'), ordered = T)


save(fc, file = 'fc corrected IEG 10.6.21 w Naive Test.R')



# DT Activity -------------------------------------------------------------

 
i.df = data.frame(score = pag$ieg.score, 
                  cl = pag$dt, rep = pag$rep_name,
                  beh = pag$beh)

beh.ieg = i.df %>% filter(cl %in% 1:19)  %>% group_by(cl,beh) %>% 
  summarise(av.ieg = mean(score), n= n())

dt.fc = beh.ieg %>% group_by() %>% mutate(sc = scales::rescale(av.ieg, c(0,5))) %>%
  group_by(cl) %>%
  mutate(cl.av.sc = mean(sc), cl.av = mean(av.ieg), nai.av = mean(sc[beh %in% c('Naive_F','Naive_M')])) %>% 
  mutate(fc = (sc/cl.av.sc)-1, diff= sc-cl.av.sc, wt_fc = (sc/cl.av.sc)-1+diff, 
         nai.diff = sc-nai.av, nai.fc = (sc/nai.av)-1, nai.wt_fc = (sc/nai.av)-1 + nai.diff)

dt.fc$beh = factor(dt.fc$beh, levels = beh.order, ordered = T)



# statistically test SMCs
cells = data.frame(cl=pag$dt[pag$stable %in% pag.cl], 
                   beh = pag$beh[pag$stable %in% pag.cl], 
                   cellname = Cells(pag)[pag$stable %in% pag.cl])

plan('multisession')

beh.t = list()
beh.n = table(data.frame(beh = pag$beh, cl = pag$dt))
scores = pag$ieg.score
for(b in unique(fc$beh)){
  message(b)
  beh.t[[b]]= sapply(as.character(1:19), FUN = function(x){
    message(x)
    if(beh.n[b,x]>0){
      t.test(
        scores[as.character(dplyr::filter(cells, beh==b & cl == x) %>%
                              dplyr::select(cellname) %>% unlist())], 
        scores[as.character(dplyr::filter(cells, beh!=b & cl == x) %>% group_by(beh) %>%
                              slice(sample(n(), min(500, n()))) %>% group_by() %>%
                              dplyr::select(cellname) %>% unlist())]) 
    }
    else{
      message(x,' too low')
      return(1)
    }
  },
  simplify = F, USE.NAMES = T)
  beh.t[[b]] = lapply(beh.t[[b]], broom::tidy) # tidy up output
  beh.t[[b]] = do.call(rbind,beh.t[[b]]) # compress into a df   
}

beh.t2 = mapply(cbind, beh.t, "Behavior"=names(beh.t), SIMPLIFY=F)
beh.t2 = sapply(beh.t2, function(x) cbind(x, Cluster=rownames(x)), simplify = F, USE.NAMES = T)

cl.pval = purrr::reduce(beh.t2, full_join)
cl.pval$fdr = p.adjust(cl.pval$p.value, method = 'fdr')

dt.fc = merge(dt.fc, cl.pval[,c(5,11:13)], by.x = c('cl','beh'), by.y = c('Cluster','Behavior'))

# make a new column for significant clusters that are both greater than ieg.score 0.5 from the cluster average and
# have a fdr corrected pvalvue of less than .05

dt.fc = dt.fc %>% group_by(cl) %>% mutate(sig = ifelse(av.ieg-cl.av >.1 & fdr < 0.051 & av.ieg>1, yes = T, no = F),
                                          sig.up = ifelse(av.ieg-cl.av >.1 & fdr < 0.051, yes = T, no = F),
                                          sig.down = ifelse(fdr < 0.051 & nai.diff<(-.1), yes = T, no = F))

# put parenting first
dt.fc$beh = factor(dt.fc$beh, levels = c('Infanticide_M','Parenting_F','Parenting_M','Virgin Parenting_F',
                                   'Aggression_M','Mating_M','Mating_F','Fear_F','Fear_M',
                                   'Exercise_F','Exercise_M','Naive_M','Naive_F'), ordered = T)

# make a heatmap of up or down activation of each dt

df = data.frame(cl = dt.fc$cl, beh = dt.fc$beh, 
                ud = 
                  ifelse(
                    dt.fc$sig.up ==T, yes = 1, 
                    no =
                      ifelse(dt.fc$sig.down == T, yes = -1,
                             no = 0)
                      )) %>% filter(! beh %in% naive) %>% mutate(ud = factor(ud))
dt.heat = function(clusters, lab = 'Weighted FC', fill = 'wt_fc', max = 7, min  = -2, mid = 1.5,
                  x.size = 10, df = fc) {
  df %>% filter(cl %in% clusters) %>% group_by() %>% mutate(cl = factor(cl, levels =clusters, ordered= T)) %>%
    ggplot(., aes(x = cl, y = beh )) +
    geom_point(aes_string(color = fill, size = 'abs.diff', alpha = 'abs.diff') ) +
    scale_color_gradient2(low = 'steelblue', mid = 'white', high = 'red' ,midpoint = mid, limits = c(min,max), oob = scales::squish, na.value = 'red') +
    ylab('Behavior') +
    xlab('Cluster') +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = x.size),
          axis.text.y = element_text(angle = 0, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank()) +
    labs(fill = paste0(lab))
}
dt.hm = dt.heat(dt.order, df = filter(dt.fc, !beh %in% naive), min = -1, max = 1, mid = 0, fill = 'wt_fc') + 
  labs(color = 'Weighted FC', size = 'Absolute\nDifference', alpha=NULL, x = '', y = '')
ggsave(dt.hm, file = 'dt.heatmap.pdf', height = 4, width = 8.5, dpi = 300, useDingbats=F)
