# delta(p4) ~ Fst 

# DEPRECATED see fst_p4_variance_explained

load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/coloc/master_coloc-1e-04.RData')
head(master_coloc)
load('/mnt/lab_data/montgomery/nicolerg/local-eqtl/admixed/annotation/fst/master_all.RData')
head(master)

# see if p4 can be explained by Fst 

# f <- as.formula(sprintf('%s ~ AFR + ASN',v))
# lm.fit <- lm(f, data=m)
# rsq <- summary(lm.fit)$r.squared
