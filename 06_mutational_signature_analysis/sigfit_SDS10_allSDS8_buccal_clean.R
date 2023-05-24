# ---
#   title: "sigfit_SDS10_buccal"
# ---
# Dec 2021


library(sigfit)


###### load data
groupname="SDS10_allSDS8_buccal"

## We know that the colony_info file and the mutcounts_matrix samples are in the same order.
# Otherwise, must correctly order
mutcounts_matrix = read.table(file="../data/SDS10_allSDS8_buccal_mutcounts_matrix.txt", header=TRUE, stringsAsFactors=FALSE, sep=" ")
colonyinfo_all = data.frame(sample=colnames(mutcounts_matrix), celltype = "SDS")

############# Make samp_type object
Nsamp = ncol(mutcounts_matrix)
samp_type = colonyinfo_all




####### testing
#mutcounts_matrix_both = mutcounts_matrix_both[,1:10]
#samp_type = samp_type[1:10,]
#samp_names = samp_names[1:10]


##### adding blood signature to cosmic signatures (they are in the same order)
mysigs = read.table("../data/mutspectrum_components.comp_distn_mean.noprior_prior.SDS10_allSDS8_buccal_divideshared.txt", stringsAsFactors = F, header=T)


##### selecting the signatures that are the union of hdp and sigprofiler results
selectsigs = t(mysigs[,2:3])
rownames(selectsigs) = c("SBSblood","SBS1")


###### fit cosmic + blood signatures
mcmc_samples_fit = fit_signatures(counts = t(mutcounts_matrix),
    signatures = selectsigs,
    iter = 100000,
    warmup = 2000,
    chains = 1,
    seed = 1114)



## refit to abundant signatures
exposures <- retrieve_pars(mcmc_samples_fit, par = "exposures", hpd_prob = 0.90)
save(samp_type, exposures, file=paste("exposures_",groupname, ".sigfit.chain1.Rdata", sep="") )
