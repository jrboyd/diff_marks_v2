#conduct a binomial test for enrichment for input gene_list within specified msigdb gene_set.  
#if no gene_set specified, entire msigdb will be checked (likely leading to over FDR correction)
#p_thresh is appplied to bonferonni corrected pvalues.
#set_regex may be used to limit msig query set further, to KEGG or CHARAFE for example.
binom_msig_enrich = function(gene_list, gene_set = NULL, p_thresh = .05, set_regex = NULL){
  if(is.null(gene_set)){
    
    names = list()
    sets = list()
    hidden = sapply(msigdb_categories, function(x){
      names <<- append(names, x$names)
    })
    hidden = sapply(msigdb_categories, function(x){
      sets <<- append(sets, x$sets)
    })
    names = unlist(names)
    #names(sets) = names
    gene_set_group = list(names = unlist(names), sets = sets)
  }else{
    if(sum(grepl(gene_set, names(msigdb_categories))) < 1){
      warning(paste(gene_set, "not found. valid gene sets are:", paste(names(msigdb_categories), collapse = ', ')))
      return(NULL)
    }
    gene_set_group = msigdb_categories[[gene_set]]
  }
  
  if(!is.null(set_regex)){
    keep = grepl(set_regex, gene_set_group$names)
    if(sum(keep) < 1){
      warning(paste("regex", set_regex, "not found in gene_set", ifelse(is.null(gene_set), character(), gene_set)))
      return(NULL)
    }
    gene_set_group$names = gene_set_group$names[keep]
    gene_set_group$sets = gene_set_group$sets[keep]
  }
  
  n_success = sapply(gene_set_group$sets, function(x){length(intersect(x, gene_list))})
  n_trials = length(gene_list)
  prob = sapply(gene_set_group$sets, function(x){length((x))}) / 25000
  all_tests = matrix(0, nrow = length(n_success), ncol = 5)
  rownames(all_tests) = gene_set_group$names
  colnames(all_tests) = c('adj_pval', 'fold-enriched', 'n_obs', 'n_exp', 'set_size')
  for(i in 1:length(n_success)){
    pval = binom.test(x = n_success[i], p = prob[i], n = n_trials)$p.value*length(n_success)
    fe = n_success[i] / n_trials / prob[i]
    
    # print(c(pval, fe))
    set_size = length(gene_set_group$sets[[i]])
    all_tests[gene_set_group$names[i],] = c(pval, fe, n_success[i], prob[i] * set_size, set_size )
  }
  keep = all_tests[,1] < p_thresh
  all_tests = all_tests[keep,]
  
  
  o = order(all_tests[,1])
  final = all_tests[o,]
  return(final)
}

binom_go_enrich = function(test_ensg){
  
  test_sym = ensg_dict[test_ensg,]$gene_name
  keep = !sapply(names(sym2gos[test_sym]), function(x)is.na(x))
  unmapped = sum(!keep)
  unmapped_sym = test_sym[!keep]
  mapped_sym = test_sym[keep]
  gos = unique(unlist(sym2gos[mapped_sym]))
  is_bp = go_table[gos,3] == "biological_process"
  gos = gos[is_bp]
  
  membership = matrix(F, nrow = length(mapped_sym), ncol = length(gos))
  rownames(membership) = mapped_sym
  colnames(membership) = gos
  
  min_size = 3
  keep = sapply(gos2sym[gos], function(x){
    length(intersect(x, mapped_sym)) >= min_size
  })
  gos = gos[keep]
  
  total_genes = length(sym2gos)
  n_trials = length(mapped_sym)
  
  tests = sapply(gos2sym[gos], function(x){
    n_success = length(intersect(x, mapped_sym))
    
    prob = length(x) / total_genes
    return(binom.test(x = n_success, n = n_trials, p = prob, alternative = "greater"))
  })
  n_exp = unlist(tests[6,]) * sapply(gos2sym[gos], length)
  
  pvals = unlist(tests[3,])
  pvals = pvals * length(gos)
  pass = pvals < .05
  pass_go = names(pass)[pass]
  n_obs = unlist(tests[1,])
  go_result = data.frame(name = go_table[pass_go,2], data.frame(adj_pval = pvals, n_obs, n_exp, fold_enriched = n_obs / n_exp)[pass_go,])
  return(go_result)
}