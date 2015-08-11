parse_gtf = function(gtf_name) {
    raw_lines = read.table(gtf_name, sep = "\n", stringsAsFactors = F)[, 1]
    get_col = function(i) {
        unlist(lapply(strsplit(raw_lines, "\t"), function(x) x[[i]]))
    }
    attribs = get_col(9)
    all_attribs = strsplit(attribs, "; ")
    
    
    get_attrib = function(key = "transcript_support_level") {
        out = lapply(all_attribs, function(x) {
            keep = grepl(key, x)
            str = "NA"
            if (sum(keep) > 0) 
                str = strsplit(x[keep], " ")[[1]][2]
            return(str)
        })
        return(unlist(out))
    }
    
    
    support = get_attrib()
    enst_id = get_attrib("transcript_id")
    ensg_id = get_attrib("gene_id")
    gene_name = get_attrib("gene_name")
    chrm = get_col(1)
    strand = get_col(7)
    start = as.numeric(get_col(4))
    end = as.numeric(get_col(5))
    
    ucsc = paste0(chrm, ":", paste0(start, "-", end))
    
    enst_dict = data.frame(chrm, start, end, strand, ucsc, support, enst_id, ensg_id, gene_name, row.names = enst_id, stringsAsFactors = F)
    
    return(enst_dict)
    # plot.pie = function(dat){ uniq = unique(dat) tmp = sapply(uniq, function(x)return(sum(x == dat))) pie(tmp) } plot.pie(support) enst2support = support
    # names(enst2support) = enst_id enst2ensg = ensg_id names(enst2ensg) = enst_id
}

if (F) {
    parse_gtf("mycounts_annotation.gff")
    save(enst_dict, file = "enst_dicts.save")
} 
