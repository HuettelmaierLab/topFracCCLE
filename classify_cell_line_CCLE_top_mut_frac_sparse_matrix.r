#####################################################################################################
#                                                                                                   #
# classify cell line based on SNPs delivered in a vcf file                                          #
#                                                                                                   #
# vcf_file: vcf format file containing sequence variations of the query sample                      #
# genome: genome version GRCh37/GRCh38 (GRCh38)                                                     #
# lines: cell lines to be tested for, if NULL all lines from data set are used (NULL)               #
# num_out: number of returned best hits (3)                                                         #
# minSNPs: minimum number of SNPs reported for a cell line from the data set to be considered (10)  #
# data_dir: directory containing SNP and meta data (data)                                           #
#####################################################################################################
classify_cell_line_CCLE_top_mut_frac = function(vcf_file, genome = "GRCh38", lines = NULL, num_out = 3, minSNPs = 10, data_dir = "data/")
{
    library(prodlim)
    library(Matrix)
    
    sep = .Platform$file.sep
    
    ## read matrix of mutations
    bin_mut_matrix = readRDS(paste0(data_dir, sep, "distinct_mutation_CCLE_full_sparse.rds"))
    
    ## read cell line metadata
    metadata = read.table(paste0(data_dir, sep, "depmap_metadata.tsv"), sep = "\t", header = T, 
                          quote = "", comment.char = "")
    
    ## remove cell lines with less than minSNPs reported SNPs
    remove = which(colSums(bin_mut_matrix$lines) < minSNPs)
    
    if(length(remove) > 0)
    {
        bin_mut_matrix$lines = bin_mut_matrix$lines[,-remove]
    }
    
    if(!is.null(lines)) ## reduce cell lines to a given subset
    {
        bin_mut_matrix$lines = bin_mut_matrix$lines[,which(colnames(bin_mut_matrix$lines) %in% lines)]
    }
    
    ## keep only the relevant genome positions
    if(genome == "GRCh37")
    {
        bin_mut_matrix$position = bin_mut_matrix$position[,-c(3,4)]
    } 
    else if(genome == "GRCh38")
    {
        bin_mut_matrix$position = bin_mut_matrix$position[,-c(1,2)]
    }
    
    m = ncol(bin_mut_matrix$lines)
    
    ############ load test data
    vcf = read.table(vcf_file, sep = "\t", header = F)
    
    ## reduce to columns present in reference
    vcf = vcf[,c(1,2,4,5)]
    
    hits = matrix(nrow = (m), ncol = 3)
    rownames(hits) = colnames(bin_mut_matrix$lines)
    colnames(hits) = c("SNPs_found", "total_SNPs", "fraction_found")
    
    pb = txtProgressBar(1, m, style = 3)
    
    ## for each cell line extract mutations and compare them with query
    for(i in 1 : m)
    {
        setTxtProgressBar(pb, i)
        
        vars = bin_mut_matrix$position[which(bin_mut_matrix$lines[,i] == 1),c(1,2,3,4)]
        
        matches = row.match(vars, vcf, nomatch = 0)
        
        hits[i, 1] = sum(matches != 0)
        hits[i, 2] = nrow(vars)
    }
    
    hits[,3] = hits[,1] / hits[,2]
    
    
    ## sort hit list
    hits = hits[order(hits[,3], decreasing = T),]
    
    hits = data.frame(line = rownames(hits), 
                      RRID = metadata[match(rownames(hits), metadata$stripped_cell_line_name), "RRID"], 
                      hits)
    
    # ToDo: add derivatives and parent cell line
    
    return(hits[1:num_out,])
}

