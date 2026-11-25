#####################################################################################################
#                                                                                                   #
# classify cell line based on SNPs delivered in a vcf file                                          #
#                                                                                                   #
# vcf_file: vcf format file containing sequence variations of the query sample                      #
# genome: genome version, either "hg19" (hg19/GRCh37), "hg38" (hg38/GRCH38) or "T2T" (T2T-CHM13v2)  #
# default "hg38"                                                                                    #
# lines: cell lines to be tested for, if NULL all lines from data set are used (NULL)               #
# num_out: number of returned best hits (3)                                                         #
# minSNPs: minimum number of SNPs reported for a cell line from the data set to be considered (10)  #
# data_dir: directory containing SNP and meta data (data)                                           #
#####################################################################################################
classify_cell_line_CCLE_top_mut_frac = function(vcf_file, genome = "hg38", lines = NULL, num_out = 3, minSNPs = 10, data_dir = "data/")
{
    library(prodlim)
    library(Matrix)
    
    sep = .Platform$file.sep
    
    ## read matrix of mutations
    bin_mut_matrix = readRDS(paste0(data_dir, sep, "distinct_mutation_CCLE_full_sparse.rds"))
    
    ## read cell line metadata
    metadata = read.table(paste0(data_dir, sep, "depmap_metadata.tsv"), sep = "\t", header = T, 
                          quote = "", comment.char = "")
    
    
    ## keep only the relevant genome positions
    if(genome == "hg19")
    {
        bin_mut_matrix$position = bin_mut_matrix$position[,-c(3,4,7,8)]
    } 
    else if(genome == "hg38")
    {
        bin_mut_matrix$position = bin_mut_matrix$position[,-c(5,6,7,8)]
    }
    else if(genome == "T2T")
    {
        bin_mut_matrix$position = bin_mut_matrix$position[,-c(3,4,5,6)]
    }
    else	## invalid genome
    {
    	stop("Invalid genome! Please choose 'hg19', 'hg38', or 'T2T'.")
    }
    
    ## remove lines with no information (NAs introduced during lift-over)
    NA_lines = which(is.na(bin_mut_matrix$position[,4]))
    
    if(length(NA_lines) > 0)
    {
        bin_mut_matrix$position = bin_mut_matrix$position[-NA_lines,]
        bin_mut_matrix$lines = bin_mut_matrix$lines[-NA_lines,]
    }
    
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
    
    m = ncol(bin_mut_matrix$lines)
    
    ############ load test data
    vcf = read.table(vcf_file, sep = "\t", header = F)
    
    ## reduce to columns present in reference
    vcf = vcf[,c(4,5,1,2)]
    
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
                      RRID = metadata[match(rownames(hits), metadata$StrippedCellLineName), "RRID"], 
                      hits)
    
    # ToDo: add derivatives and parent cell line
    
    return(hits[1:num_out,])
}

