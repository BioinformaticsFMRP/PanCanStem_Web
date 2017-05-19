replace.NA <- function(data,type.info,by = "mean"){
  if(!"group" %in% colnames(type.info)) stop("type.info must have group column")
  if(!"sample" %in% colnames(type.info)) stop("type.info must have a sample column")
  
  # Do we have NAs?
  if(!any(is.na(data))){
    message("No NAs were found")
    return(data)
  }
  # get NAs index 
  idx <- which(is.na(data) == TRUE,arr.ind=TRUE)
  count <- table(rownames(idx))
  message("======= Status Number of NA in probes ========")
  message("--------------------- Summary------------------")
  print(summary(as.numeric(count)))
  message("\n----------- Probes with more nb of NAs -----------")
  print(head(sort(count,decreasing = T)))
  message("===============================================")
                 
  idx <- cbind(idx, mean = NA, median = NA)
  
  # For each NA value calculate the mean for the same probe for the samples
  # where it belongs
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    probe <- rownames(idx)[line]
    sample <- colnames(data)[col]
    group <- type.info[type.info$sample == sample,"group"]
    samples.in.group <- type.info[type.info$group == group,]$sample
    
    # get the probe value for all samples in the same group 
    aux <- data[rownames(data) %in% probe, colnames(data) %in% samples.in.group] 
    
    idx[line,3] <- mean(as.numeric(aux),na.rm = TRUE)
    idx[line,4] <- median(as.numeric(aux),na.rm = TRUE)
  }
  # Step 2 replace
  for(line in 1:nrow(idx)){
    row <- idx[line,1]
    col <- idx[line,2]
    if(by == "mean"){
      data[idx[line,1],idx[line,2]] <- idx[line,3]  
    } else if(by == "median") { 
      data[idx[line,1],idx[line,2]] <- idx[line,4]
    }
  }
  return(data)
}



