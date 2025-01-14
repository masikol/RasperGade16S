#' @title Predict 16S rRNA GCN from sequences
#' @description A wrapper function to "one-click-and-run" the default pipeline
#' @export
#' @rdname predict_16SGCN_from_sequences
predict_16SGCN_from_sequences = function(seqs){
  if(missing(seqs)){
    seqs = system.file("extdata/Demo","demo.SILVA.fasta",package="RasperGade16S",mustWork=TRUE)
    cat(sprintf("No FASTA file supplied, running prediction using demo sequences in\n%s\n\n",seqs))
    }
  align.out = align_with_HMM_and_trim(seqs=seqs)
  epa.out = insert_query_with_EPA(seqs="RasperGade16S_align/trimmed.afa")
  pred.GCN = predict_16SGCN_from_jplace(epa.out$jplace,save2file = TRUE)
  return(pred.GCN)
}

#' @title Predict 16S rRNA GCN from sequence placements
#' @description A wrapper function to the prediction from placements
#' @export
#' @rdname predict_16SGCN_from_jplace
predict_16SGCN_from_jplace = function(jplace,numCores = 1,save2file=FALSE){
  insert.locations = parse_jplace(jplace,split = numCores)
  cat("Insert locations parsed.\n")
  if(save2file) saveRDS(insert.locations,sprintf("%s.locations.RDS",jplace))
  insert.prediction = mclapply(insert.locations,function(this.insert){
    this.res = predictHiddenStateWithPE(FMR = RasperGade16S.refdata$FMR,
                                        query.keys = this.insert$hash,laplace = FALSE)
    this.res$weight = this.insert$weight
    return(this.res)
  },mc.cores = numCores)
  cat("Copy number predicted.\n")
  insert.res = list(hsp=do.call(rbind,lapply(insert.prediction,function(x){x$hsp})),
                    error=do.call(c,lapply(insert.prediction,function(x){x$error})),
                    weight = do.call(c,lapply(insert.prediction,function(x){x$weight})))
  if(save2file) saveRDS(insert.res,sprintf("%s.prediction.RDS",jplace))
  unique.insert.res = lapply(split(x=1:length(insert.res$error),insert.res$hsp$label),function(i){
    new.error = mix_errors_by_weight(insert.res$error[i],insert.res$weight[i])
    new.stat = unname(calculate_error_mean_and_var(new.error))
    new.hsp = data.frame(node=-1,label=insert.res$hsp$label[i[1]],x=new.stat[1],var=new.stat[2])
    return(list(hsp=new.hsp,error=list(new.error)))
  })
  cat("Multiple placements merged.\n")
  insert.res = list(hsp=do.call(rbind,lapply(unique.insert.res,function(x){x$hsp})),
                    error=do.call(c,lapply(unique.insert.res,function(x){x$error})))
  insert.discrete.res = 
    do.call(rbind,
            mclapply(split(1:length(insert.res$error),
                         mod(1:length(insert.res$error),numCores)),
                   function(i){
                     discretizeResult(res = insert.res$hsp[i,],error =insert.res$error[i],laplace = FALSE)
                     },mc.cores = numCores))
  cat("Copy number discretized.\n")
  insert.GCN = insert.discrete.res$x
  names(insert.GCN) = insert.discrete.res$label
  if(save2file) saveRDS(list(discrete=insert.discrete.res,continuous=insert.res),
                        sprintf("%s.prediction.RDS",jplace))
  return(list(tab=insert.discrete.res[,-1],GCN=insert.GCN,error = insert.res$error))
}

