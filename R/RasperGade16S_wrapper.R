#' @title Predict 16S rRNA GCN from sequences
#' @description A wrapper function to "one-click-and-run" the default pipeline
#' @export
#' @rdname predict_16SGCN_from_sequences
predict_16SGCN_from_sequences = function(seqs,
                                         custom_db_dir_path=NULL,
                                         numCores=0,
                                         rmTmpFiles=FALSE){
  if(missing(seqs)) {
    seqs = system.file("extdata/Demo","demo.SILVA.fasta",package="RasperGade16S",mustWork=TRUE)
    cat(sprintf("No FASTA file supplied, running prediction using demo sequences in\n%s\n\n",seqs))
  }

  custom_db_fpaths = NULL
  if (!is.null(custom_db_dir_path)) {
    custom_db_dir_path = normalizePath(custom_db_dir_path)
    custom_db_fpaths = find_custom_db_fpaths(custom_db_dir_path)
    message(
      sprintf("Using custom database: `%s`", custom_db_dir_path)
    )
  }

  align.out = align_with_HMM_and_trim(seqs=seqs)

  if (is.null(custom_db_fpaths)) {
    epa.out = insert_query_with_EPA(
      seqs="RasperGade16S_align/trimmed.afa",
      numCores=numCores
    )
  } else {
    epa.out = insert_query_with_EPA(
      seqs="RasperGade16S_align/trimmed.afa",
      tree=custom_db_fpaths$tree,
      ref.seqs=custom_db_fpaths$ref.seqs,
      model=custom_db_fpaths$model,
      numCores=numCores
    )
  }

  pred.GCN = predict_16SGCN_from_jplace(
    jplace=epa.out$jplace,
    custom_db_fpaths=custom_db_fpaths,
    numCores=numCores,
    save2file=TRUE
  )

  if (rmTmpFiles) {
    cleanup_tmp_dirs()
  }

  return(pred.GCN)
}


#' @title Predict 16S rRNA GCN from sequence placements
#' @description A wrapper function to the prediction from placements
#' @export
#' @rdname predict_16SGCN_from_jplace
predict_16SGCN_from_jplace = function(jplace,
                                      custom_db_fpaths=NULL,
                                      numCores=1,
                                      save2file=FALSE){
  if (numCores == 0) {
    numCores = 1
  }
  insert.locations = parse_jplace(jplace, split=numCores)

  fmr = RasperGade16S.refdata$FMR
  if (!is.null(custom_db_fpaths)) {
    load(custom_db_fpaths$refdata)
    fmr = refdata$FMR
  }

  cat("Insert locations parsed.\n")
  if(save2file) saveRDS(insert.locations,sprintf("%s.locations.RDS",jplace))
  insert.prediction = mclapply(insert.locations,function(this.insert){
    this.res = predictHiddenStateWithPE(
      FMR=fmr,
      query.keys=this.insert$hash,
      laplace=FALSE
    )
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


cleanup_tmp_dirs = function() {
  tmp_dirs = c(
    file.path(getwd(), "RasperGade16S_align"),
    file.path(getwd(), "RasperGade16S_EPA")
  )
  cat("Removing temporary files...\n")
  for (dir_path in tmp_dirs) {
    cat(sprintf("  `%s`\n", dir_path))
    if (dir.exists(dir_path)) {
      unlink(dir_path, recursive=TRUE)
    }
  }
}


find_custom_db_fpaths = function(custom_db_dir_path) {

  if (!dir.exists(custom_db_dir_path)) {
    stop(
      sprintf(
        "Error. Custom dir does not exist: `%s`",
        custom_db_dir_path
      )
    )
  }

  required_file_basenames = c(
    "RasperGade16S.refdata.rda",
    "reference.tre",
    "RAxML.model.txt",
    "reference.trimmed.afa"
  )
  for (file_basename in required_file_basenames) {
    file_path = file.path(custom_db_dir_path, file_basename)
    if (!file.exists(file_path)) {
      message(
        sprintf("Error: missing required file `%s` in the custom database dir `%s`",
          file_basename,
          custom_db_dir_path
        )
      )
      message("Here is the complete list of required files:")
      for (f in required_file_basenames) {
        message(sprintf("  `%s`", f))
      }
      stop("Terminating.")
    }
  }
  return(
    list(
      refdata  = file.path(custom_db_dir_path, "RasperGade16S.refdata.rda"),
      tree     = file.path(custom_db_dir_path, "reference.tre"),
      model    = file.path(custom_db_dir_path, "RAxML.model.txt"),
      ref.seqs = file.path(custom_db_dir_path, "reference.trimmed.afa")
    )
  )
}