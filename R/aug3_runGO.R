mymxlinks = function(mx,name,sig = 5, fdr = 0.1){

  myscore = function(mat,column, sig = 1){
    df = data.frame(score = mat[,column],stringsAsFactors = F)
    df$gene = rownames(mat)
    rownames(df) = df$gene
    df$index=0
    wdf = which(df$score >= sig)
    df$index[wdf]=1
    df = df[,c('index','score','gene')]
    df
  }

  z = NULL
  for(i in 1:ncol(mx)) z[i] = list(mx[,i,drop=F])
  names(z)=name
  workdir = getwd()

  lscore = lapply(z,myscore,sig=sig)
  score3 = sapply(lscore,function(x) x = all(x[,"index"]==0))

  w = which(score3 == T)

  nam0 = names(lscore)

  if(length(w) > 0) lscore = lscore[-w]

  nam = names(lscore)

  tst = NULL

  for (i in 1:length(lscore)) tst[[i]] = myrun_go_enrichment2(fdrThresh = fdr,
    curr_exp = names(lscore)[i],score=lscore[[i]])
  dir = getwd()

  if(dir != workdir)  setwd(workdir)

  lmcl = tst

  tst2 = sapply(lmcl,function(x) x = is.null(x[2]$url_link))

  wlinks = which(tst2 == T)

  if (length(wlinks) > 0) lmcl = lmcl[-wlinks]

  m = NULL
  for (i in 1:length(lmcl)) m[[i]] = lmcl[[i]]$enrichMat_Ordered$filename[1]

  names(lmcl) = m

  nam = names(lmcl)

  lmcl

}
### i think this worked
mylinks = function(obj = out.DESeq3$results, fdr = 0.1, sig = 1){
  #source("/home/ggiaever/RProjects/2019_April_barseq_pipeline/May16_original_functionsGO.R")
  # source("/home/ggiaever/RProjects/2019_April_barseq_pipeline/2019_may20_gg_concor.r")
  #May16_original_functionsGO.R

  myscore = function(mat,column, sig = 1){
    df = data.frame(score = mat[,column],stringsAsFactors = F)
    df$gene = rownames(mat)
    rownames(df) = df$gene
    df$index=0
    wdf = which(df$score >= sig)
    df$index[wdf]=1
    df = df[,c('index','score','gene')]
    df
  }

  workdir = getwd()

  lres = lapply(obj,myres)

  nam0 = names(lres)

  lscore = lapply(lres,myscore,column = 2,sig=sig)

  score3 = sapply(lscore,function(x) x = all(x[,"index"]==0))

  w = which(score3 == T)

  if(length(w) > 0) lscore = lscore[-w]

  nam0 = names(lscore)

  tst = NULL

  for (i in 1:length(lscore)) tst[[i]] = myrun_go_enrichment2(fdrThresh = fdr,
    curr_exp = names(lscore)[i],score=lscore[[i]])
  dir = getwd()

  if(dir != workdir)  setwd(workdir)

  lmcl = tst

  tst2 = sapply(lmcl,function(x) x = is.null(x[2]$url_link))

  wlinks = which(tst2 == T)

  if (length(wlinks) > 0) lmcl = lmcl[-wlinks]

  m = NULL
  for (i in 1:length(lmcl)) m[[i]] = lmcl[[i]]$enrichMat_Ordered$filename[1]

  names(lmcl) = m

  nam = names(lmcl)

  lmcl
}

########
myhref1 = function(link,nameoflink){
  link = paste("<a href=",link,">",nameoflink,"</a>",sep="")
  #w = which(names(df)==links)
  #df = df[,-w,drop=F]
  link
}
########
mycompendium = function(lst,n = 1, fdr = 0.1){
  lenrich = lapply(lst,function(x) x = x[[1]])

  wfdr = sapply(lenrich,function(x) x = sum(x$FDR < fdr))

  lenrich2 = lapply(lenrich,myconcorqueryTocompendium, n=n)

  ayl = do.call(rbind,lenrich2)
  ayl$condition = rownames(ayl)
  ayl
}
########
########putting the two tables together
mygotable = function(lst){
  myhref1 = function(link,nameoflink){
    link = paste("<a href=",link,">",nameoflink,"</a>",sep="")
    link
  }
  links = sapply(lst,function(x) x = unlist(x[2],use.names = F))

  names(links) = names(lst)

  links = myhref1(link = links,nameoflink = names(links))

  dinks = data.frame(condition = names(lst), GO_enrichment = links,stringsAsFactors = F)

  dinks
}
########putting gotable and ayl table together works
myayl = function(gotable,concordtable,nam0){

  m = match(gotable[,"condition"],concordtable[,"condition"])

  links = cbind(gotable,concordtable[m,1:4])

  nam = setdiff(nam0,links$condition)

  lens = length(nam)

  if(length(nam) > 0){

    m = matrix(rep(links[1,],lens),byrow = T,nrow = lens)
    colnames(m)=names(links)
    dinks = data.frame(m,stringsAsFactors = F)
    dinks$condition = nam
    dinks$GO_enrichment = "no GO enrichment"
    rownames(dinks) = dinks$condition
    dinks[,3:6] = ""
    dinks = rbind(dinks,links)
    dinks = dinks[nam0,]
  }

  if(length(nam) == 0) dinks = links
    dinks
  }

mygoconcordwrapper = function(obj = out.DESeq3$results,fdr = 0.1, sig = 1){

  enrich = mylinks(obj = obj,  fdr = fdr, sig = sig)
  nam0 = names(obj)
  if(!is.null(enrich)){
  gotable = mygotable(enrich)

  concordtable = mycompendium(enrich)

  ayl = myayl(gotable = gotable,concordtable = concordtable,nam0 = nam0)

  write.table(ayl, file = "figures/AYL_concord.txt",row.names = T,sep = "\t", quote = F)
  ayl
  }

}

########
########
########

########



#### worked for all HIP data
#### Warning message:In max(nonEnrichInfo$maxOverlapGeneScore) :no non-missing arguments to max; returning -Inf
####
myconcorqueryTocompendium =
  function(enrichMat.query, useFDR=T, n = 1){
    workdir = getwd()
    psgtc = readxl::read_xlsx("/home/ggiaever/RProjects/2019_April_barseq_pipeline/2019_may19_SGTC_IDtoFILEName.xlsx")

    load("/home/ggiaever/RProjects/2019_April_barseq_pipeline/enrichSummary-GOBPv4EnrichMinOverlap1-negLog10FDR.rda")

    concor.queryToCompendium <- function(enrichMat.query, enrichSummary.comp = enrichSumm, useFDR=T) {
      if (useFDR) {
        queryVal <- enrichMat.query$FDR
      }
      else {
        # use the P values instead
        queryVal <- enrichMat.query$P
      }
      queryVal <- -log10(queryVal)
      names(queryVal) <- enrichMat.query$term
      queryVal <- queryVal[!is.na(queryVal)]

      apply(enrichSummary.comp, 2, function(otherVal) {
        uni <- intersect(names(queryVal), rownames(enrichSummary.comp)[!is.na(otherVal)])
        queryI <- match(uni, names(queryVal))
        otherI <- match(uni, rownames(enrichSummary.comp))
        concor(cbind(queryVal[queryI], otherVal[otherI])) })
    }

    tmp = concor.queryToCompendium(enrichMat.query, enrichSummary.comp = enrichSumm)

    tmp = tmp[order(tmp,decreasing = T)]
    tmp = round(tmp,2)
    tmp = tmp[1:n]

    m = match(names(tmp),psgtc$filename)

    tab = table(!is.na(m))

    if(tab) names(tmp)=psgtc$Screen.ID[m]

    link = "http://chemogenomics.pharmacy.ubc.ca/hiphop/enmaphg.php?d=MPv12.110818eP3SDfilter&exp="
    #####link to SGTC_ID.BP
    gink = paste0(link,names(tmp),".BP")

    fink = "http://chemogenomics.pharmacy.ubc.ca/hiphop/index.php?f=&q="
    #####link to SGTC_ID need just first bit
    xink = paste0(fink,names(tmp)) #,"&d=0&pv="

    df = data.frame(concord = tmp,compound = psgtc$name[m],SGTC_ID = psgtc$Screen.ID[m],
      link = xink,go = gink, present = psgtc$go[m], stringsAsFactors = F)
#df names = concord,compound, SGTC_ID,link,go,present
    df = df[,-6]
    myhref1 = function(link,nameoflink){
      link = paste0("<a href=",link)
      link = paste0(link,">")
      link2 = paste0(nameoflink,"</a>")
      link3 = paste0(link,link2)
      link3
    }
    df$SGTC_ID.BP = paste0(df$SGTC_ID,".BP")
    link = paste0("<a href=",df$go)

    link = paste0(link,">")
    link2 = paste0(df$SGTC_ID.BP,"</a>")
    link3 = paste0(link,link2)
    df$go = link3

    df$link = myhref1(link = df$link,nameoflink = df$SGTC_ID)
    #df$SGTC_ID.BP = paste0(df$SGTC_ID,".BP")



    #df$go = myhref1(link = df$go,nameoflink = df$SGTC.BP)
    w = which(names(df) %in% c("SGTC_ID","SGTC_ID.BP"))
    df = df[,-w]
    w = which(names(df) == "go")
    names(df)[w]="SGTC_ID.BP"
    w = which(names(df) == "link")
    names(df)[w]="SGTC_ID"
    #df names = concord,compound,link,go
    df
  }
########
########
myrun_go_enrichment2 =
  function (fdrThresh = 0.1, curr_exp, score){
    SupportFiles = "/home/common/labusers_enrich/SupportFiles"
    workdir = getwd()
    outdir = workdir
    Sys.setenv(PATH = "/usr/local/sbin:/usr/local/bin:/usr/bin:/usr/sbin:/sbin:/bin:/software/local/mcl/14-137/bin")
    xdir = system("echo 20`date +%y-%m-%d-%H-%M-%S`", intern = T)
    expdir = file.path(outdir, xdir)
    dir.create(expdir)
    fdrThresh = as.numeric(fdrThresh)
    TMP_ANALYSIS = expdir
    GOgeneSets_dir = SupportFiles
    setwd(expdir)
    expdir_base = basename(expdir)
    tstring = expdir_base
    bp_file = file.path(GOgeneSets_dir, "BP_geneNames_current.gmt")
    scoreMat = score
    scoreMat_stored = scoreMat
    queryGenes.mn <- sort(unique(scoreMat$gene[which(scoreMat$index >= 1)]))
    uniGenes.mn <- sort(unique(scoreMat$gene[!is.na(scoreMat$score)]))
    bp <- hiphop:::readSigs(bp_file)
    enrichMat.mn <- hiphop:::hyperG(querySet = queryGenes.mn, geneSets = bp,
      uni = uniGenes.mn, scoreMat = scoreMat, minSetSize = 5,
      maxSetSize = 300, uniSize = NA)
    queryGeneSets = list()
    queryGeneSets[[curr_exp]] = queryGenes.mn
    enrichMat.mn$filename <- curr_exp
    enrichMat_Ordered = enrichMat.mn[with(enrichMat.mn, order(FDR,
      -foldEnrichment)), ]
    scoreMat <- scoreMat[order(scoreMat$score, decreasing = T),
      ]
    scoreMat <- scoreMat[match(uniGenes.mn, scoreMat$gene), "score",
      drop = F]
    rownames(scoreMat) <- uniGenes.mn
    colnames(scoreMat) <- curr_exp
    head(scoreMat, 3)

    nonEnrichMat.mn <- hiphop:::genesNotInEnrichedTerm(queryGeneSets,
      enrichMat.mn, scoreMat, bp$NONSPECIFIC.TERMS, fdrThresh)
    bp_stored = bp
    bp <- lapply(bp, intersect, uniGenes.mn)
    lens <- sapply(bp, length)
    bp <- bp[lens >= 5 & lens <= 300]
    xgmml_file2 = paste(curr_exp, "_", tstring, "_fdrThresh_",
      as.character(fdrThresh), "_MCL_HyperG_withNonEnrich.xgmml",
      sep = "")
    xgmml_file2_stripped = gsub(pattern = ".xgmml$", replacement = "",
      x = xgmml_file2)
    q = myclusterEnrich(enrichInfo = enrichMat.mn, geneSets = bp,
      outFile = xgmml_file2, fdrThresh = fdrThresh, overlapThresh = 0.5,
      nonEnrichInfo = nonEnrichMat.mn, barModGenes = NULL,
      scoreName = "score", plotForEachEnrichedTerm = T, goTable = NULL)
    labusers_xgmml_dir = "/home/common/labusers_enrich"
    xgmml_file2_copy = file.path(labusers_xgmml_dir, xgmml_file2)
    #XML::saveXML(doc = q[["xgmml"]], file = xgmml_file2_copy)
    # url_link = paste("http://chemogenomics.pharmacy.ubc.ca/ycs/enmaphg.php?d=labusers_xgmml&exp=",
    #xgmml_file2_stripped, sep = "")

    if(!is.null(q[["xgmml"]])) XML::saveXML(doc = q[["xgmml"]], file = xgmml_file2_copy)
    url_link = paste("http://chemogenomics.pharmacy.ubc.ca/ycs/enmaphg.php?d=labusers_xgmml&exp=",
      xgmml_file2_stripped, sep = "")

    if(!is.null(q[["xgmml"]])) cat(expdir, sep = "\n")
    if(!is.null(q[["xgmml"]])) cat(url_link, sep = "\n")
    #utils::browseURL(url_link)

    #if(!is.null(q[["xgmml"]])) utils::browseURL(url_link)

    setwd(workdir)
    remove_expdir_command = paste("rm -rf", expdir)
    system(remove_expdir_command)
    ####gg added july 25
    if(!is.null(q[["xgmml"]]))
      # return(list(enrichMat_Ordered = enrichMat_Ordered, url_link = url_link,
      #   mcl_output = q))
      ####gg added august 3
      return(list(enrichMat_Ordered = enrichMat_Ordered, url_link = url_link))
  }

#####
myclusterEnrich = function (enrichInfo, geneSets, outFile, fdrThresh = 0.1, overlapThresh = 0.5,
  nonEnrichInfo = NULL, barModGenes = NULL, scoreName = "Fitness defect score",
  plotForEachEnrichedTerm = F, goTable = NULL)
{
  require(XML)
  xgmmlNS <- c(r = "http://www.cs.rpi.edu/XGMML")
  nodeSizeRange <- c(10, 40)
  prunedCol <- "#BEBEBE"
  labelWidth <- 20
  edgeWidthRange <- c(1, 5)
  overlapCoeffRange <- c(overlapThresh, 1)
  graphLabel <- strsplit(basename(outFile), "\\.")[[1]]
  graphLabel <- paste(graphLabel[-length(graphLabel)], collapse = ".")
  graphNode <- newXMLNode("graph", attrs = c(label = graphLabel),
    namespace = c(xgmml = xgmmlNS), addFinalizer = T)
  if (!is.null(nonEnrichInfo)) {
    nonEnrichInfo$maxOverlapGeneScore <- round(nonEnrichInfo$maxOverlapGeneScore,
      2)
    nonEnrichInfo$geneSetFraction <- round(nonEnrichInfo$geneSetFraction *
        100, 1)
    if (is.null(nonEnrichInfo$GOTerm.nTestedGenes)) {
      lens <- sapply(geneSets, length)
      i <- match(nonEnrichInfo$term, names(lens))
      nonEnrichInfo$GOTerm.nTestedGenes <- lens[i]
    }
    tmp <- strsplit(nonEnrichInfo$overlapGenes, "\\|")
    w <- which(is.na(tmp))
    if(length(w)>0) tmp = tmp[-w]
    if (is.null(nonEnrichInfo$unenrichedGenes)) {
      nonEnrichInfo$overlapGenes <- sapply(tmp, paste,
        collapse = "| ")
    }
    else {
      unEnriched <- strsplit(nonEnrichInfo$unenrichedGenes,
        "\\|")
      tmp.mod <- sapply(1:length(tmp), function(termI) {
        vec <- tmp[[termI]]
        i <- match(unEnriched[[termI]], tmp[[termI]])
        vec[i] <- paste("<b>", vec[i], "</b>", sep = "")
        paste(vec, collapse = "| ")
      })

      nonEnrichInfo$overlapGenes <- tmp.mod
    }
    newNode <- newXMLNode("node", attrs = c(id = "unenrichedTable",
      label = "unenrichedTable"), parent = graphNode, addFinalizer = T)
    termNode <- newXMLNode("att", attrs = c(type = "list",
      name = "terms"), parent = newNode, addFinalizer = T)
    geneNode <- newXMLNode("att", attrs = c(type = "list",
      name = "genes"), parent = newNode, addFinalizer = T)
    scoreNode <- newXMLNode("att", attrs = c(type = "list",
      name = "maxScores"), parent = newNode, addFinalizer = T)
    termSizeNode <- newXMLNode("att", attrs = c(type = "list",
      name = "termSizes"), parent = newNode, addFinalizer = T)
    termFractionNode <- newXMLNode("att", attrs = c(type = "list",
      name = "termFractions"), parent = newNode, addFinalizer = T)
    apply(nonEnrichInfo, 1, function(gene) {
      gene <- as.list(gene)
      newXMLNode("att", attrs = c(type = "string", name = "term",
        value = gene$term), parent = termNode)
      newXMLNode("att", attrs = c(type = "string", name = "geneList",
        value = gene$overlapGenes), parent = geneNode)
      newXMLNode("att", attrs = c(type = "real", name = "maxScore",
        value = gene$maxOverlapGeneScore), parent = scoreNode)
      newXMLNode("att", attrs = c(type = "string", name = "termSize",
        value = gene$GOTerm.nTestedGenes), parent = termSizeNode)
      newXMLNode("att", attrs = c(type = "string", name = "termFraction",
        value = gene$geneSetFraction), parent = termFractionNode)
    })
    rm(newNode, termNode, geneNode, scoreNode, termSizeNode,
      termFractionNode)
    if (is.null(enrichInfo)) {
      doc <- xmlDoc(graphNode)
      saveXML(doc, file = outFile)
      return()
    }
  }
  ###### key step to exit
  enrichInfo <- enrichInfo[enrichInfo$FDR <= fdrThresh, , drop = F]
  if (nrow(enrichInfo) == 0) {
    print("No enriched terms to cluster.")
    doc <- xmlDoc(graphNode)
    saveXML(doc, file = outFile)
    return()
  }
  enrichInfo$formattedLabel <- sapply(enrichInfo$term, function(curLabel) {
    curLabel <- strwrap(curLabel, labelWidth)
    paste(curLabel, collapse = "\n")
  })
  i <- match(enrichInfo$term, names(geneSets))
  if (any(is.na(i))) {
    stop("Could not find gene sets for ", sum(is.na(i)),
      " enriched terms.")
  }
  geneSets <- geneSets[i]
  if (is.null(enrichInfo$GOTerm.nTestedGenes)) {
    enrichInfo$GOTerm.nTestedGenes <- sapply(geneSets, length)
  }
  tmpSize <- -log10(enrichInfo$FDR)
  maxVal <- max(tmpSize[!is.infinite(tmpSize)])
  tmpSize[is.infinite(tmpSize)] <- maxVal + 2
  gsSizeRange <- range(tmpSize)
  if (gsSizeRange[1] == gsSizeRange[2]) {
    gsSizeRange[1] <- -log10(fdrThresh)
    gsSizeRange[2] <- gsSizeRange[2] + 1
  }
  tmpSize <- (tmpSize - gsSizeRange[1])/(gsSizeRange[2] - gsSizeRange[1])
  tmpSize <- nodeSizeRange[1] + tmpSize * (nodeSizeRange[2] -
      nodeSizeRange[1])
  enrichInfo$size <- round(tmpSize, 2)
  if (nrow(enrichInfo) == 1) {
    enrichInfo$cluster <- CLUST.COL[1]
    edgeMat <- NULL
  }
  else {
    pairI <- getUniquePairs(length(geneSets))
    distVal <- apply(pairI, 1, function(onePair) {
      hiphop:::overlapCoeff(geneSets[onePair])
    })
    distVal[distVal < overlapThresh] <- 0
    edgeMat <- data.frame(nodeA = pairI[, 1], nodeB = pairI[,
      2], coeff = distVal)
    enrichInfo$cluster <- prunedCol
    if (is.null(enrichInfo$pruneOutcome)) {
      termI <- 1:nrow(enrichInfo)
    }
    else {
      termI <- which(enrichInfo$pruneOutcome == enrichInfo$term)
    }
    if (length(termI) == 1) {
      enrichInfo$cluster[termI] <- CLUST.COL[1]
    }
    else {
      i <- which((edgeMat$nodeA %in% termI) & (edgeMat$nodeB %in%
          termI))
      clusters <- mclWrapper(edgeMat[i, , drop = F], dirname(outFile))
      mcl.in = edgeMat[i, , drop = F]
      if(is.null(mcl.in)) print("edge NULL")
      mcl.out = clusters
      clusters <- lapply(clusters, as.numeric)
      if (length(clusters) > length(CLUST.COL)) {
        stop("Need more cluster colours!")
      }
      lens <- sapply(clusters, length)
      clusters <- data.frame(id = unlist(clusters), cluster = CLUST.COL[rep(1:length(clusters),
        lens)], stringsAsFactors = F)
      enrichInfo$cluster[clusters$id] <- clusters$cluster
    }
    edgeMat <- edgeMat[edgeMat$coeff > 0, , drop = F]
    if (nrow(edgeMat) > 0) {
      edgeMat$size <- (edgeMat$coeff - overlapCoeffRange[1])/(overlapCoeffRange[2] -
          overlapCoeffRange[1])
      edgeMat$size <- edgeWidthRange[1] + edgeMat$size *
        (edgeWidthRange[2] - edgeWidthRange[1])
      edgeMat$coeff <- round(edgeMat$coeff, 2)
      edgeMat$size <- round(edgeMat$size, 2)
    }
    else {
      edgeMat <- NULL
    }
  }
  otherI <- order(enrichInfo$cluster)
  otherI <- otherI[order(enrichInfo$FDR[otherI])]
  termI <- which(enrichInfo$cluster[otherI] != prunedCol)
  if (length(termI) < length(otherI)) {
    otherI <- c(otherI[termI], otherI[-termI])
  }
  enrichInfo$id <- 1:nrow(enrichInfo)
  enrichInfo <- enrichInfo[otherI, , drop = F]
  enrichInfo$geneSetFraction <- round(enrichInfo$geneSetFraction *
      100, 1)
  enrichInfo$querySetFraction <- round(enrichInfo$querySetFraction *
      100, 1)
  allNodes <- sapply(1:nrow(enrichInfo), function(nodeI) {
    newNode <- newXMLNode("node", attrs = c(id = enrichInfo$id[nodeI],
      label = enrichInfo$term[nodeI]), parent = graphNode,
      addFinalizer = T)
    newXMLNode("att", attrs = c(type = "string", name = "formattedLabel",
      value = enrichInfo$formattedLabel[nodeI]), parent = newNode)
    newXMLNode("att", attrs = c(type = "string", name = "cluster",
      value = enrichInfo$cluster[nodeI]), parent = newNode)
    newXMLNode("att", attrs = c(type = "real", name = "FDR",
      value = enrichInfo$FDR[nodeI]), parent = newNode)
    newXMLNode("att", attrs = c(type = "real", name = "geneSetFraction",
      value = enrichInfo$geneSetFraction[nodeI]), parent = newNode)
    newXMLNode("att", attrs = c(type = "real", name = "querySetFraction",
      value = enrichInfo$querySetFraction[nodeI]), parent = newNode)
    newXMLNode("att", attrs = c(type = "integer", name = "nGenes",
      value = enrichInfo$GOTerm.nTestedGenes[nodeI]), parent = newNode)
    newNode
  })
  if (!is.null(goTable)) {
    toDoI <- match(enrichInfo$term, goTable$term)
    goid <- goTable$id[toDoI[!is.na(toDoI)]]
    toDoI <- which(!is.na(toDoI))
    if (length(toDoI) > 0) {
      sapply(1:length(toDoI), function(nodeI) {
        newXMLNode("att", attrs = c(type = "string",
          name = "GOID", value = goid[nodeI]), parent = allNodes[[toDoI[nodeI]]])
      })
    }
  }
  if (plotForEachEnrichedTerm) {
    maxScore <- max(enrichInfo$maxOverlapGeneScore,na.rm=T, na.rm = T)
    if (!is.null(nonEnrichInfo)) {
      maxScore <- max(maxScore, max(nonEnrichInfo$maxOverlapGeneScore, na.rm = T))
    }
    tmp <- sapply(enrichInfo$overlapGenes, mygenOverlapGenePlot.gChart,
      c(0, maxScore), barModGenes, scoreLabel = scoreName)
    enrichInfo$image <- tmp[3, ]
    sapply(1:nrow(enrichInfo), function(nodeI) {
      newXMLNode("att", attrs = c(type = "string", name = "image",
        value = enrichInfo$image[nodeI]), parent = allNodes[[nodeI]])
      newXMLNode("att", attrs = c(type = "real", name = "size",
        value = enrichInfo$size[nodeI]), parent = allNodes[[nodeI]])
      newXMLNode("graphics", attrs = c(type = "ELLIPSE"),
        parent = allNodes[[nodeI]])
    })
  }
  else {
    nodeInfo.plot <- t(clusterPlots(enrichInfo, scoreName,
      barModGenes))
    colnames(nodeInfo.plot) <- c("w", "h", "image")
    nodeInfo.plot <- as.data.frame(nodeInfo.plot)
    nodeInfo.plot$w <- round(as.numeric(as.character(nodeInfo.plot$w)),
      2)
    nodeInfo.plot$h <- round(as.numeric(as.character(nodeInfo.plot$h)),
      2)
    nodeInfo.plot$image <- as.character(nodeInfo.plot$image)
    nodeInfo.plot$cluster <- rownames(nodeInfo.plot)
    nodeInfo.plot$id <- substring(nodeInfo.plot$cluster,
      2)
    sapply(1:nrow(nodeInfo.plot), function(nodeI) {
      newNode <- newXMLNode("node", attrs = c(id = nodeInfo.plot$id[nodeI],
        label = nodeInfo.plot$id[nodeI]), parent = graphNode,
        addFinalizer = T)
      newXMLNode("att", attrs = c(type = "string", name = "image",
        value = nodeInfo.plot$image[nodeI]), parent = newNode)
      newXMLNode("att", attrs = c(type = "string", name = "cluster",
        value = nodeInfo.plot$cluster[nodeI]), parent = newNode)
      newXMLNode("att", attrs = c(type = "boolean", name = "isPlot",
        value = "true"), parent = newNode)
      newXMLNode("att", attrs = c(type = "real", name = "width",
        value = nodeInfo.plot$w[nodeI]), parent = newNode)
      newXMLNode("att", attrs = c(type = "real", name = "height",
        value = nodeInfo.plot$h[nodeI]), parent = newNode)
    })
    tmp <- strsplit(enrichInfo$overlapGenes, "\\|")
    enrichInfo$overlapGenes <- sapply(tmp, paste, collapse = "| ")
    sapply(1:nrow(enrichInfo), function(nodeI) {
      newXMLNode("att", attrs = c(type = "string", name = "overlapGenes",
        value = enrichInfo$overlapGenes[nodeI]), parent = allNodes[[nodeI]])
      newXMLNode("att", attrs = c(type = "real", name = "width",
        value = enrichInfo$size[nodeI]), parent = allNodes[[nodeI]])
      newXMLNode("att", attrs = c(type = "real", name = "height",
        value = enrichInfo$size[nodeI]), parent = allNodes[[nodeI]])
    })
  }
  toDoI <- which(enrichInfo$cluster != prunedCol)
  toDoI <- toDoI[enrichInfo$FDR[toDoI] == min(enrichInfo$FDR[toDoI])]
  sapply(toDoI, function(nodeI) {
    newXMLNode("att", attrs = c(type = "integer", name = "isMostSig",
      value = 5), parent = allNodes[[nodeI]])
  })
  rm(allNodes)
  if (is.null(edgeMat)) print("edgeMat is NULL")
  if (!is.null(edgeMat)) {
    apply(edgeMat, 1, function(edge) {
      edge <- as.list(edge)
      newNode <- newXMLNode("edge", attrs = c(source = edge$nodeA,
        target = edge$nodeB, label = paste(names(geneSets)[edge$nodeA],
          "(overlap)", names(geneSets)[edge$nodeB])),
        parent = graphNode, addFinalizer = T)
      newXMLNode("att", attrs = c(type = "real", name = "width",
        value = edge$size), parent = newNode)
      newXMLNode("att", attrs = c(type = "real", name = "overlapCoeff",
        value = edge$coeff), parent = newNode)
    })
  }
  doc <- xmlDoc(graphNode)
  saveXML(doc, file = outFile)
  #mau_output = list(xgmml = doc, xgmml_file = outFile, mcl.in = mcl.in,
  #mcl.out = mcl.out)
  #### added

  if(is.null(edgeMat)) mcl.out =  " "
  #### added
  mau_output = list(xgmml = doc, xgmml_file = outFile,
    mcl.out = mcl.out)
  return(mau_output)
}
# generates overlap gene (i.e. those that drive enrichment) barplots for each cluster
# enrichInfo - dataframe with info on enriched terms (one per row), with the following columns:
#              overlapGenes (genes that drive the enrichment of the term), cluster (cluster of the term)
# scoreLabel - score label to use in the barplots
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly in the
#              barplots, should they be in the top overlap genes
# scoreInfo - if provided (i.e. not NULL), a data.frame of scores with the following columns:
#              gene (gene IDs) and score (gene scores)
# maxGenes - maximum number of top overlap genes to show in a barplot
# posScores - if TRUE, enrichment is amongst genes with positive scores; if FALSE, negative scores
# plotCol - the colour of the bars, in hexidecimal format without the '#' character
# RETURNS a matrix of plot info where each column corresponds to a different cluster bar plot,
# and the rows correspond to plot width (1st), plot height (2nd) hand the google chart URL for the plot (3rd)
clusterPlots <- function(enrichInfo, scoreLabel, barModGenes=NULL, scoreInfo=NULL, maxGenes=10, posScores=T, plotCol="BEBEBE") {
  barWidth <- 10

  overlapGenes <- strsplit(enrichInfo$overlapGenes, "\\|")

  if (is.null(scoreInfo)) {
    oGenes <- unique(unlist(overlapGenes))
    genes <- strsplit(oGenes, "\\(")
    scores <- sapply(genes, function(vec) { vec[length(vec)] })
    genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
    scores <- unlist(strsplit(scores, ")"))
    scoreInfo <- data.frame(gene=genes, score=as.numeric(scores), geneStr=oGenes, stringsAsFactors=F)
  }
  else {
    scoreInfo$geneStr <- scoreInfo$gene
  }

  scoreInfo <- scoreInfo[order(scoreInfo$score, decreasing=posScores), ]
  uniGenes <- unique(scoreInfo$gene)
  i <- match(uniGenes, scoreInfo$gene)
  scoreInfo <- scoreInfo[i, ]

  clusters <- split(1:nrow(enrichInfo), enrichInfo$cluster)

  plotData <- lapply(clusters, function(clustI) {
    clusterGenes <- table(unlist(overlapGenes[clustI]))
    clusterGenes <- clusterGenes/length(clustI) * 100

    i <- match(names(clusterGenes), scoreInfo$geneStr)
    clusterGenes <- cbind(clusterGenes, scoreInfo$score[i])
    colnames(clusterGenes) <- c("% of gene sets", "score")
    rownames(clusterGenes) <- scoreInfo$gene[i]

    # sort the genes by score, then by % of gene sets they are in
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]
    orderI <- order(clusterGenes[, 1], decreasing=T)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (nrow(clusterGenes) > maxGenes) {
      clusterGenes <- clusterGenes[1:maxGenes, , drop=F]
    }

    # re-order by score
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (is.null(barModGenes)) {
      return(clusterGenes)
    }

    cbind(clusterGenes, rownames(clusterGenes) %in% barModGenes) })

  if (posScores) {
    scoreRange <- c(0, ceiling(max(scoreInfo$score, na.rm=T)))
  }
  else {
    scoreRange <- c(floor(min(scoreInfo$score, na.rm=T)), 0)
  }

  sapply(plotData, genLeadingEdgePlot.gChart, plotCol, scoreRange, barWidth, scoreLabel)
}

# generate a barplot of the common leading edge genes using google charts (can be visualized
# with the html img tag), bar length corresponds to score
# leadInfo - matrix/dataframe of leading edge genes; rownames = gene IDs, column 1 = % of gene sets, column 2 = score
#          - if a 3rd column is provided, it should provide TRUE/FALSE indicating whether or not the gene
#            should be marked
# plotCol - the colour of the bars; hexidecimal format without the '#' character
# scoreRange - a vector of the range of scores in the profile; 1st value = min,  2nd value = max
# barWidth - bar width in the plot
# scoreLabel - label for the score axis
# RETURNS a vector of plot info: plot width, plot height and the google chart URL for the plot
#
mclWrapper = function (simMat, runDir, inflation = 2)
{
  inFile <- paste(runDir, "mcl.in", sep = "/")
  write.table(simMat, inFile, row.names = F, col.names = F,
    sep = "\t", quote = F)
  outFile <- paste(runDir, "mcl.out", sep = "/")
  cmdStr <- paste("mcl", inFile, "--abc -I", inflation, "-o",
    outFile, "&>/dev/null")
  system(cmdStr)
  clusters <- scan(outFile, what = character(), sep = "\n")
  strsplit(clusters, "\t")
}

mygenLeadingEdgePlot.gChart =
  function (leadInfo, plotCol, scoreRange, barWidth, scoreLabel = "Sensitivity")
  {
    dataRange <- scoreRange[2] - scoreRange[1]
    barLens <- round((leadInfo[, 2] - scoreRange[1])/dataRange *
        100)
    if (dataRange <= 1) {
      stepSize <- 0.5
    }
    else if (dataRange <= 5) {
      stepSize <- 2
    }
    else if (dataRange <= 20) {
      stepSize <- 5
    }
    else if (dataRange <= 50) {
      stepSize <- 10
    }
    else if (dataRange <= 100) {
      stepSize <- 20
    }
    else if (dataRange <= 500) {
      stepSize <- 100
    }
    else {
      stepSize <- 250
    }
    scoreLabel <- unlist(strsplit(scoreLabel, ""))
    scoreLabel[scoreLabel == " "] <- "+"
    scoreLabel <- paste(scoreLabel, collapse = "")
    tmpStep <- 100/length(barLens)
    labelPos <- round(seq(tmpStep/2, 100, by = tmpStep))
    w <- 150
    h <- barWidth * length(barLens) + 50
    if (scoreRange[1] < 0) {
      zeroLineStr <- paste("&chp=", round(abs(scoreRange[1])/dataRange,
        2), sep = "")
    }
    else {
      zeroLineStr <- ""
    }
    if (ncol(leadInfo) > 2 && any(leadInfo[, 3] == 1)) {
      barModStr <- paste("o,000000,", which(leadInfo[, 3] ==
          1) - 1, ",-1,5", sep = "")
      barModStr <- paste("&chm=", paste(barModStr, collapse = "|"),
        sep = "")
    }
    else {
      barModStr <- ""
    }
    c(w, h, paste("http://chart.apis.google.com/chart?chxt=x,x,y&chs=",
      w, "x", h, "&cht=bhg&chd=t:", paste(barLens, collapse = "|"),
      "&chco=", plotCol, "&chxl=1:|", scoreLabel, "|2:|", paste(rev(rownames(leadInfo)),
        collapse = "|"), "&chxp=1,50|2,", paste(labelPos,
          collapse = ","), "&chxr=0,", scoreRange[1], ",",
      scoreRange[2], ",", stepSize, "&chbh=", barWidth, ",1,0",
      zeroLineStr, barModStr, sep = ""))
  }
####
##### required functions
# computes the number of unique pairs given the number of items to consider
# maxVal - the maximum number of items
# RETURNS a 2-column matrix where each row contains a different pair, specified with item indices
getUniquePairs = function (maxVal)
{
  firstI <- rep(1:(maxVal - 1), (maxVal - 1):1)
  secondI <- sapply(2:maxVal, function(x) {
    x:maxVal
  })
  cbind(firstI, unlist(secondI))
}



CLUST.COL <- c("#FF00CC","#33CCFF", "#33CC00", "#9900FF", "#FF9900", "#FFFF00", "#FFCCFF", "#FF0000", "#006600", "#009999", "#CCCC00", "#993300", "#CC99CC", "#6699CC","#CCCCFF", "#FFCC99", "#9966FF", "#CC6600", "#CCFFFF", "#99CC00", "#FF99FF", "#0066FF", "#66FFCC", "#99CCFF", "#9999CC", "#CC9900", "#CC33FF", "#006699", "#F5DF16", "#B5185E", "#99FF00", "#00FFFF", "#990000", "#CC0000", "#33CCCC", "#CC6666", "#996600", "#9999FF", "#3366FF")

# generate a barplot of the top-scoring overlap genes (driving enrichment) using google charts (can be
# visualized with the html img tag), bar length corresponds to score
# oGeneStr - |-separated string of overlap genes, and after each gene, its score is provided in parentheses,
#          genes are sorted by score in decreasing order
# scoreRange - a vector of the range of scores in the profile; 1st value = min,  2nd value = max
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly,
#          should they be in the top overlap genes
# barWidth - bar width in the plot
# maxGenes - maximum number of top overlap genes to return
# plotCol - the colour of the bars; hexidecimal format without the '#' character
# scoreLabel - label for the score axis
# RETURNS a vector of plot info: plot width, plot height and the google chart URL for the plot
mygenOverlapGenePlot.gChart <- function(oGeneStr, scoreRange, barModGenes=NULL, barWidth=10, maxGenes=10, plotCol="BEBEBE", scoreLabel="Sensitivity") {
  oGenes <- unlist(strsplit(oGeneStr, "\\|"))
  genes <- strsplit(oGenes, "\\(")
  scores <- sapply(genes, function(vec) { vec[length(vec)] })
  genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
  scores <- unlist(strsplit(scores, ")"))

  scoreMat <- data.frame(percent=rep(NA, length(genes)), score=as.numeric(scores), stringsAsFactors=F)
  rownames(scoreMat) <- genes
  if (nrow(scoreMat) > maxGenes) {
    scoreMat <- scoreMat[1:maxGenes, ]
  }

  if (!is.null(barModGenes)) {
    scoreMat$mod <- as.numeric(scoreMat$gene %in% barModGenes)
  }

  mygenLeadingEdgePlot.gChart(scoreMat, plotCol, scoreRange, barWidth, scoreLabel)
}


# generates overlap gene (i.e. those that drive enrichment) barplots for each cluster
# enrichInfo - dataframe with info on enriched terms (one per row), with the following columns:
#              overlapGenes (genes that drive the enrichment of the term), cluster (cluster of the term)
# scoreLabel - score label to use in the barplots
# barModGenes - if provided (i.e. not NULL) a vector of genes that should be marked distinctly in the
#              barplots, should they be in the top overlap genes
# scoreInfo - if provided (i.e. not NULL), a data.frame of scores with the following columns:
#              gene (gene IDs) and score (gene scores)
# maxGenes - maximum number of top overlap genes to show in a barplot
# posScores - if TRUE, enrichment is amongst genes with positive scores; if FALSE, negative scores
# plotCol - the colour of the bars, in hexidecimal format without the '#' character
# RETURNS a matrix of plot info where each column corresponds to a different cluster bar plot,
# and the rows correspond to plot width (1st), plot height (2nd) hand the google chart URL for the plot (3rd)
clusterPlots <- function(enrichInfo, scoreLabel, barModGenes=NULL, scoreInfo=NULL, maxGenes=10, posScores=T, plotCol="BEBEBE") {
  barWidth <- 10

  overlapGenes <- strsplit(enrichInfo$overlapGenes, "\\|")

  if (is.null(scoreInfo)) {
    oGenes <- unique(unlist(overlapGenes))
    genes <- strsplit(oGenes, "\\(")
    scores <- sapply(genes, function(vec) { vec[length(vec)] })
    genes <- sapply(genes, function(vec) { paste(vec[-length(vec)], collapse="(") })
    scores <- unlist(strsplit(scores, ")"))
    scoreInfo <- data.frame(gene=genes, score=as.numeric(scores), geneStr=oGenes, stringsAsFactors=F)
  }
  else {
    scoreInfo$geneStr <- scoreInfo$gene
  }

  scoreInfo <- scoreInfo[order(scoreInfo$score, decreasing=posScores), ]
  uniGenes <- unique(scoreInfo$gene)
  i <- match(uniGenes, scoreInfo$gene)
  scoreInfo <- scoreInfo[i, ]

  clusters <- split(1:nrow(enrichInfo), enrichInfo$cluster)

  plotData <- lapply(clusters, function(clustI) {
    clusterGenes <- table(unlist(overlapGenes[clustI]))
    clusterGenes <- clusterGenes/length(clustI) * 100

    i <- match(names(clusterGenes), scoreInfo$geneStr)
    clusterGenes <- cbind(clusterGenes, scoreInfo$score[i])
    colnames(clusterGenes) <- c("% of gene sets", "score")
    rownames(clusterGenes) <- scoreInfo$gene[i]

    # sort the genes by score, then by % of gene sets they are in
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]
    orderI <- order(clusterGenes[, 1], decreasing=T)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (nrow(clusterGenes) > maxGenes) {
      clusterGenes <- clusterGenes[1:maxGenes, , drop=F]
    }

    # re-order by score
    orderI <- order(clusterGenes[, 2], decreasing=posScores)
    clusterGenes <- clusterGenes[orderI, , drop=F]

    if (is.null(barModGenes)) {
      return(clusterGenes)
    }

    cbind(clusterGenes, rownames(clusterGenes) %in% barModGenes) })

  if (posScores) {
    scoreRange <- c(0, ceiling(max(scoreInfo$score, na.rm=T)))
  }
  else {
    scoreRange <- c(floor(min(scoreInfo$score, na.rm=T)), 0)
  }

  sapply(plotData, genLeadingEdgePlot.gChart, plotCol, scoreRange, barWidth, scoreLabel)
}

# generate a barplot of the common leading edge genes using google charts (can be visualized
# with the html img tag), bar length corresponds to score
# leadInfo - matrix/dataframe of leading edge genes; rownames = gene IDs, column 1 = % of gene sets, column 2 = score
#          - if a 3rd column is provided, it should provide TRUE/FALSE indicating whether or not the gene
#            should be marked
# plotCol - the colour of the bars; hexidecimal format without the '#' character
# scoreRange - a vector of the range of scores in the profile; 1st value = min,  2nd value = max
# barWidth - bar width in the plot
# scoreLabel - label for the score axis
# RETURNS a vector of plot info: plot width, plot height and the google chart URL for the plot
mygenLeadingEdgePlot.gChart <- function(leadInfo, plotCol, scoreRange, barWidth, scoreLabel="Sensitivity") {
  # express bar lengths as values in [0, 100]
  #
  # ggadded aug3
  max = max(scoreRange,na.rm = T)
  min = min(scoreRange,na.rm = T)

  dataRange <- max - min
  #dataRange <- scoreRange[2] - scoreRange[1]
  barLens <- round((leadInfo[, 2] - scoreRange[1])/dataRange * 100)

  # determine the score step size
  if (dataRange <= 1) {
    stepSize <- 0.5
  }
  else if (dataRange <= 5) {
    stepSize <- 2
  }
  else if (dataRange <= 20) {
    stepSize <- 5
  }
  else if (dataRange <= 50) {
    stepSize <- 10
  }
  else if (dataRange <= 100) {
    stepSize <- 20
  }
  else if (dataRange <= 500) {
    stepSize <- 100
  }
  else {
    stepSize <- 250
  }

  # replace any spaces in scoreLabel with a +
  scoreLabel <- unlist(strsplit(scoreLabel, ""))
  scoreLabel[scoreLabel == " "] <- "+"
  scoreLabel <- paste(scoreLabel, collapse="")

  # determine the positions of the gene labels
  tmpStep <- 100/length(barLens)
  labelPos <- round(seq(tmpStep/2, 100, by=tmpStep))

  # compute plot size
  w <- 150
  h <- barWidth*length(barLens) + 50

  # specify the zero line if have negative values
  if (scoreRange[1] < 0) {
    zeroLineStr <- paste("&chp=", round(abs(scoreRange[1])/dataRange, 2), sep="")
  }
  else {
    zeroLineStr <- ""
  }

  # if a 3rd column is in leadInfo, use it to determine which gene bars to mark
  if (ncol(leadInfo) > 2 && any(leadInfo[, 3] == 1)) {
    barModStr <- paste("o,000000,", which(leadInfo[,3]==1) - 1, ",-1,5", sep="")
    barModStr <- paste("&chm=", paste(barModStr, collapse="|"), sep="")
  }
  else {
    barModStr <- ""
  }

  c(w, h, paste("http://chart.apis.google.com/chart?chxt=x,x,y&chs=", w, "x", h,
    "&cht=bhg&chd=t:", paste(barLens, collapse="|"),
    "&chco=", plotCol,
    "&chxl=1:|", scoreLabel, "|2:|", paste(rev(rownames(leadInfo)), collapse="|"),
    "&chxp=1,50|2,", paste(labelPos, collapse=","),
    "&chxr=0,", scoreRange[1], ",", scoreRange[2], ",", stepSize,
    "&chbh=", barWidth, ",1,0",
    zeroLineStr, barModStr, sep=""))
}

