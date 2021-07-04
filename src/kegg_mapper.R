#################################################
#                                               #
#        KEGG MAPPER UTILITY FUNCTIONS          #
#                                               #
#################################################

### MAPPING DATA LOCATIONS (PART OF REPO, SHOULD UPDATE REGULARLY OR CREATE PIPELINE TO GENERATE AUTOMATICALLY)
KEGGGENE_TO_KEGGPATHWAY_MAP='etc/kegg_data/kegg_human_genes_to_pathways.csv'
GENESYMBOL_TO_KEGGGENE_MAP='etc/ncbi_data/gene_to_entrez_and_kegg.csv'

### OTHER PROTIGY NAMES THAT SHOULD BE GLOBALIZED
COMPARISON_LIST_FIELD_NAME = 'grp.comp'
FOLD_CHANGE_NAME_BASE = 'logFC'
SIGNIFICANCE_NAME_BASE = 'adj.P.Val'
DATA_OBJECT_NAME = 'data'
DATA_OUTPUT_OBJECT_NAME = 'output'
GENE_SYMBOL_NAME = 'id.mapped'  # POTENTIALLY TO BE REPLACED BY ENTREZ GENE ID LATER
ID_CONCAT_NAME = 'id.concat'

### OTHER MISC GLOBALS
JOINER_CHARACTER = '.'
KEGG_GENE_ACCESSION_FIELDNAME = 'keggGene'
DIRECTION_NAME_BASE = 'direction'
GLOBAL_DIRECTION_NAME_BASE = 'global_direction'
defaultKEGGorganismtoken = "hsa"
kurl<-"https://www.genome.jp/kegg-bin/mcolor_pathway"
kurlbase<-'http://www.genome.jp';

map_kegg_pathways <- function(protigy_parameters, protigy_results, output='data/output.html') {
  #######################################################################################
  #                                                                                     #
  #          MAIN KEGG MAPPING WORKFLOW                                                 #
  #                                                                                     #
  #                                                                                     #
  # param: protigy_parameters - global.param.imp from protigy data structure            #
  #                           - must contain 'grp.comp.all' list                        #
  #                                                                                     #
  # param: protigy_results    - global.results.imp from protigy data structure          #
  #                           - must contain 'data' list with 'output' dataframe        #
  #                                                                                     #
  #                                                                                     #
  #######################################################################################

  # Gather potential comparisons
  comparison_table <- gather_comparisons_and_make_column_names(protigy_parameters = protigy_parameters)
  
  # Map gene symbols to KEGG ids, replace original data table
  protigy_results <- annotate_genes_with_kegg_ids(protigy_results = protigy_results)
  
  # Determine number of comparisons and initialize html generation object
  n_comparisons <- dim(comparison_table)[1]
  html_generation_object<-list()

  # Load cached genes to pathway data
  kegg_genes_to_pathways <- read.csv(KEGGGENE_TO_KEGGPATHWAY_MAP, stringsAsFactors = FALSE, row.names = 1)
    
  # Main loop for each comparison in Protigy results
  for (j in 1:n_comparisons) {
    
    # Mark genes up/down/both based on statistics (set the filtering level!)
    protigy_results <- prepareKEGGhelperFunction(
      protigy_results_with_kegg = protigy_results,
      fcfield = comparison_table[['foldchange_name']][j],
      sigfield = comparison_table[['significance_name']][j],
      site_dirfield = comparison_table[['direction_name']][j],
      global_dirfield = comparison_table[['global_direction_name']][j]
    )
  
    # Select and merge these results with pw data
    required_data_for_pathway_computation <- 
      select_minimal_data_and_merge(
        protigy_results_with_directions = protigy_results,
        kegg_pathway_data = kegg_genes_to_pathways,
        fcfield = comparison_table[['foldchange_name']][j],
        sigfield = comparison_table[['significance_name']][j],
        site_dirfield = comparison_table[['direction_name']][j],
        global_dirfield = comparison_table[['global_direction_name']][j]
    )
    
    # Compute pathway enrichments
    pathway_enrichment_table<-findKEGGsetSizes(
      DF = required_data_for_pathway_computation,
      keggGeneField=KEGG_GENE_ACCESSION_FIELDNAME,
      AdjPValField=comparison_table[['significance_name']][j],
      AdjPValCut=0.1,
      alternative='greater',
      useBest=FALSE
    )

    # Assign to the html generation object
    html_generation_object$conditions[[j]]<-comparison_table[['comparison_name']][j]
    html_generation_object$datasets[[j]]<-required_data_for_pathway_computation
    html_generation_object$pwenrichments[[j]]<-pathway_enrichment_table
  }  

  # Generate HTML output document
  html<-generateHTMLdocFromExtendedKEGGplotObject(html_generation_object,output,sigcut=1)

}


gather_comparisons_and_make_column_names <- function(protigy_parameters) {
  #######################################################################################
  #                                                                                     #
  # function: gather_comparisons_and_make_column_names                                  #
  #   makes a dataframe containing all relevant comparisons and names the columns for   #
  #   use in determining significance and directionality                                #
  #                                                                                     #
  # param: protigy_parameters - global.param.imp from protigy data structure            #
  #                           - must contain 'grp.comp.all' list                        #
  #                                                                                     #
  # returns: dataframe with the following fields:                                       #
  #    comparison_name                                                                  #
  #    significance_name                                                                #
  #    foldchange_name                                                                  #
  #    direction_name                                                                   #
  #                                                                                     #
  #######################################################################################
  
  comparison_names = protigy_parameters[[COMPARISON_LIST_FIELD_NAME]]
  
  comparison_nomenclature_df = data.frame(
    comparison_name = comparison_names,
    significance_name = paste(SIGNIFICANCE_NAME_BASE, comparison_names, sep=JOINER_CHARACTER),
    foldchange_name = paste(FOLD_CHANGE_NAME_BASE, comparison_names, sep=JOINER_CHARACTER),
    direction_name = paste(DIRECTION_NAME_BASE, comparison_names, sep=JOINER_CHARACTER),
    global_direction_name = paste(GLOBAL_DIRECTION_NAME_BASE, comparison_names, sep=JOINER_CHARACTER),
    row.names = comparison_names,
    stringsAsFactors = FALSE
  )

  return(comparison_nomenclature_df)    
}


annotate_genes_with_kegg_ids <- function(protigy_results){

  
  # Load the gene symbol to kegg accession mapping
  gene_to_kegg_table = read.csv(GENESYMBOL_TO_KEGGGENE_MAP, stringsAsFactors = FALSE, row.names = 1)
  
  # Extract the whole protigy result table
  protigy_result_table = protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]

  # Extract just the gene symbols and row ids from the protigy results, note uppercasing of gene symbol names
  temporary_data_frame_for_merging = data.frame(
    id = rownames(protigy_result_table),
    id.concat = protigy_result_table[[ID_CONCAT_NAME]],
    id.mapped = toupper(protigy_result_table[[GENE_SYMBOL_NAME]]),
    row.names = rownames(protigy_result_table),
    stringsAsFactors = FALSE 
  )
  
  # Perform the merge
  temp_merge_result = merge(
    x = temporary_data_frame_for_merging,
    y = gene_to_kegg_table,
    by.x = 'id.mapped',
    by.y = 'Approved.symbol',
    all.x = TRUE,
    sort = FALSE
  )

  # Reorder to match original order and select the kegg accessions
  rownames(temp_merge_result) <- temp_merge_result[['id']]    
  ordered_kegg_accessions <- temp_merge_result[rownames(protigy_result_table),KEGG_GENE_ACCESSION_FIELDNAME]
  
  # Rebind to original result table
  protigy_result_table[[KEGG_GENE_ACCESSION_FIELDNAME]] <- ordered_kegg_accessions
  
  # Replace the original data table object in the larger result object
  protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]<-protigy_result_table
  
  # returns a modified form of the input object - should be reassigned to the input parameter
  return(protigy_results)
}




prepareKEGGhelperFunction <- function(protigy_results_with_kegg, fcfield = 'logFC', 
                                      sigfield = 'adj.P.Val', 
                                      site_dirfield = 'direction_site',
                                      global_dirfield = ' direction_global',
                                      sig_threshold = 0.1) 
  {
  # NOTE: protigy results must have already gone through annotate_genes_with_kegg_ids
  
  protigy_result_table = protigy_results_with_kegg[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]
  
  # Mark each row in the table as up, down, both or unch[anged]  
  protigy_result_table[[site_dirfield]] <-apply(
    X = protigy_result_table,
    MARGIN = 1,
    FUN = function (x) { 
      if (as.numeric(x[[sigfield]]) > sig_threshold 
                         | 
          is.na(as.numeric(x[[sigfield]])) 
                         | 
          is.na(as.numeric(x[[fcfield]]))
      ) 
      {
        return("unch") 
      } else { 
        return(if (as.numeric(x[[fcfield]]) < 0) "down" else "up") 
      }
    }
  )
  
  # Find unique combinations for each gene
  combinations<-(unique(protigy_result_table[,c(KEGG_GENE_ACCESSION_FIELDNAME,site_dirfield)]))
  
  #initialize null result
  collapsed_directions<-"unch"
  
  # consolidate multiple observations of the same gene to a single direction
  if (all(c('up','down') %in% colnames(table(combinations)))) 
  {
    collapsed_directions<-apply(
      X = table(combinations),
      MARGIN = 1,
      FUN = function (x) { if (x['up']==1 & x['down']==1) 
        { 
          return("both") 
        } else { 
          if (x['up']==1) 
          {
            return ("up")
          } else if (x['down']==1) 
          { 
            return ("down") 
          } 
          else 
          {
            return("unch")
          } 
        }
      }
    )
  } else {
    if (any(c('up','down') %in% colnames(table(combinations)))) { warning("Only one type of directionality (up or down) is present. Color plots are wrong.")}
  }
  
  # prepare for merge back to result table
  collapsed_direction_df<-as.data.frame(collapsed_directions)
  colnames(collapsed_direction_df)[1] <- global_dirfield
  
  # perform the merge
  temporary_merge_result<-merge(protigy_result_table,collapsed_direction_df,by.x=KEGG_GENE_ACCESSION_FIELDNAME,by.y=0,all.x=TRUE,sort=FALSE)

  # reassign the new results back into the overall data object
  protigy_results_with_kegg[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]<-temporary_merge_result
  
  return(protigy_results_with_kegg)
  
}


select_minimal_data_and_merge <- function(protigy_results_with_directions, kegg_pathway_data, fcfield = 'logFC', 
                                          sigfield = 'adj.P.Val', 
                                          site_dirfield = 'direction_site',
                                          global_dirfield = ' direction_global') {
  # keggGene - joining key, from both
  # keggPathway - kegg_genes_to_pathways
  # keggPWname - kegg_genes_to_pathways
  # keggNumber - kegg_genes_to_pathways
  # logFC field - from protigy results
  # global direction - from protigy results
  
  protigy_data_table <- protigy_results_with_directions[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]

  fields_needed <- c(KEGG_GENE_ACCESSION_FIELDNAME,
                     global_dirfield,
                     sigfield
                    )
  
  mimimal_protigy_data_to_merge_with_pathway_data <- protigy_data_table[fields_needed]

  merged_protigy_and_pathway_data <- merge(x = mimimal_protigy_data_to_merge_with_pathway_data, 
                                           y = kegg_pathway_data, 
                                           by = KEGG_GENE_ACCESSION_FIELDNAME,
                                           sort = FALSE) 

  return(merged_protigy_and_pathway_data)
   
}


###

findKEGGsetSizes<-function(DF, keggGeneField='keggGene', AdjPValField='adj.P.Val', AdjPValCut=0.1, alternative='greater', useBest=FALSE) {

  ### TODO: NEEDS BETTER DOCS, MAYBE REFACTORING?
    
  DiffExpressedKEGGids<-DF[keggGeneField][DF[AdjPValField] <= AdjPValCut]
  

  totalKeggGeneIDs<-length(unique(unlist(DF[keggGeneField])));
  totalDiffExpressedKEGGids<-length(unique(DiffExpressedKEGGids));
  
  #iterate through all pathways
  
  uniquePWcombos<-unique(DF[,c('keggPathway','keggPWname')])
  uniquePathways = as.character(unlist(uniquePWcombos$keggPathway))
  uniquePathwaysDesc = as.character(unlist(uniquePWcombos$keggPWname))
  outputData<-data.frame(pathwayCode=NULL,pathwayDescription=NULL,fisherPval=NULL,nObs=NULL, nSig=NULL, sigLevel = NULL, stringsAsFactors=FALSE);
  
  for (n in 1:length(uniquePathways)) {
    #for (n in 1:10) {
    allObservedGenesInPathway = unique(subset(DF, keggPathway == uniquePathways[n], select=keggGeneField))
    diffObservedInPathway = unique(subset(allObservedGenesInPathway, keggGene %in% DiffExpressedKEGGids, select=keggGeneField))
    InCategoryDiffExpressed <- length(diffObservedInPathway$keggGene);
    OutCategoryDiffExpressed <- totalDiffExpressedKEGGids - InCategoryDiffExpressed;
    InCategoryNotDiffExpressed <- length(allObservedGenesInPathway$keggGene) - InCategoryDiffExpressed;
    OutCategoryNotDiffExpressed <- totalKeggGeneIDs - InCategoryDiffExpressed - OutCategoryDiffExpressed - InCategoryNotDiffExpressed;
    # print(paste("Pathway",uniquePathways[n],"T:", length(allObservedGenesInPathway$keggGene), "ID:",  InCategoryDiffExpressed, "IN:", InCategoryNotDiffExpressed, "OD:",OutCategoryDiffExpressed,"ON:",OutCategoryNotDiffExpressed,sep=" "))
    fpv<-fisher.test(matrix(c(InCategoryDiffExpressed,InCategoryNotDiffExpressed,OutCategoryDiffExpressed,OutCategoryNotDiffExpressed),nrow=2,ncol=2),alternative=alternative);
    outputData<-rbind(outputData,data.frame(pathwayCode=uniquePathways[n],pathwayDescription=uniquePathwaysDesc[n],fisherPval=fpv$p.value, nObs=length(allObservedGenesInPathway$keggGene), nSig=InCategoryDiffExpressed, sigLevel = AdjPValCut, stringsAsFactors=FALSE))
  }
  return(outputData)
}


################################################################
# HTML GENERATION FUNCTIONS BELOW - PERHAPS MOVE TO SUBLIBRARY #
# ALSO MAYBE BETTER DOCS, REFACTORIZATION                      #
################################################################

writeFormElementToHTML<-function(url,params,connection) {
  cat(paste0("<form action=\"",url,"\" method=\"post\" enctype=\"multipart/form-data\" target=\"_blank\">\n"),file=connection)
  cat(paste0("  <input readonly=\"readonly\" type=\"hidden\" name=\"map\" value=\"",params$map,"\">\n"),file=connection)
  cat(paste0("  <input readonly=\"readonly\" type=\"hidden\" name=\"mode\" value=\"",params$mode,"\">\n"),file=connection)
  cat(paste0("  <input readonly=\"readonly\" type=\"hidden\" name=\"reference\" value=\"",params$reference,"\">\n"),file=connection)
  cat("  <textarea hidden readonly=\"readonly\" id=\"s_q\" name=\"unclassified\" cols=\"1\" rows=\"1\">\n",file=connection)
  cat(gsub("\t"," ",params$unclassified),file=connection);
  cat("  </textarea>\n",file=connection)
  cat("  <input type=\"submit\" value=\"KEGG\">\n",file=connection)
  cat("</form>\n",file=connection)
}

initiateHTMLdoc<-function(connection) {
  
  cat("<!DOCTYPE html>\n",file=connection)
  cat("<html>\n",file=connection)
  
  cat("<head>
      <meta charset=\"utf-8\"><meta http-equiv=\"X-UA-Compatible\" content=\"IE=edge\">
      <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
      <meta name=\"google-site-verification\" content=\"qpiygib7oAZsD4rfeIa2NrdI0CyHVN0MlNteNW4cmaI\">
      <!-- Favicon--><link href=\"//assets.clue.io/clue/public/img/favicon.ico\" rel=\"shortcut icon\" type=\"image/x-icon\">
      <link href=\"https://fonts.googleapis.com/css?family=Roboto:400,100,300,500,700|Roboto+Condensed:400,300,700|Roboto+Slab:100,300,400,700|Istok+Web:400,700|Open+Sans:300,400,400i,600,600i,700,700i,800|Source+Serif+Pro:400,600|Roboto+Mono:400,500\" rel=\"stylesheet\" type=\"text/css\">
      <link href=\"https:/clue.io//public/css/dist/clue.all.min.css?rel=af4c7c85c2\" rel=\"stylesheet\">
      <link href=\"https:/clue.io//public/css/dist/clue.css?rel=f02bebaa6f\" rel=\"stylesheet\">
      <style>
      .tableheader {
      font-weight: bold;
      }
      .tableheader td {
      background-color: #ffffff;
      }
      .spannedtableheading {
      text-align: center;
      }
      .tablesubheader {
      font-style: italic;
      }
      .tablesubheader td {
      background-color: #ffffff;
      }
      .keggtable {
      border-spacing: 2px;
      width: 1024px;
      }
      .keggtable tr:nth-child(odd) {
      background-color: #f2f2f2;
      }
      .keggtable tr:hover {
      background-color: #fbb4b9;
      }
      .keggcontainer {
      width: 1200px;
      padding: 25px 25px 25px 25px;
      }
      </style>
      <title id=\"cluecmap\">KEGG Pathway Mapper</title>
      </head>
      
      ",file=connection)
  
  cat("<body>
      <div id=\"main\" class=\"keggcontainer\">\n",file=connection)
  
  return()
}

writeHeadingAndInitiateTable<-function(connection, heading, conditions) {
  cat("<h1>",heading,"</h1>
      ",file=connection)
  cat("<table class=\"keggtable\">
      ",file=connection)
  cat("<tr class=\"tableheader\">
      <td>Pathway</td>
      ",file=connection)
  for (j in 1:length(conditions)) {
    cat("  <td colspan=\"3\" class=\"spannedtableheading\">",file=connection)
    cat(conditions[j], file=connection)
    cat("</td>\n",file=connection)
  }
  cat("<tr>",file=connection)
  
  cat("<tr class=\"tablesubheader\">
      <td></td>
      ",file=connection)
  for (j in 1:length(conditions)) {
    cat("  <td>",file=connection)
    cat("p.val</td><td>nSig</td><td>nObs</td>", file=connection)
  }
  cat("<tr>",file=connection)
  
  
}

finishTableAndWritePostTableElements<-function(connection, trailer=NULL) {
  cat("</table>
      ",file=connection)
}


finalizeHTMLdoc<-function(connection) {
  
  cat("</div>\n",file=connection)
  cat("</body>\n",file=connection)
  cat("</html>\n",file=connection)
  
  return()
}

generateHTMLdocFromExtendedKEGGplotObject<-function(EKPobject, file, sigcut=0.01, indexCondition=1, sort = TRUE) {
  fcon<-file(file,"w")
  initiateHTMLdoc(fcon);
  
  
  nConditions<-length(EKPobject$conditions)
  conditions<-EKPobject$conditions
  headingText<-"KEGG Pathway Enrichments"
  
  writeHeadingAndInitiateTable(fcon,headingText, conditions);
  
  #clean out any NA rows and assign rownames
  for (k in 1:nConditions) {
    EKPobject$pwenrichments[[k]]<-EKPobject$pwenrichments[[k]][!(is.na(EKPobject$pwenrichments[[k]]$pathwayCode)),]
    rownames(EKPobject$pwenrichments[[k]])<-EKPobject$pwenrichments[[k]]$pathwayCode
  }
  
  #align all pathways enrichments
  pathwayCodeInFirstCondition<-EKPobject$pwenrichments[[1]]$pathwayCode
  for (k in 1:nConditions) {
    EKPobject$pwenrichments[[k]]<-EKPobject$pwenrichments[[k]][pathwayCodeInFirstCondition,]
  }
  
  if (sort) {
    pathwayEnrichmentOrder<-order(EKPobject$pwenrichments[[indexCondition]]$fisherPval)
    for (j in 1:length(EKPobject$pwenrichments)) {
      EKPobject$pwenrichments[[j]]<-EKPobject$pwenrichments[[j]][pathwayEnrichmentOrder,]
    }
  }
  
  
  filteredDataIndices<-which(EKPobject$pwenrichments[[indexCondition]]$fisherPval <= sigcut)
  for (j in 1:length(filteredDataIndices)) {
    currentPathwayName<-EKPobject$pwenrichments[[indexCondition]]$pathwayDescription[j]
    currentPathwayCode<-EKPobject$pwenrichments[[indexCondition]]$pathwayCode[j]
    currentPathwayEnrichmentValues<-double(nConditions);
    currentPathwayNObs<-integer(nConditions);
    currentPathwayNSig<-integer(nConditions);
    for (k in 1:nConditions) {
      currentPathwayEnrichmentValues[k]<-EKPobject$pwenrichments[[k]]$fisherPval[j]
      currentPathwayNObs[k]<-EKPobject$pwenrichments[[k]]$nObs[j]
      currentPathwayNSig[k]<-EKPobject$pwenrichments[[k]]$nSig[j]
    }
    if (!(is.na(currentPathwayName))) { 
      writeKEGGPathwayTableRow(fcon, EKPobject,currentPathwayName,currentPathwayCode,currentPathwayEnrichmentValues,currentPathwayNObs,currentPathwayNSig)
    }
  }
  
  finishTableAndWritePostTableElements(fcon)
  finalizeHTMLdoc(fcon);
  close(fcon)
  return(EKPobject)
}

writeKEGGPathwayTableRow<-function(connection, EKPobject, pathwayName, pathwayCode, enrichmentValues, nObs, nSig) {
  
  textAreaText<-getTextAreaText(EKPobject,pathwayName);
  formParameters<-list(map=pathwayCode,mode="color",reference="white",unclassified=textAreaText);
  
  cat("<tr>
      ", file=connection)
  
  cat("  <td>",file=connection)
  cat(paste0("    ",pathwayName), file=connection)
  cat("</td>\n",file=connection)
  
  for (j in 1:length(enrichmentValues)) {
    cat("  <td>",file=connection)
    cat(formatC(enrichmentValues[j],format = "e", digits=2),file=connection)
    cat("</td>\n",file=connection)
    
    cat("  <td>",file=connection)
    cat(nSig[j],file=connection)
    cat("</td>\n",file=connection)
    
    cat("  <td>",file=connection)
    cat(nObs[j],file=connection)
    cat("</td>\n",file=connection)
    
  }
  
  cat("  <td>",file=connection)
  #cat(paste0("Placeholder"), file=connection)
  writeFormElementToHTML(kurl,formParameters,connection)
  cat("</td>\n",file=connection)
  
  cat("</tr>\n\n", file=connection)
  
  return();
}

getTextAreaText<-function(keggdataListObject, pathwayName, organismCode=defaultKEGGorganismtoken, plotColor="#ad2e91,#eeeeee",diagnostics=FALSE) {
  
  conditions<-keggdataListObject$conditions
  nConditions<-length(conditions)
  kText <- paste0("#",organismCode," ",paste(conditions,collapse=" "),"\n");
  
  datasubsets<-list();
  
  for (n in 1:nConditions) {
    datasubsets[[n]]<-subset(keggdataListObject$datasets[[n]],keggPWname == pathwayName, select = c("keggPathway","keggNumber",grep('global_direction',colnames(keggdataListObject$datasets[[n]]), value=TRUE)));
    colnames(datasubsets[[n]])[1]<-paste0("keggPathway.",n);
    colnames(datasubsets[[n]])[3]<-paste0("direction.",n);
  }
  
  datamerge<-unique(datasubsets[[1]]);
  if (nConditions > 1) {
    for (n in 2:nConditions) {
      datamerge<-merge(datamerge,unique(datasubsets[[n]]),by='keggNumber',all=TRUE)
    }
  }
  
  return.stuff<-datamerge
  
  
  kPath = unique(datasubsets[[1]]$keggPathway);
  
  for (j in 1:length(datamerge$keggNumber)) {
    localPlotColor=plotColor; #default is both, purple
    colorSpecString<-'';
    for (n in 1:nConditions) {
      directionName<-paste0('direction.',n);
      if (is.na(datamerge[j,directionName])) {
        colorSpecString=paste0(colorSpecString,'');
      } else {
        if (datamerge[j,directionName] == "up") {
          colorSpecString=paste0(colorSpecString,"red,#eeeeee");
        }
        if (datamerge[j,directionName] == "down") {
          colorSpecString=paste0(colorSpecString,"blue,#eeeeee");
        }
        if (datamerge[j,directionName] == "both") {
          colorSpecString=paste0(colorSpecString,localPlotColor);
        }
        if (datamerge[j,directionName] == "unch") {
          colorSpecString=paste0(colorSpecString,"#333333,#eeeeee");
        }
      }
      if (n < nConditions) {
        colorSpecString=paste0(colorSpecString,'\t');
      }
    }
    kText <- paste0(kText,datamerge$keggNumber[j]," ",colorSpecString,"\n")
  }
  return.stuff<-(kText);
  
  #h=basicHeaderGatherer()
  #webresult<-postForm(kurl, map = kPath, unclassified = kText, mode = "color", submit = "Exec", reference = "white", .opts= curlOptions(headerfunction = h$update))
  #newurl<-paste0(kurlbase,as.character(unlist((h$value())["Location"])));
  return(return.stuff)
}

#' Generate 1 KEGG map file per Protigy comparison in dataset.
#'
#' Generates a single .html file for each comparison indicated in the Protigy parameters that links to KEGG.
#'
#' This function is useful for large Protigy analyses (e.g. > ~2 comparisons) that would generate unweildy output tables.
#' Instead, it will take the 3 main Protigy objects and generate a single HTML file for each comparison indicated in the
#' \code{global.imp.params$grp.comp} list
#' 
#' @param protigy_input The \code{global.imp.input} object from Protigy
#' @param protigy_parameters The \code{global.imp.parameters} object from Protigy
#' @param protigy_results The \code{global.imp.results} object from Protigy
#' @param skip A list (of strings) of comparisons that should be explicitly skipped. Default is \code{c()}
#' @return None
#' 
#' Note: the function will check to the existence of the output file before proceeding.  Presence of the file will result
#' in skipping of the comparison.  Erase the file to get it to regenerate.
#' 
#' @examples 
#' one_at_a_time_helper(global.imp.input, global.imp.parameters, global.imp.results)
one_at_a_time_helper <- function (protigy_input, protigy_parameters, protigy_results, skip=c()) {

  # CURRENTLY THIS WORKFLOW IS SPECIFIC TO DIA DATA SUPPLIED FROM SEBASTIAN VACA
  
  # There are issues with how the gene symbol is determined or parsed
  # Sometime it appears in the 'id.mapped' field of global.imp.results$data$output
  # Sometimes it is in the geneSymbol field of global.imp.input$table.anno
  # It seems like perhaps it appears in id.mapped if the user-specified "id" field is the uniprot accession number
  # This probably causes Protigy to do the mapping
  # I'm not sure where why is showed up in geneSymbol in Sebastian's data
  if ( all(
    make.names(rownames(protigy_input$table.anno))
    ==
    make.names(rownames(protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]))
  ) ) 
  {
    protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]][[GENE_SYMBOL_NAME]] <- 
    vapply(strsplit(protigy_input$table.anno$geneSymbol,"\\|"), 
           `[`, 1, FUN.VALUE=character(1))
  } else {
    stop("Annotation row names do not match data output rownames.")
  }

  for (j in 1:length(protigy_parameters[[COMPARISON_LIST_FIELD_NAME]])) {

    temp_param = list()
    temp_result = list()
    
    temp_param[[COMPARISON_LIST_FIELD_NAME]]<-protigy_parameters[[COMPARISON_LIST_FIELD_NAME]][j]

    output_filename <- paste0('data/',temp_param[[COMPARISON_LIST_FIELD_NAME]],'.html')
    
    if (file.exists(output_filename)) {
      print(paste0("Exists: ", protigy_parameters[[COMPARISON_LIST_FIELD_NAME]][j]))
      next
    }
    
    if ( protigy_parameters[[COMPARISON_LIST_FIELD_NAME]][j] %in% skip ) {
      print(paste0("Skipping: ", protigy_parameters[[COMPARISON_LIST_FIELD_NAME]][j]))
      next
    } else {
      print(paste0("Working: ", protigy_parameters[[COMPARISON_LIST_FIELD_NAME]][j]))
    }
    
    grep_string = paste(temp_param[[COMPARISON_LIST_FIELD_NAME]],GENE_SYMBOL_NAME,ID_CONCAT_NAME,sep="|")
    grep_result <- grep(grep_string, colnames(protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]]))    
    temp_result[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]] <- 
      protigy_results[[DATA_OBJECT_NAME]][[DATA_OUTPUT_OBJECT_NAME]][,grep_result]
    
    result_object <- map_kegg_pathways(protigy_parameters = temp_param, protigy_results = temp_result, output = output_filename)

  }
  
}
