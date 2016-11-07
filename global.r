################################################################################################################
## Filename: global.r
## Created: October 09, 2015
## Author(s): Karsten Krug
##
## Purpose: Shiny-app to perform differential expression analysis, primarily on proteomics data, to perform
##          simple data QC, to interactively browse through the results and to download high-quality result
##          figures.
##
## This file defines global parameters, loads all required R-packages and defines the actual functions to perform data filtering,
## data normalization, the moderated test statistics, and visualization. The code for moderated t-tests,
## two-component normalization and the reproducibility filter has been written by Mani DR and
## adopted by me for intergration into a Shiny-Server environment.
##
##
## required packages:
##
## cran.pckg <- c('pheatmap', 'RColorBrewer', 'hexbin', 'Hmisc', 'grid', 'scatterplot3d', 'plotly', 'WriteXLS', 'reshape','nlme', 'BlandAltmanLeh', 'mice','mixtools', 'mclust')
## bioc.pgkg <- c( 'preprocessCore', 'limma')
##
## changelog: 20160614 included 'na' to indicate missing values
##                     outsourced Mani's code to a separate file 'modT.r'
################################################################################################################

source('modT.r')
source('pheatmap.r')

#################################################################
## global parameters
#################################################################
## version number
VER="0.5.0"
## maximal filesize for upload
MAXSIZEMB <<- 400
## list of strings indicating missing data
NASTRINGS <<- c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?", "na", "#VALUE!")
## speparator tested in the uploaded file
SEPARATOR <<- c('\t', ',', ';')
## Colors used throughout the app to color the defined groups
GRPCOLORS <<- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Set2"))
## number of characters to display in plots/tables for column names
STRLENGTH <<- 20
## operating system
OS <<- Sys.info()['sysname']
## temp directory to write the Excel file
TMPDIR <<- ifelse(OS=='Windows', "./", "/tmp/")
## app name
APPNAME <<- sub('.*/','',getwd())
## directory to store data files
DATADIR <<- ifelse(OS=='Windows', "./", "/local/shiny-data/")

#################################################################
## load required packages
#################################################################
library(shiny)
## heatmap
##library(pheatmap)
library(scales)
library(gtable)
## moderated t-test
library(limma)
## colors
library (RColorBrewer)
## multiscatter
library(hexbin)
library(Hmisc)
library(grid)
## pca
library(scatterplot3d)
library(plotly)
## export
library(WriteXLS)
## reproducibility filter
library(reshape)
library(nlme)
library(BlandAltmanLeh)
## normalization Quantile
library(preprocessCore)
## normalization 2-component
library (mice)
library (mixtools)
library (mclust)


#################################################################################
##     Heatmap of expression values combining all of the results from all tests
##
##
##
#################################################################################
plotHM <- function(res,
                   grp,
                   grp.col,
                   grp.col.legend,
                   hm.clust,
                   hm.title,
                   hm.scale,
                   style,
                   filename=NA, cellwidth=NA, cellheight=NA, max.val=NA, fontsize_col, fontsize_row, ...){

    ## convert to data matrix
    res <- data.matrix(res)

    #########################################
    ## different 'styles' for different tests
    ## - reorder columns
    ## - gaps between experiments
    if(style == 'One-sample mod T'){
        res <- res[, names(grp[order(grp)])]
        gaps_col=cumsum(table(grp[order(grp)]))
        gapsize_col=20
    }
    if(style == 'Two-sample mod T'){
        res <- res[, names(grp[order(grp)])]
        gaps_col=NULL
        gapsize_col=0
    }
    if(style == 'mod F' | style == 'none'){
        res <- res[, names(grp[order(grp)])]
        gaps_col=NULL
        gapsize_col=0
    }
    #########################################
    ## scaling
    if(hm.scale == 'row')
        res <- t(apply(res, 1, function(x)(x-mean(x, na.rm=T))/sd(x, na.rm=T)))
    if(hm.scale == 'column')
        res <- apply(res, 2, function(x)(x-mean(x, na.rm=T))/sd(x, na.rm=T))

    ##########################################
    ##          cluster
    ## 20160309 NA handling
    ##
    ##########################################
    na.idx.row <- na.idx.col <- NULL

    ## column clustering
    if(hm.clust == 'column'){
        Rowv=FALSE
        colv.dist = dist(t(res), method='euclidean', diag=T, upper=T)
        na.idx.col <- which(apply(as.matrix(colv.dist), 1, function(x) sum(is.na(x))) > 0)
        if(length(na.idx.col)> 0){
            colv.dist <- colv.dist[-na.idx.col, ]
            colv.dist <- colv.dist[, -na.idx.col]
        }
        Colv=hclust(as.dist(colv.dist), method='complete')
    ## row clustering
    } else if( hm.clust == 'row'){
        rowv.dist <- as.matrix(dist(res, method='euclidean', diag=T, upper=T))
        na.idx.row <- which(apply(as.matrix(rowv.dist), 1, function(x) sum(is.na(x))) > 0)
        if(length(na.idx.row)> 0){
            rowv.dist <- rowv.dist[-na.idx.row, ]
            rowv.dist <- rowv.dist[, -na.idx.row]
        }
        Rowv=hclust(as.dist(rowv.dist), method='complete')
        Colv=FALSE

    ## row and column clustering
    } else if(hm.clust == 'both'){

        ## row clustering
        rowv.dist <- as.matrix(dist(res, method='euclidean', diag=T, upper=T))
        na.idx.row <- which(apply(as.matrix(rowv.dist), 1, function(x) sum(is.na(x))) > 0)
        if(length(na.idx.row)> 0){
            rowv.dist <- rowv.dist[-na.idx.row, ]
            rowv.dist <- rowv.dist[, -na.idx.row]
        }
        Rowv=hclust(as.dist(rowv.dist), method='complete')

        ## column clustering
        colv.dist = dist(t(res), method='euclidean', diag=T, upper=T)
        na.idx.col <- which(apply(as.matrix(colv.dist), 1, function(x) sum(is.na(x))) > 0)
        if(length(na.idx.col)> 0){
            colv.dist <- colv.dist[-na.idx.col, ]
            colv.dist <- colv.dist[, -na.idx.col]
        }
        Colv=hclust(as.dist(colv.dist), method='complete')
    } else {
        Rowv=Colv=FALSE
    }

    #########################################
    ## capping
    if(!is.na(max.val)){
        res[ res < -max.val ] <- -max.val
        res[ res > max.val ] <- max.val
    }
    #########################################
    ## min/max value
    max.val = ceiling( max( abs(res), na.rm=T) )
    min.val = -max.val

    ##########################################
    ## colors
    color.breaks = seq( min.val, max.val, length.out=12 )
    color.hm = rev(brewer.pal (length(color.breaks)-1, "RdBu"))
    ##color.hm = rev(brewer.pal (length(color.breaks)-1, "Spectral"))

    ##############################################
    ## annotation of columns
    anno.col=data.frame(Group=grp)
    anno.col.color=list(Group=grp.col.legend)

    ##############################################
    ## heatmap title
    ##hm.title = paste(hm.title, '\nsig / total: ', nrow(res), ' / ', ,sep='')
    if(!is.null(na.idx.row) | !is.null(na.idx.col))
        hm.title = paste(hm.title, '\nremoved rows / columns: ', length(na.idx.row), ' / ' , length(na.idx.col), sep='')

    ############################################
    ## plot the heatmap
    pheatmap(res, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
             cluster_rows=Rowv, cluster_cols=Colv, border_col=NA, col=color.hm, filename=filename, main=hm.title, annotation_col=anno.col, annotation_colors=anno.col.color, labels_col=chopString(colnames(res), STRLENGTH), breaks=color.breaks, scale='none', cellwidth=cellwidth, cellheight=cellheight, gaps_col=gaps_col, gapsize_col=gapsize_col, labels_row=chopString(rownames(res), STRLENGTH), na_col='black')
}



#################################################################################################
##                     multiscatterplot using hexagonal binning
## - mat    numerical matrix of expression values, rows are features, columns are samples
##
## changelog: 2015116 implementation
#################################################################################################
my.multiscatter <- function(mat, hexbin=30, hexcut=5, cor=c('pearson', 'spearman', 'kendall'), repro.filt=NULL, grp, grp.col.legend, define.max=F, max.val=3, min.val=-3){

    ## cor method
    corm = match.arg(cor)
    ## correlation
    cm = cor(mat, use='pairwise.complete', method=corm)

    ## number of samples to compare
    N = ncol(mat)

    ## define limits
    if(define.max){
        lim=c(min.val, max.val)

    } else{
        lim=max( abs( mat ), na.rm=T )
        lim=c(-lim, lim)
    }

    ###########################################################################
    ## help function to set up the viewports
    ## original code from:  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
    multiplot <- function(plots, cols=1) {
        ## Make a list from the ... arguments and plotlist
        ##plots <- c(list(...))
        ## number of plots
        numPlots = length(plots)
        ## layout matrix
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
        ## Set up the page
        grid.newpage()
        ## grid layout
        la <-  grid.layout(nrow(layout), ncol(layout))
        pushViewport(viewport(layout = la))
        ## Make each plot, in the correct location
        for (i in numPlots:1) {
            ## Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col)

            ## textplot: correlation coefficient
            if(matchidx$row < matchidx$col){
                numb = plots[[i]]
                col='black'
                size = min(max(abs(90*as.numeric(numb)), 25), 50)
                grid.rect(width=unit(.85, 'npc'), height=unit(.85, 'npc'), vp=vp, gp=gpar(fill='grey95', col='transparent'))
                grid.text(numb, vp=vp, gp=gpar(fontsize=size, col=col))
            } else {
                print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                layout.pos.col = matchidx$col))
            }
        }
    } ## end function 'multiplot'
    #########################################################################
    ##
    ##                     do the actual plotting
    ##
    #########################################################################

    ## list to store the plots
    plotList=vector('list', N*N)
    count=1
    for(i in 1:N)
        for(j in 1:N){

            ## extract pairwise data
            dat <- data.frame(x=mat[,i], y=mat[,j])
            rownames(dat) <- rownames(mat)

            ## filter according to xlim/ylim
            dat$x[ which(dat$x < lim[1] | dat$x > lim[2]) ] <- NA
            dat$y[ which(dat$y < lim[1] | dat$y > lim[2]) ] <- NA

            ## extract groups
            current.group <- unique(grp[names(grp)[c(i,j)]])

            ##cat(current.group, '\n')
            ##str(add.points)

            ###########################
            ## lower triangle
            if(i < j){

                ## hexbin
                ##hex <- hexbin(dat$x, dat$y, hexbin, xbnds=range(dat$x, na.rm=T), ybnds=range(dat$y, na.rm=T) )
                hex <- hexbin(dat$x, dat$y, hexbin, xbnds=lim, ybnds=lim )
                gghex <- data.frame(hcell2xy(hex), c = cut2(hex@count, g = hexcut))
                p <- ggplot(gghex) + geom_hex(aes(x = x, y = y, fill = c) ,stat = "identity") + guides(fill=FALSE) + theme( plot.margin=unit(rep(0, 4), 'cm')) + xlab('') + ylab('') + xlim(lim[1], lim[2]) + ylim(lim[1], lim[2])

                ##if(length(current.group) == 1)
                ##    p <- p + scale_fill_manual( values=paste(rep( grp.col.legend[current.group], hexcut) ))
                ##else

                p <- p + scale_fill_manual( values=paste('grey', ceiling(seq(70, 20, length.out=hexcut)), sep=''))

                ## add filtered values
                if(!is.null(repro.filt) & length(current.group) == 1){
                    not.valid.idx <- repro.filt[[current.group]]
                    dat.repro <- dat[not.valid.idx, ]
                    ##cat(not.valid.idx)
                    ##cat(dim(dat.repro))
                    p = p + geom_point( aes(x=x, y=y ), data=dat.repro, colour=my.col2rgb('red', 50), size=1)
                }
            }
            ###########################
            ## diagonal
            if(i == j){
                p = ggplot(dat, aes(x=x)) + geom_histogram(fill=grp.col.legend[current.group], colour=grp.col.legend[current.group], binwidth=sum(abs(range(dat$x, na.rm=T)))/50) + ggtitle(colnames(mat)[i]) + theme(plot.title=element_text(size=9)) + theme( panel.background = element_blank(), plot.margin=unit(rep(0, 4), 'cm')) + xlab(paste('N',sum(!is.na(dat$x)), sep='=')) + ylab('') + xlim(lim[1], lim[2]) ##+ annotate('text', label=sum(!is.na(dat$x)), x=unit(0, 'npc'), y=unit(0, 'npc'))


            }
            ###########################
            ## upper triangle
            if(i > j){
                cortmp = cm[i,j]

                p=paste(round(cortmp,2))

            }

            plotList[[count]] <- p
            count=count+1
        }

    multiplot( plotList, cols=N)
}


##########################################################################################################
##                     translate a color name into rgb space
##
## changelog:  20100929 implementation
##########################################################################################################
my.col2rgb <- function(color, alpha=80, maxColorValue=255){

    out <- vector( "character", length(color) )

    for(col in 1:length(color)){

        col.rgb <- col2rgb(color[col])

        out[col] <- rgb(col.rgb[1], col.rgb[2], col.rgb[3], alpha=alpha, maxColorValue=maxColorValue)

    }
    return(out)
}


#############################################################################################
##
##              different normalization methods for expression data
##
## 20160235
#############################################################################################
normalize.data <- function(data, id.col, method=c('Median', 'Quantile', 'Median-MAD', '2-component')){
    cat('\n\n-- normalize data --\n\n')

    method = match.arg(method)

    ids = data[, id.col]
    data = data[ , -grep(paste('^', id.col, '$', sep=''), colnames(data))]

    data <- data.matrix(data)

    ## quantile
    if(method == 'Quantile'){
        require("preprocessCore")
        data.norm <- normalize.quantiles(data)
        rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data))


        ## shift median to zero
        data.norm <- apply(data.norm, 2, function(x) x - median(x, na.rm=T))
    }
    ## median only
    if(method == 'Median'){
        data.norm <- apply(data, 2, function(x) x - median(x, na.rm=T))
        ##rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data), sep='.')
    }
    ## median & MAD
    if(method == 'Median-MAD'){
        data.norm <- apply(data, 2, function(x) (x - median(x, na.rm=T))/mad(x, na.rm=T) )
        ##rownames(data.norm) <- rownames(data)
        colnames(data.norm) <- paste( colnames(data), sep='.')
    }
    ## 2-component normalization
    if(method == '2-component'){
        data.norm.list = apply(data, 2, two.comp.normalize, type="unimodal")
        ##cat('\n\n here 1\n\nlength list:', length(data.norm.list), '\n')
        ##save(data.norm.list, file='test.RData')
        ## check if successful
        for(i in 1:length(data.norm.list)){
            ##cat('\ni=', i, '\n')
            if(length(data.norm.list[[i]]) == 1){
                if(data.norm.list[[i]] == 'No_success'){
                    ##cat('\n\nno success\n\n')
                    return(paste( colnames(data)[i] ))
                }
            }
            ##cat('\n length:', length(data.norm.list[[i]]), '\n')
        }
        ##cat('\n\n here \n\n')
        data.norm = matrix( unlist(lapply(data.norm.list, function(x)x$norm.sample)), ncol=length(data.norm.list), dimnames=list(rownames(data), names(data.norm.list)) )
    }
    ## add id column
    data.norm <- data.frame(ids, data.norm)
    colnames(data.norm)[1] <- id.col
    return(data.norm)
}


##############################################################################
##
##  - perform principle component analysis
##  - calculate variances explained by components
##  - plot the results
##
## changelog: 20131001 implementation
##            20131007 3D plot now plot pc1 vs. pc3 vs. pc2
##                     instead of pc1 vs. pc2 vs. pc3
##            20151208 legend
##            20161103 check number of rows (N), 2D plot only
##############################################################################
my.prcomp <- function(x, col=NULL, cor=T, plot=T, rgl=F, scale=T, pch=20, cex.points=3, rgl.point.size=30, main="PCA", leg.vec=NULL, leg.col=NULL, ...){

    cex.font = 1.8

    ##View(x)

    ## number of data columns, N=2 -> 2D plot only
    N <- nrow(x)

    ## color
    if( is.null(col) ) col="black"

    ## perform pca
    pca <- prcomp(x, scale=scale)

    ##View(pca$x)

    ## calculate variance
    comp.var <- eigen(cov(pca$x))$values

    ## extract the principle components
    pc1=pca$x[,1]
    pc2=pca$x[,2]
    if(N>2)
        pc3=pca$x[,3]

    ##############
    # rgl plot
    ##############
    if(rgl & N > 2){
        require(rgl)
        plot3d(pc1, pc2, pc3, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), type="s", col=col, expand=1.2, size=rgl.point.size)
    }

    ########################################
    # scatterplot 2D/3D
    ########################################
    if( plot){

         require(scatterplot3d)

        if(N > 2)
            par(mfrow=c(1,3), mar=c(7,7,3,1))
        if(N <= 2)
            par(mfrow=c(1,2), mar=c(7,7,3,1))

         ## PC 1-2
         plot(pc1, pc2, xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), pch=pch, main=main, col=col, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:2]/sum(comp.var)),1),"%", sep=""), cex=cex.points, ylim=c( min(pc2),  max(pc2)+.15*max(pc2)), cex.axis=cex.font, cex.lab=cex.font, cex.sub=cex.font )


        if(N > 2) {
            ## PC 1-3
            scatterplot3d( pca$x[,1], pca$x[,3], pca$x[,2], xlab=paste("PC 1 (", round(100*comp.var[1]/sum(comp.var),1),"%)", sep=""), zlab=paste("PC 2 (",round(100*comp.var[2]/sum(comp.var),1),"%)", sep=""), ylab=paste("PC 3 (", round(100*comp.var[3]/sum(comp.var),1),"%)", sep=""), color=col,  cex.symbols=cex.points, pch=pch, main=main, sub=paste("Cumulative variance = ", round(100*sum(comp.var[1:3]/sum(comp.var)),1),"%", sep=""), type="h" )

        }
        ## legend
         plot.new()
         plot.window(xlim=c(0,1), ylim=c(0, 1))
         if(!is.null(leg.vec) & !is.null(leg.col))
             legend('topleft', legend=leg.vec, col=leg.col, pch=pch, pt.cex=max(1, cex.points-1.5), ncol=ifelse( length(leg.vec)> 10, 2, 1), bty='n', cex=2 )
       par(mfrow=c(1,1))
    }

    return(pca)
}

#####################################################
## color ramp
##
## ToDo: opacity!
####################################################
myColorRamp <- function(colors, values, range=NULL) {

    if(is.null(range))
        v <- (values - min(values))/diff(range(values))
    else
        v <- (values - min(values, na.rm=T))/diff( range )

    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

#################################################
##   Given a string and a number of characters
##   the function chops the string to the
##   specified number of characters and adds
##   '...' to the end.
## parameter
##   string     - character
##   nChar      - numeric
## value
##   string of 'nChar' characters followed
##     by '...'
##################################################
chopString <- function(string, nChar=10, add.dots=T)
{

    string.trim <- strtrim(string, nChar)

    if(add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ], '...')
    if(!add.dots)
        string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ])

    return(string.trim)

}


########################################################################
## 20160224
##                   reproducibility filter
##
## n=2: Bland-Altman
## n>2: lmm-model written by Mani DR
##
## - replaces not reprodicibly measuered values in 'tab' with 'NA'
##
########################################################################
my.reproducibility.filter <- function(tab, grp.vec, id.col='id', alpha=0.05){

    ## extract groups
    groups = unique(grp.vec)

    ## list to store index of filtered values per group
    values.filt <- vector('list', length(groups))
    names(values.filt) <- groups

    ## add rownames to tab
    ##rownames(tab) <- tab[, id.col]

    tab.repro.filter <- tab
   ## View(tab.repro.filter)

    ############################################
    ## loop over replicate groups
    for(gg in groups){

        gg.idx = names(grp.vec)[ which(grp.vec == gg) ]

        ########################################
        ## if there are more than 2 replicates
        ## use the Mani's lmm model
        if( length(gg.idx) > 2 ){
            repro.idx <- reproducibility.filter( tab[, c(id.col, gg.idx)], id.col=id.col, alpha=alpha)

            if(length(repro.idx) != nrow(tab)) stop('Reproducibility vector not of same length as matrix!\n')

            not.repro.idx <- which(!repro.idx)

            if(length(not.repro.idx) > 0)
                tab[not.repro.idx, gg.idx] <- NA

            values.filt[[gg]] <- not.repro.idx
        }
        ########################################
        ## if there are two replicates use
        ## Blandt-Altmann filter
        ## R-package 'BlandAltmanLeh'
        if( length(gg.idx) == 2 ){

            ## Bland-Altman
            ba <-  bland.altman.stats(as.numeric( as.character( tab[, gg.idx[1] ]) ), as.numeric( as.character( tab[,  gg.idx[2] ] )), two=3.290527 )
            ## calculate diffs on my own..
            my.diffs <- tab[, gg.idx[1]] - tab[, gg.idx[2]]
            ## index of outliers
            ##not.repro.idx <- which( ba$diffs < ba$lower.limit | ba$diffs > ba$upper.limit)
            not.repro.idx <- which( my.diffs < ba$lower.limit | my.diffs > ba$upper.limit)

            ## set values of outliers to NA
            if(length(not.repro.idx) > 0)
                tab[not.repro.idx, gg.idx] <- NA

            ## store the results
            values.filt[[gg]] <- rownames(tab)[ not.repro.idx ]
            rm(not.repro.idx)
        }

    }
    return(list(table=tab, values.filtered=values.filt))
}


##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n){
    if( n < 50)
        height=500
    if( n >= 50 & n <= 100)
        height=800
    if(n >=100)
        height=800+n
    return(height)
}

###################################################################
##
##       generate the boxplots under the 'QC' tab
##
###################################################################
makeBoxplot <- function(tab, id.col, grp, grp.col, grp.col.leg, legend=T, cex.lab=1.5, mar=c(4,15,2,6)){

	 cat('\n-- makeBoxplot --\n')

    ## table
    tab <- tab[, setdiff(colnames(tab), id.col)]

	##	cat(id.col)

    ##########################################
    ## order after groups
    ord.idx <- order(grp)
    grp <- grp[ord.idx]
    tab <- tab[, ord.idx]
    grp.col <- grp.col[ ord.idx]


    at.vec=1:ncol(tab)
    ##########################################
    ## plot
    par(mar=mar)
    boxplot(tab, pch=20, col='white', outline=T, horizontal=T, las=2, xlab=expression(log[2](ratio)), border=grp.col, at=at.vec, axes=F, main='', cex=2, xlim=c(0, ifelse(legend, ncol(tab)+2, ncol(tab)) ))
    ##legend('top', legend=names(grp.col.leg), ncol=2, bty='n', border = names(grp.col.leg), fill='white', cex=1.5)
    if(legend)
        legend('top', legend=names(grp.col.leg), ncol=length(grp.col.leg), bty='n', border = grp.col.leg, fill=grp.col.leg, cex=cex.lab)
    ##legend('top', legend=c(input$label.g1, input$label.g2), ncol=2, bty='n', border = c('grey10', 'darkblue'), fill='white', cex=1.5, lwd=3)
    mtext( paste('N=',unlist(apply(tab,2, function(x)sum(!is.na(x)))), sep=''), at=at.vec, side=4, las=2, adj=0, cex.lab=cex.lab)
    axis(1)
    axis(2, at=at.vec, labels=chopString(colnames(tab), STRLENGTH), las=2, cex=cex.lab)

}

###################################################################
##
##       generate the profile plots under the 'QC' tab
##
###################################################################
makeProfileplot <- function(tab, id.col, grp, grp.col, grp.col.leg, legend=T, cex.lab=1.5, mar=c(5,5,3,1), ... ){

    cat('\n-- makeProfileplot --\n')

    ## table
    tab <- tab[, setdiff(colnames(tab), id.col)]

    xlim=max(abs(tab), na.rm=T)

    ## caclulate densities
    dens <- apply(tab, 2, density, na.rm=T)

    ## ylim
    ylim <- max(unlist(lapply(dens, function(x) max(x$y))))

    ##########################################
    ## plot
    par(mar=mar)
    for(i in 1:ncol(tab)){
        if(i == 1)
            plot(dens[[i]], xlab='expression', xlim=c(-xlim, xlim), ylim=c(0, ylim), col=my.col2rgb(grp.col[i], alpha=100), lwd=3, cex.axis=2, cex.lab=2, cex.main=1.5, ...)
        else
            lines(dens[[i]], col=my.col2rgb(grp.col[i], alpha=100), lwd=3)

        ## divide legend if there are too many experiments
        N.exp <- length(names(grp.col.leg))
        if( N.exp > 15){
            legend('topright', legend=names(grp.col.leg)[1:floor(N.exp/2)], col=grp.col.leg[1:floor(N.exp/2)], lty='solid', bty='n', cex=1.5, lwd=3)
            legend('topleft', legend=names(grp.col.leg)[ceiling(N.exp/2):N.exp], col=grp.col.leg[ceiling(N.exp/2):N.exp], lty='solid', bty='n', cex=1.5, lwd=3)

        } else
            legend('topright', legend=names(grp.col.leg), col=grp.col.leg, lty='solid', bty='n', cex=1.5, lwd=3)
    }

    cat('\n-- makeProfileplot exit --\n')
}
