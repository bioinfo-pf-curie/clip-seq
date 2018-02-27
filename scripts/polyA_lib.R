#PeakIntron_LastPeakSum<- function(x){
#    y <- x[grep("^RNA|\\|RNA",x$status),]
#      if( nrow(y) > 0){
#                CLIP <- x[grep("^RNA|\\|RNA",x$status,invert=TRUE), ]
#                CLIP <- CLIP[grep("CLIP",CLIP$status),c(1:3,6:(ncol(x)-1))]
#          if(nrow(CLIP) > 0){
#		RNA <- colSums(y[,c(9:(ncol(x)-1)) ])
#                names(RNA) <- paste0("RNA.",names(RNA))
#                return(cbind(CLIP, as.list(RNA)))
#              }else return(NULL)
#          } else return(NULL)
#}

PeakIntron_LastPeakSum<- function(x){
    y <- x[grep("^RNA|\\|RNA",x$status),]
      if( nrow(y) > 0){
                NLE <- x[grep("^RNA|\\|RNA",x$status,invert=TRUE), ]
                NLE <- NLE[grep("CLIP",NLE$status),c(1:3,6:(ncol(x)-1))]
          if(nrow(NLE) > 0){
                LE <- colSums(y[,c(9:(ncol(x)-1)) ])
                names(LE) <- paste0("RNA.",names(LE))
                return(cbind(NLE, as.list(LE)))
              }else return(NULL)
          } else return(NULL)
}
