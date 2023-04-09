
#' A function apply marginal test on the MRGN inferred models
#'
#' Since MRGN only applies conditional test with the inferred models, it does not account for the marginal test that could influence the results. So, we use this function to update the inferred models.
#'
#' @param res A dataframe of dimension 1 x 14 (generated from analyzeTrios()):
#'
#' @return The input dataframe with the updated inferred models column
#' @export
#'
#' @seealso [analyzeTrios()] to generate the input data (causal inference models data) for this function
#'
#' @examples #Use the function to perform the MRGN inference with confounding variables
#' @examples final.result = analyzeTrios(meth, gene, cna, final.trios.df, pc.meth.pos[[1]], pc.gene.pos[[1]], pc.meth.pos[[2]]$sig.asso.covs, pc.gene.pos[[2]]$sig.asso.covs, clinical.pos, pc.meth.pos[[3]], pc.gene.pos[[3]], 2, 3, 5, 3)
#' @examples final.result
#'
#' @examples df = postFilter(final.result)
#' @examples df
#'


postFilter <- function(res){

    na_row = which(is.na(res))

        if(length(na_row) > 0){

              res <- res[-which(is.na(res), arr.ind=TRUE),]
        }else{

              res = res
        }

    #create a new column for new model name after filtering
    res$Inferred.Model2[1:nrow(res)] <- rep(NA)

    #loop through all the rows in the data
    for(i in 1:nrow(res)){

      print(i)

      #if classification is either M0, M1, or M2
      if(res$Inferred.Model[i] == "M0.1" | res$Inferred.Model[i] == "M0.2" | res$Inferred.Model[i] == "M2.1" | res$Inferred.Model[i] == "M2.2"){

        #we check the marginal p values and if they are less than 0.01, we re classify as "Other"
        if(res$`pV1:T2`[i] < 0.01 & res$`pV1:T1`[i] < 0.01){

          res$Inferred.Model2[i] <- "Other"

        }else{

          res$Inferred.Model2[i] <- res$Inferred.Model[i]

        }

      }else{

        res$Inferred.Model2[i] <- res$Inferred.Model[i]

      }

      #for all the models, we check the marginal pvalues, if they are greather than 0.05, we re classify as "Other"
      if(res$`pV1:T2`[i] > 0.05 & res$`pV1:T1`[i] > 0.05){

        res$Inferred.Model2[i] <- "Other"

      }else{

        res$Inferred.Model2[i] <- res$Inferred.Model[i]

      }

    }


    #if classification is M1.1
    if(res$Inferred.Model[i] == "M1.1"){

      #we check the marginal p values between V1 and T2. If they are greather than 0.05, we re classify as "Other"
      if(res$`pV1:T2`[i] > 0.05){

        res$Inferred.Model2[i] <- "Other"

      }
    }else{

      res$Inferred.Model2[i] <- res$Inferred.Model[i]

    }

    #if classification is M1.2
    if(res$Inferred.Model[i] == "M1.2"){

      #we check the marginal p values between V1 and T1. If they are greather than 0.05, we re classify as "Other"
      if(res$`pV1:T1`[i] > 0.05){

        res$Inferred.Model2[i] <- "Other"

      }
    }else{

      res$Inferred.Model2[i] <- res$Inferred.Model[i]

    }

    return(res)

}
