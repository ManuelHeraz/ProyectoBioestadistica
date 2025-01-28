library("sybil")
svr <- "glpk" #can be cplex or glpk

if (svr == "cplex") {
  library(cplexAPI)
  SYBIL_SETTINGS("SOLVER", "cplexAPI")
}

bac_models <- readRDS("metamouse-20230510.RDS")

flux_thr <- 1e-06 #flux threshold accounting for potential numerical inaccuracies of the solver's results

i <- 0
stats_fva <- vector()
fva_reactions <- vector(mode = "list", length = length(bac_models)) #initializing the vector that will include the lists of active reactions

for (model in bac_models){
  i <- i + 1
  print(i)
  fba <- optimizeProb(model, algorithm = "fba")
  if (mod_obj(fba)>= flux_thr){ #checks if the model can grow
    opt <- fluxVar(model, percentage = 99, verbose = 0) #minimum percentage of the objective function's optimum value is set to 99 to allow for slightly suboptimal, yet feasible, fluxes distributions
    st <- checkOptSol(opt)@status_code
    np <- num_of_prob(opt)
    r_ids <- react_id(react(opt))
    minfl <- lp_obj(opt)[1:(np/2)] # select FVA minimums
    maxfl <- lp_obj(opt)[(np/2+1):np] # select FVA maximums
    r_bounds <- lapply(1:(np/2), function(r_num){ # for each reaction, selects the respective minimum and maximum flux obtained with FVA
      r_id <- r_ids[r_num]
      lb <- minfl[r_num]
      ub <- maxfl[r_num]
      res <- c(lb, ub)
    })
    names(r_bounds) <- r_ids
    
    active <- sapply(r_bounds, function(reac){ #select all active reactions based on reaction bounds given the diet constraints
      (abs(reac[1]) >= flux_thr) || (abs(reac[2]) >= flux_thr)
    })

    active <- names(active)[active]

    res <- list (active)
    names(res) <- c("active")
    stats_fva[[i]] <- st
    names(stats_fva)[[i]] <- model@mod_name
    fva_reactions[[i]] <- res
    names(fva_reactions)[[i]] <- model@mod_name
  }
  else {
    fva_reactions[[i]] <- "no growth"
    names(fva_reactions)[[i]] <- model@mod_name
    stats_fva[[i]] <- "no growth"
    names(stats_fva)[[i]] <- model@mod_name
  }
}

for (name in names(stats_fva)){
  if(svr == "glpk"){
    if(stats_fva[name]!=5){ #5 is the "OK" status code of the glpk solver. If the status code is different, the results from the respective model need to be checked
      print(name)
      print(stats_fva[name])
    }
  } else if (svr == "cplex"){ #1 is the "OK" status code of the cplex solver. If the status code is different, the results from the respective model need to be checked
    if(stats_fva[name]!=1){
      print(name)
      print(stats_fva[name])
    }
  }
}

saveRDS(fva_reactions,"fva_99_reactions_thr06.RDS")


