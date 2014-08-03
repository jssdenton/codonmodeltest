## codonmodeltest.R (2013-09-30)

## Copyright 2013 John S. S. Denton

cat("\n\n-------------------------------------------\n")
cat("CodonModelTest, v0.22 by John S. S. Denton\n")
cat("     initial release November, 2013        \n")
cat("           jdenton@amnh.org               \n")
cat("-------------------------------------------\n\n")

require(ape, quiet = TRUE)
require(stringr, quiet = TRUE)
source("convert.R")

os <- Sys.info()[1]

# Number of cores:

if (os == "Windows") {
	
	if (as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) > 1) {
    			cat("Detected", Sys.getenv('NUMBER_OF_PROCESSORS'), "cores. To use multicore codonPhyML via OpenMP", 
    					"ensure the appropriate configuration option is set when compiling.\n\n") 
    			}
	if (as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) == 1) {
 				cat("Detected 1 core. It is recommended that codon model selection be run", 
 						"on multicore machines to improve performance on large datasets.\n\n") 
 				}

		} else { 
	
		n.cores <- system("/usr/sbin/sysctl -n hw.ncpu 2>/dev/null", intern = TRUE)
		
		if (n.cores > 1) {
    			cat("Detected", n.cores, "cores. To use multicore codonPhyML via OpenMP", 
    					"ensure the appropriate configuration option is set when compiling.\n\n") 
    			}
		if (n.cores == 1) {
 				cat("Detected 1 core. It is recommended that codon model selection be run", 
 						"on multicore machines to improve performance on large datasets.\n\n") 
 				}
			}
	
cat(paste("You are running", R.Version()$version.string, sep = " "))

cat("\n To run CodonModelTest, ensure that your working directory is set to the folder containing the script. \n")

cat("\n\nCall function as codonmodeltest('FILENAME', genetic.code, method, model.set, refine.m3)\n")

codonmodeltest <- function(seqfile, format = "interleaved", itree = NULL,
                      exclude = NULL, execname = NULL, append = TRUE, genetic.code = NULL, method = NULL, model.set, refine.m3 = TRUE)

{
  	
## Check sequence file format:

seq.check(seqfile)
	
  	if (is.null(genetic.code)) stop("NULL GCODE: You must specify a genetic code. Options are ", 
  				"STANDARD, TVMC, TYMC, THMPCMCMSC, THIMC, THCDHNC, THEFMC, THENC, THBAPPC, THAYNC, THAMC", 
  				"THAFMC, BLNC, CHMC, TRMC, SCOMC, or THMC. See NCBI Genetic Codes page for explanation.")
       
    
## Note: fix the exe search to begin from the home directory. Currently, needs to have exe in same path as script and sequence file.

    if (is.null(execname)) {
    	message("Locating codonPHYML executable...")
        if (os == "Linux") stop("COMPATIBILITY: Check codonphyml website for Linux executable.")
        if (os == "Darwin")  execname <- list.files(pattern = "codonPhyML_macosx", recursive = TRUE, full.names = TRUE)        
        if (os == "Windows") execname <- list.files(pattern = "codonPhyML_windows_binary", recursive = TRUE, full.names = TRUE)
    }
 
# Clear previous version of likelihood output file:

	if (file.exists(paste(seqfile, "_codonphyml_stats.txt", sep = ""))) file.remove(paste(seqfile, "_codonphyml_stats.txt", sep = ""))

    if (is.null(execname))
        stop("MISSING EXECUTABLE: You must give an executable file name for codonPhyML.")
   
   if (model.set == "phylo") {
    	cat("\n\nFitting Goldman and Yang (GY94) model variants for downstream phylogenetic inference.\n\n")
	   codonmodeltest.model <-
    		c("GY94+M0+F1xCodon", "GY94+M0+F1x4", "GY94+M0+F3x4",
      		  "GY94+M3+F1xCodon","GY94+M3+F1x4","GY94+M3+F3x4", 
      		  "GY94+M5+F1xCodon","GY94+M5+F1x4","GY94+M5+F3x4","GY94+M0+Gamma[3]+F1xCodon","GY94+M0+Gamma[3]+F1x4",
      		  "GY94+M0+Gamma[3]+F3x4")
		
		N <- length(codonmodeltest.model)

# Correct the number of free parameters depending on the genetic code selected (only applies to CODON-based frequencies):
# Note -- MG + CODON variants are not included, so no correction is needed for them in model.set = "all", below.
	
	corr.vec.gy <- c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0)
	codonmodeltest.nfp <- c(62, 5, 11, 66, 9, 15, 63, 6, 12, 64, 7, 13)
	
	if (toupper(genetic.code) == "STANDARD") {
		codonmodeltest.nfp <- codonmodeltest.nfp - 0 * corr.vec.gy 
		} else if (toupper(genetic.code) == "TVMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp - corr.vec.gy 
		} else if (toupper(genetic.code) == "TYMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "THMPCMCMSC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "THIMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "THCDHNC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + 2 * corr.vec.gy 
		} else if (toupper(genetic.code) == "THEFMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "THENC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "THBAPPC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + 0 * corr.vec.gy 
		} else if (toupper(genetic.code) == "THAYNC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + 0 * corr.vec.gy 
		} else if (toupper(genetic.code) == "THAMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "THAFMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + 2 * corr.vec.gy 
		} else if (toupper(genetic.code) == "BLNC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "CHMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "TRMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.gy 
		} else if (toupper(genetic.code) == "SCOMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp + 0 * corr.vec.gy 
		} else if (toupper(genetic.code) == "THMC") {
			codonmodeltest.nfp <- codonmodeltest.nfp - corr.vec.gy 
		}

		}
   
   if (model.set == "all")  {
   		cat("Fitting both Goldman and Yang (GY94) and Muse and Gaut (MG94) model variants.\n\n") 

	   codonmodeltest.model <-
    		c("GY94+M0+F1xCodon", "GY94+M0+F1x4", "GY94+M0+F3x4",
      		  "GY94+M3+F1xCodon","GY94+M3+F1x4","GY94+M3+F3x4", 
      		  "GY94+M5+F1xCodon","GY94+M5+F1x4","GY94+M5+F3x4","GY94+M0+Gamma[3]+F1xCodon","GY94+M0+Gamma[3]+F1x4",
      		  "GY94+M0+Gamma[3]+F3x4","MG94+M0+F1x4","MG94+M0+F3x4","MG94+M3+F1x4","MG94+M3+F3x4",
			  "MG94+M5+F1x4","MG94+M5+F3x4","MG94+M0+Gamma[3]+F1x4","MG94+M0+Gamma[3]+F3x4")

       N <- length(codonmodeltest.model) 
 
 	corr.vec.all <- c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    		
    codonmodeltest.nfp <-
		    c(62, 5, 11, 66, 9, 15, 63, 6, 12, 64, 7, 13, 5, 11, 9, 15,
      			6, 12, 7, 13)

	if (toupper(genetic.code) == "STANDARD") {
		codonmodeltest.nfp <- codonmodeltest.nfp - 0 * corr.vec.all } 
	if (toupper(genetic.code) == "TVMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp - corr.vec.all }
	if (toupper(genetic.code) == "TYMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "THMPCMCMSC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "THIMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "THCDHNC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + 2 * corr.vec.all }
	if (toupper(genetic.code) == "THEFMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "THENC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "THBAPPC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + 0 * corr.vec.all }
	if (toupper(genetic.code) == "THAYNC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + 0 * corr.vec.all }
	if (toupper(genetic.code) == "THAMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "THAFMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + 2 * corr.vec.all }
	if (toupper(genetic.code) == "BLNC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "CHMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "TRMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + corr.vec.all }
	if (toupper(genetic.code) == "SCOMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp + 0 * corr.vec.all }
	if (toupper(genetic.code) == "THMC") {
		codonmodeltest.nfp <- codonmodeltest.nfp - corr.vec.all }
 
		}
    
    format <- match.arg(format, c("interleaved", "sequential"))
    fmt <- rep("", N)
    if (format != "interleaved") fmt[] <- "-q"
  
# Options for genetic code: Standard is universal gcode, and can build initial topology using KOSI07. 
# Anything else requires either JC distances, or a specified empirically-derived matrix as input.
# Default, here, is to use JC distances if gcode dne Standard.


if (toupper(genetic.code) == "STANDARD")  {
			gcode <- rep("-g STANDARD", N)
			i.tree.dist <- "--dist_tree_model KOSI07"
				} 	else { gcode <- rep(paste("-g", toupper(genetic.code), sep=" "), N)
			i.tree.dist <- "--dist_tree_model JC69"
		}


    boot <- rep("-b 0", N)
    
  if (model.set == "phylo") {
  		mdl <- paste("-m", rep("GY", 15))
  		wparam <- paste("-w", rep(c("DM0","DMODEL","DGAMMA","DM0"), c(3,3,3,3)))
  		ratevar <- rep(c("-c 1", "-c 3"), c(9,3))
  		freqmodel <- paste("--fmodel", c(rep(c("F1XCODONS", "F1X4", "F3X4"), 4)))
  		}
  		
  if (model.set == "all")  {
  		mdl <- paste("-m", rep(c("GY","MG"), c(12,8)))
    	wparam <- paste("-w", rep(c("DM0","DMODEL","DGAMMA","DM0","DMODEL","DGAMMA","DM0"), c(3,3,3,5,2,2,2)))
   		ratevar <- rep(c("-c 1", "-c 3", "-c 1", "-c 1", "-c 3"), c(9,3,4,2,2))
    	freqmodel <- paste("--fmodel", c(rep(c("F1XCODONS", "F1X4", "F3X4"), 4), rep(c("F1X4", "F3x4"), 4)))
    	}
	iters <- rep("--optHeuristic 4,1", N)	## These values are for a VERY quick-and-dirty model testing.
    tstv <- rep("-t e", N)
    wclasses <- rep("--wclasses 3", N) ## ignored for M0 model.
    alpha <- rep("-a e", N)    
    freqopt <- rep("-f e", N)   ## Use empirical estimates to calculate frequencies. 
    							## This will be switched to "-f o" in later, faster versions. 
      
    cmd <- paste(execname, "-i", seqfile, fmt, gcode, boot, tstv, mdl, wparam, wclasses, 
    				alpha, ratevar, freqmodel, freqopt, i.tree.dist, iters, "-o lr", "--append")
    	outfile <- paste(seqfile, "_codonphyml_stats.txt", sep = "")
    if (!append) {
        unlink(outfile)
        unlink(paste(seqfile, "_codonphyml_tree.txt", sep = ""))
    }
    
    imod <- 1:N
    
    for (i in imod) { 
    		message(paste("Testing model", i, "of", N, ":", codonmodeltest.model[i], sep = " "))
			
		if (os == "Darwin")	{
				system(cmd[i], intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE) 
							} else {
							system(cmd[i], intern = FALSE, invisible = TRUE, show.output.on.console = FALSE)  ## Flip to invisible=FALSE, show...=TRUE to debug.
							}
							
			}

    l <- readLines(outfile, warn = FALSE)
    l <- grep("\\.\\s*Log-likelihood:", l, value = TRUE)
	
# Now, rename the LL-containing file generated by the first pass, to free up name for M3 refinement. 
	
	file.rename(outfile, paste(seqfile, "_model_fit_first_pass.txt", sep = ""))
    
# in case there were already some results in the output file:

     if (dd <- length(l) - length(imod)) l <- l[-(1:dd)]
 	 loglik <- as.numeric(sub(". Log-likelihood:", "", l))

ic.table <- matrix(nrow = length(codonmodeltest.model), ncol = length(c("k", "Log-Lik", "AIC", "AICc", "BIC")))

	nsites <- as.numeric(str_split(str_trim(readLines(seqfile, warn = FALSE)[[1]][1], side="both"), "\\s")[[1]][2]) / 3
		
    ic.table[,1] <- codonmodeltest.nfp    ## N parameters
    ic.table[,2] <- loglik
    ic.table[,3] <- 2 * (codonmodeltest.nfp - loglik)      ## AIC
    ic.table[,4] <- 2 * (codonmodeltest.nfp - loglik) + 2 * ((codonmodeltest.nfp * (codonmodeltest.nfp + 1))) / (nsites - codonmodeltest.nfp - 1)     ## AICc
    ic.table[,5] <- -2 * loglik + (codonmodeltest.nfp * log(nsites))    ## BIC
    ic.table <- as.data.frame(ic.table, row.names = codonmodeltest.model)
	names(ic.table) <- c("k", "Log-Lik", "AIC", "AICc", "BIC")

message ("\n\nCreating and populating Results folder...\n\n")

cmt <- list.files(pattern = "codonmodeltest\\.", recursive = FALSE, full.names = TRUE)
pdir <- dirname(cmt)
results.dir <- paste(strsplit(basename(seqfile), "\\.")[[1]][1], "_codonmodeltest", sep="")
	
dir.create(results.dir)

if (is.null(method)) {

		if (os == "Darwin") {

		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)

			} else {

		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
			}

		}

if (toupper(method) == "AIC") {
			delta.aic <- ic.table$AIC - min(ic.table$AIC)
			weights <- round(exp(-0.5*delta.aic) / sum(exp(-0.5*delta.aic)), digits=10)
			ic.table <- cbind(ic.table, delta.aic, weights)
			ic.table <- ic.table[with(ic.table, order(-weights)), ]
			ic.table <- ic.table[ic.table$weights!=0, ]
			
			if (length(ic.table$weights) == 1) {
				message("A single model has weight 1.0. Check parameter estimates for possible convergence issues.")
					} else if (length(ic.table$weights > 1 && abs(ic.table$delta.aic[1] - ic.table$delta.aic[2])) <= 4) {
				warning("MODEL: Supports separating the top two models do not exceed the recommended threshold.") 
				} 
						
	if (os =="Darwin") {
	
		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
		}


		}
		
if (toupper(method) == "AICC")
		{	delta.aicc <- ic.table$AICc - min(ic.table$AICc)
			weights <- round(exp(-0.5*delta.aicc) / sum(exp(-0.5*delta.aicc)), digits=10)
			ic.table <- cbind(ic.table, delta.aicc, weights)
			ic.table <- ic.table[with(ic.table, order(-weights)), ]
			ic.table <- ic.table[ic.table$weights!=0, ]

			if (length(ic.table$weights) == 1) {
				message("A single model has weight 1.0. Check parameter estimates for possible convergence issues.")
					} else if (length(ic.table$weights > 1 && abs(ic.table$delta.aicc[1] - ic.table$delta.aicc[2])) <= 4) {
				warning("MODEL: Supports separating the top two models do not exceed the recommended threshold.") 
				} 
			

	if (os =="Darwin") {
	
		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
		}


		}

if (toupper(method) == "BIC")
		{	delta.bic <- ic.table$BIC - min(ic.table$BIC)
			weights <- round(exp(-0.5*delta.bic) / sum(exp(-0.5*delta.bic)), digits=10)
			ic.table <- cbind(ic.table, delta.bic, weights)
			ic.table <- ic.table[with(ic.table, order(-weights)), ]
			ic.table <- ic.table[ic.table$weights!=0, ]
			
			if (length(ic.table$weights) == 1) {
				message("A single model has weight 1.0. Check parameter estimates for possible convergence issues.")
					} else if (length(ic.table$weights > 1 && abs(ic.table$delta.bic[1] - ic.table$delta.bic[2])) <= 4) {
				warning("MODEL: Supports separating the top two models do not exceed the recommended threshold.") 
				} 
			

	if (os =="Darwin") {
	
		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(ic.table, file = paste(getwd(), "/", results.dir, "/", "codonmodeltest.out", sep=""), sep = "\t", append = FALSE)
		
		}


		}

best <- row.names(ic.table[1,])

cat("------------------------------------\n")
cat("      CodonModelTest Results        \n")
cat("------------------------------------\n")
cat(paste("Best-fit model is:", best, "\n\n", sep=" ")) 

elements <- str_split(best, "\\+")

if (length(elements[[1]]) < 4) { 
	gamma.para <- "Gamma[1]"
	freq.para <- elements[[1]][3] } else { 
	gamma.para <- elements[[1]][3]
	freq.para <- elements[[1]][4] }

g.para <- elements[[1]][1]
omega.para <- elements[[1]][2]

if (omega.para == "M3" && model.set == "phylo" && refine.m3 == "TRUE") {

{
		m3.model <-
    		c("GY94+M3[2]+F1xCodon", "GY94+M3[3]+F1xCodon", "GY94+M3[4]+F1xCodon", "GY94+M3[5]+F1xCodon", 
    		  "GY94+M3[2]+F1x4", "GY94+M3[3]+F1x4", "GY94+M3[4]+F1x4", "GY94+M3[5]+F1x4", 
    		  "GY94+M3[2]+F3x4", "GY94+M3[3]+F3x4", "GY94+M3[4]+F3x4", "GY94+M3[5]+F3x4")
		
		m3.N <- length(m3.model)
		
		m3.nfp <-
    		c(64, 66, 68, 70, 7, 9, 11, 13, 13, 15, 17, 19)

		corr.m3 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

	if (toupper(genetic.code) == "STANDARD") {
		m3.nfp <- m3.nfp - 0 * corr.m3 } 
	if (toupper(genetic.code) == "TVMC") {
		m3.nfp <- m3.nfp - corr.m3 }
	if (toupper(genetic.code) == "TYMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THMPCMCMSC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THIMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THCDHNC") {
		m3.nfp <- m3.nfp + 2 * corr.m3 }
	if (toupper(genetic.code) == "THEFMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THENC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THBAPPC") {
		m3.nfp <- m3.nfp + 0 * corr.m3 }
	if (toupper(genetic.code) == "THAYNC") {
		m3.nfp <- m3.nfp + 0 * corr.m3 }
	if (toupper(genetic.code) == "THAMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THAFMC") {
		m3.nfp <- m3.nfp + 2 * corr.m3 }
	if (toupper(genetic.code) == "BLNC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "CHMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "TRMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "SCOMC") {
		m3.nfp <- m3.nfp + 0 * corr.m3 }
	if (toupper(genetic.code) == "THMC") {
		m3.nfp <- m3.nfp - corr.m3 }

  		mdl.m3 <- paste("-m", rep("GY", 12))
  		wparam.m3 <- paste("-w", rep("DMODEL", 12))
  		ratevar.m3 <- rep(c("-c 1"), 12)
  		freqmodel.m3 <- paste("--fmodel", c(rep(c("F1XCODONS", "F1X4", "F3X4"), c(4,4,4))))
  		wclasses.m3 <- paste("--wclasses", c(rep(c(2, 3, 4, 5), 3)))
  		}
  		
if (omega.para == "M3" && model.set == "all" && refine.m3 == "TRUE")  {

		m3.model <-
    		c("GY94+M3[2]+F1xCodon", "GY94+M3[3]+F1xCodon", "GY94+M3[4]+F1xCodon", "GY94+M3[5]+F1xCodon", 
    		  "GY94+M3[2]+F1x4", "GY94+M3[3]+F1x4", "GY94+M3[4]+F1x4", "GY94+M3[5]+F1x4", 
    		  "GY94+M3[2]+F3x4", "GY94+M3[3]+F3x4", "GY94+M3[4]+F3x4", "GY94+M3[5]+F3x4",
    		  "MG+M3[2]+F1x4", "MG+M3[3]+F1x4", "MG+M3[4]+F1x4", "MG+M3[5]+F1x4",
    		  "MG+M3[2]+F3x4", "MG+M3[3]+F3x4", "MG+M3[4]+F3x4", "MG+M3[5]+F3x4")
		
		m3.N <- length(m3.model)

	m3.nfp <-
    		c(64, 66, 68, 70, 7, 9, 11, 13, 13, 15, 17, 19, 7, 9, 11, 13, 13, 15, 19, 21)

	corr.m3 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
		
	if (toupper(genetic.code) == "STANDARD") {
		m3.nfp <- m3.nfp - 0 * corr.m3 } 
	if (toupper(genetic.code) == "TVMC") {
		m3.nfp <- m3.nfp - corr.m3 }
	if (toupper(genetic.code) == "TYMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THMPCMCMSC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THIMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THCDHNC") {
		m3.nfp <- m3.nfp + 2 * corr.m3 }
	if (toupper(genetic.code) == "THEFMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THENC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THBAPPC") {
		m3.nfp <- m3.nfp + 0 * corr.m3 }
	if (toupper(genetic.code) == "THAYNC") {
		m3.nfp <- m3.nfp + 0 * corr.m3 }
	if (toupper(genetic.code) == "THAMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "THAFMC") {
		m3.nfp <- m3.nfp + 2 * corr.m3 }
	if (toupper(genetic.code) == "BLNC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "CHMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "TRMC") {
		m3.nfp <- m3.nfp + corr.m3 }
	if (toupper(genetic.code) == "SCOMC") {
		m3.nfp <- m3.nfp + 0 * corr.m3 }
	if (toupper(genetic.code) == "THMC") {
		m3.nfp <- m3.nfp - corr.m3 }
	
  		mdl.m3 <- paste("-m", rep(c("GY","MG"), c(12,8)))
    	wparam.m3 <- paste("-w", rep("DMODEL", 20))
   		ratevar.m3 <- rep(c("-c 1"), 20)
    	freqmodel.m3 <- paste("--fmodel", c(rep(c("F1XCODONS", "F1X4", "F3X4"), 4), rep(c("F1X4", "F3x4"), 4)))
  		wclasses.m3 <- paste("--wclasses", c(rep(c(2, 3, 4, 5), 5)))

    	}
	iters <- rep("--optHeuristic 4,1", N)	## These values are for a quick-and-dirty model testing.
	tstv <- rep("-t e", N)
    alpha <- rep("-a e", N)    
    freqopt <- rep("-f e", N)
      
    cmd.m3 <- paste(execname, "-i", seqfile, fmt, gcode, boot, tstv, mdl.m3, wparam.m3, wclasses.m3, 
    				alpha, ratevar.m3, freqmodel.m3, freqopt, i.tree.dist, iters, "-o lr", "--append")
    outfile.m3 <- paste(seqfile, "_codonphyml_stats.txt", sep = "")
    if (!append) {
        unlink(outfile.m3)
        unlink(paste(seqfile, "_codonphyml_tree.txt", sep = ""))
    }
    
    imod.m3 <- 1 : m3.N
    
	message("\n\n Fitting K = 2, 3, 4, 5 for M3 model.\n")
	
    for (i in imod.m3) { 
    		message(paste("Testing M3 subvariant", i, "of", m3.N, sep = " "))

			if (os == "Darwin")	{
				system(cmd[i], intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE) 
							} else {
							system(cmd[i], intern = FALSE, invisible = TRUE, show.output.on.console = FALSE) }
			}

    l.m3 <- readLines(outfile.m3, warn = FALSE)
    l.m3 <- grep("\\.\\s*Log-likelihood:", l, value = TRUE)
    
      if (dd <- length(l.m3) - length(imod.m3)) l.m3 <- l[-(1:dd)]
 	 loglik.m3 <- as.numeric(sub(". Log-likelihood:", "", l.m3))

# Now, rename the LL-containing file generated by M3 refinement to clean up file names. 
	
	file.rename(outfile.m3, paste(seqfile, "_model_fit_m3.txt", sep = ""))
	 
m3.table <- matrix(nrow = length(m3.model), ncol = length(c("k", "Log-Lik", "AIC", "AICc", "BIC")))

    m3.table[,1] <- m3.nfp    ## N parameters
    m3.table[,2] <- loglik.m3
    m3.table[,3] <- 2 * (m3.nfp - loglik.m3)      ## AIC
    m3.table[,4] <- 2 * (m3.nfp - loglik.m3) + 2 * ((m3.nfp * (m3.nfp + 1))) / (nsites - m3.nfp - 1)     ## AICc
    m3.table[,5] <- -2 * loglik.m3 + (m3.nfp * log(nsites))    ## BIC
    m3.table <- as.data.frame(m3.table, row.names = m3.model)
	names(m3.table) <- c("k", "Log-Lik", "AIC", "AICc", "BIC")

if(is.null(method)) {

	if (os =="Darwin") {
	
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		}
}


if (toupper(method) == "AIC") 
		{	delta.aic.m3 <- m3.table$AIC - min(m3.table$AIC)
			weights.m3 <- round(exp(-0.5*delta.aic.m3) / sum(exp(-0.5*delta.aic.m3)), digits=10)
			m3.table <- cbind(m3.table, delta.aic.m3, weights.m3)
			m3.table <- m3.table[with(m3.table, order(-weights.m3)), ]
			m3.table <- m3.table[m3.table$weights.m3!=0, ]
			
			if (length(m3.table$weights.m3) == 1) {
				message("A single model has weight 1.0. Check parameter estimates for possible convergence issues.")
					} else if (length(m3.table$weights.m3 > 1 && abs(m3.table$delta.aic.m3[1] - m3.table$delta.aic.m3[2])) <= 4) {
				warning("MODEL: Supports separating the top two models do not exceed the recommended threshold.") 
				} 
			
			
	if (os =="Darwin") {
	
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		}


		}
		
if (toupper(method) == "AICC")
		{	delta.aicc.m3 <- m3.table$AICc - min(m3.table$AICc)
			weights.m3 <- round(exp(-0.5*delta.aicc.m3) / sum(exp(-0.5*delta.aicc.m3)), digits=10)
			m3.table <- cbind(m3.table, delta.aicc.m3, weights.m3)
			m3.table <- m3.table[with(m3.table, order(-weights.m3)), ]
			m3.table <- m3.table[m3.table$weights.m3!=0, ]
			
			if (length(m3.table$weights.m3) == 1) {
				message("A single model has weight 1.0. Check parameter estimates for possible convergence issues.\n")
					} else if (length(m3.table$weights.m3 > 1 && abs(m3.table$delta.aicc.m3[1] - m3.table$delta.aicc.m3[2])) <= 4) {
				warning("\nMODEL: Supports separating the top two models do not exceed the recommended threshold.") 
				} 
			
	if (os =="Darwin") {
	
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		}


		}

if (toupper(method) == "BIC")
		{	delta.bic.m3 <- m3.table$BIC - min(m3.table$BIC)
			weights.m3 <- round(exp(-0.5*delta.bic.m3) / sum(exp(-0.5*delta.bic.m3)), digits=10)
			m3.table <- cbind(m3.table, delta.bic.m3, weights.m3)
			m3.table <- m3.table[with(m3.table, order(-weights.m3)), ]
			m3.table <- m3.table[m3.table$weights.m3!=0, ]
			
			if (length(m3.table$weights.m3) == 1) {
				message("A single model has weight 1.0. Check parameter estimates for possible convergence issues.")
					} else if (length(m3.table$weights.m3 > 1 && abs(m3.table$delta.bic.m3[1] - m3.table$delta.bic.m3[2])) <= 4) {
				warning("MODEL: Supports separating the top two models do not exceed the recommended threshold.") 
				} 
			
	if (os =="Darwin") {
	
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		} else {
		
		write.table(m3.table, file = paste(getwd(), "/", results.dir, "/", "m3refined.out", sep=""), sep = "\t", append = FALSE)
		
		}


		}
  
best.m3 <- row.names(m3.table[1,])

cat("------------------------------------\n")
cat("      CodonModelTest Results        \n")
cat("     (M3 subvariant optimized)      \n")
cat("------------------------------------\n\n")
cat(paste("Best-fit model is an M3 subvariant:", best.m3, "\n\n", sep=" ")) 

elements.m3 <- str_split(best.m3, "\\+")

if (length(elements.m3[[1]]) < 4) { 
	gamma.para.m3 <- "Gamma[1]"
	freq.para.m3 <- elements.m3[[1]][3] } else { 
	gamma.para.m3 <- elements.m3[[1]][3]
	freq.para.m3 <- elements.m3[[1]][4] }

omega.para.m3 <- elements.m3[[1]][2]
   
}
 

if (g.para == "GY94") {

# Define best-fit model elements (as a string to parse):

# Set up software-specific parsing of genetic codes and parameters:

# Assign rj-MCMC procedure in MrBayes 3.2.x to nst, if deltaIC < 4:
	
	if (length(ic.table[,6] == 1)) states <- 1					## Regular GY has nst=1
	else if (abs(ic.table[,6][1] - ic.table[,6][2]) <= 4)	{ 
				states <- "mixed" 
				} else { states <- 2 }

if (toupper(genetic.code) == "STANDARD") gcode.mb <- "Universal"
else if (toupper(genetic.code) == "TVMC") gcode.mb <- "Vertmt"
else if (toupper(genetic.code) == "THMPCMCMSC") gcode.mb <- "Mycoplasma"
else if (toupper(genetic.code) == "THIMC") gcode.mb <- "Metmt"
else if (toupper(genetic.code) == "THCDHNC") gcode.mb <- "Ciliates"
else if (toupper(genetic.code) == "TYMC") gcode.mb <- "Yeast"

else { gcode.mb <- toupper(genetic.code)
	warning("GCODE: Genetic code provided is more specific than codes allowed in MrBayes.", 
				"Gcode has been directly written to MrBayes block, but an approximation may be", 
				"required. Check translation products are correct if using an approximation.")
}

if (omega.para == "M5") stop("M5 model not (yet) implemented in phylogenetic search programs. No model blocks written.")
if (omega.para == "M3" && refine.m3 == "FALSE") message("M3 model comparison conducted here assumes K = 3 as a coarse approximation. ", 
									"Check a variety of K values to further optimize model selection within M3.")

if (freq.para == "F1x4" || freq.para == "F3x4") {		## For F1x4, or F3x4, use empirical frequencies.
	mb.freq <- "fixed(empirical)" } else {
	mb.freq <- "dirichlet(1.0)"							## Otherwise, is F1x61, and frequencies are estimated.
	}
									
if (omega.para == "M0") omega.mb <- "equal"
if (omega.para == "M3") omega.mb <- "M3"

if (gamma.para == "Gamma[1]") {
	rates.mb <- "equal" } else { rates.mb <- "gamma" }

if (omega.para == "M0") {
		omega.garli <- "none"
		nrates.garli <- 1 
			} else {
				omega.garli <- "nonsynonymous"
				nrates.garli <- as.numeric(str_split(omega.para.m3, "")[[1]][5])
		}
		
if (toupper(genetic.code) == "STANDARD") gcode.garli <- "standard"
else if (toupper(genetic.code) == "TVMC") gcode.garli <- "vertmito"
else if (toupper(genetic.code) == "THIMC") gcode.garli <- "invertmito"
else { gcode.garli <- toupper(genetic.code)
warning("GCODE: Genetic code is more specific than codes allowed by GARLI.", 
			"Gcode has been directly written to GARLI conf, but an approximation may be required.", 
			"Check translation products are correct if using approximation.")
}


# Write MrBayes and Garli blocks, using best-fit model. Blocks written 
# assuming user implements BEAGLE library, and that the user will be running
# analysis on GPU:

# MRBAYES:

	if (os =="Darwin") {
	
		sink(file = paste(getwd(), "/", results.dir, "/", "codonmodel.mb", sep=""), append = FALSE, split = FALSE)
		
		} else {
		
		sink(file = paste(getwd(), "/", results.dir, "/", "codonmodel.mb", sep=""), append = FALSE, split = FALSE)
		
		}

cat("begin mrbayes;\n\n")
cat("set usebeagle = yes beagledevice = gpu beaglescaling = dynamic beaglesse = yes autoclose = no nowarn = no;\n")
cat(paste("lset nucmodel = codon nst =", states, "omegavar =", omega.mb, "rates=", rates.mb, "Ngammacat=", 
			as.numeric(str_extract(gamma.para, "[0-9]")), "code=", gcode.mb, ";\n"), sep="")
cat("unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) omega=(all) ratemultiplier=(all) topology=(all) brlens=(all);\n") 
cat(paste("prset statefreqpr=", mb.freq, "applyto=(all);\n"), sep="")
cat("mcmcp ngen= 100000000 relburnin=yes burninfrac=0.5  printfreq=10 samplefreq=1000 nchains=4 savebrlens=yes;\n")
cat("mcmc;\n")

sink()

# GARLI:

	if (os =="Darwin") {
	
		sink(file = paste(getwd(), "/", results.dir, "/", "codonmodel.conf", sep=""), append = FALSE, split = FALSE)
		
		} else {
		
		sink(file = paste(getwd(), "/", results.dir, "/", "codonmodel.conf", sep=""), append = FALSE, split = FALSE)
		
		}

cat(paste("[", strsplit(basename(seqfile), "\\.")[[1]][1], " codon model]\n", sep=""))
cat("datatype = codon\n")
cat("ratematrix = 1rate\n")
cat(paste("statefrequencies =", freq.para, "\n"))
cat(paste("ratehetmodel =", omega.garli, "\n"))
cat("numratecats =", nrates.garli, "\n")
cat("invariantsites = none\n")
cat(paste("geneticcode =", gcode.garli, "\n"))

sink()
flush.console()

} else { 

warning("MODELTYPE: Best-fit model is an MG94 variant. MG parameterization is not (yet) implemented in traditional ", 
			"topology search programs. No model blocks written.") 

}

message(paste("Output written to ", getwd(), "/", results.dir, " folder. \n\n", sep = "")) 

message("F3x4 calculations have been conducted by the standard method to maintain compatibility with GARLI. ", 
		"Check model optimality using corrected frequencies (CF3x4) in CodonPhyML or HyPhy.")

rm(list=ls())
		
}
