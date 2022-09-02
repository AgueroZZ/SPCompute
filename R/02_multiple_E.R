#' Compute the required power for binary response, a SNP G and two continuous covariates that are conditionally independent given G, using the Semi-Sim method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample to approximate the fisher information matrix, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE, muE, sigmaE and gammaG have to be vectors of length 2.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param seed An integer number that indicates the seed used for the simulation to compute the approximate fisher information matrix, by default is 123.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param LargePowerApproxi TRUE or FALSE indicates whether to use the large power approximation formula.
#' @return The power that can be achieved at the given sample size (using semi-sim method).
#' @noRd
Compute_Power_Sim_BCC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8, LargePowerApproxi = FALSE){
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (2*pG*qG + 2*pG^2)
    gamma02 <- muE2 - gammaG2 * (2*pG*qG + 2*pG^2)


    betaG <- parameters$betaG
    varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2

    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}


    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)


    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)


    ### Simulate for SE: by averaging B times
    set.seed(seed)
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = 1, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError2)

      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
      weight <- stats::dlogis(eta)
      I <- I + weight* X %*% t(X)
    }
    I <- I/B
    if(LargePowerApproxi){
      SE <- sqrt((solve(I)[2,2]))/sqrt(n)
      return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
    }
    else{
      compute_power <- function(n){
        ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
        SE = sqrt((solve(I)[2,2]))/sqrt(n)
        ### Once know this SE of betaG hat, compute its power at this given sample size n:
        Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
        Power
      }
      compute_power(n)
    }
  }
  else if(mode == "dominant"){
    preva <- parameters$preva
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)
    betaG <- parameters$betaG

    varG <- pG*qG
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}

    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)

    solveForbeta0_dom_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_dom_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    ### Simulate for SE: by averaging B times
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    set.seed(seed)

    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError2)
      X <- matrix(c(1,G,E1,E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
      weight <- stats::dlogis(eta)
      I <- I + weight* X %*% t(X)
    }
    I <- I/B
    if(LargePowerApproxi){
      SE <- sqrt((solve(I)[2,2]))/sqrt(n)
      return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
    }
    else{
      compute_power <- function(n){
        ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
        SE = sqrt((solve(I)[2,2]))/sqrt(n)
        ### Once know this SE of betaG hat, compute its power at this given sample size n:
        Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
        Power
      }
      compute_power(n)
    }
  }
  else if(mode == "recessive") {
    preva <- parameters$preva
    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    pG <- parameters$pG^2
    qG <- 1 - pG

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)

    betaG <- parameters$betaG

    if((sigmaE1^2) <= (gammaG1^2) * qG*pG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * qG*pG)
    if((sigmaE2^2) <= (gammaG2^2) * qG*pG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * qG*pG)


    solveForbeta0_rec_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_rec_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    ### Simulate for SE: by averaging B times
    I <- matrix(data = 0, nrow = 4, ncol = 4)
    set.seed(seed)

    for (i in 1:B) {
      G <- sample(c(0,1), size = 1, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(1,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(1,sd = sigmaError2)

      X <- matrix(c(1,G,E1, E2), ncol = 1)
      eta <- beta0 + betaG*G + betaE1*E1 + betaE2*E2
      weight <- stats::dlogis(eta)
      I <- I + weight* X %*% t(X)
    }
    I <- I/B
    if(LargePowerApproxi){
      SE <- sqrt((solve(I)[2,2]))/sqrt(n)
      return(stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE))
    }
    else{
      compute_power <- function(n){
        ### Once know this (expected) single I, scale by dividing sqrt(n): n is sample size
        SE = sqrt((solve(I)[2,2]))/sqrt(n)
        ### Once know this SE of betaG hat, compute its power at this given sample size n:
        Power = stats::pnorm(-stats::qnorm(1-(alpha/2)) + betaG/SE ) + stats::pnorm(-stats::qnorm(1-(alpha/2)) - betaG/SE)
        Power
      }
      compute_power(n)
    }
  }
}




#' Compute the required power for binary response, a SNP G and two continuous covariates that are conditionally independent given G, using the empirical method.
#'
#' @param n An integer number that indicates the sample size.
#' @param B An integer number that indicates the number of simulated sample, by default is 10000.
#' @param parameters Refer to SPCompute::Compute_Power_Sim; Except betaE, muE, sigmaE and gammaG have to be vectors of length 2.
#' @param mode A string of either "additive", "dominant" or "recessive", indicating the genetic mode, by default is "additive".
#' @param alpha A numeric value that denotes the significance level used in the study, by default is 0.05.
#' @param searchSizeBeta0 The interval radius for the numerical search of beta0, by default is 8. Setting to higher values may solve some numerical problems at the cost of longer runtime.
#' @param seed An integer number that indicates the seed used for the simulation, by default is 123.
#' @return The power that can be achieved at the given sample size (computed from empirical power).
#' @noRd
Compute_Power_Emp_BCC <- function(n, B = 10000, parameters, mode = "additive", alpha = 0.05, seed = 123, searchSizeBeta0 = 8){
  if(mode == "additive"){
    preva <- parameters$preva
    pG <- parameters$pG
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (2*pG*qG + 2*pG^2)
    gamma02 <- muE2 - gammaG2 * (2*pG*qG + 2*pG^2)

    betaG <- parameters$betaG
    varG <- (2*pG*qG + 4*pG^2) - (2*pG*qG + 2*pG^2)^2

    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}

    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)

    solveForbeta0_add_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1,2), size = B, replace = TRUE, prob = c(qG^2, 2*pG*qG, pG^2))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_add_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    set.seed(seed)
    correct <- c()
    for (i in 1:B) {
      G <- sample(c(0,1,2), size = n, replace = TRUE, prob = c(qG^2,2*pG*qG, pG^2))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError2)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(glm(y~G + E1 + E2, family = binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  else if(mode == "dominant"){
    preva <- parameters$preva
    pG <- 2*parameters$pG*(1-parameters$pG) + parameters$pG^2
    qG <- 1 - pG

    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)
    betaG <- parameters$betaG
    varG <- qG*pG
    if((sigmaE1^2) <= (gammaG1^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}
    if((sigmaE2^2) <= (gammaG2^2) * varG){return(message("Error: SigmaE must be larger to be compatible with other parameters"))}

    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * varG)
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * varG)

    solveForbeta0_dom_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_dom_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    set.seed(seed)
    correct <- c()
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError2)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(glm(y~G + E1 + E2, family = binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B

  }
  else if(mode == "recessive") {
    preva <- parameters$preva
    gammaG1 <- parameters$gammaG[1]
    muE1 <- parameters$muE[1]
    sigmaE1 <- parameters$sigmaE[1]
    betaE1 <- parameters$betaE[1]

    gammaG2 <- parameters$gammaG[2]
    muE2 <- parameters$muE[2]
    sigmaE2 <- parameters$sigmaE[2]
    betaE2 <- parameters$betaE[2]

    pG <- parameters$pG^2
    qG <- 1 - pG

    gamma01 <- muE1 - gammaG1 * (pG)
    gamma02 <- muE2 - gammaG2 * (pG)

    betaG <- parameters$betaG

    if((sigmaE1^2) <= (gammaG1^2) * qG*pG){return(message("Error: SigmaE[1] must be larger to be compatible with other parameters"))}
    sigmaError1 <- sqrt(sigmaE1^2 - (gammaG1^2) * qG*pG)
    if((sigmaE2^2) <= (gammaG2^2) * qG*pG){return(message("Error: SigmaE[2] must be larger to be compatible with other parameters"))}
    sigmaError2 <- sqrt(sigmaE2^2 - (gammaG2^2) * qG*pG)


    solveForbeta0_rec_con <- function(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2){
      qG <- 1 - pG
      ComputeP <- function(beta0){
        set.seed(seed)
        G <- sample(c(0,1), size = B, replace = TRUE, prob = c(qG,pG))
        E1 <- gamma01 + gammaG1 * G + stats::rnorm(B, sd = sigmaError1)
        E2 <- gamma02 + gammaG2 * G + stats::rnorm(B, sd = sigmaError2)
        y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(B)
        P <- mean(ifelse(y > 0, 1, 0))
        P - preva
      }
      stats::uniroot(ComputeP, c(-searchSizeBeta0, searchSizeBeta0))$root
    }
    beta0 <- solveForbeta0_rec_con(preva, betaG, betaE1, betaE2, pG, gammaG1, gammaG2)

    set.seed(seed)
    correct <- c()
    for (i in 1:B) {
      G <- sample(c(0,1), size = n, replace = TRUE, prob = c(qG, pG))
      E1 <- gamma01 + gammaG1*G + stats::rnorm(n,sd = sigmaError1)
      E2 <- gamma02 + gammaG2*G + stats::rnorm(n,sd = sigmaError2)
      y <- beta0 + betaG * G + betaE1 * E1 + betaE2 * E2 + stats::rlogis(n)
      y <- ifelse(y > 0, 1, 0)
      correct[i] <- summary(glm(y~G + E1 + E2, family = binomial("logit")))$coefficients[2,4] <= alpha
    }
    Power <- sum(correct)/B
  }
  Power
}







#
# Compute_Power_Sim_BCC(n = 500, B = 10000,
#                       parameters = list(preva = 0.2, betaG = 0.6, betaE = c(0.9,0.5), gammaG = c(0.2,0.1), muE = c(0,0.1), sigmaE = c(1,2), pG = 0.3),
#                       mode = "additive")
#
#
#
# Compute_Power_Emp_BCC(n = 500, B = 1000,
#                       parameters = list(preva = 0.2, betaG = 0.6, betaE = c(0.9,0.5), gammaG = c(0.2,0.1), muE = c(0,0.1), sigmaE = c(1,2), pG = 0.3),
#                       mode = "additive")




