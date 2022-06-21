dea <- function(base = NULL, frontier = NULL,
                noutput = 1, orientation=1, rts = 1, onlytheta =
                    FALSE) {

  if(is.null(frontier))
    frontier <- base

  if(!is.null(base) & !is.null(frontier)){
    if(is.null(nrow(base)))
      base <- matrix(base, nrow = 1)
    else
      base <- as.matrix(base)

    frontier <- as.matrix(frontier)
  }

  if(ncol(base) != ncol(frontier))
    stop("Number of columns in base matrix and frontier matrix should be the same!")

  s <- noutput
  m <- ncol(base) - s
  n <- nrow(base)
  nf <- nrow(frontier)

  Y.f <- t(frontier[, 1:s])
  X.f <- t(frontier[, (s+1):(s+m)])
  if(n == 1){
    Y.b <- matrix(base[, 1:s], ncol = 1)
    X.b <- matrix(base[, (s+1):(s+m)], ncol = 1)
  }
  else{
    Y.b <- t(base[, 1:s])
    X.b <- t(base[, (s+1):(s+m)])
  }

  re <- list()

  first.obj <- c(1, rep(0, nf))
  first.dir <- rep(">=", m + s)
  for(i in 1:n){
    if(onlytheta == TRUE){
      if(rts == 1){
        if(orientation == 1){
          ## input-oriented CRS
          obj <- first.obj
          dir <- first.dir
          rhs <- c(rep(0, m), Y.b[, i])
          con1 <- cbind(X.b[, i], -X.f)
          con2 <- cbind(rep(0, s), Y.f)
          con <- rbind(con1, con2)
          re[[i]] <- lp("min", obj, con, dir, rhs)$solution[1]
        }
        else{
          ## output-oriented CRS
          obj <- first.obj
          dir <- first.dir
          rhs <- c(-X.b[, i], rep(0, s))
          con1 <- cbind(rep(0, m), -X.f)
          con2 <- cbind(-Y.b[, i], Y.f)
          con <- rbind(con1, con2)
          re[[i]] <- lp("max", obj, con, dir, rhs)$solution[1]
        }
      }
      else{
        if(orientation == 1){
          ## input-oriented VRS
          obj <- first.obj
          dir <- c(first.dir, "==")
          rhs <- c(rep(0, m), Y.b[, i], 1)
          con1 <- cbind(X.b[, i], -X.f)
          con2 <- cbind(rep(0, s), Y.f)
          con3 <- c(0, rep(1, nf))
          con <- rbind(con1, con2, con3)
          re[[i]] <- lp("min", obj, con, dir, rhs)$solution[1]
        }
        else{
          ## output-oriented VRS
          obj <- first.obj
          dir <- c(first.dir, "==")
          rhs <- c(-X.b[, i], rep(0, s), 1)
          con1 <- cbind(rep(0, m), -X.f)
          con2 <- cbind(-Y.b[, i], Y.f)
          con3 <- c(0, rep(1, nf))
          con <- rbind(con1, con2, con3)
          re[[i]] <- lp("max", obj, con, dir, rhs)$solution[1]
        }
      }
    }

    else{ ## with lambda and slack
      if(rts == 1){
        if(orientation == 1){
          ## input-oriented CRS
          obj <- c(first.obj, rep(0, m + s))
          dir <- rep("==", m + s)
          rhs <- c(rep(0, m), Y.b[, i])
          con1 <- cbind(X.b[, i], -X.f, -diag(m), matrix(0, m, s))
          con2 <- cbind(rep(0, s), Y.f, matrix(0, s, m), -diag(s))
          con <- rbind(con1, con2)
          re[[i]] <- lp("min", obj, con, dir, rhs)$solution
        }
        else{
          ## output-oriented CRS
          obj <- c(first.obj, rep(0, m + s))
          dir <- rep("==", m + s)
          rhs <- c(-X.b[, i], rep(0, s))
          con1 <- cbind(rep(0, m), -X.f, -diag(m), matrix(0, m, s))
          con2 <- cbind(Y.b[, i], -Y.f, matrix(0, s, m), diag(s))
          con <- rbind(con1, con2)
          re[[i]] <- lp("max", obj, con, dir, rhs)$solution
        }
      }
      else{
        if(orientation == 1){
          ## input-oriented VRS
          obj <- c(first.obj, rep(0, m + s))
          dir <- rep("==", m + s + 1)
          rhs <- c(rep(0, m), Y.b[, i], 1)
          con1 <- cbind(X.b[, i], -X.f, -diag(m), matrix(0, m, s))
          con2 <- cbind(rep(0, s), Y.f, matrix(0, s, m), -diag(s))
          con3 <- c(0, rep(1, nf), rep(0, m + s))
          con <- rbind(con1, con2, con3)
          re[[i]] <- lp("min", obj, con, dir, rhs)$solution
        }
        else{
          ## output-oriented VRS
          obj <- c(first.obj, rep(0, m + s))
          dir <- rep("==", m + s + 1)
          rhs <- c(-X.b[, i], rep(0, s), 1)
          con1 <- cbind(rep(0, m), -X.f, -diag(m), matrix(0, m, s))
          con2 <- cbind(Y.b[, i], -Y.f, matrix(0, s, m), diag(s))
          con3 <- c(0, rep(1, nf), rep(0, m + s))
          con <- rbind(con1, con2, con3)
          re[[i]] <- lp("max", obj, con, dir, rhs)$solution
        }
      }
    }
  }
  re <- data.frame(do.call(rbind, re))
  
  if(onlytheta == TRUE){
    names(re) <- "eff"
  }
  else{
    names(re) <-
      c("eff", paste("lambda.", 1:nf, sep = ""),
        paste("slack.x", 1:m, sep = ""), paste("slack.y", 1:s, sep = ""))
  }
  return(re)
}
