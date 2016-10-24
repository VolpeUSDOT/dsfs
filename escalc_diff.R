function (measure, formula, ai, bi, ci, di, n1i, n2i, x1i, x2i, 
          t1i, t2i, m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, ni, 
          data, slab, subset, add = 1/2, to = "only0", drop00 = FALSE, 
          vtype = "LS", var.names = c("yi", "vi"), add.measure = FALSE, 
          append = TRUE, replace = TRUE, digits = 4, ...) 


  if (is.element(measure, c("MD", "SMD", "SMDH", "ROM", "RPB", 
                            "RBIS", "D2OR", "D2ORN", "D2ORL"))) {
    mf.m1i <- mf[[match("m1i", names(mf))]]
    mf.m2i <- mf[[match("m2i", names(mf))]]
    mf.sd1i <- mf[[match("sd1i", names(mf))]]
    mf.sd2i <- mf[[match("sd2i", names(mf))]]
    mf.n1i <- mf[[match("n1i", names(mf))]]
    mf.n2i <- mf[[match("n2i", names(mf))]]
    m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
    m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
    sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
    sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
    n1i <- eval(mf.n1i, data, enclos = sys.frame(sys.parent()))
    n2i <- eval(mf.n2i, data, enclos = sys.frame(sys.parent()))
    if (!is.null(subset)) {
      m1i <- m1i[subset]
      m2i <- m2i[subset]
      sd1i <- sd1i[subset]
      sd2i <- sd2i[subset]
      n1i <- n1i[subset]
      n2i <- n2i[subset]
    }
    if (length(m1i) == 0L || length(m2i) == 0L || length(sd1i) == 
        0L || length(sd2i) == 0L || length(n1i) == 0L || 
        length(n2i) == 0L) 
      stop("Cannot compute outcomes. Check that all of the required \n  information is specified via the appropriate arguments.")
    if (!all(length(m1i) == c(length(m1i), length(m2i), length(sd1i), 
                              length(sd2i), length(n1i), length(n2i)))) 
      stop("Supplied data vectors are not all of the same length.")
    if (any(c(sd1i, sd2i) < 0, na.rm = TRUE)) 
      stop("One or more standard deviations are negative.")
    if (any(c(n1i, n2i) < 0, na.rm = TRUE)) 
      stop("One or more sample sizes are negative.")
    ni.u <- n1i + n2i
    k <- length(m1i)
    ni <- ni.u
    mi <- ni - 2
    sdpi <- sqrt(((n1i - 1) * sd1i^2 + (n2i - 1) * sd2i^2)/mi)
    di <- (m1i - m2i)/sdpi
    if (measure == "MD") {
      yi <- m1i - m2i
      if (length(vtype) == 1L) 
        vtype <- rep(vtype, k)
      vi <- rep(NA_real_, k)
      for (i in seq_len(k)) {
        if (vtype[i] == "UB" || vtype[i] == "LS") 
          vi[i] <- sd1i[i]^2/n1i[i] + sd2i[i]^2/n2i[i]
        if (vtype[i] == "HO") 
          vi[i] <- sdpi[i]^2 * (1/n1i[i] + 1/n2i[i])
      }
    }
    if (measure == "SMD") {
      warn.before <- getOption("warn")
      options(warn = -1)
      cmi <- .cmicalc(mi)
      options(warn = warn.before)
      yi <- cmi * di
      if (length(vtype) == 1L) 
        vtype <- rep(vtype, k)
      vi <- rep(NA_real_, k)
      mnwyi <- sum(ni * yi, na.rm = TRUE)/sum(ni, na.rm = TRUE)
      for (i in seq_len(k)) {
        if (vtype[i] == "UB") 
          vi[i] <- 1/n1i[i] + 1/n2i[i] + (1 - (mi[i] - 
                                                 2)/(mi[i] * cmi[i]^2)) * yi[i]^2
        if (vtype[i] == "LS") 
          vi[i] <- 1/n1i[i] + 1/n2i[i] + yi[i]^2/(2 * 
                                                    ni[i])
        if (vtype[i] == "HO") 
          vi[i] <- 1/n1i[i] + 1/n2i[i] + mnwyi^2/(2 * 
                                                    ni[i])
      }
    }
    if (measure == "SMDH") {
      warn.before <- getOption("warn")
      options(warn = -1)
      cmi <- .cmicalc(mi)
      options(warn = warn.before)
      si <- sqrt((sd1i^2 + sd2i^2)/2)
      yi <- cmi * (m1i - m2i)/si
      vi <- yi^2 * (sd1i^4/(n1i - 1) + sd2i^4/(n2i - 1))/(2 * 
                                                            (sd1i^2 + sd2i^2)^2) + (sd1i^2/(n1i - 1) + sd2i^2/(n2i - 
                                                                                                                 1))/((sd1i^2 + sd2i^2)/2)
      vi <- cmi^2 * vi
    }
    if (measure == "ROM") {
      yi <- log(m1i/m2i)
      if (length(vtype) == 1L) 
        vtype <- rep(vtype, k)
      vi <- rep(NA_real_, k)
      for (i in seq_len(k)) {
        if (vtype[i] == "LS") 
          vi[i] <- sd1i[i]^2/(n1i[i] * m1i[i]^2) + sd2i[i]^2/(n2i[i] * 
                                                                m2i[i]^2)
        if (vtype[i] == "HO") 
          vi[i] <- sdpi[i]^2/(n1i[i] * m1i[i]^2) + sdpi[i]^2/(n2i[i] * 
                                                                m2i[i]^2)
      }
    }
    if (is.element(measure, c("RPB", "RBIS"))) {
      hi <- mi/n1i + mi/n2i
      yi <- di/sqrt(di^2 + hi)
      if (measure == "RPB") {
        if (length(vtype) == 1L) 
          vtype <- rep(vtype, k)
        vi <- rep(NA_real_, k)
        for (i in seq_len(k)) {
          if (vtype[i] == "ST" || vtype[i] == "LS") 
            vi[i] <- hi[i]^2/(hi[i] + di[i]^2)^3 * (1/n1i[i] + 
                                                      1/n2i[i] + di[i]^2/(2 * ni[i]))
          if (vtype[i] == "CS") 
            vi[i] <- (1 - yi[i]^2)^2 * (ni[i] * yi[i]^2/(4 * 
                                                           n1i[i] * n2i[i]) + (2 - 3 * yi[i]^2)/(2 * 
                                                                                                   ni[i]))
        }
      }
    }
    if (measure == "RBIS") {
      p1i <- n1i/ni
      p2i <- n2i/ni
      zi <- qnorm(p1i, lower.tail = FALSE)
      fzi <- dnorm(zi)
      yi <- sqrt(p1i * p2i)/fzi * yi
      yi.t <- ifelse(abs(yi) > 1, sign(yi), yi)
      vi <- 1/(ni - 1) * (p1i * p2i/fzi^2 - (3/2 + (1 - 
                                                      p1i * zi/fzi) * (1 + p2i * zi/fzi)) * yi.t^2 + 
                            yi.t^4)
    }
    if (is.element(measure, c("D2OR", "D2ORL"))) {
      yi <- pi/sqrt(3) * di
      vi <- pi^2/3 * (1/n1i + 1/n2i + di^2/(2 * ni))
    }
    if (measure == "D2ORN") {
      yi <- 1.65 * di
      vi <- 1.65^2 * (1/n1i + 1/n2i + di^2/(2 * ni))
    }
  }
 


  if (is.element(measure, c("MC", "SMCC", "SMCR", "SMCRH", 
                            "ROMC"))) {
    mf.m1i <- mf[[match("m1i", names(mf))]]
    mf.m2i <- mf[[match("m2i", names(mf))]]
    mf.sd1i <- mf[[match("sd1i", names(mf))]]
    mf.sd2i <- mf[[match("sd2i", names(mf))]]
    mf.ni <- mf[[match("ni", names(mf))]]
    mf.ri <- mf[[match("ri", names(mf))]]
    m1i <- eval(mf.m1i, data, enclos = sys.frame(sys.parent()))
    m2i <- eval(mf.m2i, data, enclos = sys.frame(sys.parent()))
    sd1i <- eval(mf.sd1i, data, enclos = sys.frame(sys.parent()))
    sd2i <- eval(mf.sd2i, data, enclos = sys.frame(sys.parent()))
    ni <- eval(mf.ni, data, enclos = sys.frame(sys.parent()))
    ri <- eval(mf.ri, data, enclos = sys.frame(sys.parent()))
    if (!is.null(subset)) {
      m1i <- m1i[subset]
      m2i <- m2i[subset]
      sd1i <- sd1i[subset]
      sd2i <- sd2i[subset]
      ni <- ni[subset]
      ri <- ri[subset]
    }
   


     
  class(dat) <- c("escalc", "data.frame")
  return(dat)
}