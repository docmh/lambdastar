# Source this file, then run lambdastar_tutorial(). No analysis runs on source().
# See system.file("doc", "LEARN_LAMBDASTAR.md", package="lambdastar").

.ls_seed <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(code)
}

.ls_integer <- function(x, name, minimum) {
  if (length(x)!=1L || !is.numeric(x) || !is.finite(x) ||
      x!=floor(x) || x<minimum || x>.Machine$integer.max-10)
    stop(name, " must be an integer >= ", minimum, call.=FALSE)
}

# A controlled teaching design, not a fitted empirical measurement model.
# Four orthogonal profiles have sample variance 1. The fourth is unobserved.
# Only the measurement errors are redrawn after this design has been generated.
lambdastar_tutorial_data <- function(n=160L, items=4L, noise_sd=1.2, seed=17L) {
  .ls_integer(n,"n",30); .ls_integer(items,"items",2); .ls_integer(seed,"seed",0)
  if (length(noise_sd)!=1L || !is.numeric(noise_sd) || !is.finite(noise_sd) || noise_sd<=0)
    stop("noise_sd must be positive; the T=0 example is separate.",call.=FALSE)
  .ls_seed(seed, {
    raw <- matrix(stats::rnorm(n*4), n, 4)
    profiles <- qr.Q(qr(scale(raw, center=TRUE, scale=FALSE)))*sqrt(n-1)
    d <- data.frame(id=seq_len(n), workload=profiles[,1], resources=profiles[,2],
                    unused=profiles[,3])
    signal <- -.9*d$workload + .7*d$resources + .5*profiles[,4]
    errors <- matrix(stats::rnorm(n*items), n, items)
    indicators <- paste0("item",seq_len(items))
    d[indicators] <- signal + noise_sd*errors
    d$score <- rowMeans(d[indicators])
    list(data=d, indicators=indicators, signal=signal,
      settings=list(n=n, items=items, noise_sd=noise_sd, seed=seed),
      truth=c(temperature=noise_sd^2/items, variance_signal=1.55,
              kappa_workload=.81/1.55, kappa_full=1.30/1.55))
  })
}

# Oracle values use the known latent signal and generating noise variance.
# They are not estimates computed from observed items.
.ls_oracle <- function(example, variables) {
  d <- example$data; s <- example$signal-mean(example$signal)
  X <- cbind(1, as.matrix(d[variables]))
  fitted <- stats::lm.fit(X,s)$fitted.values
  total <- sum(s^2); captured <- sum((fitted-mean(fitted))^2)
  T <- unname(example$truth["temperature"])
  c(temperature=T, lambda_star=total/T, lambda_h=captured/T,
    lambda_perp=(total-captured)/T, kappa=captured/total)
}

.ls_markdown_table <- function(x) {
  cells <- lapply(x,function(v) if(is.numeric(v)) format(signif(v,6),trim=TRUE) else as.character(v))
  body <- apply(as.data.frame(cells,stringsAsFactors=FALSE),1,paste,collapse=" | ")
  c(paste0("| ",paste(names(x),collapse=" | ")," |"),
    paste0("| ",paste(rep("---",ncol(x)),collapse=" | ")," |"),paste0("| ",body," |"))
}

# Default bootstrap settings are 5000 draws and 95%; type7 is explicit so the
# lesson can verify the endpoints using a plain quantile() call.
lambdastar_tutorial <- function(output_dir="lambda-star-tutorial", n=160L,
    items=4L, noise_sd=1.2, seed=17L, B=5000L, conf.level=.95, T_q=1/12,
    pause=interactive(), show_plots=interactive()) {
  if (!requireNamespace("lambdastar",quietly=TRUE)) stop("Install lambdastar first.")
  .ls_integer(B,"B",20)
  cfg <- lambdastar::lambda_bootstrap_control(B=B,conf.level=conf.level,seed=seed+1L,
    min_valid=min(B,200L),keep_draws=TRUE,keep_indices=TRUE,percentile_interpolation="type7")
  example <- lambdastar_tutorial_data(n,items,noise_sd,seed)
  if (length(T_q)!=1L || !is.numeric(T_q) || !is.finite(T_q) || T_q<=0)
    stop("T_q must be a positive resolution variance on the mean-score scale.")
  for (flag in list(pause,show_plots))
    if (!is.logical(flag) || length(flag)!=1L || is.na(flag)) stop("pause/show_plots must be TRUE or FALSE.")
  if (!is.character(output_dir) || length(output_dir)!=1L || is.na(output_dir) || !nzchar(output_dir))
    stop("output_dir must name a directory.")
  # Never overwrite a previous lesson run or the user's files.
  if (dir.exists(output_dir) && length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
    stop("output_dir is not empty. Choose a new folder for this experiment.",call.=FALSE)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  figures <- file.path(output_dir,"figures");dir.create(figures)
  d <- example$data; indicators <- example$indicators
  blue <- "#286090"; gold <- "#B17A16"; ink <- "#30343B"; grey <- "#9AA1AA"
  report <- c("# Learning Lambda Star: your worked example", "",
    "Synthetic data; Lambda Star estimates use the installed package. Known values come from the explicit generating model.",
    sprintf("N=%d; K=%d; item-error SD=%.3g; seed=%d; B=%d; confidence level=%.3g; T_q=%.6g.",
      n,items,noise_sd,seed,B,conf.level,T_q),
    "Estimator: paper. Resampling: whole persons with replacement. Percentiles: type 7.","")
  add <- function(...) {report <<- c(report,...,"");cat(paste(c(...),collapse="\n"),"\n\n")}
  step <- function(title, question) {
    add(paste0("## ",title),"",paste0("**Predict before running:** ",question))
    if (pause && interactive()) readline("Make your prediction, then press Enter: ")
  }
  code <- function(...) add("```r",...,"```")
  table <- function(name,x) {
    utils::write.csv(x,file.path(output_dir,paste0(name,".csv")),row.names=FALSE)
    add(.ls_markdown_table(x))
  }
  plot_file <- function(name,draw) {
    render <- function() {
      graphics::par(mar=c(5,4.7,4,1.4),mgp=c(2.8,.7,0),tcl=-.25,
        col=ink,col.axis=ink,col.lab=ink,col.main=ink,family="sans",las=1)
      draw()
    }
    grDevices::png(file.path(figures,paste0(name,".png")),width=1600,height=1000,res=150)
    device <- grDevices::dev.cur()
    tryCatch(render(),finally=grDevices::dev.off(device))
    if (show_plots) {
      old_par <- graphics::par(no.readonly=TRUE)
      tryCatch(render(),finally=graphics::par(old_par))
    }
    add(sprintf("![%s](figures/%s.png)",gsub("-"," ",name),name))
  }
  utils::write.csv(cbind(d,known_signal=example$signal),file.path(output_dir,"synthetic-data.csv"),row.names=FALSE)
  step("1. Repeated measurements", "Which differences persist across items, and which fluctuate?")
  add("Imagine repeated parallel measurements of wellbeing. Scores are continuous arbitrary units, not a validated questionnaire.",
    "The known stable signal is -0.9*workload + 0.7*resources + 0.5*an unobserved profile. The unused predictor contributes no stable signal.",
    "The profiles are orthogonal with sample variance 1: this controlled construction makes the true contributions inspectable. Each item adds independent normal error.")
  code("d <- lesson$example$data", "indicators <- lesson$example$indicators", "head(d)")
  plot_file("01-measurements",function() {
    graphics::par(mfrow=c(1,2))
    selected <- order(example$signal)[round(seq(1,n,length.out=12))]
    Y <- as.matrix(d[selected,indicators]);xx <- seq_along(selected)
    graphics::plot(xx,example$signal[selected],type="n",ylim=range(Y,example$signal[selected]),xaxt="n",
      xlab="People ordered by known signal (IDs)",ylab="Score units",main="Repeated measurements")
    graphics::mtext(sprintf("12 of %d people; %d parallel items",n,items),3,.5,cex=.8)
    graphics::axis(1,at=xx,labels=d$id[selected],cex.axis=.7)
    for(j in seq_len(items)) graphics::points(xx+(j-(items+1)/2)*.5/items,Y[,j],pch=1,col=grey,cex=.7)
    graphics::points(xx,example$signal[selected],pch=18,col=blue,cex=1.2)
    graphics::points(xx,d$score[selected],pch=17,col=gold,cex=1)
    graphics::legend("topleft",c("Individual items","Known signal","Item mean"),pch=c(1,18,17),
      col=c(grey,blue,gold),bty="n",cex=.7)
    lim <- range(example$signal,d$score)
    graphics::plot(example$signal,d$score,xlim=lim,ylim=lim,pch=16,col=grDevices::adjustcolor(blue,.6),
      xlab="Known stable signal",ylab="Observed item mean",main="Mean score and stable signal")
    graphics::mtext(sprintf("All %d people; dashed line: exact agreement",n),3,.5,cex=.8)
    graphics::abline(0,1,lty=2,col=ink)
  })
  add("The mean is less noisy than an individual item, but it still contains measurement error. The latent signal is available here only because we generated the data.")

  step("2. Estimate the background and total structure", "Will all between-person variation in the item mean count as stable structure?")
  measurement <- lambdastar::lambda_measure(d,indicators)
  oracle <- .ls_oracle(example,c("workload","resources"))
  observed <- sum((d$score-mean(d$score))^2)
  budget <- measurement$d*measurement$temperature
  energy <- data.frame(component=c("Observed energy","Noise budget d*T","Raw signal energy"),
    value=c(observed,budget,observed-budget))
  code("measurement <- lambdastar::lambda_measure(d, indicators)","measurement$estimates")
  table("measurement",data.frame(parameter=c("T (mean-score variance)","lambda_star"),
    known=unname(oracle[1:2]),estimated=unname(measurement$estimates)))
  table("measurement-energy",energy)
  add("On the mean scale, T estimates item-error variance divided by K. Here its generating value is noise_sd^2/K.",
    "With d=N-1, raw signal energy = observed centered energy - d*T, and raw lambda_star = raw signal energy/T. The reported signal estimate is constrained to be nonnegative.",
    "Lambda_star is total stable energy relative to background fluctuation. It is not a percentage or a proportion of explained variance.")
  plot_file("02-measurement",function() {
    graphics::barplot(energy$value,names.arg=c("Observed energy","Noise budget: d*T","Raw signal estimate"),
      col=c(grey,gold,blue),border=ink,ylim=range(0,energy$value)*1.15,
      ylab="Squared centered score units",main="Observed energy and measurement correction")
    graphics::mtext(sprintf("N=%d, d=%d; subtraction uses estimated mean-scale T",n,measurement$d),3,.5,cex=.85)
    graphics::abline(h=0,col=ink)
  })

  step("3. Localize structure in a hypothesis", "How much should workload capture? What remains after adding resources?")
  formulas <- list(H0=~1,workload=~workload,full=~workload+resources,extra=~workload+resources+unused)
  variables <- list(character(),"workload",c("workload","resources"),c("workload","resources","unused"))
  fits <- lapply(formulas,function(f) lambdastar::lambda_model(d,f,indicators,T_q))
  truths <- lapply(variables,function(v) .ls_oracle(example,v))
  localization <- do.call(rbind,lapply(seq_along(fits),function(i) data.frame(model=names(fits)[i],
    p=fits[[i]]$df_h,lambda_star=fits[[i]]$lambda_star,lambda_h=fits[[i]]$lambda_h,
    lambda_perp=fits[[i]]$lambda_perp,kappa=fits[[i]]$kappa,true_kappa=unname(truths[[i]]["kappa"]))))
  rownames(localization) <- NULL
  code("T_q <- lesson$settings$T_q",
    "full <- lambdastar::lambda_model(d, ~workload + resources, indicators, T_q)","full$estimates")
  table("localization",localization)
  add("X_H contains the actual encoded predictor directions. Projection separates observed energy into captured and residual parts.",
    "The paper estimator subtracts p*T from captured energy and (d-p)*T from residual energy. With positive estimated total signal, kappa_H = lambda_H/lambda_star.",
    sprintf("Known capture is %.1f%% for workload and %.1f%% for workload plus resources. The remaining hidden profile is stable structure, not item measurement noise.",
      100*example$truth["kappa_workload"],100*example$truth["kappa_full"]),
    "H0 contains no hypothesis direction after the intercept is removed. An unused predictor can capture some measurement fluctuation in a particular dataset.",
    "Estimated lambda_H can decrease after adding a direction: the extra noise correction can exceed its gain in observed captured energy.")
  plot_file("03-localization",function() {
    graphics::par(mar=c(5,6,4,1.4))
    values <- rbind(localization$lambda_h,localization$lambda_perp)
    graphics::barplot(values,horiz=TRUE,names.arg=localization$model,col=c(blue,gold),border=ink,
      xlab="Estimated noncentrality (common mean-scale T)",main="Captured and remaining structure",
      xlim=c(0,max(1,colSums(values))*1.12),ylim=c(0,nrow(localization)*1.2+.8))
    graphics::mtext(sprintf("Same %d people and %d items for all four models",n,items),3,.5,cex=.85)
    graphics::legend("topright",c("Captured lambda_H","Remaining lambda_perp"),fill=c(blue,gold),bty="n",cex=.8)
  })

  step("4. Inspect the ordinary percentile bootstrap", "Which rows travel together when a person is sampled twice?")
  bootstrap <- lambdastar::lambda_model(d,formulas$full,indicators,T_q,bootstrap=cfg)
  ci <- bootstrap$bootstrap$ci
  code("ix <- lesson$bootstrap$bootstrap$indices[1, ]", "d[ix[1:10], c('id', 'workload', 'resources', indicators)]",
    "draws <- lesson$bootstrap$bootstrap$draws", "confint(lesson$bootstrap, parm='kappa')",
    "# When all kappa draws are defined (check the diagnostics first):",
    "alpha <- (1-lesson$settings$conf.level)/2",
    "if (all(is.finite(draws[, 'kappa']))) quantile(draws[, 'kappa'], c(alpha, 1-alpha), type=7)")
  table("intervals",ci[c("parameter","estimate","lower","upper","n_valid","n_failed","status")])
  utils::write.csv(bootstrap$bootstrap$draws,file.path(output_dir,"bootstrap-draws.csv"),row.names=FALSE)
  add(sprintf("All %d draws sample N complete rows with replacement: items and predictors stay together. The selected paper estimator is recalculated each time.",B),
    sprintf("The central %.1f%% is bounded by the %.2f%% and %.2f%% percentiles. This lesson explicitly chooses type-7 interpolation; the package also supports the usual boot interpolation.",
      100*conf.level,50*(1-conf.level),50*(1+conf.level)),
    "No endpoints are widened or replaced. Undefined draws are counted, and an affected regular CI is unavailable; zeros in a defined statistic remain real observations.",
    "A percentile CI describes uncertainty under the chosen row-resampling procedure. One interval containing the known signal is not a coverage validation.")
  plot_file("04-bootstrap",function() {
    graphics::par(mfrow=c(1,2))
    for (key in c("temperature","kappa")) {
      values <- bootstrap$bootstrap$draws[,key];valid <- is.finite(values)
      row <- ci[ci$parameter==key,]
      if (!any(valid)) {
        graphics::plot.new();graphics::title(main=key);graphics::text(.5,.5,"No finite draws; no regular CI")
      } else {
        h <- graphics::hist(values[valid],breaks="FD",plot=FALSE)
        limits <- c(row$lower,row$upper)
        graphics::plot(h,col=blue,border="white",ylim=c(0,max(h$counts)*1.3),
          xlim=range(h$breaks,unname(oracle[key]),limits,na.rm=TRUE),
          main=if(key=="kappa") "Bootstrap capture fraction" else "Bootstrap temperature",
          xlab=if(key=="kappa") "kappa_H" else "T: mean-score variance",ylab="Number of draws")
        graphics::mtext(sprintf("Finite draws: %d/%d; CI status: %s",sum(valid),B,row$status),3,.5,cex=.8)
        graphics::segments(unname(oracle[key]),0,unname(oracle[key]),max(h$counts)*1.02,lty=2,lwd=2,col=ink)
        has_ci <- all(is.finite(limits))
        if (has_ci) graphics::segments(limits,0,limits,max(h$counts)*1.02,col=gold,lwd=2)
        graphics::legend("topright",c("Known generating value",if(has_ci) "Actual CI endpoints"),
          lty=c(2,if(has_ci)1),col=c(ink,if(has_ci)gold),lwd=2,bty="n",cex=.7)
      }
    }
  })

  step("5. Compare explicitly coded models", "Does reducing residual error necessarily reduce F_H?")
  comparison <- list(useful=lambdastar::lambda_compare(d,formulas$workload,formulas$full,indicators,T_q),
    unused=lambdastar::lambda_compare(d,formulas$full,formulas$extra,indicators,T_q))
  parsimony <- do.call(rbind,lapply(names(fits),function(name) {
    p <- fits[[name]]$parsimony
    data.frame(model=name,RSS=p$rss,dimension=p$dimension_cost,geometry=p$geometry_cost,F_H=p$free_energy)
  }))
  increments <- do.call(rbind,lapply(names(comparison),function(name) {
    c <- comparison[[name]];a <- c$model_a$parsimony;b <- c$model_b$parsimony
    data.frame(addition=name,delta_RSS=b$rss-a$rss,delta_dimension=b$dimension_cost-a$dimension_cost,
      delta_geometry=b$geometry_cost-a$geometry_cost,delta_F=c$delta_F)
  }))
  code("lambdastar::lambda_compare(d, ~workload, ~workload + resources, indicators, T_q)")
  table("parsimony",parsimony);table("model-increments",increments)
  add("F_H = RSS + p*Xi_eta(T,T_q) - Xi_a(T,T_q)*log(V_H). The package evaluates Xi_eta numerically; V_H comes from the declared coding of X_H.",
    sprintf("The same mean-scale resolution variance T_q=%.6g is used for both models in every comparison. It is an explicit illustrative assumption, not estimated from these continuous data.",T_q),
    sprintf("Adding resources changes F_H by %.6g; adding the unused predictor changes it by %.6g. A negative change favors the larger model by this functional; a positive change favors the smaller one.",
      comparison$useful$delta_F,comparison$unused$delta_F),
    "Extra directions can reduce RSS while adding dimension or geometry costs. An irrelevant predictor can occasionally reduce F_H through sampling fluctuation; a single example does not guarantee model recovery.")
  plot_file("05-parsimony",function() {
    graphics::par(mfrow=c(1,2))
    graphics::barplot(parsimony$F_H,names.arg=parsimony$model,col=blue,border=ink,
      ylab="F_H (squared score units)",main="Parsimony functional by model",ylim=c(0,max(parsimony$F_H)*1.1))
    graphics::mtext(sprintf("Common data and T_q=%.4g",T_q),3,.5,cex=.8)
    v <- unlist(increments[increments$addition=="unused",-1],use.names=FALSE)
    graphics::barplot(v,names.arg=c("RSS","Dimension","Geometry","Total"),col=c(grey,gold,gold,blue),border=ink,
      ylim=range(0,v)*1.3,ylab="Change in squared score units",main="Adding the unused predictor")
    graphics::mtext("Component changes sum to delta F; separate vertical scale",3,.5,cex=.75)
    graphics::abline(h=0,col=ink)
  })

  step("6. Change the measurement conditions", "What happens to lambda_star and kappa when noise, K, or N changes?")
  grids <- list(noise_sd=sort(unique(c(.4,.7,1,1.3,1.6,2,noise_sd))),
    items=sort(unique(c(2,3,4,6,8,12,items))),n=sort(unique(c(40,80,120,160,240,320,n))))
  experiments <- do.call(rbind,lapply(names(grids),function(knob) do.call(rbind,lapply(grids[[knob]],function(value) {
    settings <- list(n=n,items=items,noise_sd=noise_sd,seed=seed);settings[[knob]] <- value
    e <- do.call(lambdastar_tutorial_data,settings)
    f <- lambdastar::lambda_model(e$data,formulas$full,e$indicators,T_q)
    t <- .ls_oracle(e,c("workload","resources"))
    data.frame(experiment=knob,value=value,n=settings$n,items=settings$items,noise_sd=settings$noise_sd,
      T=f$temperature,true_T=t["temperature"],lambda_star=f$lambda_star,true_lambda_star=t["lambda_star"],
      kappa=f$kappa,true_kappa=t["kappa"])
  }))));rownames(experiments) <- NULL
  table("experiments",experiments)
  add("More item noise raises T and reduces the true noncentrality. More independent parallel items lower the noise of their mean and raise noncentrality.",
    "Here Var(S)=1.55 is held exactly fixed by construction. Thus true lambda_star=(N-1)*1.55/T grows with N. The true capture fraction stays 1.30/1.55.",
    "Each point is one dataset, not a simulation average or a confidence band. At fixed N, the seed retains the same profiles and reuses item errors when noise or K changes. Changing N generates a new standardized design.")
  plot_file("06-experiments",function() {
    graphics::par(mfrow=c(2,3),mar=c(4.5,4.5,3.4,1),cex=.9)
    for (metric in c("lambda_star","kappa")) for (knob in names(grids)) {
      x <- experiments[experiments$experiment==knob,]
      truth <- x[[paste0("true_",metric)]];estimate <- x[[metric]]
      ylim <- if(metric=="kappa")c(0,1) else range(0,truth,estimate,na.rm=TRUE)*1.08
      graphics::plot(x$value,truth,type="b",pch=5,lty=2,col=ink,ylim=ylim,
        xlab=switch(knob,noise_sd="Item-error SD",items="Number of items K",n="Number of people N"),
        ylab=metric,main=paste(metric,"and",switch(knob,noise_sd="noise",items="K",n="N")))
      graphics::lines(x$value,estimate,type="b",pch=16,col=blue)
      graphics::mtext("One dataset per setting; known vs estimated",3,.35,cex=.68)
      if(metric=="lambda_star" && knob=="noise_sd") graphics::legend("topright",c("Known","Estimated"),
        col=c(ink,blue),pch=c(5,16),lty=c(2,1),bty="n",cex=.75)
    }
  })

  step("7. Change the coding, preserve the space", "Will an invertible recoding change captured structure, geometry, or both?")
  recoded <- d;recoded$resources_tilted <- recoded$resources+2*recoded$workload
  coding <- lambdastar::lambda_compare(recoded,~workload+resources,~workload+resources_tilted,indicators,T_q)
  coding_table <- do.call(rbind,lapply(c("model_a","model_b"),function(name) {
    f <- coding[[name]];data.frame(coding=name,lambda_h=f$lambda_h,kappa=f$kappa,
      RSS=f$parsimony$rss,volume=f$parsimony$volume,geometry=f$parsimony$geometry_cost,F_H=f$parsimony$free_energy)
  }))
  code("d$resources_tilted <- d$resources + 2*d$workload",
    "lambdastar::lambda_compare(d, ~workload+resources, ~workload+resources_tilted, indicators, T_q)")
  table("coding",coding_table)
  add("Both codings span the same hypothesis space: their projections, lambda_H and kappa agree. The normalized column volume changes, so the geometry term and F_H can change.",
    "Choose a substantively meaningful coding and report it. A convenient QR basis for projection does not replace the declared model coding in F_H.")
  plot_file("07-coding",function() {
    graphics::par(mfrow=c(1,2))
    normalized <- scale(as.matrix(recoded[c("workload","resources","resources_tilted")]))
    limits <- range(normalized)
    for (j in 1:2) {
      z <- normalized[,c(1,j+1)]
      graphics::plot(z[,1],z[,2],pch=16,col=grDevices::adjustcolor(blue,.65),
        xlim=limits,ylim=limits,asp=1,
        xlab="Workload (standardized)",ylab=if(j==1)"Resources (standardized)" else "Recoded resources (standardized)",
        main=if(j==1)"Original coding" else "Invertible recoding")
      graphics::mtext(sprintf("N=%d; V_H=%.3f; same hypothesis space",n,coding_table$volume[j]),3,.5,cex=.8)
    }
  })
  add("## Your next experiment", "Change one argument and choose a new output folder. Predict the direction first, then inspect the new tables and plots.",
    "For example: lambdastar_tutorial(output_dir='more-noise', noise_sd=2). Try items=8, n=320, or T_q=1/48 next.",
    "Changing T_q leaves lambda/T/kappa unchanged for the same data, but can change the costs in F_H. Every model comparison still uses one common data matrix and resolution.",
    "The separate ADVANCED_EXAMPLES.md introduces population estimation, common nuisance adjustment and cluster resampling.")
  result <- list(example=example,measurement=measurement,oracle=oracle,fits=fits,
    bootstrap=bootstrap,comparison=comparison,coding=coding,
    tables=list(energy=energy,localization=localization,parsimony=parsimony,increments=increments,
      experiments=experiments,coding=coding_table),output_dir=output_dir,
    settings=list(B=B,conf.level=conf.level,T_q=T_q,version=as.character(utils::packageVersion("lambdastar"))))
  saveRDS(result,file.path(output_dir,"lesson.rds"))
  writeLines(report,file.path(output_dir,"WALKTHROUGH.md"))
  utils::capture.output(utils::sessionInfo(),file=file.path(output_dir,"session-info.txt"))
  cat("Saved lesson:",file.path(output_dir,"WALKTHROUGH.md"),"\n")
  invisible(result)
}
