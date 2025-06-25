library(mlbench)
library(splines)
library(caret)


kl_divergence <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  p <- ifelse(p == 0, 1e-10, p) 
  q <- ifelse(q == 0, 1e-10, q)
  sum(p * log(p / q))
}

# Sigmoidal
sigmoid_kernel <- function(x, kappa = 1, theta = 0) {
  x <- as.matrix(x)
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("Input 'x' for sigmoid_kernel must be a finite numeric matrix without NAs or non-finite values.")
  }
  return(tanh(kappa * (x %*% t(x)) + theta))
}

evaluate_kappa_kl <- function(data, labels, kappa_values, theta = 0) {
  results <- data.frame(kappa = kappa_values, KL = rep(NA, length(kappa_values)))
  
  for (k_idx in seq_along(kappa_values)) {
    kappa <- kappa_values[k_idx]
    
    K <- tryCatch({
      sigmoid_kernel(data, kappa = kappa, theta = theta)
    }, error = function(e) {
      message(paste("Error calculating Sigmoid kernel for kappa =", kappa, ":", e$message))
      return(NULL)
    })
    
    if (is.null(K) || anyNA(K) || any(!is.finite(K)) || all(K == 0) || all(K == 1)) {
      results$KL[k_idx] <- NA
      next
    }
    
    classes <- unique(labels)
    class_dists <- list()
    
    for (cl in classes) {
      idx <- which(labels == cl)
      if (length(idx) < 2) {
        warning(paste("Class '", cl, "' has less than 2 samples. Skipping KL calculation for this class.", sep=""))
        next
      }
      class_dist <- rowMeans(K[idx, ])
      
      sum_cd <- sum(class_dist, na.rm = TRUE)
      if (is.na(sum_cd) || sum_cd == 0 || !is.finite(sum_cd)) {
        class_dists[[as.character(cl)]] <- rep(1 / length(class_dist), length(class_dist))
      } else {
        class_dists[[as.character(cl)]] <- class_dist / sum_cd
      }
      if (anyNA(class_dists[[as.character(cl)]]) || any(!is.finite(class_dists[[as.character(cl)]])) || sum(class_dists[[as.character(cl)]], na.rm=TRUE) == 0) {
        class_dists[[as.character(cl)]] <- rep(1 / length(class_dist), length(class_dist))
      }
    }
    
    kl_values <- c()
    cl_names <- names(class_dists)
    
    if (length(cl_names) >= 2) {
      for (i in 1:(length(cl_names) - 1)) {
        for (j in (i + 1):length(cl_names)) {
          p <- class_dists[[cl_names[i]]]
          q <- class_dists[[cl_names[j]]]
          
          if (length(p) == 0 || length(q) == 0 || anyNA(p) || anyNA(q) || any(!is.finite(p)) || any(!is.finite(q))) {
            warning(paste("Skipping KL for invalid distributions between", cl_names[i], "and", cl_names[j]))
            next
          }
          
          kl_val <- kl_divergence(p, q)
          if (!is.na(kl_val) && is.finite(kl_val)) {
            kl_values <- c(kl_values, kl_val)
          }
        }
      }
      if (length(kl_values) > 0) {
        results$KL[k_idx] <- mean(kl_values)
      } else {
        results$KL[k_idx] <- NA
      }
    } else {
      results$KL[k_idx] <- NA
    }
  }
  return(results)
}

selecionar_kappa_ideal_spline_max <- function(kl_results) {
  valid_points <- !is.na(kl_results$KL) & !is.infinite(kl_results$KL)
  
  if (sum(valid_points) < 4) {
    warning("Not enough valid points to fit a spline. Returning kappa for max observed KL or NA.")
    if (any(valid_points)) {
      return(kl_results$kappa[which.max(kl_results$KL[valid_points])])
    }
    return(NA)
  }
  
  kappa_valid <- kl_results$kappa[valid_points]
  kl_valid <- kl_results$KL[valid_points]
  
  if (length(unique(kappa_valid)) < 2) { 
    warning("Not enough unique kappa values for spline. Returning max observed or NA.")
    if (any(valid_points)) {
      return(kl_results$kappa[which.max(kl_results$KL[valid_points])])
    }
    return(NA)
  }
  
  spline_fit <- tryCatch({
    smooth.spline(kappa_valid, kl_valid, spar = 0.6)
  }, error = function(e) {
    warning(paste("Error fitting spline:", e$message, ". Returning max observed KL or NA."))
    if (any(valid_points)) {
      return(kl_results$kappa[which.max(kl_results$KL[valid_points])])
    }
    return(NA)
  })
  
  if (is.null(spline_fit)) {
    if (any(valid_points)) {
      return(kl_results$kappa[which.max(kl_results$KL[valid_points])])
    }
    return(NA)
  }
  
  finer_kappa_grid <- seq(min(kappa_valid), max(kappa_valid), length.out = 500)
  predicted_kl <- predict(spline_fit, x = finer_kappa_grid)$y
  
  if (all(is.na(predicted_kl)) || all(!is.finite(predicted_kl))) {
    warning("Spline prediction resulted in all NA/non-finite values. Returning max observed KL or NA.")
    if (any(valid_points)) {
      return(kl_results$kappa[which.max(kl_results$KL[valid_points])])
    }
    return(NA)
  }
  
  max_kl_idx <- which.max(predicted_kl)
  if (length(max_kl_idx) == 0) {
    warning("Could not find a maximum on the spline. Returning max observed KL or NA.")
    if (any(valid_points)) {
      return(kl_results$kappa[which.max(kl_results$KL[valid_points])])
    }
    return(NA)
  }
  
  ideal_kappa <- finer_kappa_grid[max_kl_idx]
  
  if (ideal_kappa == min(finer_kappa_grid) || ideal_kappa == max(finer_kappa_grid)) {
    warning("Ideal kappa found at the edge of the search range. Consider expanding 'kappa_seq'.")
  }
  
  return(ideal_kappa)
}


avaliar_dataset_sigmoid <- function(nome_dataset, kappa_range_default = c(0.0001, 0.1), num_kappas_default = 100, theta_param = 1) { # Range padrão diminuído
  cat("\n✨ Dataset:", nome_dataset, "\n")
  data(list = nome_dataset, package = "mlbench")
  dados <- get(nome_dataset)
  
  dados <- na.omit(as.data.frame(dados))
  
  y <- dados[[ncol(dados)]]
  X <- dados[, -ncol(dados)]
  
  X_clean <- as.data.frame(sapply(X, function(col) {
    if (is.factor(col)) {
      as.numeric(as.factor(col))
    } else if (is.character(col)) { 
      as.numeric(as.factor(col))
    } else {
      as.numeric(col) 
    }
  }))
  
  X_clean[is.na(X_clean)] <- 0
  X_clean[!is.finite(as.matrix(X_clean))] <- 0
  
  cols_to_scale <- sapply(X_clean, function(col) var(col, na.rm = TRUE) > 1e-9)
  
  if (any(cols_to_scale)) {
    X_scaled <- as.data.frame(scale(X_clean[, cols_to_scale, drop = FALSE]))
    names(X_scaled) <- names(X_clean)[cols_to_scale]
    
    if (sum(!cols_to_scale) > 0) {
      X_final <- cbind(X_scaled, X_clean[, !cols_to_scale, drop = FALSE])
      X_final <- X_final[, names(X_clean)]
    } else {
      X_final <- X_scaled
    }
  } else {
    warning("No variable columns to scale in dataset ", nome_dataset, ". Using unscaled (but numeric) data.")
    X_final <- X_clean
  }
  
  set.seed(123)
  if (!is.factor(y)) {
    y <- as.factor(y)
  }
  
  if (any(table(y) < 2)) {
    warning(paste("Dataset '", nome_dataset, "' has classes with fewer than 2 samples. createDataPartition might fail or result in empty training classes.", sep=""))
  }
  
  idx <- createDataPartition(y, p = 0.7, list = FALSE)
  
  if (length(idx) == 0 || nrow(X_final[idx, , drop = FALSE]) == 0) {
    warning(paste("Conjunto de treino vazio para o dataset '", nome_dataset, "' após o particionamento. Pulando avaliação.", sep=""))
    return(list(kappa_ideal = NA, kl_results = data.frame(kappa=numeric(0), KL=numeric(0)), dataset_name = nome_dataset))
  }
  
  X_train <- X_final[idx, , drop = FALSE]
  y_train <- y[idx]
  
  y_train <- droplevels(y_train)
  if (length(levels(y_train)) < 2) {
    warning(paste("y_train para '", nome_dataset, "' tem menos de 2 classes distintas após particionamento. Impossível calcular KL.", sep=""))
    return(list(kappa_ideal = NA, kl_results = data.frame(kappa=numeric(0), KL=numeric(0)), dataset_name = nome_dataset))
  }
  

  current_kappa_range <- kappa_range_default
  num_kappas <- num_kappas_default
  current_theta <- theta_param 
  
  if (nome_dataset == "PimaIndiansDiabetes") {
    current_kappa_range <- c(0.0001, 0.1) # Ajuste focado
    num_kappas <- 200 # Mais pontos
    current_theta <- 0.05 # Theta mais sutil
  } else if (nome_dataset == "BreastCancer") {
    current_kappa_range <- c(0.0001, 0.05) # Ajuste focado
    num_kappas <- 200 # Mais pontos
    current_theta <- 0.01 # Theta bem pequeno
  } else if (nome_dataset == "Vehicle") {
    current_kappa_range <- c(0.0001, 0.1) # Ajuste focado
    num_kappas <- 200 # Mais pontos
    current_theta <- 0.05
  } else if (nome_dataset == "Ozone") { 
    current_kappa_range <- c(0.0001, 0.01) # Ajuste focado para Ozone
    num_kappas <- 200 # Mais pontos
    current_theta <- 0.01 
    if (!is.factor(y_train)) {
      breaks_ozone <- unique(quantile(y_train, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE))
      if (length(breaks_ozone) < 2) {
        breaks_ozone <- c(min(y_train, na.rm=TRUE) - 1e-6, max(y_train, na.rm=TRUE) + 1e-6)
      }
      y_train_binned <- cut(y_train, breaks = breaks_ozone, include.lowest = TRUE, ordered_result = TRUE)
      valid_levels <- names(table(y_train_binned)[table(y_train_binned) > 0])
      y_train_binned <- factor(y_train_binned, levels = valid_levels)
      
      if (length(levels(y_train_binned)) < 2) {
        warning("Ozone target binning resulted in fewer than 2 distinct classes. Cannot compute KL divergence.")
        return(list(kappa_ideal = NA, kl_results = data.frame(kappa=numeric(0), KL=numeric(0)), dataset_name = nome_dataset))
      }
      y_train <- y_train_binned
    }
  } else if (nome_dataset == "Ionosphere") {
    current_kappa_range <- c(0.0001, 0.05) # Ajuste focado
    num_kappas <- 200 # Mais pontos
    current_theta <- 0.01
  }
  
  
  kappa_seq <- seq(current_kappa_range[1], current_kappa_range[2], length.out = num_kappas)
  
  kl_results <- evaluate_kappa_kl(X_train, y_train, kappa_seq, theta = current_theta) 
  
  if (sum(!is.na(kl_results$KL) & is.finite(kl_results$KL)) < 1) {
    warning(paste("No valid finite KL results for dataset", nome_dataset, ". Cannot find ideal kappa reliably."))
    kappa_ideal <- NA
  } else {
    kappa_ideal <- selecionar_kappa_ideal_spline_max(kl_results)
  }
  
  cat("Kappa ideal:", round(kappa_ideal, 4), "\n")
  
  return(list(kappa_ideal = kappa_ideal, kl_results = kl_results, dataset_name = nome_dataset))
}


gera_graficos_sigmoid <- function(kappa_ideal, kl_results, dataset_name) {
  valid_points <- !is.na(kl_results$KL) & !is.infinite(kl_results$KL)
  
  main_title <- paste("KL vs Kappa (Sigmoidal) - ", dataset_name)
  
  if (sum(valid_points) == 0) { 
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1), 
         main=paste(main_title, "\n(Sem Dados Válidos para Plotar)"))
    text(0.5, 0.5, "Nenhum dado KL válido para plotar.", cex=1.2)
    warning(paste("Nenhum dado KL válido para plotar para ", dataset_name, "."))
    return(invisible(NULL))
  }
  
  kl_values_for_ylim <- kl_results$KL[valid_points]
  ylim_min <- min(kl_values_for_ylim) * 0.95
  ylim_max <- max(kl_values_for_ylim) * 1.05
  if (ylim_min == ylim_max) {
    ylim_min <- ylim_min - 0.1
    ylim_max <- ylim_max + 0.1
  }
  
  if (sum(valid_points) < 4) {
    plot(kl_results$kappa, kl_results$KL, type = "p", main = paste(main_title, " - Sem Ajuste Spline"),
         xlab = "Kappa", ylab = "Divergência KL", col = "black", pch = 16, ylim = c(ylim_min, ylim_max))
    if (!is.na(kappa_ideal) && is.finite(kappa_ideal)) {
      abline(v = kappa_ideal, col = "red", lty = 2, lwd = 1.5)
      text(kappa_ideal, max(kl_results$KL, na.rm = TRUE),
           paste("Kappa ideal:", round(kappa_ideal, 3)), pos = 4, col = "black")
    }
    warning(paste("Pontos válidos insuficientes para gerar spline para ", dataset_name, ". Mostrando apenas pontos brutos."))
    return(invisible(NULL))
  }
  
  spline_fit <- smooth.spline(kl_results$kappa[valid_points], kl_results$KL[valid_points], spar = 0.6)
  plot(kl_results$kappa, kl_results$KL, type = "p", main = main_title,
       xlab = "Kappa", ylab = "Divergência KL", col = "black", pch = 16, ylim = c(ylim_min, ylim_max))
  lines(spline_fit, col = "blue", lwd = 1.5)
  
  if (!is.na(kappa_ideal) && is.finite(kappa_ideal)) {
    abline(v = kappa_ideal, col = "red", lty = 2, lwd = 1.5)
    predicted_kl_at_ideal <- predict(spline_fit, x = kappa_ideal)$y
    text_y_pos <- if (!is.na(predicted_kl_at_ideal) && is.finite(predicted_kl_at_ideal)) predicted_kl_at_ideal else max(kl_values_for_ylim)
    
    if (text_y_pos < ylim_min * 1.05) text_y_pos <- ylim_min * 1.05
    if (text_y_pos > ylim_max * 0.95) text_y_pos <- ylim_max * 0.95
    
    text(kappa_ideal, text_y_pos, round(kappa_ideal, 3), pos = 4, col = "black")
  }
  
  legend("bottomright", legend = c("Ajuste Spline", "Kappa Ideal"),
         col = c("blue", "red"), lty = c(1, 2), bty = "n")
}


datasets_to_evaluate <- c("PimaIndiansDiabetes", "BreastCancer", "Vehicle", "Ozone", "Ionosphere")
results_list <- list()


for (ds_name in datasets_to_evaluate) {

  dataset_result <- avaliar_dataset_sigmoid(ds_name, kappa_range_default = c(0.0001, 0.1), num_kappas_default = 100, theta_param = 0.01) 
  results_list[[ds_name]] <- dataset_result
  gera_graficos_sigmoid(dataset_result$kappa_ideal, dataset_result$kl_results, dataset_result$dataset_name)

}


print("Valores ideais de Kappa para Kernel Sigmoidal (encontrados por máximo da spline):")
print(sapply(results_list, `[[`, "kappa_ideal"))