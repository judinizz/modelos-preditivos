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

rbf_kernel <- function(x, sigma = 1) {
  x <- as.matrix(x)
  n <- nrow(x)
  K <- matrix(0, n, n)
  dist_sq <- as.matrix(dist(x))^2
  K <- exp(-sigma * dist_sq)
  return(K)
}

# RBF
evaluate_sigma_kl <- function(data, labels, sigma_values) {
  results <- data.frame(sigma = sigma_values, KL = NA)
  
  numeric_cols <- sapply(data, is.numeric)
  if (any(numeric_cols)) {
    data[, numeric_cols] <- scale(data[, numeric_cols])
  } else {
    warning("No numeric columns found for scaling in evaluate_sigma_kl.")
  }

  for (s_idx in seq_along(sigma_values)) {
    sigma <- sigma_values[s_idx]
    
    K <- tryCatch({
      rbf_kernel(data, sigma = sigma)
    }, error = function(e) {
      message(paste("Error calculating RBF kernel for sigma =", sigma, ":", e$message))
      return(NULL)
    })
    
    if (is.null(K) || any(is.na(K)) || all(K == 0)) {
      results$KL[s_idx] <- NA 
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
      class_dists[[as.character(cl)]] <- class_dist / sum(class_dist)
    }
    
    kl_values <- c()
    cl_names <- names(class_dists)
    

    if (length(cl_names) >= 2) {
      for (i in 1:(length(cl_names) - 1)) {
        for (j in (i + 1):length(cl_names)) {
          p <- class_dists[[cl_names[i]]]
          q <- class_dists[[cl_names[j]]]
          p <- p + 1e-10
          q <- q + 1e-10
          kl_values <- c(kl_values, kl_divergence(p, q))
        }
      }
      if (length(kl_values) > 0) {
        results$KL[s_idx] <- mean(kl_values)
      } else {
        results$KL[s_idx] <- NA 
      }
    } else {
      results$KL[s_idx] <- NA
    }
  }
  return(results)
}


selecionar_sigma_ideal_spline_max <- function(kl_results) {
  valid_points <- !is.na(kl_results$KL) & !is.infinite(kl_results$KL)
  
  if (sum(valid_points) < 4) {
    warning("Not enough valid points to fit a spline. Returning sigma for max observed KL.")
    return(kl_results$sigma[which.max(kl_results$KL)])
  }
  

  spline_fit <- tryCatch({
    smooth.spline(kl_results$sigma[valid_points], kl_results$KL[valid_points], spar = 0.6)
  }, error = function(e) {
    warning(paste("Error fitting spline:", e$message, ". Returning sigma for max observed KL."))
    return(NULL)
  })

  if (is.null(spline_fit)) {
    return(kl_results$sigma[which.max(kl_results$KL)])
  }
  
  finer_sigma_grid <- seq(min(kl_results$sigma[valid_points]), max(kl_results$sigma[valid_points]), length.out = 500)
  predicted_kl <- predict(spline_fit, x = finer_sigma_grid)$y
  

  max_kl_idx <- which.max(predicted_kl)
  if (length(max_kl_idx) == 0) {
    warning("Could not find a maximum on the spline. Returning sigma for max observed KL.")
    return(kl_results$sigma[which.max(kl_results$KL)])
  }
  
  ideal_sigma <- finer_sigma_grid[max_kl_idx]
  

  if (ideal_sigma == min(finer_sigma_grid) || ideal_sigma == max(finer_sigma_grid)) {
      warning("Ideal sigma found at the edge of the search range. Consider expanding 'sigma_seq'.")
  }
  
  return(ideal_sigma)
}

avaliar_dataset_rbf <- function(nome_dataset, sigma_range = c(0.01, 5), num_sigmas = 30) {
  cat("\n✨ Dataset:", nome_dataset, "\n")
  data(list = nome_dataset, package = "mlbench")
  dados <- get(nome_dataset)
  dados <- na.omit(as.data.frame(dados)) 
  
  y <- dados[[ncol(dados)]]
  X <- dados[, -ncol(dados)]
  
  X <- data.frame(lapply(X, function(col) {
    if (is.factor(col)) {
      as.numeric(as.factor(col))
    } else if (is.character(col)) {
      as.numeric(as.factor(col)) 
    } else {
      col 
    }
  }))
  

  set.seed(123)

  if (!is.factor(y)) {
    y <- as.factor(y)
  }
  idx <- createDataPartition(y, p = 0.7, list = FALSE)
  X_train <- X[idx, ]
  y_train <- y[idx]
  
  sigma_seq <- seq(sigma_range[1], sigma_range[2], length.out = num_sigmas)
  
  kl_results <- evaluate_sigma_kl(X_train, y_train, sigma_seq)
  
  sigma_ideal <- selecionar_sigma_ideal_spline_max(kl_results)
  cat("Sigma ideal:", round(sigma_ideal, 4), "\n")
  
  return(list(sigma_ideal = sigma_ideal, kl_results = kl_results))
}

gera_graficos_rbf <- function(sigma_ideal, kl_results) {
  valid_points <- !is.na(kl_results$KL) & !is.infinite(kl_results$KL)
  
  if (sum(valid_points) < 4) {
    plot(kl_results$sigma, kl_results$KL, type = "p", main = "KL vs Sigma (RBF) - No Spline Fit",
         xlab = "Sigma", ylab = "Divergência KL", col = "black", pch = 16)
    if (!is.na(sigma_ideal)) {
      abline(v = sigma_ideal, col = "red", lty = 2, lwd = 1.5)
      text(sigma_ideal, max(kl_results$KL, na.rm = TRUE), 
           paste("Sigma ideal:", round(sigma_ideal, 3)), pos = 4, col = "black")
    }
    warning("Not enough valid points to generate a spline for plotting. Showing raw points only.")
    return(invisible(NULL))
  }
  
  spline_fit <- smooth.spline(kl_results$sigma[valid_points], kl_results$KL[valid_points], spar = 0.6)
  
  plot(kl_results$sigma, kl_results$KL, type = "p", main = "KL vs Sigma (RBF)",
       xlab = "Sigma", ylab = "Divergência KL", col = "black", pch = 16)
  lines(spline_fit, col = "blue", lwd = 1.5)
  
  if (!is.na(sigma_ideal)) {
    abline(v = sigma_ideal, col = "red", lty = 2, lwd = 1.5)
    predicted_kl_at_ideal <- predict(spline_fit, x = sigma_ideal)$y
    text(sigma_ideal, predicted_kl_at_ideal, 
         paste("Sigma Ideal:", round(sigma_ideal, 3)), pos = 4, col = "black")
  }
  
  legend("bottomright", legend = c("Spline Fit", "Sigma Ideal"),
         col = c("blue", "red"), lty = c(1, 2), bty = "n")
}


datasets_to_evaluate <- c("PimaIndiansDiabetes", "BreastCancer", "Vehicle", "Ozone", "Ionosphere")
results_list <- list()

for (ds_name in datasets_to_evaluate) {

  current_sigma_range <- c(0.01, 5) 
  if (ds_name == "iris") {
    current_sigma_range <- c(0.01, 2) 
  }
  
  dataset_result <- avaliar_dataset_rbf(ds_name, sigma_range = current_sigma_range, num_sigmas = 50)
  results_list[[ds_name]] <- dataset_result
  
  gera_graficos_rbf(dataset_result$sigma_ideal, dataset_result$kl_results)
}

print("Ideal Sigma values for RBF Kernel (found by spline maximum):")
sapply(results_list, `[[`, "sigma_ideal")