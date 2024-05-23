library(Matrix)
library(geoR)
library(MASS)
library(corpcor)
library(mvnfast)
library(ggplot2)
library(reshape2)

# Define RMSE and SI functions
rmse <- function(o, p) sqrt(mean((o - p)^2))
si <- function(o, p) sqrt(mean(((o - mean(o)) - (p - mean(p)))^2)) / mean(o^2)


# Simulation parameters
lon <- seq(1, 15, 1)
lat <- seq(1, 15, 1)
lonlat <- expand.grid(lon, lat)
D <- as.matrix(dist(lonlat))
par_x <- c(6, 2, 3/2)
sigma_x <- par_x[1] * matern_cov(par_x[2:3], D)

# Simulate data
n <- nrow(lonlat)
x <-  MASS::mvrnorm(n, mu = rep(0, nrow(lonlat)), Sigma = sigma_x)
par_beta <- c(1, 1, 3/2)
beta <- simulate_beta(D, par_beta, mu = rep(0, nrow(lonlat)), cov.model = "matern")
s <- 1
y <- x %*% beta + rnorm(n, sd = s)

# Fit the model
result <- em_gr(x = x, y = y, lonlat = lonlat, cov.model = "matern")


# Evaluate the model
predicted_y <- x %*% result$mu
rmse_val <- rmse(y, predicted_y)
si_val <- si(y, predicted_y)

# Print results
cat("RMSE:", rmse_val, "\n")
cat("SI:", si_val, "\n")

# Plot the simulated and predicted beta spatial fields
predicted_beta <- result$mu

data <- data.frame(lonlat, beta = beta, predicted_beta = predicted_beta)
data_melt <- melt(data, id.vars = c("Var1", "Var2"), variable.name = "type", value.name = "value")

ggplot(data_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_wrap(~ type, scales = "free") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Simulated vs Predicted Beta Spatial Field", x = "Longitude", y = "Latitude") +
  theme_minimal()


# Fit the EM algorithm for Generalized Ridge regression with spatial covariates
result_matern <- em_gr(x = x, y = y, lonlat = lonlat, cov.model = "matern")
result_car <- em_gr(x = x, y = y, lonlat = lonlat, cov.model = "CAR")
result_ridge <- em_gr(x = x, y = y, lonlat = lonlat, cov.model = "ridge")

# Prepare data for plotting
data <- data.frame(lonlat, beta = beta, 
                   beta_matern = result_matern$mu, 
                   beta_car = result_car$mu, 
                   beta_ridge = result_ridge$mu)
data_melt <- melt(data, id.vars = c("Var1", "Var2"), variable.name = "type", value.name = "value")

# Plot the simulated and predicted beta spatial fields
ggplot(data_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_wrap(~ type, scales = "free", labeller = as_labeller(c(
    beta = "True Beta",
    beta_matern = "Estimated Beta (Matern)",
    beta_car = "Estimated Beta (CAR)",
    beta_ridge = "Estimated Beta (Ridge)"
  ))) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Simulated vs Estimated Beta Spatial Fields", x = "Longitude", y = "Latitude") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
