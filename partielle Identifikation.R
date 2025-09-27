library(MASS) #für Multinomialverteilung
library(ivreg) #eingebauter TSLS

set.seed(123) # zur Replizierung der Ergebnisse

S <- 1000   # Anzahl Simulationen
N <- 50   # Größe der Stichproben

# korrekte Werte für die Koeffizienten:
beta_0 <- 2
beta_1 <- 1

# korrekter Wert der Missspezifikation
Modellfehler <- 0.8


#Funktion "data" erstellt alle Variablen für einen Simulationsdurchgang
data <- function(){
  # cov_xz bezeichnet jeweilige Kovarianzen (hier von Regressor x und IV z), e ist ein Teil des Fehlerterms der Regression
  cov_xz <- 0.5           # Kovarianz beschreibt die Stärke des Instruments
  cov_xe <- 0             # Exogenitätsannahme wird aufgrund des später hinzugefügten Messfehlers ungültig
  cov_ze <- Modellfehler  # falls Modell missspezifiziert sein soll, gilt Cov(Z_i, e_i) =/= 0
  
  daten <- mvrnorm(N, c(1, 1, 0), matrix(c(1, cov_xz, cov_xe, # Regressor x, IV Z und Fehlerterm e wird erstellt
                                           cov_xz, 1, cov_ze,
                                           cov_xe, cov_ze, 1), 
                                         nrow=3, ncol=3, byrow=TRUE))
  
  x_true <- daten[,1]
  Z <- daten[,2]  # Instrument, korreliert mit x und e
  e <- daten[,3]  # Fehlerterm
  
  
  Messfehler <- rnorm(N, 0, 1)
  x_beob <- x_true + Messfehler # x mit Messfehler, die einzigen Daten, die im Algorithmus weiter genutzt werden
  # Kovarianz von x_beob und eps* ist also ungleich 0 in der später folgenden Regression
  # eps* = eps - beta_1 * Messfehler
  
  nu <- rnorm(N, 0, 1)
  y <- beta_0 + beta_1 * x_true + e/sqrt(N) + nu  #IV ist also korreliert mit e
  
  return(list(x_messfehl = x_beob, IV = Z, y = y))
}


#Two-stage least squares wird genutzt
TSLS <- function(x, y, IV){
  model_TSLS <- ivreg(y ~ x | IV)
  coeff_beta0 <- model_TSLS$coefficients["(Intercept)"]
  coeff_beta1 <- model_TSLS$coefficients["x"]
  
  conf_int <- confint(model_TSLS)[2,] # Konfidenzintervall mit Niveau 95%
  
  return(list(beta0=coeff_beta0, beta1=coeff_beta1, ci=conf_int))
}


# Berechnung von M_min (im part. Ident. Algorithmus) --> nur in einem überidentifizierten Modell brauchbar
M_min <- function(g_hat, W_hat, Sigma_hat, Gamma_hat, c){
  # Berechne die untere Grenze für M, siehe hierzu Appendix B von Armstrong und Kolesár
  
  J <- N * t(g_hat) %*% W_hat %*% g_hat # J-Statistik
  
  # nun wie in Appendix B beschrieben, berechne den Wert, gegen den die J-Statistik konvergiert: J_crit= c' * Sigma^(-0.5) * R * Sigma^(-0.5) * c
  
  # Sig_minus_half beschreibt Sigma^(-0.5)
  ev <- eigen(Sigma_hat, symmetric = TRUE)
  Dmh <- diag(1 / sqrt(pmax(ev$values, 1e-12)))
  Sig_minus_half <- ev$vectors %*% Dmh %*% t(ev$vectors)
  
  R_hat <- diag(2) - Sig_minus_half %*% Gamma_hat %*% solve(t(Gamma_hat) %*% solve(Sigma_hat) %*% Gamma_hat) %*% t(Gamma_hat) %*% Sig_minus_half  # für Definition von R: siehe Appendix B
  nc_par <- as.numeric(t(c) %*% Sig_minus_half %*% R_hat %*% Sig_minus_half %*% c) # dies ist der nichtzentrale Parameter, gegen den die J-Statistik konvergiert
  
  matrix_A <- Sig_minus_half %*% R_hat %*% Sig_minus_half
  # Kleinster Eigenwert von Matrix_A^{-1}
  evals <- eigen(ginv(matrix_A), symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(evals)
  
  J_crit <- sqrt(qchisq(1- 0.05, df=1, ncp=nc_par^2))
  M_min <- sqrt( max(0, (J - J_crit)/ (N* lambda_min) ) )  # berechne untere Grenze für M
  
  return(M_min)
}

# Algorithmus partielle Identifikation
partIdent <- function(x, y, IV){
  # schreibe die Variablen in Matrixform:
  Z <- cbind(1, IV) # Instrument
  X <- cbind(1, x)  # Regressor
  
  # 1) Berechne den initialen Schätzer von beta über TSLS-Algorithmus
  model2 <- TSLS(x, y, IV)
  b_0 <- model2$beta0 # TSLS-Schätzer für beta_0
  b_1 <- model2$beta1 # TSLS-Schätzer für beta_1
  beta_iv <- as.numeric(c(b_0, b_1))  # initialer Schätzer von beta
  
  # -----------------------
  # 2) Momentenbedingung g(beta) = mean(Z_i * (Y_i - X_i' beta_hat))
  eps_hat <- as.vector(y - X %*% beta_iv)
  g_i <- Z * eps_hat
  g_hat <- colMeans(g_i)    # g_hat = (mean(eps], mean(Z*eps)), Vektor der Länge 2
  
  # dazu Schätzung von Sigma = Var(g(beta)), wird für CI später gebraucht
  Sigma_hat <- cov(g_i)       # (2*2)-Matrix
  
  # -----------------------
  # 3) Berechne Schätzer von Sensitivität k für H=(0,1), also für beta_1
  
  # Schätzung von Gamma = - E[Z_i X_i']
  Gamma_hat <- - crossprod(Z, X) / N   # (2*2)-Matrix
  
  # Schätzung von Gewichtungsmatrix W (verwende Sigma^{-1})
  W_hat <- tryCatch(solve(Sigma_hat), error = function(e) MASS::ginv(Sigma_hat))
  
  # nun zur Schätzung von k: k' = - H (Gamma' W Gamma)^{-1} Gamma' W
  H <- matrix(c(0, 1), nrow = 1)
  
  A <- t(Gamma_hat) %*% W_hat %*% Gamma_hat
  A_inv <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  k_t <- - H %*% A_inv %*% t(Gamma_hat) %*% W_hat  # Vektor mit Länge 2
  k <- as.numeric(t(k_t))  # k wird zu einem vektor definiert
  
  #------------------------
  # 4) Berechne den initialen Schätzer für beta_1 = beta_{1,iv} + k' * g(beta)
  h_hat <- beta_iv[2] + t(k) %*% g_hat
  
  # -----------------------
  # 5) Berechne nun den Worst-case bias = M * sup_{c aus C(1)} |k * c|
  #definiere hierzu Menge C(1):
  B <- matrix(c(mean(Z[,2]), mean(Z[,2]^2)), nrow=2)  # Missspezifizierung betrifft nur Instrument Z_{i1}=Z[,2], nur diese Werte werden in B berücksichtigt
  c <- B # Supremum aus C(1)
  
  # Berechne nun mögliche Werte für M:
  #M_min <- M_min(g_hat, W_hat, Sigma_hat, Gamma_hat, c) # wenn d_Z = d_b 
  # --> Berechnung von M_min im genau identifizierten Fall rechnerisch immer ca. 0, die J-Statistik ist hier leider nicht sinnvoll
  
  M <- Modellfehler/sqrt(N)  # setze M = 0 im korrekten Fall, wähle M aus (0, unendlich) im missspezierten Fall
  
  # worst case-bias:
  bias_hat <- M * k %*% c
  
  # -----------------------
  # Nun geht es an die Berechnung des Konfidenzintervalls:
  # 6) Berechne Varianz von k
  var_k <- as.numeric(t(k) %*% Sigma_hat %*% k) 
  se_hat <- sqrt(var_k / N)
  
  # -----------------------
  # 7) Kritischer Wert cv_alpha(t) (für zweiseitiges 1-alpha Konfidenzintervall)
  t_val <- if (var_k > 0) bias_hat / sqrt(var_k) else 0   # t ist der "nichtzentrale Parameter" der Chi^2 Verteilung, t = M*bias_{C(1)} / sqrt(k * Sigma * k)
  alpha <- 0.05
  cv <- sqrt(qchisq(1-alpha, df=1, ncp=t_val^2))
  
  # -----------------------
  # 8) Konfidenzintervall für beta_1
  beta1_hat <- h_hat
  ci_lower <- beta1_hat - cv * se_hat
  ci_upper <- beta1_hat + cv * se_hat
  
  confint_pI <- c(ci_lower, ci_upper)
  
  return(list(h_hat=beta1_hat, conf_int = confint_pI))
}


# In die Listen werden die jeweiligen Schätzer (für jede Methode) für beta_1 gespeichert.
coeff_TSLS <- numeric()
coeff_pI <- numeric()
CI_pI <- matrix(nrow=S, ncol=2)
CI_TSLS <- matrix(nrow=S, ncol=2)

Zähler_pI <- 0
Zähler_TSLS <- 0

# Hier werden die Simulationsdurchgänge gestartet
for (s in 1:S){
  # Hier werden die Daten für den Durchgang generiert und in die Variablen abgespeichert
  variablen <- data() 
  x <- variablen$x_messfehl
  y <- variablen$y
  IV <- variablen$IV
  
  # Two stage least squares
  TSLS_variablen <- TSLS(x, y, IV)
  coeff_TSLS[s] <- TSLS_variablen$beta1
  CI_TSLS[s,] <- TSLS_variablen$ci
  
  # Partielle Identifikation
  partIdent_variablen <- partIdent(x, y, IV)
  coeff_pI[s] <- partIdent_variablen$h_hat
  CI_pI[s,] <- partIdent_variablen$conf_int
  
  # Prüfen der Konfidenzniveaus der Intervalle
  if (beta_1 >= CI_pI[s,1] & beta_1 <= CI_pI[s,2]) {Zähler_pI <- Zähler_pI + 1}
  if (beta_1 >= CI_TSLS[s,1] & beta_1 <= CI_TSLS[s,2] ) {Zähler_TSLS <- Zähler_TSLS + 1}
}

#-----------------
# Hier wird das Histogramm von Schätzwerten für beta_1 erstellt
Grenze_plot <- 0.5

# partielle Identifikation:
hist(coeff_pI, xlab="Wert für Schätzer von beta_1", ylab="Häufigkeit", 
     main="Vergleich der Schätzer von partieller Identifikation und TSLS", 
     col = rgb(1,0,0,0.2), xlim = c(beta_1-Grenze_plot, beta_1+Grenze_plot), breaks="FD") 

# TSLS
hist(coeff_TSLS, col = rgb(0,1,0,0.2), add=TRUE, breaks="FD")

abline(v=beta_1, col="red")  #Linie an der Stelle, wo der wahre Wert von beta_1 liegt
legend("topleft", c("TSLS", "partielle Identifikation"), fill=c(rgb(0,1,0,0.2), rgb(1,0,0,0.2)))


#------------------
# Konfidenzniveaus
prob_pI <- Zähler_pI/S
prob_TSLS <- Zähler_TSLS/S


# Werte für TSLS:
mean(coeff_TSLS); var(coeff_TSLS) 
mean(CI_TSLS[,2]-CI_TSLS[,1]); median(CI_TSLS[,2]-CI_TSLS[,1])

# Werte für partielle Identifikation:
mean(coeff_pI); var(coeff_pI) 
mean(CI_pI[,2]-CI_pI[,1]); median(CI_pI[,2]-CI_pI[,1])

# Differenz der Schätzer
mean(coeff_TSLS- coeff_pI)


