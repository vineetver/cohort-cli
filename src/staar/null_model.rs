use faer::prelude::*;
use faer::Mat;

/// Fitted null model — immutable after construction.
///
/// Stores components of the projection matrix P, NOT the full n×n matrix.
/// P*G is computed on-the-fly as G - X*(X'X)^-1*X'*G which is O(n*k*m)
/// instead of O(n^2*m) and avoids materializing a 50K×50K matrix.
pub struct NullModel {
    pub residuals: Mat<f64>,
    pub x_matrix: Mat<f64>,
    pub xtx_inv: Mat<f64>,
    pub sigma2: f64,
    pub n_samples: usize,
    /// Fitted probabilities μ_i = P(Y_i = 1) under the null. Binary traits only.
    pub fitted_values: Option<Vec<f64>>,
}

impl NullModel {
    /// Compute (I - H) * G: project G into the residual space of X.
    ///
    /// Returns the HAT-complement projection, NOT scaled by 1/σ².
    /// Score tests handle the σ² scaling explicitly to keep the math transparent
    /// and match the STAARpipeline/GMMAT formulation exactly.
    ///
    /// (I-H)G = G - X(X'X)^-1 X'G
    ///
    /// Complexity: O(n * k * m) where k = covariates, m = variants in gene.
    pub fn project(&self, g: &Mat<f64>) -> Mat<f64> {
        let n = g.nrows();
        let m = g.ncols();
        assert_eq!(n, self.n_samples, "G rows must match n_samples");

        // X'G  (k × m)
        let xt_g = self.x_matrix.transpose() * g;
        // (X'X)^-1 X'G  (k × m)
        let xtx_inv_xt_g = &self.xtx_inv * &xt_g;
        // H*G = X(X'X)^-1 X'G  (n × m)
        let mut h_g = &self.x_matrix * &xtx_inv_xt_g;

        // (I - H) * G — reuse h_g allocation
        for j in 0..m {
            for i in 0..n {
                h_g[(i, j)] = g[(i, j)] - h_g[(i, j)];
            }
        }
        h_g
    }
}

/// Fit GLM null model for continuous trait via QR decomposition.
///
/// y = X * beta + epsilon
///
/// Caller must prepend an intercept column to X if desired.
pub fn fit_glm(y: &Mat<f64>, x: &Mat<f64>) -> NullModel {
    let n = y.nrows();
    let k = x.ncols();
    assert_eq!(x.nrows(), n, "X rows must match y rows");
    assert!(n > k, "Need more samples than covariates");

    // Solve via normal equations: beta = (X'X)^-1 X'y
    // For n >> k (typical: 50K samples, 10 covariates), this is efficient and stable.
    let xtx = x.transpose() * x;           // k × k
    let xty = x.transpose() * y;           // k × 1
    let eye_k = Mat::<f64>::identity(k, k);
    let xtx_inv = xtx.col_piv_qr().solve(&eye_k); // (X'X)^-1
    let beta = &xtx_inv * &xty;            // beta = (X'X)^-1 X'y

    // Residuals: r = y - X * beta
    let y_hat = x * &beta;
    let mut residuals = Mat::zeros(n, 1);
    for i in 0..n {
        residuals[(i, 0)] = y[(i, 0)] - y_hat[(i, 0)];
    }

    // Residual variance: sigma2 = ||r||^2 / (n - k)
    let rss: f64 = (0..n).map(|i| residuals[(i, 0)].powi(2)).sum();
    let sigma2 = rss / (n - k) as f64;

    NullModel {
        residuals,
        x_matrix: x.to_owned(),
        xtx_inv,
        sigma2,
        n_samples: n,
        fitted_values: None,
    }
}

/// Fit logistic null model for binary trait via IRLS.
///
/// y = 0/1, logit(P(y=1)) = X * beta
pub fn fit_logistic(y: &Mat<f64>, x: &Mat<f64>, max_iter: usize) -> NullModel {
    let n = y.nrows();
    let k = x.ncols();
    let max_iter = if max_iter == 0 { 25 } else { max_iter };

    // Initialize beta to zeros
    let mut beta = Mat::zeros(k, 1);

    for _ in 0..max_iter {
        // mu = sigmoid(X * beta)
        let eta = x * &beta;
        let mut mu = Mat::zeros(n, 1);
        let mut w_diag = Mat::zeros(n, 1); // diagonal of W
        for i in 0..n {
            let eta_i = eta[(i, 0)].clamp(-500.0, 500.0); // prevent exp overflow
            let p = 1.0 / (1.0 + (-eta_i).exp());
            mu[(i, 0)] = p;
            w_diag[(i, 0)] = p * (1.0 - p);
        }

        // Working response: z = eta + (y - mu) / w
        // Weighted least squares: X' * W * X * delta = X' * W * z
        let mut xtwx: Mat<f64> = Mat::zeros(k, k);
        let mut xtwy: Mat<f64> = Mat::zeros(k, 1);

        for i in 0..n {
            let wi = w_diag[(i, 0)].max(1e-10);
            let zi = eta[(i, 0)] + (y[(i, 0)] - mu[(i, 0)]) / wi;
            for j in 0..k {
                xtwy[(j, 0)] += x[(i, j)] * wi * zi;
                for l in 0..k {
                    xtwx[(j, l)] += x[(i, j)] * wi * x[(i, l)];
                }
            }
        }

        let new_beta: Mat<f64> = xtwx.col_piv_qr().solve(&xtwy);

        // Check convergence
        let mut max_diff = 0.0f64;
        for j in 0..k {
            let diff: f64 = new_beta[(j, 0)] - beta[(j, 0)];
            max_diff = max_diff.max(diff.abs());
        }
        beta = new_beta;
        if max_diff < 1e-8 {
            break;
        }
    }

    // Final residuals and fitted probabilities.
    // Score test uses U = G'r where r = Y - μ (raw score residuals, NOT working).
    // SPA needs the fitted μ_i to compute the exact CGF of the score distribution.
    let eta = x * &beta;
    let mut residuals = Mat::zeros(n, 1);
    let mut fitted = Vec::with_capacity(n);
    for i in 0..n {
        let eta_i = eta[(i, 0)].clamp(-500.0, 500.0);
        let mu = 1.0 / (1.0 + (-eta_i).exp());
        fitted.push(mu);
        residuals[(i, 0)] = y[(i, 0)] - mu;
    }

    let xtx = x.transpose() * x;
    let eye_k = Mat::<f64>::identity(k, k);
    let xtx_inv = xtx.col_piv_qr().solve(&eye_k);

    NullModel {
        residuals,
        x_matrix: x.to_owned(),
        xtx_inv,
        sigma2: 1.0,
        n_samples: n,
        fitted_values: Some(fitted),
    }
}
