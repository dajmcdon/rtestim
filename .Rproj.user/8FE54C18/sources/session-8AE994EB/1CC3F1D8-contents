library(rlang)

generate_I <- function(n) {
    stopifnot(n > 0, is_integerish(n))
    Matrix::bandSparse(n,
        k = 0,
        diagonals = list(rep(1, n))
    )
}

# divided difference matrices of dim (n-2)*n:
generate_D <- function(n) {
    stopifnot(n > 0, is_integerish(n))
    Matrix::bandSparse(n,
        k = c(0, 1),
        diagonals = list(rep(-1, n), rep(1, n - 1))
    )[-n, ]
}

generate_Dk <- function(n, k = 1) {
    stopifnot(n > 0, n > k + 1, k >= 0)
    stopifnot(is_integerish(n), is_integerish(k))
    D0 <- generate_D(n)
    if (k > 0) D0 <- Matrix::diff(D0, differences = k)
    D0
}

generate_D2 <- function(n, type = c("2", "a", "b")) {
    type <- match.arg(type)
    switch(type,
        `2` = generate_Dk(n, 1),
        a = generate_D(n)[-(n - 1), ],
        b = cbind(0, generate_D(n - 1))
    )
}
