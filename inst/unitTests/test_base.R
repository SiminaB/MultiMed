##test that results are the same whether E and Y are vectors or matrices
##and whether M is a data frame or a matrix
##for the single mediator case, also tests whether the results are the same
##whether M is a vector or a matrix
test_vectMatDataFrame <- function() {
    set.seed(20183)
    E <- rnorm(100)
    ##multiple mediators
    M <- matrix(rnorm(100*4), nrow=100, ncol=4)
    Y <- rnorm(100)
    set.seed(100)
    medTest1 <- medTest(E, M, Y, nperm=10)
    set.seed(100)
    medTest2 <- medTest(matrix(E, ncol=1), M, matrix(Y, ncol=1), nperm=10)
    set.seed(100)
    medTest3 <- medTest(matrix(E, ncol=1), data.frame(M),
                        matrix(Y, ncol=1), nperm=10)
    set.seed(100)
    medTest4 <- medTest(E, data.frame(M),
                        matrix(Y, ncol=1), nperm=10)

    ##single mediator
    set.seed(100)
    medTest11 <- medTest(E, M[,1], Y, nperm=10)
    set.seed(100)
    medTest12 <- medTest(E, M[,1,drop=FALSE], Y, nperm=10)

    checkEquals(medTest1, medTest2)
    checkEquals(medTest1, medTest3)
    checkEquals(medTest1, medTest4)
    checkEquals(medTest11, medTest12)
}

##test that the largest statistic has the smallest p-value for multiple
##mediators
test_largestSsmallestP <- function()
{
    ##do this for 4 and 9 mediators, to check different scenarios
    set.seed(20183)
    E <- rnorm(100)
    ##multiple mediators
    M <- matrix(rnorm(100*4), nrow=100, ncol=4)
    Y <- rnorm(100)
    set.seed(100)
    medTest1 <- medTest(E, M, Y, nperm=10)
    argMaxS1 <- which.max(medTest1[,"S"])
    argMinP1 <- which.min(medTest1[,"p"])

    set.seed(20183)
    E <- rnorm(100)
    ##multiple mediators
    M <- matrix(rnorm(100*9), nrow=100, ncol=9)
    Y <- rnorm(100)
    set.seed(100)
    medTest2 <- medTest(E, M, Y, nperm=10)
    argMaxS2 <- which.max(medTest2[,"S"])
    argMinP2 <- which.min(medTest2[,"p"])

    ##first check that the statistics, p-values are not the same
    ##for all the potential mediators
    checkTrue(length(unique(medTest1[,"S"])) > 1)
    checkTrue(length(unique(medTest1[,"p"])) > 1)
    checkTrue(length(unique(medTest2[,"S"])) > 1)
    checkTrue(length(unique(medTest2[,"p"])) > 1)
    ##now check that largest test statistic corresponds to smallest
    ##p-value
    checkEquals(argMaxS1, argMinP1)
    checkEquals(argMaxS2, argMinP2)
}

##test that equal values of w give the same result, whether it's a
##single number or a vector of the same length as E
test_equalW <- function()
{
    set.seed(20183)
    E <- rnorm(100)
    M <- matrix(rnorm(100*4), nrow=100, ncol=4)
    Y <- rnorm(100)
    set.seed(100)
    medTest1 <- medTest(E, M, Y, nperm=10, w=1)
    set.seed(100)
    medTest2 <- medTest(E, M, Y, nperm=10, w=5)
    set.seed(100)
    medTest3 <- medTest(E, M, Y, nperm=10, w=rep(1,100))
    set.seed(100)
    medTest4 <- medTest(E, M, Y, nperm=10, w=rep(5,100))
    set.seed(100)
    medTest5 <- medTest(E, M, Y, nperm=10, w=rep(1/100,100))

    checkEquals(medTest1, medTest2)
    checkEquals(medTest1, medTest3)
    checkEquals(medTest1, medTest4)
    checkEquals(medTest1, medTest5)
}

