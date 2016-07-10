u <- UGP$new(package='GPfit',X=matrix(runif(10),5,2),Z=runif(5))
u$predict(matrix(runif(8),4,2))
u$update()
u$delete()
