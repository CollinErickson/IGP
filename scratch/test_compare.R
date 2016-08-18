compare.UGP(packages=c("laGP"), func=TestFunctions::banana, D=2, N=50, Npred = 1e3, reps = 5)
compare.UGP(packages=c("laGP", "mlegp", "GauPro"), func=TestFunctions::banana, D=2, N=50, Npred = 1e3, reps = 5)
compare.UGP(packages=c("GauPro", "GauPro","laGP", "GPy","sklearn"), func=TestFunctions::RFF_get(D=4), D=4, N=200, Npred = 1e3, reps = 10)
compare.UGP(packages=c("GauPro", "GauPro","laGP", "mlegp", "GPy","sklearn"), func=TestFunctions::RFF_get(D=2), D=2, N=40, Npred = 1e3, reps = 5)
