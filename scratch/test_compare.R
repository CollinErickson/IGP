compare.UGP(packages=c("laGP"), func=TestFunctions::banana, D=2, N=50, Npred = 1e3, reps = 5)
compare.UGP(packages=c("laGP", "mlegp", "GauPro"), func=TestFunctions::banana, D=2, N=50, Npred = 1e3, reps = 5)
compare.UGP(packages=c("GauPro-Par", "GauPro-NP","laGP", "GPy"), func=TestFunctions::RFF_get(D=2), D=2, N=50, Npred = 1e3, reps = 4,
            , init_list = list('1'=list(parallel=T), '2'=list(parallel=F)))
compare.UGP(packages=c("GauPro-Dev", "GauPro-LLH","laGP", "GPy", "sklearn"), func=TestFunctions::RFF_get(D=2), D=2, N=50, Npred = 1e3, reps = 4,
            , init_list = list('1'=list(useLLH=F), '2'=list(useLLH=T)))
compare.UGP(packages=c("GauPro", "GauPro","laGP", "mlegp", "GPy","sklearn"), func=TestFunctions::RFF_get(D=2), D=2, N=40, Npred = 1e3, reps = 5)
