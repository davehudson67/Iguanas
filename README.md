# Iguanas

So have initially tried to fit a mortliaty model to the data... just so we can get an underlying shape for the mortality which we then can carry forward to add covariates etc.

The Fit... R code files do this for each of exponential, gompertz, gompertz-makeham and siler models, the known age one was just to see if I could get it working. Then move onto using all the data (unknown births and deaths) with the double censored files.

Then the InitialModel Comparisons should compare the model fits and create boostrapped CI etc but thats where its throwing a few wierd things and hoping Tj may be able to help... he hasnt come back to me yet.


I have today been trying to fit growth models to the data - which hopefully can inform the model and help narrow down the ranges for the unknown births etc. But havent managed to get it to work properly.

Aadpting code from this paper:
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13842
