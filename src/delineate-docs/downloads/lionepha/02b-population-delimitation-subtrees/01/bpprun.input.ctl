
          seed =  -1

       seqfile = bpprun.input.chars.txt
      Imapfile = bpprun.input.imap.txt
       outfile = results.out.txt
      mcmcfile = results.mcmc.txt

  speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 0        * species tree NNI/SPR
   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 9  L_casta_AK_Prince_of_Wales_Island L_casta_AK_Ketchikan L_casta_CA_West_Branch_Mill_Creek L_casta_CA_Mt_Tamalpais L_casta_CA_Soda_Creek L_casta_OR_Marys_Peak L_casta_OR_School_Creek L_casta_OR_Lost_Prairie L_casta_WA_Taneum_Creek
                    1 1 1 1 1 1 1 1 1
                    (((L_casta_OR_Lost_Prairie,(L_casta_OR_School_Creek,(L_casta_WA_Taneum_Creek,L_casta_OR_Marys_Peak))),(L_casta_CA_West_Branch_Mill_Creek,(L_casta_CA_Mt_Tamalpais,L_casta_CA_Soda_Creek))),(L_casta_AK_Ketchikan,L_casta_AK_Prince_of_Wales_Island));

        diploid =   1 1 1 1 1 1 1 1 1 

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 8  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.004   # invgamma(a, b) for theta
      tauprior = 3 0.002    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 10
       nsample = 1000
