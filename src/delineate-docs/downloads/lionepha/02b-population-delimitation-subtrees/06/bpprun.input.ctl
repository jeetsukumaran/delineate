
          seed =  -1

       seqfile = bpprun.input.chars.txt
      Imapfile = bpprun.input.imap.txt
       outfile = results.out.txt
      mcmcfile = results.mcmc.txt

  speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 0        * species tree NNI/SPR
   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 15  L_osculans_CA_Cold_Creek L_osculans_CA_Rainbow L_osculans_CA_Stanislaus_NF L_osculans_CA_Warner_Range L_osculans_CA_Los_Padres_NF L_osculans_CA_Strawberry_Creek L_osculans_CA_Carson_Spur L_osculans_CA_Nanny_Creek L_osculans_CA_Pike_County_Gulch L_osculans_OR_School_Creek L_osculans_OR_Middle_Fork_Berry_Creek L_osculans_OR_Goodman_Creek L_osculans_OR_Eugene L_osculans_OR_Walton_Lake L_osculans_OR_Little_Philips_Creek
                    2 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                    (((((L_osculans_CA_Stanislaus_NF,L_osculans_CA_Rainbow),(L_osculans_CA_Los_Padres_NF,L_osculans_CA_Warner_Range)),(L_osculans_CA_Cold_Creek,L_osculans_CA_Strawberry_Creek)),(L_osculans_CA_Carson_Spur,(L_osculans_CA_Nanny_Creek,L_osculans_CA_Pike_County_Gulch))),(L_osculans_OR_Little_Philips_Creek,((L_osculans_OR_Walton_Lake,L_osculans_OR_Middle_Fork_Berry_Creek),(L_osculans_OR_School_Creek,(L_osculans_OR_Eugene,L_osculans_OR_Goodman_Creek)))));

        diploid =   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 

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
