
          seed =  -1

       seqfile = bpprun.input.chars.txt
      Imapfile = bpprun.input.imap.txt
       outfile = results.out.txt
      mcmcfile = results.mcmc.txt

  speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 0        * species tree NNI/SPR
   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 8  L_lindrothi_CA_South_Fork_Bishop_Creek L_lindrothi_CA_Emerald_Lake L_lindrothi_CA_East_Fork_Kaweah_River L_lindrothi_CA_Kaiser_Pass L_lindrothi_CA_Tioga_Lake L_lindrothi_CA_Sonora_Pass L_lindrothi_CA_Deadman_Creek L_lindrothi_CA_Long_Valley_Creek
                    1 3 1 1 1 1 1 1
                    ((L_lindrothi_CA_Tioga_Lake,(L_lindrothi_CA_Long_Valley_Creek,(L_lindrothi_CA_Deadman_Creek,L_lindrothi_CA_Kaiser_Pass))),(L_lindrothi_CA_South_Fork_Bishop_Creek,(L_lindrothi_CA_Emerald_Lake,(L_lindrothi_CA_East_Fork_Kaweah_River,L_lindrothi_CA_Sonora_Pass))));

        diploid =   1 1 1 1 1 1 1 1 

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
