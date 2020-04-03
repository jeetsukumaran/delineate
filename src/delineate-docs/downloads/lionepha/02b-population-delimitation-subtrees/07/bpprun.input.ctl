
          seed =  -1

       seqfile = bpprun.input.chars.txt
      Imapfile = bpprun.input.imap.txt
       outfile = results.out.txt
      mcmcfile = results.mcmc.txt

  speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 0        * species tree NNI/SPR
   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 30  L_probata_BC_Summit_Creek L_probata_CA_Middle_Martis_Creek L_probata_CA_South_Fork_Bishop_Creek L_probata_CA_Sherman_Pass L_probata_CA_Strawberry_Creek L_probata_CA_White_Mountains L_probata_CA_Algoma_Camp L_probata_CA_Warner_Range L_probata_CA_Nanny_Creek L_probata_CA_Tamarack_Lake L_probata_CA_Ellery_Lake L_probata_CA_Squaw_Valley_Resort L_probata_ID_Galena_Summit L_probata_ID_Baker_Creek L_probata_ID_Park_Creek L_probata_MT_Prospect_Creek L_probata_MT_Mill_Creek L_probata_NV_Ruby_Mountains L_probata_OR_Steens_Mountains L_probata_OR_Mt_Ashland L_probata_OR_Lost_Prairie L_probata_OR_Odell_Creek L_probata_OR_Lonesome_Spring L_probata_OR_Lostine_River_Valley L_probata_OR_Little_Philips_Creek L_probata_UT_Stansbury_Mtns L_probata_UT_Shingle_Creek L_probata_UT_Tushar_Mountains L_probata_WA_Taneum_Creek L_probata_WA_Blue_Mountains
                    1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 2 1 1 1 1
                    ((((L_probata_OR_Odell_Creek,(L_probata_WA_Taneum_Creek,(L_probata_OR_Lost_Prairie,(L_probata_CA_Algoma_Camp,L_probata_CA_Middle_Martis_Creek)))),(((L_probata_CA_Strawberry_Creek,L_probata_CA_Nanny_Creek),(L_probata_CA_Ellery_Lake,L_probata_OR_Mt_Ashland)),(L_probata_CA_Sherman_Pass,(L_probata_CA_Tamarack_Lake,L_probata_CA_Squaw_Valley_Resort)))),(((L_probata_ID_Park_Creek,L_probata_MT_Prospect_Creek),((L_probata_WA_Blue_Mountains,(L_probata_CA_Warner_Range,(L_probata_CA_White_Mountains,L_probata_CA_South_Fork_Bishop_Creek))),(L_probata_OR_Lostine_River_Valley,(L_probata_BC_Summit_Creek,L_probata_OR_Steens_Mountains)))),(L_probata_MT_Mill_Creek,((L_probata_OR_Lonesome_Spring,L_probata_ID_Baker_Creek),(L_probata_OR_Little_Philips_Creek,(L_probata_NV_Ruby_Mountains,L_probata_ID_Galena_Summit)))))),(L_probata_UT_Shingle_Creek,(L_probata_UT_Stansbury_Mtns,L_probata_UT_Tushar_Mountains)));

        diploid =   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 

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
