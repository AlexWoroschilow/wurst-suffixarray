
/*
 *  wurstimbiss.rescue.c
 *  snippets 
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 04/03/07 15:59:33 CEST
 *  
 */



		 /*
	  if(verbose_flag) {
	  	printf("lambda: %19.16e\n", info->lambda);
	  	printf("check: %19.16e\n", 
			checklambda(info->scr, alphabet->domainsize, info->sf, 
			  				avgsum, info->lambda));
	   }
	  */
	  /*info->H= relentropy(info->sortind, info->scr, alphabet->domainsize,
				info->scrf, info->lambda);
	  info->K=entropyK(info->sortind, info->scr, alphabet->domainsize, 
		  		avgsum, info->scrf, info->lambda, info->H);
	  printf("K:%19.16e\n", info->K);
	  printf("inputscr: %f\n", inputscr);
	  printf("exp: %19.16e\n", exp(-info->lambda*inputscr));
	  printf("p:%19.16e\n", info->K*exp(-info->lambda*inputscr));
	*/



	/*blastlike statistics
	imbiss->df=dbfreq(space, sequences, noofvecs, alphabet, 1);
	imbiss->sf=seqfreq(space, input, alphabet);
	
  
	  scr = logoddscr(space, imbiss->df, imbiss->sf, alphabet);
	  imbiss->score = scr;
		for(i=0; i < alphabet->domainsize; i++) {
	
			avgsum += imbiss->df[i]*imbiss->score[i];
	  }
	
	  for (i=0; i< input->length; i++) {
	  	inputscr += scr[input->sequence[i]]; 
	  }
	  
	  printf("\nBLAST statistics:\n-------------------\n");
	  printf("E(score): %f\n", avgsum);
	
	  printf("inputscr: %f\n", inputscr);
	  imbiss->scrf = scorefreq(space, scr, alphabet->domainsize, alphabet, 
		  				imbiss->sf, imbiss->df);
	  writescores(scr, alphabet->domainsize, "score.csv");	
      writescores(imbiss->scrf, alphabet->domainsize, "scoref.csv");
	  writescores(imbiss->df, alphabet->domainsize, "freq.csv");
	  writescores(imbiss->sf, alphabet->domainsize, "seqf.csv");	   
	  imbiss->sortind = quickSort(space, scr, alphabet->domainsize, 
		  				cmp_dbl, NULL);
	
	  evd=ALLOCMEMORY(space, NULL, sizeof(evdparam), 1);
	  evd->noofscores = alphabet->domainsize;
	  evd->probs = imbiss->df;
	  evd->scores = imbiss->score;

	  imbiss->lambda = uniroot(0, 1, score_evd, 0.0000001, evd); 
	  printf("lambda: %19.16e\n", imbiss->lambda); 
	 		 	
	  printf("check: %19.16e\n", 
			checklambda(imbiss->score, alphabet->domainsize, imbiss->df, 
			  				avgsum, imbiss->lambda)); 

	  imbiss->H= relentropy(imbiss->sortind, imbiss->score, alphabet->domainsize,
				imbiss->df, imbiss->lambda);
	  	
	  printf("rel. entropy H: %19.16e\n", imbiss->H);
	
	  imbiss->K=entropyK(imbiss->sortind, imbiss->score, alphabet->domainsize, 
		  		avgsum, imbiss->df, imbiss->lambda, imbiss->H);
	 */
