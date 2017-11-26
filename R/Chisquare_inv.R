###############################################################################################################
#
# CHISQUARE_INV  Inverse of chi-square cumulative distribution function (cdf).
#
#   X = chisquare_inv(P,V) returns the inverse of chi-square cdf with V
#   degrees of freedom at fraction P.
#   This means that P*100 percent of the distribution lies between 0 and X.
#
#   To check, the answer should satisfy:   P==gammainc(X/2,V/2)
#
#   Uses FMIN and CHISQUARE_SOLVE
#
#   Written January 1998 by C. Torrence
#-------------------------------------------------------------------------------------------------------------

chisquare_inv<-function(P,V){

	nargin=length(as.list(match.call()))-1	

	if(nargin<2){stop('Must input both P and V')}
	if((1-P)< 1E-4){stop('P must be < 0.9999')}
	
	if((P==0.95) & (V==2)){
	X = 5.9915					# this is a no-brainer
	return(X)
	}
		
	MINN=0.01        				# hopefully this is small enough
	MAXX=1            				# actually starts at 10 (see while loop below)
	X = 1
	TOLERANCE=1E-4   				# this should be accurate enough

	while((X+TOLERANCE)>=MAXX){  			# should only need to loop thru once
		MAXX=MAXX*10

	# this calculates value for X, NORMALIZED by V

            X=as.numeric(fminbnd(chisquare_solve,a=MINN,b=MAXX,tol=TOLERANCE,P=P,V=V)[1])

		MINN=MAXX
	  }


	X=X*V  						# put back in the goofy V factor

	return(X)
}