
   
# takes a matrix object corresponding to a triad and 
# outputs what type of triad it is.


itpeval = function(tijk){
	#for debugging:
	#ijk = c(2,3,4)
	#tijk = submat
	
	#tijk is a vector 
	pos = c(0,0,0) # Creates a vector for output
	stijk = sum(tijk)
	if(stijk==0){ # triad empty: should never happen #?
		pos=c(1,1,1)
	}
	
	if (stijk==1){ #one arc in triad: likely most cases
		sndr = which(tijk==1, arr.ind=TRUE)[1]
		rcvr = which(tijk==1, arr.ind=TRUE)[2]
		# the 'which(condition, arr.ind=TRUE)' function returns the row x col 
		# position in the matrix which meets the condition specified
		# the sender is on the row and the reciever is on the column
		
		#islt=setdif(1:3,sndr||rcvr); is equivalent to:
		islt = c(1:3)[!(1:3 %in% c(sndr,rcvr))]
		pos[sndr]=2
		pos[rcvr]=3
		pos[islt]=4
	}
	
	if(stijk==2){
		m=(.5)*(sum(diag(tijk%*%tijk))) # take the trace of the matrix * matrix
		od = apply(tijk, 1, sum) #row sum
		.in = apply(tijk, 2, sum) #col sum
		if (m==0){
			if ((max(od)==2) & (max(.in)==1)){ #o21d 
				sndr = which(od==2)
				rcvr = c(1:3)[!(1:3 %in% sndr)] # note that setdiff(1:3,c(sndr,rcvr)) -- 2 ff's in setdiff! -- 
				# would have also worked
				pos[sndr]=7
				pos[rcvr[1]]=8
				pos[rcvr[2]]=8
			} else if((max(od)==1) & (max(.in)==2)){ #o21u 
				sndr=which(od==1)
				rcvr = c(1:3)[!(1:3 %in% sndr)]
				pos[sndr[1]]=9
				pos[sndr[2]]=9
				pos[rcvr]=10
			} else{ #o21c 
				dif = od-.in
				sndr = which(dif==1);
				brdg = which(dif==0);
				endr = which(dif==-1);
				pos[sndr]=11
				pos[brdg]=12
				pos[endr]=13
				#cat(tijk); #cat(ijk); #cat('021C')
			} 
		} else{ #must be 102
			dyad=which(od==1)
			islt=c(1:3)[!(1:3 %in% dyad)]
			pos[dyad[1]]=5
			pos[dyad[2]]=5
			pos[islt]=6
		} 
	} #stijk==2
	
	if(stijk==3){
		m=(.5)*(sum(diag(tijk%*%tijk))) 
		od = apply(tijk, 1, sum) #row sum
		
		if(m==1){
			if(max(od)==1){
				.in = apply(tijk, 2, sum) #col sum
#				/* #cat '111D' i j k tijk; */
				brdg=which(.in==2)
				sndr=which(.in==0)
				endr=which(.in==1)
				pos[sndr]=14
				pos[brdg]=15
				pos[endr]=16
			} else{
#				/* #cat '111U' i j k tijk; */
				brdg=which(od==2)
				sndr=which(od==1)
				endr=which(od==0)
				pos[sndr]=17
				pos[brdg]=18
				pos[endr]=19
			}#111 loop
		} else{ #must be 030
			if (min(od)==0){#030T
				brdg=which(od==1)
				sndr=which(od==2)
				endr=which(od==0);
				pos[sndr]=20;
				pos[brdg]=21;
				pos[endr]=22;
				#cat(tijk); #cat(ijk); #cat('030T')
			}else{ #030C
				pos=c(23,23, 23)
			} 		
		} 	
	} #stijk==3 
	
	if(stijk==4) {
		m=(.5)*(sum(diag(tijk%*%tijk))) 
		od = apply(tijk, 1, sum) #row sum
		
		if(m==1){
			.in = apply(tijk, 2, sum) #col sum
			if(min(od)==1 & min(.in)==0){
#				/* #cat '120D' i j k tijk; */
				sndr=which(.in==0)
				rcvr=c(1:3)[!(1:3 %in% sndr)]
				pos[sndr]=26
				pos[rcvr[1]]=27
				pos[rcvr[2]]=27
			} else if(min(od)==0 & min(.in)==1){
#				/* #cat '120U' i j k tijk; */
				rcvr=which(od==0)
				sndr=c(1:3)[!(1:3 %in% rcvr)]
				pos[rcvr]=28
				pos[sndr[1]]=29
				pos[sndr[2]]=29
			}else{ #120c
				dif=od-.in;
				sndr=which(dif==1)
				brdg=which(dif==0)
				rcvr=which(dif==-1)
				pos[sndr]=30
				pos[brdg]=31
				pos[rcvr]=32 
			} 			
			# man=120
		}else{ #must be 201
			brdg=which(od==2)
			sndr=c(1:3)[!(1:3 %in% brdg)]
			pos[brdg]=25
			pos[sndr[1]]=24
			pos[sndr[2]]=24
		} 
	}#stijk==4 
	
	if(stijk==5){ #5 arcs in triad
		od = apply(tijk, 1, sum) #row sum
		.in = apply(tijk, 2, sum) #row sum
		dif=od-.in
		sndr=which(dif==1)
		brdg=which(dif==-1)
		endr=which(dif==0)
		pos[sndr]=33
		pos[brdg]=34
		pos[endr]=35
	} 
	if(stijk==6){ #must be complete triad
		pos=c(36, 36, 36)
	}		
	return(pos)
}#itpeval


# itpcen takes an igraph object, and then iteratively calls itpeval over
# every triad in the graph

triadcensus = function(g){
	if(length(V(g)) < 3){
		cat(paste("Network needs at least 3 nodes!"))
		itpm = matrix( NA, nrow = max( 1, length(V(g)$name) ), ncol = 36 )
		rownames(itpm) = V(g)$name
		colnames(itpm) = c("003", "012_S", "012_E", "012_I", "102_D", "102_I",
				"021D_S", "021D_E", "021U_S", "021U_E", 
				"021C_S", "021C_B", "021C_E", 
				"111D_S", "111D_B", "111D_E",
				"111U_S", "111U_B", "111U_E",
				"030T_S", "030T_B", "030T_E",
				"030C", "201_S", "201_B", 
				"120D_S", "120D_E", "120U_E", "120U_S",
				"120C_S", "120C_B", "120C_E",
				"210_S", "210_B", "210_E", "300")
		
		return(itpm)
	}
	adjl = get.adjlist(g, mode = "all")
	
	#adjl=choose(adjl=.,0,adjl); is equivalent to:
	adjl=ifelse(is.na(adjl),0,adjl)
	# this deals with missing values, in SAS ".", which in R is "NA" 
	
	# This matrix will hold the positions of each node in the graph
	# we set the number of rows in the matrix to the length of the 
	# adjlist, which should correspond to nodes in the graph
	itpm=matrix(0, nrow = length(adjl), ncol = 36)
	
	.nr=length(adjl)-2
	for(i in 1:(length(adjl)-2)){
		if(100%%i==0){
			cat(paste("Case: ", i, " of: ", .nr, "\n"))
		}
		icont=adjl[[i]] # nodes that i sends to
		icont = icont +1 #this corrects igraph's indexing of verteces from 0
		
		for(j in (i+1):(length(adjl)-1)){
			#cat(ijk)
			if(j %in% icont){  # if ego sends to alter
				for(k in (j+1):length(adjl)){
					ijk= c(i,j,k) # i||j||k in original SAS code
					#ijk = c(3,4,5) # for debugging
					#cat(ijk)
					#submat=inmat[ijk,ijk]; # original SAS code
					submat = get.adjacency(g)[ijk,ijk]
					# we return the adjacency matrix of the graph, and then 
					# subset the matrix by rows i,j,k and columns i,j,k.
					itp=itpeval(submat) #call the itpeval function defined above
					# what we've done is to create a subgraph consisting of 
					# each triad in the graph, and then return the triad
					# position of each node.
					#cat(ijk);##cat(itp);
					itpm[ijk[1],itp[1]]=itpm[ijk[1],itp[1]]+1;
					itpm[ijk[2],itp[2]]=itpm[ijk[2],itp[2]]+1;
					itpm[ijk[3],itp[3]]=itpm[ijk[3],itp[3]]+1;
					## This adds 1 to every row corresponding to ijk 
					# and every column corresponding to the triad 
					# position of each node (which we got by calling itpeval).
				}	
			}else{ #ego does not send to alter */
				jcont = adjl[[j]] # nodes that j sends/recieves?
				jcont = jcont +1 #this corrects igraph's indexing of verteces from 0 
				jbiger = ifelse(jcont>j, 1, 0)
				jcont = jcont*jbiger # this gives us only nodes w/ index >j
				#cat(jcont)
				icnt_gtj = ifelse(icont>j, 1, 0) #nodes i sends to that are greater than j
				icnt_gtj = icont*icnt_gtj
				jcont = union(icnt_gtj,jcont)
				jcont=setdiff(jcont,0)
				
				if(length(jcont)>0){ #if jcont is not empty
					#cat(c(i, j, jcont))
					for(k in 1:length(jcont)){
						ijk=c(i,j,jcont[k]) # this is causing an error
						#cat(ijk)
						
						submat = get.adjacency(g)[ijk,ijk]
						# we return the adjacency matrix of the graph, and then 
						# subset the matrix by rows i,j,k and columns i,j,k.
						
						itp=itpeval(submat) #call the itpeval function defined above
						# what we've done is to create a subgraph consisting of 
						# each triad in the graph, and then return the triad
						# position of each node.
						#cat( ijk);##cat(itp) 
						itpm[ijk[1],itp[1]]=itpm[ijk[1],itp[1]]+1;
						itpm[ijk[2],itp[2]]=itpm[ijk[2],itp[2]]+1;
						itpm[ijk[3],itp[3]]=itpm[ijk[3],itp[3]]+1;
						
					}
				}	
			} 
		} 		
	}
	n=nrow(itpm)
	tottri=((n-1)*(n-2))/2
#	itpm[,1]=tottri-itpm[,+];
	itpm[,1]=tottri-apply(itpm, 1,sum)
	
	#define the column names:
	colnames(itpm) = c("003", "012_S", "012_E", "012_I", "102_D", "102_I",
			"021D_S", "021D_E", "021U_S", "021U_E", 
			"021C_S", "021C_B", "021C_E", 
			"111D_S", "111D_B", "111D_E",
			"111U_S", "111U_B", "111U_E",
			"030T_S", "030T_B", "030T_E",
			"030C", "201_S", "201_B", 
			"120D_S", "120D_E", "120U_E", "120U_S",
			"120C_S", "120C_B", "120C_E",
			"210_S", "210_B", "210_E", "300")
	rownames(itpm) = V(g)$name
	return(itpm)
}

