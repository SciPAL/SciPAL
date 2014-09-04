mergeRes_list=[]
mergeRes_list.append("\n\n#\tn_rows\tCPU_factorization_pure\tGPU_fac\tGPU_fac_incl_data_transfer\tdiag_update\tfactorize_diagonal_block\tlo_update\tstrip_update\tDevice\tshared_memory\tdouble_precision\n")
for DEVICE in range(0,2):
	for SHARED_MEMORY in range(0,2):
		for DOUBLE_PRECISION in range(0,2):
			erg = open("chol_fac_times_"+str(DEVICE)+"_"+str(SHARED_MEMORY)+"_"+str(DOUBLE_PRECISION)+".dat")
			erg_list = erg.readlines()
			erg.close()
			i=0
			for line in erg_list:
				if(str.find(line,"n rows CPU factorization") == -1):
					mergeRes_list.append(erg_list[i][0:len(erg_list[i])-1] + "\t"+str(DEVICE)+"\t"+str(SHARED_MEMORY)+"\t"+str(DOUBLE_PRECISION)+"\n")
				i+=1
			mergeRes_list.append("\n\n")
				
mergeRes = open("mergedResults.dat","w");
mergeRes.writelines(mergeRes_list)
mergeRes.close()
