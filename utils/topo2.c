// final revision 19/11/12
//----------------------------------------------------------------------
#define cmax_length 1000
#define species_length 10
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
int main(int argc, char *argv[])
{
	FILE *fp;
	
	// fp related
	char current_folder[cmax_length],file_path[cmax_length],word[cmax_length];
	
	// dummies
	int int_buffer;
	char buffer[cmax_length];
	char *fcheck;
	
	//
	int i,j,k,l;
	int atoms,bonds,angles;
	int *atomic_ID;
	double *x,*y,*z;
	char **species,**unique_species;		

	int flag,temp;
	int element,T1,T2,T3,dihedrals,init_dihedrals;
	int max_neighbors,*neighbors;						
	
	int *init_B1,*init_B2,*B1,*B2;
	int **bonds_array;
	
	int *A1,*A2,*A3;
	
	int *init_D1,*init_D2,*init_D3,*init_D4;
	int *D1,*D2,*D3,*D4;
	
	int **bonds_ID;
	int **bonds_ID_final,bonds_ID_counter;
	int *bond_ID_export;
	
	int **angles_ID;
	int **angles_ID_final,angles_ID_counter;
	int *angles_ID_export;
	
	int **dihedrals_ID;
	int **dihedrals_ID_final,dihedrals_ID_counter;
	int *dihedrals_ID_export;
	
	//------------------------------------------------------------------
	// print exe info
	if(argc==1)
	{
		printf(
		"\n*Reads mol2 files and generates topology files containing the full topological information (bonds, angles, dihedrals)\n\n"
		"./topo <filename>.mol2\n\n"
		);
		exit(-1);
	}
	
	// read the current directory
	fcheck=getcwd(current_folder,cmax_length);
	
	// import data from .mol2 file - the filename is the first exe argument (argv[1])
	sprintf(file_path,"%s/%s",current_folder,argv[1]);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
	
	for(i=0;i<2;++i){fcheck=fgets(buffer,cmax_length,fp);}	// ignore first and second line; read from the third
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms,&bonds);
	//for(i=0;i<4;++i){fcheck=fgets(buffer,cmax_length,fp);}	// ignore lines...
	while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
	
	// preallocations	
	x=(double*)malloc(atoms*sizeof(double));		// coordinates
	y=(double*)malloc(atoms*sizeof(double));
	z=(double*)malloc(atoms*sizeof(double));
	init_B1=(int*)malloc(bonds*sizeof(int));		// bonding info
	init_B2=(int*)malloc(bonds*sizeof(int));
	species=(char**)malloc(atoms*sizeof(char*));	// species
	for (i=0;i<atoms;++i){species[i]=(char*)malloc(species_length*sizeof(char));}
	
	// read species and coordinates
	for(i=0;i<atoms;++i)
	{
		fcheck=fgets(buffer,cmax_length,fp);
        // based on the atom name column of the mol2 file (second column in @<TRIPOS>ATOM)
		sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf",&int_buffer,species[i],&x[i],&y[i],&z[i]);
        // based on the atom type column of the mol2 file (sixth column in @<TRIPOS>ATOM)
        //sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s",&int_buffer,word,&x[i],&y[i],&z[i],species[i]);
	}
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	// read bonding info
	for(i=0;i<bonds;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d",&int_buffer,&init_B1[i],&init_B2[i]);}
	
	fclose(fp);
	
	// preallocate the arrays for the full bonding info
	B1=(int*)malloc(2*bonds*sizeof(int));
	B2=(int*)malloc(2*bonds*sizeof(int));
	// populate
	for(i=0;i<bonds;++i)
	{
		B1[i]=init_B1[i];
		B2[i]=init_B2[i];
		B1[i+bonds]=init_B2[i];
		B2[i+bonds]=init_B1[i];
	}
		
	// bubble sort with respect to the B1 array
	flag=1;
	while(flag==1)
	{
		flag=0;
		for(i=1;i<2*bonds;++i)
		{
			if(B1[i]<B1[i-1])
			{
				temp=B1[i-1];
				B1[i-1]=B1[i];
				B1[i]=temp;
				
				temp=B2[i-1];
				B2[i-1]=B2[i];
				B2[i]=temp;
				
				flag=1;
			}
		}	
	}
	
	// array to store the number of atoms bonded to every single atom
	neighbors=(int*)malloc(atoms*sizeof(int));
	// initialize
	for(i=0;i<atoms;++i){neighbors[i]=0;}
	// populate array
	for(i=0;i<2*bonds;++i){neighbors[B1[i]-1]=neighbors[B1[i]-1]+1;}
	// calculate the maximum number of bonded atoms - needed for the following preallocation
	max_neighbors=0;
	for(i=0;i<atoms;++i){if(neighbors[i]>=max_neighbors){max_neighbors=neighbors[i];}}
	
	// bonding info array preallocation and initialization
	bonds_array=(int**)malloc(atoms*sizeof(int*));
	for(i=0;i<atoms;++i){bonds_array[i]=(int*)malloc(max_neighbors*sizeof(int));}
	for(i=0;i<atoms;++i){for(j=0;j<max_neighbors;++j){bonds_array[i][j]=0;}}
	
	// populate bonds array -- the following routine works correctly provided that the B2 array is sorted with respect to the B1 array!
	k=-1;
	for(i=0;i<atoms;++i)
	{
		for(j=0;j<neighbors[i];++j)
		{
			k=k+1;
			bonds_array[i][j]=B2[k];
		}
	}
	
	// calculate the total number of angles
	angles=0;
	for(i=0;i<atoms;++i){angles=angles+(neighbors[i]*(neighbors[i]-1))/2;}
	
	A1=(int*)malloc(angles*sizeof(int));
	A2=(int*)malloc(angles*sizeof(int));
	A3=(int*)malloc(angles*sizeof(int));
	
	// populate angles array
	l=-1;
	for(i=0;i<atoms;++i)
	{
		// the combinations double loop
		for(j=0;j<neighbors[i]-1;++j)
		{
			for(k=j+1;k<neighbors[i];++k)
			{
				l=l+1;
				A1[l]=bonds_array[i][j];
				A2[l]=i+1;
				A3[l]=bonds_array[i][k];
			}
		}
	}

	// calculate the initial number of dihedrals examining the branching info
	init_dihedrals=0;
	for(i=0;i<angles;++i)
	{
		// check the branching starting from the left side
		element=A1[i];
		T1=A2[i];
		T2=A3[i];
		for(j=0;j<neighbors[element-1];++j)
		{
			T3=bonds_array[element-1][j];
			if(T3!=T1 && T3!=T2)
			{
				init_dihedrals=init_dihedrals+1;
			}
		}
		// check the branching starting from the right side
		element=A3[i];
		T1=A1[i];
		T2=A2[i];
		for(j=0;j<neighbors[element-1];++j)
		{
			T3=bonds_array[element-1][j];
			if(T3!=T1 && T3!=T2)
			{
				init_dihedrals=init_dihedrals+1;
			}
		}
	}
	
	// preallocate the array for all dihedrals
	init_D1=(int*)malloc(init_dihedrals*sizeof(int));
	init_D2=(int*)malloc(init_dihedrals*sizeof(int));
	init_D3=(int*)malloc(init_dihedrals*sizeof(int));
	init_D4=(int*)malloc(init_dihedrals*sizeof(int));
	
	// populate
	k=0;
	for(i=0;i<angles;++i)
	{
		// left
		element=A1[i];
		T1=A2[i];
		T2=A3[i];
		for(j=0;j<neighbors[element-1];++j)
		{
			T3=bonds_array[element-1][j];
			if(T3!=T1 && T3!=T2)
			{
				if(T3<T2)	// (1)<(4)
				{
					init_D1[k]=T3;
					init_D2[k]=element;
					init_D3[k]=T1;
					init_D4[k]=T2;
				}
				else
				{
					init_D1[k]=T2;
					init_D2[k]=T1;
					init_D3[k]=element;
					init_D4[k]=T3;
				}
				k=k+1;
			}
		}
		// right
		element=A3[i];
		T1=A1[i];
		T2=A2[i];
		for(j=0;j<neighbors[element-1];++j)
		{
			T3=bonds_array[element-1][j];
			if(T3!=T1 && T3!=T2)
			{
				if(T1<T3)	// (1)<(4)
				{
					init_D1[k]=T1;
					init_D2[k]=T2;
					init_D3[k]=element;
					init_D4[k]=T3;
				}
				else
				{
					init_D1[k]=T3;
					init_D2[k]=element;
					init_D3[k]=T2;
					init_D4[k]=T1;
				}
				k=k+1;
			}
		}
	}

	// erase multiple entries
	for(i=0;i<init_dihedrals;++i)
	{
		for(j=0;j<i-1;++j)
		{
			if(init_D1[j]==init_D1[i] && init_D2[j]==init_D2[i] && init_D3[j]==init_D3[i] && init_D4[j]==init_D4[i])
			{
				init_D1[j]=0;
				init_D2[j]=0;
				init_D3[j]=0;
				init_D4[j]=0;
			}
		}
		for(j=i+1;j<init_dihedrals;++j)
		{
			if(init_D1[j]==init_D1[i] && init_D2[j]==init_D2[i] && init_D3[j]==init_D3[i] && init_D4[j]==init_D4[i])
			{
				init_D1[j]=0;
				init_D2[j]=0;
				init_D3[j]=0;
				init_D4[j]=0;
			}
		}
	}

	// bubble sort
	flag=1;
	while(flag==1)
	{
		flag=0;
		for(i=1;i<init_dihedrals;++i)
		{
			if(init_D1[i]<init_D1[i-1])
			{
				temp=init_D1[i-1];init_D1[i-1]=init_D1[i];init_D1[i]=temp;
				temp=init_D2[i-1];init_D2[i-1]=init_D2[i];init_D2[i]=temp;
				temp=init_D3[i-1];init_D3[i-1]=init_D3[i];init_D3[i]=temp;
				temp=init_D4[i-1];init_D4[i-1]=init_D4[i];init_D4[i]=temp;
				
				flag=1;
			}
		}	
	}

	// count the final number of dihedrals based on non-zero entries
	dihedrals=0;
	for(i=0;i<init_dihedrals;++i){if(init_D1[i]!=0){dihedrals=dihedrals+1;}}
	
	// preallocate and populate final dihedrals array
	D1=(int*)malloc(dihedrals*sizeof(int));
	D2=(int*)malloc(dihedrals*sizeof(int));
	D3=(int*)malloc(dihedrals*sizeof(int));
	D4=(int*)malloc(dihedrals*sizeof(int));
	
	k=-1;
	for(i=0;i<init_dihedrals;++i)
	{
		if(init_D1[i]!=0)
		{
			k=k+1;
			D1[k]=init_D1[i];
			D2[k]=init_D2[i];
			D3[k]=init_D3[i];
			D4[k]=init_D4[i];
		}
	}

	// output
	printf("Number of bonds: %d\n",bonds);
	printf("Number of angles: %d\n",angles);
	printf("Number of dihedrals: %d\n",dihedrals);
	
	//------------------------------------------------------------------
	
	// assign an integer counter starting from 1 for every different species present in the dat file
	atomic_ID=(int*)malloc(atoms*sizeof(int));			// the array for the new species representation based on integers
	unique_species=(char**)malloc(atoms*sizeof(char*));	// the array to store the unique species
	for (i=0;i<atoms;++i){unique_species[i]=(char*)malloc(species_length*sizeof(char));}
	
	// identify unique secies
	sprintf(unique_species[0],"%s",species[0]);
	k=1;
	// check the rest species
	for(i=1;i<atoms;++i)
	{
		flag=0;
		for(j=0;j<k;++j)
		{
			if(strcmp(species[i],unique_species[j])==0){flag=1;}
		}
		if(flag==0)
		{
			sprintf(unique_species[k],"%s",species[i]);
			k=k+1;
		}
	}
	// output
	printf("Found %d different species:\n",k);
	for(i=0;i<k;++i){printf("%s\t%d\n",unique_species[i],i+1);}
	
	// populate the atomic_ID array
	for(i=0;i<atoms;++i)
	{
		for(j=0;j<k;++j)
		{
			if(strcmp(species[i],unique_species[j])==0)
			{
				atomic_ID[i]=j+1;
			}
		}
	}
	
	//------------------------------------------------------------------
	// write topo.species
	
	sprintf(file_path,"%s/topo.species",current_folder);
	fp=fopen(file_path,"w+");
	
	for(i=0;i<k;++i){fprintf(fp,"%s\n",unique_species[i]);}
	
	fclose(fp);
	printf("\nGenerated %s\n",file_path);
	
	//------------------------------------------------------------------
	// Bonds
	
	bonds_ID=(int**)malloc(bonds*sizeof(int*));		// 2D array to store the species counters
	for(i=0;i<bonds;++i){bonds_ID[i]=(int*)malloc(2*sizeof(int));}
	// populate
	for(i=0;i<bonds;++i)
	{
		bonds_ID[i][0]=atomic_ID[init_B1[i]-1];
		bonds_ID[i][1]=atomic_ID[init_B2[i]-1];
	}
	
	// erase multiple entries
	for(i=0;i<bonds;++i)
	{
		for(j=0;j<i-1;++j)
		{
			if(bonds_ID[j][0]==bonds_ID[i][0]&&bonds_ID[j][1]==bonds_ID[i][1])
			{
				bonds_ID[j][0]=0;
				bonds_ID[j][1]=0;
			}
			else if(bonds_ID[j][1]==bonds_ID[i][0]&&bonds_ID[j][0]==bonds_ID[i][1])
			{
				bonds_ID[j][0]=0;
				bonds_ID[j][1]=0;
			}
		}
		for(j=i+1;j<bonds;++j)
		{
			if(bonds_ID[j][0]==bonds_ID[i][0]&&bonds_ID[j][1]==bonds_ID[i][1])
			{
				bonds_ID[j][0]=0;
				bonds_ID[j][1]=0;
			}
			else if(bonds_ID[j][1]==bonds_ID[i][0]&&bonds_ID[j][0]==bonds_ID[i][1])
			{
				bonds_ID[j][0]=0;
				bonds_ID[j][1]=0;
			}
		}
	}
	
	// calculate number of unique bonds	
	bonds_ID_counter=0;
	for(i=0;i<bonds;++i)
	{
		if(bonds_ID[i][0]!=0)
		{
			bonds_ID_counter=bonds_ID_counter+1;
		}
	}
	
	// save the unique bonding info
	bonds_ID_final=(int**)malloc(bonds_ID_counter*sizeof(int*));	// preallocate
	for(i=0;i<bonds_ID_counter;++i){bonds_ID_final[i]=(int*)malloc(2*sizeof(int));}
	// populate
	j=-1;
	for(i=0;i<bonds;++i)
	{
		if(bonds_ID[i][0]!=0)
		{
			j=j+1;
			bonds_ID_final[j][0]=bonds_ID[i][0];
			bonds_ID_final[j][1]=bonds_ID[i][1];
		}
	}
	// bond type counter
	bond_ID_export=(int*)malloc(bonds*sizeof(int));
	// populate based on comparison
	for(i=0;i<bonds;++i)
	{
		for(j=0;j<bonds_ID_counter;++j)
		{
			if(atomic_ID[init_B1[i]-1]==bonds_ID_final[j][0] && atomic_ID[init_B2[i]-1]==bonds_ID_final[j][1])
			{
				bond_ID_export[i]=j+1;
			}
			else if(atomic_ID[init_B1[i]-1]==bonds_ID_final[j][1] && atomic_ID[init_B2[i]-1]==bonds_ID_final[j][0])
			{
				bond_ID_export[i]=j+1;
			}
		}		
	}
	// output
	printf("Bonds\n%d types\n",bonds_ID_counter);
	for(i=0;i<bonds;++i)
	{
		if(bonds_ID[i][0]!=0)
		{
			printf("[%d]\t%s\t%s\n",bond_ID_export[i],unique_species[bonds_ID[i][0]-1],unique_species[bonds_ID[i][1]-1]);
		}
	}
	
	//------------------------------------------------------------------
	// Angles
		
	angles_ID=(int**)malloc(angles*sizeof(int*));	// array to store the species ID
	for(i=0;i<angles;++i){angles_ID[i]=(int*)malloc(3*sizeof(int));}
	// populate
	for(i=0;i<angles;++i)
	{
		angles_ID[i][0]=atomic_ID[A1[i]-1];
		angles_ID[i][1]=atomic_ID[A2[i]-1];
		angles_ID[i][2]=atomic_ID[A3[i]-1];
	}
	// erase multiple entries
	for(i=0;i<angles;++i)
	{
		for(j=0;j<i-1;++j)
		{
			if(angles_ID[j][0]==angles_ID[i][0]&&angles_ID[j][1]==angles_ID[i][1]&&angles_ID[j][2]==angles_ID[i][2])
			{
				angles_ID[j][0]=0;
				angles_ID[j][1]=0;
				angles_ID[j][2]=0;
			}
			else if(angles_ID[j][2]==angles_ID[i][0]&&angles_ID[j][1]==angles_ID[i][1]&&angles_ID[j][0]==angles_ID[i][2])
			{
				angles_ID[j][0]=0;
				angles_ID[j][1]=0;
				angles_ID[j][2]=0;
			}
		}
		for(j=i+1;j<angles;++j)
		{
			if(angles_ID[j][0]==angles_ID[i][0]&&angles_ID[j][1]==angles_ID[i][1]&&angles_ID[j][2]==angles_ID[i][2])
			{
				angles_ID[j][0]=0;
				angles_ID[j][1]=0;
				angles_ID[j][2]=0;
			}
			else if(angles_ID[j][2]==angles_ID[i][0]&&angles_ID[j][1]==angles_ID[i][1]&&angles_ID[j][0]==angles_ID[i][2])
			{
				angles_ID[j][0]=0;
				angles_ID[j][1]=0;
				angles_ID[j][2]=0;
			}
		}
	}
	// unique angles counter
	angles_ID_counter=0;
	for(i=0;i<angles;++i)
	{
		if(angles_ID[i][0]!=0)
		{
			angles_ID_counter=angles_ID_counter+1;
		}
	}
	// store unique angles info
	angles_ID_final=(int**)malloc(angles_ID_counter*sizeof(int*));	// preallocate
	for(i=0;i<angles_ID_counter;++i){angles_ID_final[i]=(int*)malloc(3*sizeof(int));}
	// populate
	j=-1;
	for(i=0;i<angles;++i)
	{
		if(angles_ID[i][0]!=0)
		{
			j=j+1;
			angles_ID_final[j][0]=angles_ID[i][0];
			angles_ID_final[j][1]=angles_ID[i][1];
			angles_ID_final[j][2]=angles_ID[i][2];
		}
	}
	// angle type counter
	angles_ID_export=(int*)malloc(angles*sizeof(int));
	// populate
	for(i=0;i<angles;++i)
	{
		for(j=0;j<angles_ID_counter;++j)
		{
			if(atomic_ID[A1[i]-1]==angles_ID_final[j][0] && atomic_ID[A2[i]-1]==angles_ID_final[j][1] && atomic_ID[A3[i]-1]==angles_ID_final[j][2])
			{
				angles_ID_export[i]=j+1;
			}
			else if(atomic_ID[A1[i]-1]==angles_ID_final[j][2] && atomic_ID[A2[i]-1]==angles_ID_final[j][1] && atomic_ID[A3[i]-1]==angles_ID_final[j][0])
			{
				angles_ID_export[i]=j+1;
			}
		}
			
	}
	// output
	printf("Angles\n%d types\n",angles_ID_counter);
	for(i=0;i<angles;++i)
	{
		if(angles_ID[i][0]!=0)
		{
			printf("[%d]\t%s\t%s\t%s\n",angles_ID_export[i],unique_species[angles_ID[i][0]-1],unique_species[angles_ID[i][1]-1],unique_species[angles_ID[i][2]-1]);
		}
	}
	
	//------------------------------------------------------------------
	// Dihedrals
		
	dihedrals_ID=(int**)malloc(dihedrals*sizeof(int*));	// array to store the species ID
	for(i=0;i<dihedrals;++i){dihedrals_ID[i]=(int*)malloc(4*sizeof(int));}
	// populate
	for(i=0;i<dihedrals;++i)
	{
		dihedrals_ID[i][0]=atomic_ID[D1[i]-1];
		dihedrals_ID[i][1]=atomic_ID[D2[i]-1];
		dihedrals_ID[i][2]=atomic_ID[D3[i]-1];
		dihedrals_ID[i][3]=atomic_ID[D4[i]-1];
	}
	// erase multiple entries
	for(i=0;i<dihedrals;++i)
	{
		for(j=0;j<i-1;++j)
		{
			if(dihedrals_ID[j][0]==dihedrals_ID[i][0]&&dihedrals_ID[j][1]==dihedrals_ID[i][1]&&dihedrals_ID[j][2]==dihedrals_ID[i][2]&&dihedrals_ID[j][3]==dihedrals_ID[i][3])
			{
				dihedrals_ID[j][0]=0;
				dihedrals_ID[j][1]=0;
				dihedrals_ID[j][2]=0;
				dihedrals_ID[j][3]=0;
			}
			
			else if(dihedrals_ID[j][3]==dihedrals_ID[i][0]&&dihedrals_ID[j][2]==dihedrals_ID[i][1]&&dihedrals_ID[j][1]==dihedrals_ID[i][2]&&dihedrals_ID[j][0]==dihedrals_ID[i][3])
			{
				dihedrals_ID[j][0]=0;
				dihedrals_ID[j][1]=0;
				dihedrals_ID[j][2]=0;
				dihedrals_ID[j][3]=0;
			}
			
		}
		for(j=i+1;j<dihedrals;++j)
		{
			if(dihedrals_ID[j][0]==dihedrals_ID[i][0]&&dihedrals_ID[j][1]==dihedrals_ID[i][1]&&dihedrals_ID[j][2]==dihedrals_ID[i][2]&&dihedrals_ID[j][3]==dihedrals_ID[i][3])
			{
				dihedrals_ID[j][0]=0;
				dihedrals_ID[j][1]=0;
				dihedrals_ID[j][2]=0;
				dihedrals_ID[j][3]=0;
			}
			
			else if(dihedrals_ID[j][3]==dihedrals_ID[i][0]&&dihedrals_ID[j][2]==dihedrals_ID[i][1]&&dihedrals_ID[j][1]==dihedrals_ID[i][2]&&dihedrals_ID[j][0]==dihedrals_ID[i][3])
			{
				dihedrals_ID[j][0]=0;
				dihedrals_ID[j][1]=0;
				dihedrals_ID[j][2]=0;
				dihedrals_ID[j][3]=0;
			}
			
		}
	}
	// unique dihedrals counter
	dihedrals_ID_counter=0;
	for(i=0;i<dihedrals;++i)
	{
		if(dihedrals_ID[i][0]!=0)
		{
			dihedrals_ID_counter=dihedrals_ID_counter+1;
		}
	}
	// store unique dihedrals info
	dihedrals_ID_final=(int**)malloc(dihedrals_ID_counter*sizeof(int*));	// preallocate
	for(i=0;i<dihedrals_ID_counter;++i){dihedrals_ID_final[i]=(int*)malloc(4*sizeof(int));}
	// populate
	j=-1;
	for(i=0;i<dihedrals;++i)
	{
		if(dihedrals_ID[i][0]!=0)
		{
			j=j+1;
			dihedrals_ID_final[j][0]=dihedrals_ID[i][0];
			dihedrals_ID_final[j][1]=dihedrals_ID[i][1];
			dihedrals_ID_final[j][2]=dihedrals_ID[i][2];
			dihedrals_ID_final[j][3]=dihedrals_ID[i][3];
		}
	}
	// dihedrals type counter
	dihedrals_ID_export=(int*)malloc(dihedrals*sizeof(int));
	// populate
	for(i=0;i<dihedrals;++i)
	{
		for(j=0;j<dihedrals_ID_counter;++j)
		{
			if(atomic_ID[D1[i]-1]==dihedrals_ID_final[j][0] && atomic_ID[D2[i]-1]==dihedrals_ID_final[j][1] && atomic_ID[D3[i]-1]==dihedrals_ID_final[j][2] && atomic_ID[D4[i]-1]==dihedrals_ID_final[j][3])
			{
				dihedrals_ID_export[i]=j+1;
			}
			if(atomic_ID[D1[i]-1]==dihedrals_ID_final[j][3] && atomic_ID[D2[i]-1]==dihedrals_ID_final[j][2] && atomic_ID[D3[i]-1]==dihedrals_ID_final[j][1] && atomic_ID[D4[i]-1]==dihedrals_ID_final[j][0])
			{
				dihedrals_ID_export[i]=j+1;
			}
		}
			
	}
	// output
	printf("Dihedrals\n%d types\n",dihedrals_ID_counter);
	for(i=0;i<dihedrals;++i)
	{
		if(dihedrals_ID[i][0]!=0)
		{
			printf("[%d]\t%s\t%s\t%s\t%s\n",dihedrals_ID_export[i],unique_species[dihedrals_ID[i][0]-1],unique_species[dihedrals_ID[i][1]-1],unique_species[dihedrals_ID[i][2]-1],unique_species[dihedrals_ID[i][3]-1]);
		}
	}
	
	//------------------------------------------------------------------
	// write topo.log
	sprintf(file_path,"%s/topo.log",current_folder);
	fp=fopen(file_path,"w+");
	
	fprintf(fp,"Bonds\n%d(%d) types\n",bonds_ID_counter,bonds_ID_counter);
	for(i=0;i<bonds;++i)
	{
		if(bonds_ID[i][0]!=0)
		{
			fprintf(fp,"[%d]\t%s\t%s\n",bond_ID_export[i],unique_species[bonds_ID[i][0]-1],unique_species[bonds_ID[i][1]-1]);
		}
	}
	fprintf(fp,"Angles\n%d(%d) types\n",angles_ID_counter,angles_ID_counter);
	for(i=0;i<angles;++i)
	{
		if(angles_ID[i][0]!=0)
		{
			fprintf(fp,"[%d]\t%s\t%s\t%s\n",angles_ID_export[i],unique_species[angles_ID[i][0]-1],unique_species[angles_ID[i][1]-1],unique_species[angles_ID[i][2]-1]);
		}
	}
	fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals_ID_counter,dihedrals_ID_counter);
	for(i=0;i<dihedrals;++i)
	{
		if(dihedrals_ID[i][0]!=0)
		{
			fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",dihedrals_ID_export[i],unique_species[dihedrals_ID[i][0]-1],unique_species[dihedrals_ID[i][1]-1],unique_species[dihedrals_ID[i][2]-1],unique_species[dihedrals_ID[i][3]-1]);
		}
	}
	
	fclose(fp);
	printf("\nGenerated %s\n",file_path);
	
	//------------------------------------------------------------------
	// write topo.out
	sprintf(file_path,"%s/topo.out",current_folder);
	fp=fopen(file_path,"w+");
	
	fprintf(fp,"%d bonds\n%d angles\n%d dihedrals\n\n",bonds,angles,dihedrals);
	fprintf(fp,"%d bond types\n%d angle types\n%d dihedral types\n\n",bonds_ID_counter,angles_ID_counter,dihedrals_ID_counter);
	fprintf(fp,"Bonds\n\n");
	for(i=0;i<bonds;++i){fprintf(fp,"%d\t%d\t%d\t%d\t--\t%s-%s\n",i+1,bond_ID_export[i],init_B1[i],init_B2[i],species[init_B1[i]-1],species[init_B2[i]-1]);}
	fprintf(fp,"\nAngles\n\n");
	for(i=0;i<angles;++i){fprintf(fp,"%d\t%d\t%d\t%d\t%d\t--\t%s-%s-%s\n",i+1,angles_ID_export[i],A1[i],A2[i],A3[i],species[A1[i]-1],species[A2[i]-1],species[A3[i]-1]);}
	fprintf(fp,"\nDihedrals\n\n");
	for(i=0;i<dihedrals;++i){fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t--\t%s-%s-%s-%s\n",i+1,dihedrals_ID_export[i],D1[i],D2[i],D3[i],D4[i],species[D1[i]-1],species[D2[i]-1],species[D3[i]-1],species[D4[i]-1]);}
	
	fclose(fp);
	printf("\nGenerated %s\n",file_path);
	
	//------------------------------------------------------------------
	// write topo.neigh (for the impropers identification)
	
	sprintf(file_path,"%s/topo.neigh",current_folder);
	fp=fopen(file_path,"w+");
	
	fprintf(fp,"%d atoms\n",atoms);
	for(i=0;i<atoms;++i){fprintf(fp,"%d\t%lf\t%lf\t%lf\n",i+1,x[i],y[i],z[i]);}
	for(i=0;i<atoms;++i)
	{
		fprintf(fp,"%d\t%s\t",i+1,unique_species[atomic_ID[i]-1]);
		
		for(j=0;j<max_neighbors;++j)
		{
			if(bonds_array[i][j]!=0){fprintf(fp,"%d\t%s\t",bonds_array[i][j],unique_species[atomic_ID[bonds_array[i][j]-1]-1]);}
		}
		fprintf(fp,"\n");
		
	}
	
	fclose(fp);
	printf("\nGenerated %s\n",file_path);
	
	//------------------------------------------------------------------
	// free memory
	
	for (i=0;i<dihedrals_ID_counter;++i){free(dihedrals_ID_final[i]);}free(dihedrals_ID_final);
	for (i=0;i<dihedrals;++i){free(dihedrals_ID[i]);}free(dihedrals_ID);
	for (i=0;i<angles_ID_counter;++i){free(angles_ID_final[i]);}free(angles_ID_final);
	for (i=0;i<angles;++i){free(angles_ID[i]);}free(angles_ID);
	for (i=0;i<bonds_ID_counter;++i){free(bonds_ID_final[i]);}free(bonds_ID_final);
	for (i=0;i<bonds;++i){free(bonds_ID[i]);}free(bonds_ID);
	for (i=0;i<atoms;++i){free(bonds_array[i]);}free(bonds_array);
	for (i=0;i<atoms;++i){free(species[i]);}free(species);
	for (i=0;i<atoms;++i){free(unique_species[i]);}free(unique_species);
	free(neighbors);
	free(bond_ID_export);free(angles_ID_export);free(dihedrals_ID_export);
	free(x);free(y);free(z);
	free(atomic_ID);
	free(init_B1);free(init_B2);free(B1);free(B2);
	free(A1);free(A2);free(A3);
	free(init_D1);free(init_D2);free(init_D3);free(init_D4);free(D1);free(D2);free(D3);free(D4);
	return 0;
}

