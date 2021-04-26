#define cmax_length     1000
#define sub_length      10
#define equal_diff      1e-8
#define rcut            12.0
#define LJ14            0.5
#define MASS_ID         "MASS"
#define BONDS_ID        "BOND"
#define ANGLES_ID       "ANGLE"
#define DIHEDRALS_ID    "DIHE"
#define LJ_ID           "NONB"
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<string.h>
// tokenize functions
void tokenize2(char *buffer,int *elements_out,int *depth,char ***res);
void clear_tok(int elements,char ***res);
int resolve_elements(FILE *fp,char *match);
int resolve_elements(FILE *fp,char *match)
{
    char buffer[cmax_length],word[cmax_length];
    int lines=0,element=0,start=0,stop=0,found=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        ++lines;sprintf(word,"%s","");sscanf(buffer,"%s",word);
        if(strcmp(word,match)==0){start=lines;found=1;}
        if(strcmp(buffer,"\n")==0 && found==1){stop=lines;break;}
    }
    element=stop-start-1;
    rewind(fp);
    return element;
}
int equal(double A,double B);
int equal(double A,double B)
{
    int value;
    if(fabs(A-B)<equal_diff){value=1;}else{value=0;}
    return value;
}
void find_unique_d(int N, int M, double **array,int **out);
void find_unique_d(int N, int M, double **array,int **out)
{
    int i,j,k,sum;
    int *init,*final;
    
    init=(int*)malloc(N*sizeof(int));
    final=(int*)malloc(N*sizeof(int));
    
    for(i=0;i<N;++i)
    {
        init[i]=i+1;
        final[i]=i+1;
    }
    
    for(i=0;i<N-1;++i)
    {
        for(j=i+1;j<N;++j)
        {
            sum=0;
            for(k=0;k<M;++k)sum=sum+equal(array[i][k],array[j][k]);
            if(sum==M)
            {
                final[j]=-1;
            }
        }
    }
    
    j=0;
    for(i=0;i<N;++i)
    {
        if(final[i]!=-1)
        {
            j=j+1;
            final[i]=j;
        }
    }
    
    for(i=0;i<N;++i)
        out[i][2]=final[i];
    
    for(i=0;i<N;++i)
    {
        if(final[i]==-1)
        {
            for(j=0;j<N;++j)
            {
                if(final[j]!=-1)
                {
                    sum=0;
                    for(k=0;k<M;++k)sum=sum+equal(array[i][k],array[j][k]);
                    if(sum==M)
                    {
                        final[i]=final[j];
                        break;
                    }
                }
            }
        }
    }
    
    for(i=0;i<N;++i)
    {
        out[i][0]=init[i];
        out[i][1]=final[i];
    }
    
    free(init);free(final);
    
}
int main(int argc,char **argv)
{
    FILE *fp;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],word[cmax_length];
    int atom_types,bond_types,angle_types,dihedral_types,i,j,int_buffer;
    char **B1,**B2,**A1,**A2,**A3,**D1,**D2,**D3,**D4;
    double *r0,*Kb,*theta0,*Ka,*Kd,*phi0,*nd;
    double mass;
    
    double **fu_in;
    int **fu_out,cols,rows,alterB=0,alterA=0,alterD=0;
    
    char lstring[cmax_length],rstring[cmax_length];
    int divide,k;
    
    getcwd(current_folder,cmax_length);
    
    if(argc==2)
    {
        
        sprintf(file_path,"%s/%s",current_folder,argv[1]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}

        /*
        
         Force field parameters from prmtop
         MASS
         ca    12.010
         ha     1.008
         
         BOND
         ca-ca    461.100   1.398
         ca-ha    345.800   1.086
         
         ANGLE
         ca-ca-ca     66.600  120.020
         ha-ca-ca     48.200  119.880
         
         DIHE
         ca-ca-ca-ca    1     3.62500000  180.000   2.0    SCEE=1.2 SCNB=2.0
         ca-ca-ca-ha    1     3.62500000  180.000   2.0    SCEE=1.2 SCNB=2.0
         ha-ca-ca-ha    1     3.62500000  180.000   2.0    SCEE=1.2 SCNB=2.0
         
         IMPROPER
         ca-ca-ca-ca     1.10000000  180.000   2.0
         ca-ca-ca-ha     1.10000000  180.000   2.0
         
         NONB
         ca    1.90800000   0.08600000
         ha    1.45900000   0.01500000
         
         
        */
        
        // types
        atom_types=resolve_elements(fp,MASS_ID);//printf("atom types %d\n",atom_types);
        bond_types=resolve_elements(fp,BONDS_ID);//printf("bond types %d\n",bond_types);
        angle_types=resolve_elements(fp,ANGLES_ID);//printf("angle types %d\n",angle_types);
        dihedral_types=resolve_elements(fp,DIHEDRALS_ID);//printf("dihedral types %d\n",dihedral_types);
        
        // masses
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,MASS_ID)==0)break;}
        for(i=0;i<atom_types;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf",word,&mass);printf("mass\t%s\t%lf\n",word,mass);}
        
        rewind(fp);
        
        // bonds
        divide=5;
        rows=bond_types;
        cols=2;
        B1=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)B1[i]=(char*)malloc(sub_length*sizeof(char));
        B2=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)B2[i]=(char*)malloc(sub_length*sizeof(char));
        r0=(double*)malloc(rows*sizeof(double));
        Kb=(double*)malloc(rows*sizeof(double));
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,BONDS_ID)==0)break;}
        for(i=0;i<rows;++i)
        {
            /*
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf",word,&Kb[i],&r0[i]);
            for(j=0;j<cmax_length;++j){if(word[j]=='\0')break;if(word[j]=='-')word[j]='\t';}
            sscanf(word,"%s\t%s",B1[i],B2[i]);
            */
            fgets(buffer,cmax_length,fp);
            k=-1;
            for(j=0;j<divide;++j)
            {
                k=k+1;
                lstring[k]=buffer[j];
            }
            lstring[k+1]='\0';
            k=-1;
            for(j=divide;j<cmax_length;++j)
            {
                if(buffer[j]=='\0')break;
                k=k+1;
                rstring[k]=buffer[j];
            }
            for(j=0;j<cmax_length;++j){if(lstring[j]=='\0')break;if(lstring[j]=='-')lstring[j]='\t';}
            sscanf(lstring,"%s\t%s",B1[i],B2[i]);
            sscanf(rstring,"%lf\t%lf",&Kb[i],&r0[i]);
            
        }
        
        fu_in=(double**)malloc(rows*sizeof(double*));for(i=0;i<rows;++i)fu_in[i]=(double*)malloc(cols*sizeof(double));
        fu_out=(int**)malloc(rows*sizeof(int*));for(i=0;i<rows;++i)fu_out[i]=(int*)malloc(3*sizeof(double)); // always 3!
        //
        for(i=0;i<rows;++i){
            fu_in[i][0]=Kb[i];
            fu_in[i][1]=r0[i];
        }
        //
        find_unique_d(rows,cols,fu_in,fu_out);
        printf("bond_style\t%s\n","harmonic");
        for(i=0;i<rows;++i)
        {
            printf("bond\t%d\t%s\t%s\t%s\t%lf\t%lf\n",fu_out[i][1],B1[i],B2[i],"harmonic",fu_in[i][1],fu_in[i][0]);
        }
        for(i=0;i<rows;++i)if(fu_out[i][2]==-1){alterB=1;break;}
        for(i=0;i<rows;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        
        // angles
        divide=8;
        rows=angle_types;
        cols=2;
        A1=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)A1[i]=(char*)malloc(sub_length*sizeof(char));
        A2=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)A2[i]=(char*)malloc(sub_length*sizeof(char));
        A3=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)A3[i]=(char*)malloc(sub_length*sizeof(char));
        theta0=(double*)malloc(rows*sizeof(double));
        Ka=(double*)malloc(rows*sizeof(double));
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,ANGLES_ID)==0)break;}
        for(i=0;i<rows;++i)
        {
            /*
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf",word,&Ka[i],&theta0[i]);
            for(j=0;j<cmax_length;++j){if(word[j]=='\0')break;if(word[j]=='-')word[j]='\t';}
            sscanf(word,"%s\t%s\t%s",A1[i],A2[i],A3[i]);
            */
            fgets(buffer,cmax_length,fp);
            k=-1;
            for(j=0;j<divide;++j)
            {
                k=k+1;
                lstring[k]=buffer[j];
            }
            lstring[k+1]='\0';
            k=-1;
            for(j=divide;j<cmax_length;++j)
            {
                if(buffer[j]=='\0')break;
                k=k+1;
                rstring[k]=buffer[j];
            }
            for(j=0;j<cmax_length;++j){if(lstring[j]=='\0')break;if(lstring[j]=='-')lstring[j]='\t';}
            sscanf(lstring,"%s\t%s\t%s",A1[i],A2[i],A3[i]);
            sscanf(rstring,"%lf\t%lf",&Ka[i],&theta0[i]);
        }
        
        fu_in=(double**)malloc(rows*sizeof(double*));for(i=0;i<rows;++i)fu_in[i]=(double*)malloc(cols*sizeof(double));
        fu_out=(int**)malloc(rows*sizeof(int*));for(i=0;i<rows;++i)fu_out[i]=(int*)malloc(3*sizeof(double)); // always 3!
        //
        for(i=0;i<rows;++i){
            fu_in[i][0]=Ka[i];
            fu_in[i][1]=theta0[i];
        }
        //
        find_unique_d(rows,cols,fu_in,fu_out);
        printf("angle_style\t%s\n","harmonic");
        for(i=0;i<rows;++i)
        {
            printf("angle\t%d\t%s\t%s\t%s\t%s\t%lf\t%lf\n",fu_out[i][1],A1[i],A2[i],A3[i],"harmonic",fu_in[i][1],fu_in[i][0]);
        }
        for(i=0;i<rows;++i)if(fu_out[i][2]==-1){alterA=1;break;}
        for(i=0;i<rows;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        
        //
        
        // dihedrals
        
        // read
        divide=11;
        rows=dihedral_types;
        D1=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D1[i]=(char*)malloc(sub_length*sizeof(char));
        D2=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D2[i]=(char*)malloc(sub_length*sizeof(char));
        D3=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D3[i]=(char*)malloc(sub_length*sizeof(char));
        D4=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D4[i]=(char*)malloc(sub_length*sizeof(char));
        phi0=(double*)malloc(rows*sizeof(double));
        Kd=(double*)malloc(rows*sizeof(double));
        nd=(double*)malloc(rows*sizeof(double));
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,DIHEDRALS_ID)==0)break;}
        for(i=0;i<rows;++i)
        {
            /*
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d\t%lf\t%lf\t%lf",word,&int_buffer,&Kd[i],&phi0[i],&nd[i]);
            for(j=0;j<cmax_length;++j){if(word[j]=='\0')break;if(word[j]=='-')word[j]='\t';}
            sscanf(word,"%s\t%s\t%s\t%s",D1[i],D2[i],D3[i],D4[i]);
            */
            fgets(buffer,cmax_length,fp);
            k=-1;
            for(j=0;j<divide;++j)
            {
                k=k+1;
                lstring[k]=buffer[j];
            }
            lstring[k+1]='\0';
            k=-1;
            for(j=divide;j<cmax_length;++j)
            {
                if(buffer[j]=='\0')break;
                k=k+1;
                rstring[k]=buffer[j];
            }
            for(j=0;j<cmax_length;++j){if(lstring[j]=='\0')break;if(lstring[j]=='-')lstring[j]='\t';}
            sscanf(lstring,"%s\t%s\t%s\t%s",D1[i],D2[i],D3[i],D4[i]);
            sscanf(rstring,"%d\t%lf\t%lf\t%lf",&int_buffer,&Kd[i],&phi0[i],&nd[i]);
        }
    //    for(i=0;i<dihedral_types;++i)printf("%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\n",D1[i],D2[i],D3[i],D4[i],Kd[i],phi0[i],nd[i]);
        
        char **D1b,**D2b,**D3b,**D4b;
        D1b=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D1b[i]=(char*)malloc(sub_length*sizeof(char));
        D2b=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D2b[i]=(char*)malloc(sub_length*sizeof(char));
        D3b=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D3b[i]=(char*)malloc(sub_length*sizeof(char));
        D4b=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D4b[i]=(char*)malloc(sub_length*sizeof(char));
        for(i=0;i<dihedral_types;++i){sprintf(D1b[i],"%s",D1[i]);sprintf(D2b[i],"%s",D2[i]);sprintf(D3b[i],"%s",D3[i]);sprintf(D4b[i],"%s",D4[i]);}
        for(i=0;i<dihedral_types;++i)
        {
            for(j=i+1;j<dihedral_types;++j)
            {
                if((strcmp(D1b[i],D1b[j])==0&&strcmp(D2b[i],D2b[j])==0&&strcmp(D3b[i],D3b[j])==0&&strcmp(D4b[i],D4b[j])==0)
                   ||
                   (strcmp(D1b[i],D4b[j])==0&&strcmp(D2b[i],D3b[j])==0&&strcmp(D3b[i],D2b[j])==0&&strcmp(D4b[i],D1b[j])==0))
                {
                    sprintf(D1b[j],"%d",0);
                    sprintf(D2b[j],"%d",0);
                    sprintf(D3b[j],"%d",0);
                    sprintf(D4b[j],"%d",0);
                }
            }
        }
        //for(i=0;i<dihedral_types;++i)printf("%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\n",D1b[i],D2b[i],D3b[i],D4b[i],Kd[i],phi0[i],nd[i]);
        
        int unique_D=dihedral_types;
        for(i=0;i<dihedral_types;++i)if(strcmp(D1b[i],"0")==0)unique_D=unique_D-1;
        //printf("%d\n",unique_D);
        int max_occur=1,occur;
        for(i=0;i<dihedral_types;++i)
        {
            if(strcmp(D1b[i],"0")!=0)
            {
                occur=0;
                for(j=0;j<dihedral_types;++j)
                {
                    if((strcmp(D1b[i],D1[j])==0&&strcmp(D2b[i],D2[j])==0&&strcmp(D3b[i],D3[j])==0&&strcmp(D4b[i],D4[j])==0)
                       ||
                       (strcmp(D1b[i],D4[j])==0&&strcmp(D2b[i],D3[j])==0&&strcmp(D3b[i],D2[j])==0&&strcmp(D4b[i],D1[j])==0))
                    {
                        ++occur;
                    }
                }
                if(occur>max_occur)max_occur=occur;
            }
        }
        //printf("%d\n",max_occur);
        rows=unique_D;
        cols=max_occur*3;
        fu_in=(double**)malloc(rows*sizeof(double*));for(i=0;i<rows;++i)fu_in[i]=(double*)malloc(cols*sizeof(double));
        fu_out=(int**)malloc(rows*sizeof(int*));for(i=0;i<rows;++i)fu_out[i]=(int*)malloc(3*sizeof(double)); // always 3!
        //
        for(i=0;i<rows;++i)
            for(j=0;j<cols;++j)
                fu_in[i][j]=0.0;
        
        int *fourier_count;
        fourier_count=(int*)malloc(unique_D*sizeof(int));
        
        int k;
        k=-1;
        for(i=0;i<dihedral_types;++i)
        {
            if(strcmp(D1b[i],"0")!=0)
            {
                ++k;
                occur=0;
                for(j=0;j<dihedral_types;++j)
                {
                    if((strcmp(D1b[i],D1[j])==0&&strcmp(D2b[i],D2[j])==0&&strcmp(D3b[i],D3[j])==0&&strcmp(D4b[i],D4[j])==0)
                       ||
                       (strcmp(D1b[i],D4[j])==0&&strcmp(D2b[i],D3[j])==0&&strcmp(D3b[i],D2[j])==0&&strcmp(D4b[i],D1[j])==0))
                    {
                        ++occur;
                        fu_in[k][0+3*(occur-1)]=Kd[i];
                        fu_in[k][1+3*(occur-1)]=phi0[i];
                        fu_in[k][2+3*(occur-1)]=nd[i];
                    }
                }
                fourier_count[k]=occur;
            }
        }
        /*
        for(i=0;i<rows;++i)
        {
            printf("%d\n",fourier_count[i]);
            for(j=0;j<cols;++j)printf("%lf\t",fu_in[i][j]);
            printf("\n");
        }
        */
        find_unique_d(rows,cols,fu_in,fu_out);
        
        j=-1;
        printf("dihedral_style\tfourier\n");
        for(i=0;i<dihedral_types;++i)
        {
            if(strcmp(D1b[i],"0")!=0)
            {
                ++j;
                printf("dihedral\t%d\t%s\t%s\t%s\t%s\tfourier\t%d",fu_out[j][1],D1b[i],D2b[i],D3b[i],D4b[i],fourier_count[j]);
                for(k=0;k<fourier_count[j];++k)printf("\t%lf\t%lf\t%d",fu_in[j][0+3*k],fu_in[j][1+3*k],(int)fu_in[j][2+3*k]);
                printf("\n");
            }
        
        }
        
        for(i=0;i<rows;++i)if(fu_out[i][2]==-1){alterD=1;break;}
        
        for(i=0;i<rows;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        
        for(i=0;i<dihedral_types;++i){free(D1b[i]);free(D2b[i]);free(D3b[i]);free(D4b[i]);}free(D1b);free(D2b);free(D3b);free(D4b);
        
        free(fourier_count);
        
        //
        
        printf("groupFF\n");
        if(alterB==1)printf("alterB\n");else{printf("#alterB\n");}
        if(alterA==1)printf("alterA\n");else{printf("#alterA\n");}
        if(alterD==1)printf("alterD\n");else{printf("#alterD\n");}
        printf("mapFF\n");
        printf("LJ_amber\n");
        printf("pair_style\tlj/cut\t%lf\n",rcut);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,LJ_ID)==0)break;}
        for(i=0;i<atom_types;++i)
        {
            fgets(buffer,cmax_length,fp);printf("pair_coeff_amber\t%s",buffer);
        }
        printf("special\nspecial_bonds\tlj\t0.0\t0.0\t%lf\n",LJ14);
        
        printf("#\tB\tA\tD\tI\tLJ\n");
        printf("#\t%d\t%d\t%d\t%d\t%d\n",bond_types,angle_types,unique_D,0,atom_types);
            
        fclose(fp);
        
        free(r0);free(Kb);for(i=0;i<bond_types;++i){free(B1[i]);free(B2[i]);}free(B1);free(B2);
        free(theta0);free(Ka);for(i=0;i<angle_types;++i){free(A1[i]);free(A2[i]);free(A3[i]);}free(A1);free(A2);free(A3);
        free(phi0);free(Kd);free(nd);for(i=0;i<dihedral_types;++i){free(D1[i]);free(D2[i]);free(D3[i]);free(D4[i]);}free(D1);free(D2);free(D3);free(D4);
        
        //----------------------------------------------------------------------

    }

    if(argc>2)
    {
        int atoms,bonds,angles,dihedrals,k,l;
        char **neutral_types,**charged_types;
        char **B1N,**B2N,**B1C,**B2C;
        char **A1N,**A2N,**A3N,**A1C,**A2C,**A3C;
        char **D1N,**D2N,**D3N,**D4N,**D1C,**D2C,**D3C,**D4C;
        char actual[cmax_length];
        int bond_types_charged,angle_types_charged,dihedral_types_charged;
        
        char E1[sub_length],E2[sub_length],E3[sub_length],E4[sub_length];
        double d1,d2,max_fourier;
        int *fourier_count;
        
        int atom_types_charged;
        
        // tokenize variables
        int elements_out,depth;
        char **values;
    
        // read neutral mol2
        sprintf(file_path,"%s/%s",current_folder,argv[1]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);
        neutral_types=(char**)malloc(atoms*sizeof(char*));for(i=0;i<atoms;++i)neutral_types[i]=(char*)malloc(sub_length*sizeof(char));
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
        for(i=0;i<atoms;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%s",&int_buffer,neutral_types[i]);
        }
        fclose(fp);
        
        // read charged mol2
        sprintf(file_path,"%s/%s",current_folder,argv[2]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        charged_types=(char**)malloc(atoms*sizeof(char*));for(i=0;i<atoms;++i)charged_types[i]=(char*)malloc(sub_length*sizeof(char));
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
        for(i=0;i<atoms;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%s",&int_buffer,charged_types[i]);
        }
        fclose(fp);
        
        // atoms types
        for(i=0;i<atoms-1;++i)for(j=i+1;j<atoms;++j)if(strcmp(charged_types[i],charged_types[j])==0)sprintf(charged_types[j],"%s","-");
        atom_types_charged=atoms;
        for(i=0;i<atoms;++i)if(strcmp(charged_types[i],"-")==0)atom_types_charged=atom_types_charged-1;
        //for(i=0;i<atoms;++i)printf("%s\t%s\n",neutral_types[i],charged_types[i]);
        
        // read ff file
        // go to bottom and read number of types
        sprintf(file_path,"%s/%s",current_folder,argv[3]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL){}
        sscanf(buffer,"%s\t%d\t%d\t%d\t%d\t%d",word,&bond_types,&angle_types,&dihedral_types,&int_buffer,&atom_types);
        rewind(fp);
        //printf("%d\t%d\t%d\t%d\n",bond_types,angle_types,dihedral_types,atom_types);
        
        // masses
        
        for(i=0;i<atoms;++i)
        {
            if(strcmp(charged_types[i],"-")!=0)
            {
                while(fgets(buffer,cmax_length,fp)!=NULL)
                {
                    sprintf(word,"%s","");sscanf(buffer,"%s",word);
                    if(strcmp(word,"mass")==0)
                    {
                        sscanf(buffer,"%s\t%s\t%lf",word,E1,&d1);
                        if(strcmp(E1,neutral_types[i])==0)printf("mass\t%s\t%lf\n",charged_types[i],d1);
                    }
                }
                rewind(fp);
            }
        }
        
        
        
        fclose(fp);
        
        // read NEUTRAL topo
        sprintf(file_path,"%s/%s",current_folder,argv[4]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // bonds, angles, dihedrals - NO IMPROPERS!!!
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
        B1N=(char**)malloc(bonds*sizeof(char*));for(i=0;i<bonds;++i)B1N[i]=(char*)malloc(sub_length*sizeof(char));
        B2N=(char**)malloc(bonds*sizeof(char*));for(i=0;i<bonds;++i)B2N[i]=(char*)malloc(sub_length*sizeof(char));
        A1N=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)A1N[i]=(char*)malloc(sub_length*sizeof(char));
        A2N=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)A2N[i]=(char*)malloc(sub_length*sizeof(char));
        A3N=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)A3N[i]=(char*)malloc(sub_length*sizeof(char));
        D1N=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D1N[i]=(char*)malloc(sub_length*sizeof(char));
        D2N=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D2N[i]=(char*)malloc(sub_length*sizeof(char));
        D3N=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D3N[i]=(char*)malloc(sub_length*sizeof(char));
        D4N=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D4N[i]=(char*)malloc(sub_length*sizeof(char));
        // locate bonds
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Bonds\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<bonds;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%s\t%s",&int_buffer,&int_buffer,&int_buffer,&int_buffer,word,actual);
            for(j=0;j<cmax_length;++j){if(actual[j]=='\0')break;if(actual[j]=='-')actual[j]='\t';}sscanf(actual,"%s\t%s",B1N[i],B2N[i]);}
        //for(i=0;i<bonds;++i)printf("%s\t%s\n",B1N[i],B2N[i]);
        // locate angles
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Angles\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<angles;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%s\t%s",&int_buffer,&int_buffer,&int_buffer,&int_buffer,&int_buffer,word,actual);
            for(j=0;j<cmax_length;++j){if(actual[j]=='\0')break;if(actual[j]=='-')actual[j]='\t';}sscanf(actual,"%s\t%s\t%s",A1N[i],A2N[i],A3N[i]);}
        // locate dihedrals
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Dihedrals\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<dihedrals;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s",&int_buffer,&int_buffer,&int_buffer,&int_buffer,&int_buffer,&int_buffer,word,actual);
            for(j=0;j<cmax_length;++j){if(actual[j]=='\0')break;if(actual[j]=='-')actual[j]='\t';}sscanf(actual,"%s\t%s\t%s\t%s",D1N[i],D2N[i],D3N[i],D4N[i]);}
        fclose(fp);
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // charged bonds
        B1C=(char**)malloc(bonds*sizeof(char*));for(i=0;i<bonds;++i)B1C[i]=(char*)malloc(sub_length*sizeof(char));
        B2C=(char**)malloc(bonds*sizeof(char*));for(i=0;i<bonds;++i)B2C[i]=(char*)malloc(sub_length*sizeof(char));
        sprintf(file_path,"%s/%s",current_folder,argv[5]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // charged types
        for(i=0;i<4;++i)fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bond_types_charged);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angle_types_charged);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedral_types_charged);
        // locate bonds
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Bonds\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<bonds;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%s\t%s",&int_buffer,&int_buffer,&int_buffer,&int_buffer,word,actual);
            for(j=0;j<cmax_length;++j){if(actual[j]=='\0')break;if(actual[j]=='-')actual[j]='\t';}sscanf(actual,"%s\t%s",B1C[i],B2C[i]);}
        fclose(fp);
        // refine bond types: CHARGED
        for(i=0;i<bonds-1;++i)
        {
            for(j=i+1;j<bonds;++j)
            {
                if((strcmp(B1C[i],B1C[j])==0&&strcmp(B2C[i],B2C[j])==0)||(strcmp(B1C[i],B2C[j])==0&&strcmp(B2C[i],B1C[j])==0))
                {
                    sprintf(B1C[j],"%s","-");sprintf(B2C[j],"%s","-");
                }
            }
        }
        //for(i=0;i<bonds;++i)printf("%s\t%s\t|\t%s\t%s\n",B1N[i],B2N[i],B1C[i],B2C[i]);
        // open ff file
        sprintf(file_path,"%s/%s",current_folder,argv[3]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // bonds
        rows=bond_types_charged;
        cols=2;
        B1=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)B1[i]=(char*)malloc(sub_length*sizeof(char));
        B2=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)B2[i]=(char*)malloc(sub_length*sizeof(char));
        r0=(double*)malloc(rows*sizeof(double));
        Kb=(double*)malloc(rows*sizeof(double));
        j=-1;
        for(i=0;i<bonds;++i)
        {
            if(strcmp(B1C[i],"-")!=0)
            {
                ++j;
                while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"bond_style")==0)break;}
                for(k=0;k<bond_types;++k)
                {
                    fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%lf\t%lf",word,&int_buffer,E1,E2,word,&d1,&d2);
                    if((strcmp(E1,B1N[i])==0&&strcmp(E2,B2N[i])==0)||(strcmp(E1,B2N[i])==0&&strcmp(E2,B1N[i])==0))
                    {
                        r0[j]=d1;
                        Kb[j]=d2;
                        sprintf(B1[j],"%s",B1C[i]);
                        sprintf(B2[j],"%s",B2C[i]);
                        //printf("%s\t%s\t%lf\t%lf\n",B1C[i],B2C[i],d1,d2);
                    }
                }
                rewind(fp);
            }
        }
        //
        fu_in=(double**)malloc(rows*sizeof(double*));for(i=0;i<rows;++i)fu_in[i]=(double*)malloc(cols*sizeof(double));
        fu_out=(int**)malloc(rows*sizeof(int*));for(i=0;i<rows;++i)fu_out[i]=(int*)malloc(3*sizeof(double)); // always 3!
        //
        for(i=0;i<rows;++i){
            fu_in[i][0]=Kb[i];
            fu_in[i][1]=r0[i];
        }
        //
        find_unique_d(rows,cols,fu_in,fu_out);
        printf("bond_style\t%s\n","harmonic");
        for(i=0;i<rows;++i)
        {
            printf("bond\t%d\t%s\t%s\t%s\t%lf\t%lf\n",fu_out[i][1],B1[i],B2[i],"harmonic",fu_in[i][1],fu_in[i][0]);
        }
        for(i=0;i<rows;++i)if(fu_out[i][2]==-1){alterB=1;break;}
        for(i=0;i<rows;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        fclose(fp);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // charged angles
        A1C=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)A1C[i]=(char*)malloc(sub_length*sizeof(char));
        A2C=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)A2C[i]=(char*)malloc(sub_length*sizeof(char));
        A3C=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)A3C[i]=(char*)malloc(sub_length*sizeof(char));
        sprintf(file_path,"%s/%s",current_folder,argv[5]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // locate angles
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Angles\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<angles;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%s\t%s",&int_buffer,&int_buffer,&int_buffer,&int_buffer,&int_buffer,word,actual);
            for(j=0;j<cmax_length;++j){if(actual[j]=='\0')break;if(actual[j]=='-')actual[j]='\t';}sscanf(actual,"%s\t%s\t%s",A1C[i],A2C[i],A3C[i]);}
        fclose(fp);
        // refine angle types: CHARGED
        for(i=0;i<angles-1;++i)
        {
            for(j=i+1;j<angles;++j)
            {
                if((strcmp(A1C[i],A1C[j])==0&&strcmp(A2C[i],A2C[j])==0&&strcmp(A3C[i],A3C[j])==0)
                   ||
                   (strcmp(A1C[i],A3C[j])==0&&strcmp(A2C[i],A2C[j])==0&&strcmp(A3C[i],A1C[j])==0))
                {
                    sprintf(A1C[j],"%s","-");sprintf(A2C[j],"%s","-");sprintf(A3C[j],"%s","-");
                }
            }
        }
        // open ff file
        sprintf(file_path,"%s/%s",current_folder,argv[3]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // angles
        rows=angle_types_charged;
        cols=2;
        A1=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)A1[i]=(char*)malloc(sub_length*sizeof(char));
        A2=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)A2[i]=(char*)malloc(sub_length*sizeof(char));
        A3=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)A3[i]=(char*)malloc(sub_length*sizeof(char));
        theta0=(double*)malloc(rows*sizeof(double));
        Ka=(double*)malloc(rows*sizeof(double));
        j=-1;
        for(i=0;i<angles;++i)
        {
            if(strcmp(A1C[i],"-")!=0)
            {
                ++j;
                while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"angle_style")==0)break;}
                for(k=0;k<angle_types;++k)
                {
                    fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%s\t%lf\t%lf",word,&int_buffer,E1,E2,E3,word,&d1,&d2);
                    if((strcmp(E1,A1N[i])==0&&strcmp(E2,A2N[i])==0&&strcmp(E3,A3N[i])==0)||(strcmp(E1,A3N[i])==0&&strcmp(E2,A2N[i])==0&&strcmp(E3,A1N[i])==0))
                    {
                        theta0[j]=d1;
                        Ka[j]=d2;
                        sprintf(A1[j],"%s",A1C[i]);
                        sprintf(A2[j],"%s",A2C[i]);
                        sprintf(A3[j],"%s",A3C[i]);
                    }
                }
                rewind(fp);
            }
        }
        //
        fu_in=(double**)malloc(rows*sizeof(double*));for(i=0;i<rows;++i)fu_in[i]=(double*)malloc(cols*sizeof(double));
        fu_out=(int**)malloc(rows*sizeof(int*));for(i=0;i<rows;++i)fu_out[i]=(int*)malloc(3*sizeof(double)); // always 3!
        //
        for(i=0;i<rows;++i){
            fu_in[i][0]=Ka[i];
            fu_in[i][1]=theta0[i];
        }
        //
        find_unique_d(rows,cols,fu_in,fu_out);
        printf("angle_style\t%s\n","harmonic");
        for(i=0;i<rows;++i)
        {
            printf("angle\t%d\t%s\t%s\t%s\t%s\t%lf\t%lf\n",fu_out[i][1],A1[i],A2[i],A3[i],"harmonic",fu_in[i][1],fu_in[i][0]);
        }
        for(i=0;i<rows;++i)if(fu_out[i][2]==-1){alterA=1;break;}
        for(i=0;i<rows;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        fclose(fp);
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // charged dihedrals
        D1C=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D1C[i]=(char*)malloc(sub_length*sizeof(char));
        D2C=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D2C[i]=(char*)malloc(sub_length*sizeof(char));
        D3C=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D3C[i]=(char*)malloc(sub_length*sizeof(char));
        D4C=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)D4C[i]=(char*)malloc(sub_length*sizeof(char));
        sprintf(file_path,"%s/%s",current_folder,argv[5]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // locate dihedrals
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Dihedrals\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<dihedrals;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s",&int_buffer,&int_buffer,&int_buffer,&int_buffer,&int_buffer,&int_buffer,word,actual);
            for(j=0;j<cmax_length;++j){if(actual[j]=='\0')break;if(actual[j]=='-')actual[j]='\t';}sscanf(actual,"%s\t%s\t%s\t%s",D1C[i],D2C[i],D3C[i],D4C[i]);}
        fclose(fp);
        // refine dihedral types: CHARGED
        for(i=0;i<dihedrals-1;++i)
        {
            for(j=i+1;j<dihedrals;++j)
            {
                if((strcmp(D1C[i],D1C[j])==0&&strcmp(D2C[i],D2C[j])==0&&strcmp(D3C[i],D3C[j])==0&&strcmp(D4C[i],D4C[j])==0)
                   ||
                   (strcmp(D1C[i],D4C[j])==0&&strcmp(D2C[i],D3C[j])==0&&strcmp(D3C[i],D2C[j])==0&&strcmp(D4C[i],D1C[j])==0)
                   )
                {
                    sprintf(D1C[j],"%s","-");sprintf(D2C[j],"%s","-");sprintf(D3C[j],"%s","-");sprintf(D4C[j],"%s","-");
                }
            }
        }
        
        //for(i=0;i<dihedrals;++i)printf("%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\n",D1N[i],D2N[i],D3N[i],D4N[i],D1C[i],D2C[i],D3C[i],D4C[i]);
        
        // open ff file
        sprintf(file_path,"%s/%s",current_folder,argv[3]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"dihedral_style")==0)break;}
        max_fourier=1;
        for(i=0;i<dihedral_types;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d",word,&int_buffer,word,word,word,word,word,&j);if(j>max_fourier)max_fourier=j;
        }
        rewind(fp);
        
        // dihedrals
        rows=dihedral_types_charged;
        cols=3*max_fourier;
        D1=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D1[i]=(char*)malloc(sub_length*sizeof(char));
        D2=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D2[i]=(char*)malloc(sub_length*sizeof(char));
        D3=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D3[i]=(char*)malloc(sub_length*sizeof(char));
        D4=(char**)malloc(rows*sizeof(char*));for(i=0;i<rows;++i)D4[i]=(char*)malloc(sub_length*sizeof(char));
        
        fu_in=(double**)malloc(rows*sizeof(double*));for(i=0;i<rows;++i)fu_in[i]=(double*)malloc(cols*sizeof(double));
        fu_out=(int**)malloc(rows*sizeof(int*));for(i=0;i<rows;++i)fu_out[i]=(int*)malloc(3*sizeof(double)); // always 3!
        
        fourier_count=(int*)malloc(dihedral_types_charged*sizeof(int));
        
        for(i=0;i<rows;++i)
            for(j=0;j<cols;++j)
                fu_in[i][j]=0.0;
        
        j=-1;
        for(i=0;i<dihedrals;++i)
        {
            if(strcmp(D1C[i],"-")!=0)
            {
                ++j;
                while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"dihedral_style")==0)break;}
                for(k=0;k<dihedral_types;++k)
                {
                    fgets(buffer,cmax_length,fp);
                    
                    tokenize2(buffer,&elements_out,&depth,&values);
                    // 0        1   2   3   4   5   6       7   8           9           10
                    // dihedral	1	ca	ca	ca	ca	fourier	1	3.625000	180.000000	2
                    sprintf(E1,"%s",values[2]);
                    sprintf(E2,"%s",values[3]);
                    sprintf(E3,"%s",values[4]);
                    sprintf(E4,"%s",values[5]);
                    
                    
                    if((strcmp(E1,D1N[i])==0&&strcmp(E2,D2N[i])==0&&strcmp(E3,D3N[i])==0&&strcmp(E4,D4N[i])==0)
                       ||
                       (strcmp(E1,D4N[i])==0&&strcmp(E2,D3N[i])==0&&strcmp(E3,D2N[i])==0&&strcmp(E4,D1N[i])==0))
                    {

                        sprintf(D1[j],"%s",D1C[i]);
                        sprintf(D2[j],"%s",D2C[i]);
                        sprintf(D3[j],"%s",D3C[i]);
                        sprintf(D4[j],"%s",D4C[i]);
                        
                        for(l=0;l<atoi(values[7]);++l)
                        {
                            fu_in[j][0+l*3]=atof(values[8 +l*3]);
                            fu_in[j][1+l*3]=atof(values[9 +l*3]);
                            fu_in[j][2+l*3]=atof(values[10+l*3]);
                        }
                        fourier_count[j]=atoi(values[7]);
                    }
                    
                    clear_tok(elements_out,&values);
                     
                }
                rewind(fp);
            }
        }
        
        
        //for(i=0;i<rows;++i)printf("- %s\t%s\t%s\t%s\n",D1[i],D2[i],D3[i],D4[i]);
        
        
        find_unique_d(rows,cols,fu_in,fu_out);
        
        printf("dihedral_style\tfourier\n");
        for(i=0;i<rows;++i)
        {

                printf("dihedral\t%d\t%s\t%s\t%s\t%s\tfourier\t%d",fu_out[i][1],D1[i],D2[i],D3[i],D4[i],fourier_count[i]);
                for(k=0;k<fourier_count[i];++k)printf("\t%lf\t%lf\t%d",fu_in[i][0+3*k],fu_in[i][1+3*k],(int)fu_in[i][2+3*k]);
                printf("\n");
            
            
        }
        
        for(i=0;i<rows;++i)if(fu_out[i][2]==-1){alterD=1;break;}
        
        for(i=0;i<rows;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        
        free(fourier_count);
        
        fclose(fp);

        
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        printf("groupFF\n");
        if(alterB==1)printf("alterB\n");else{printf("#alterB\n");}
        if(alterA==1)printf("alterA\n");else{printf("#alterA\n");}
        if(alterD==1)printf("alterD\n");else{printf("#alterD\n");}
        printf("mapFF\n");
        printf("LJ_amber\n");
        printf("pair_style\tlj/cut\t%lf\n",rcut);
        
        // read ff file
        sprintf(file_path,"%s/%s",current_folder,argv[3]);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        // LJ
        for(i=0;i<atoms;++i)
        {
            if(strcmp(charged_types[i],"-")!=0)
            {
                while(fgets(buffer,cmax_length,fp)!=NULL)
                {
                    sprintf(word,"%s","");sscanf(buffer,"%s",word);
                    if(strcmp(word,"pair_coeff_amber")==0)
                    {
                        sscanf(buffer,"%s\t%s\t%lf\t%lf",word,E1,&d1,&d2);
                        if(strcmp(E1,neutral_types[i])==0)printf("pair_coeff_amber\t%s\t%lf\t%lf\n",charged_types[i],d1,d2);
                    }
                }
                rewind(fp);
            }
        }
        fclose(fp);
        
        printf("special\nspecial_bonds\tlj\t0.0\t0.0\t%lf\n",LJ14);
        
        printf("#\tB\tA\tD\tI\tLJ\n");
        printf("#\t%d\t%d\t%d\t%d\t%d\n",bond_types_charged,angle_types_charged,dihedral_types_charged,0,atom_types_charged);
        
        
        //
        
        
        free(r0);free(Kb);for(i=0;i<bond_types_charged;++i){free(B1[i]);free(B2[i]);}free(B1);free(B2);
        free(theta0);free(Ka);for(i=0;i<angle_types_charged;++i){free(A1[i]);free(A2[i]);free(A3[i]);}free(A1);free(A2);free(A3);
        
        for(i=0;i<atoms;++i)free(neutral_types[i]);free(neutral_types);
        for(i=0;i<atoms;++i)free(charged_types[i]);free(charged_types);
        for(i=0;i<bonds;++i){free(B1N[i]);free(B2N[i]);free(B1C[i]);free(B2C[i]);}free(B1N);free(B2N);free(B1C);free(B2C);
        for(i=0;i<angles;++i){free(A1N[i]);free(A2N[i]);free(A3N[i]);free(A1C[i]);free(A2C[i]);free(A3C[i]);}free(A1N);free(A2N);free(A3N);free(A1C);free(A2C);free(A3C);
        for(i=0;i<dihedrals;++i){free(D1N[i]);free(D2N[i]);free(D3N[i]);free(D4N[i]);free(D1C[i]);free(D2C[i]);free(D3C[i]);free(D4C[i]);}
        free(D1N);free(D2N);free(D3N);free(D4N);
        free(D1C);free(D2C);free(D3C);free(D4C);
        
        
    }

    
    //
    
    
    
    return 0;
}

//

void tokenize2(char *buffer,int *elements_out,int *depth,char ***res)
{
    //
    char second_buffer[cmax_length];
    int j,k,l;
    int end;
    int elements;
    
    int i,*positions;
    
    // loop over buffer and tokenize spaces and tabs
    end=-1;
    for(j=0;j<cmax_length;++j)
    {
        end=end+1;
        if(buffer[j]=='\0')break;
        if(buffer[j]==' ' || buffer[j]=='\t' || buffer[j]=='\n')buffer[j]='$';
    }
    // erase multiple dollars
    for(j=0;j<end-1;++j)    // end-1 because we check pairs j,j+1
    {
        if(buffer[j]=='$' && buffer[j+1]=='$')  // if you find two consecutive $s
        {
            l=-1;                               // index for writing to second_buffer
            for(k=0;k<=j;++k)                   // copy the left segment (including the first $) to memory (in second_buffer)
            {
                l=l+1;
                second_buffer[l]=buffer[k];
            }
            for(k=j+2;k<=end;++k)               // copy the remaining segment, without the second $
            {
                l=l+1;
                second_buffer[l]=buffer[k];
            }
            sprintf(buffer,"%s",second_buffer); // overwrite buffer
            j=j-1;                              // retrack to account for multiple $s
            end=end-1;                          // shrink dimension by one position
        }
    }
    // strip leading token
    if(buffer[0]=='$')
    {
        sprintf(second_buffer,"%s",buffer);
        for(j=0;j<cmax_length;++j)
        {
            if(second_buffer[j]=='\0')break;
            buffer[j]=second_buffer[j+1];
        }
    }
    // count the number of final tokens
    elements=0;
    for(j=0;j<cmax_length;++j)
    {
        if(buffer[j]=='\0')break;
        if(buffer[j]=='$')elements=elements+1;
    }
    // replace tokens
    /*
     for(j=0;j<cmax_length;++j)
     {
     if(buffer[j]=='\0')break;
     if(buffer[j]=='$')buffer[j]='\t';
     }
     */
    *depth=-1;
    if(elements>0){
        // resolve depth
        positions=(int*)malloc(elements*sizeof(int));
        k=-1;
        for(j=0;j<cmax_length;++j)
        {
            if(buffer[j]=='\0')break;
            if(buffer[j]=='$'){++k;positions[k]=j;}
        }
        
        *depth=positions[0];for(j=1;j<elements;++j)if(positions[j]-positions[j-1]>*depth)*depth=positions[j]-positions[j-1];
        *depth=*depth+1;
        // allocate results array
        *res=(char**)malloc(elements*sizeof(char*));for(j=0;j<elements;++j)(*res)[j]=(char*)malloc(*depth*sizeof(char));
        // populate results array
        i=-1; // column index
        k=0;  // row index
        for(j=0;j<cmax_length;++j)
        {
            if(buffer[j]=='\0')break;
            if(buffer[j]!='$'){++i;(*res)[k][i]=buffer[j];}
            if(buffer[j]=='$'){(*res)[k][i+1]='\0';i=-1;++k;}
        }
        
        //for(i=0;i<elements;++i)printf("%s\n",(*res)[i]);
        
        
        free(positions);
    }
    *elements_out=elements;
    //
}

void clear_tok(int elements,char ***res)
{
    int i;
    if(elements>0){for(i=0;i<elements;++i)free((*res)[i]);free(*res);}
}

