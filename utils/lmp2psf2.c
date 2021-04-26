
// sed -i -e 's/ \+ /\t/g' instruct.txt ; sed -i -e 's/ /\t/g' instruct.txt

#define cmax_length 1000
#define sub_length 10
#define Hmass 1.00794
#define Bmass 10.811
#define Cmass 12.0107
#define Nmass 14.0067
#define Omass 15.999
#define Smass 32.065
#define Fmass 18.9984
#define Clmass 35.453
#define Brmass 79.904
#define Imass 126.90447
#define pi 3.14159265359
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<math.h>

void tokenize(char *buffer,int *ignore_out,int *elements_out);

void resolve_atomic_number(int *atomic_Z,char *atom_type);
void resolve_atomic_number(int *atomic_Z,char *atom_type)
{
    *atomic_Z=0;
    // H
    if(strcmp(atom_type,"H")==0 || strcmp(atom_type,"H.spc")==0 || strcmp(atom_type,"H.t3p")==0){*atomic_Z=1;}
    // B
    if(strcmp(atom_type,"B")==0){*atomic_Z=5;}
    // C
    if(strcmp(atom_type,"C.3")==0 || strcmp(atom_type,"C.2")==0 || strcmp(atom_type,"C.1")==0 || strcmp(atom_type,"C.ar")==0 || strcmp(atom_type,"C.cat")==0){*atomic_Z=6;}
    // N
    if(strcmp(atom_type,"N.3")==0 || strcmp(atom_type,"N.2")==0 || strcmp(atom_type,"N.1")==0 || strcmp(atom_type,"N.ar")==0 || strcmp(atom_type,"N.am")==0 || strcmp(atom_type,"N.pl3")==0 || strcmp(atom_type,"N.4")==0){*atomic_Z=7;}
    // O
    if(strcmp(atom_type,"O.3")==0 || strcmp(atom_type,"O.2")==0 || strcmp(atom_type,"O.co2")==0 || strcmp(atom_type,"O.spc")==0 || strcmp(atom_type,"O.t3p")==0){*atomic_Z=8;}
    // S
    if(strcmp(atom_type,"S.3")==0 || strcmp(atom_type,"S.2")==0 || strcmp(atom_type,"S.O")==0 || strcmp(atom_type,"S.O2")==0){*atomic_Z=16;}
    // F
    if(strcmp(atom_type,"F")==0){*atomic_Z=9;}
    // Cl
    if(strcmp(atom_type,"Cl")==0){*atomic_Z=17;}
    // Br
    if(strcmp(atom_type,"Br")==0){*atomic_Z=35;}
    // I
    if(strcmp(atom_type,"I")==0){*atomic_Z=53;}
    // Du - dummy
    if(strcmp(atom_type,"Du")==0){*atomic_Z=999;}
    if(*atomic_Z==0){printf("! ... Unsupported atom type %s ... Exiting ... !\n",atom_type);exit(-1);}
}

void resolve_species(int atomic_Z,char *species);
void resolve_species(int atomic_Z,char *species)
{
    if(atomic_Z==1)sprintf(species,"%s","H");
    if(atomic_Z==5)sprintf(species,"%s","B");
    if(atomic_Z==6)sprintf(species,"%s","C");
    if(atomic_Z==7)sprintf(species,"%s","N");
    if(atomic_Z==8)sprintf(species,"%s","O");
    if(atomic_Z==9)sprintf(species,"%s","F");
    if(atomic_Z==16)sprintf(species,"%s","S");
    if(atomic_Z==17)sprintf(species,"%s","Cl");
    if(atomic_Z==35)sprintf(species,"%s","Br");
    if(atomic_Z==53)sprintf(species,"%s","I");
    if(atomic_Z==999)sprintf(species,"%s","Du");
}

double resolve_mass(int atomic_Z);
double resolve_mass(int atomic_Z)
{
    double mass;
    if(atomic_Z==1)mass=Hmass;
    if(atomic_Z==5)mass=Bmass;
    if(atomic_Z==6)mass=Cmass;
    if(atomic_Z==7)mass=Nmass;
    if(atomic_Z==8)mass=Omass;
    if(atomic_Z==9)mass=Fmass;
    if(atomic_Z==16)mass=Smass;
    if(atomic_Z==17)mass=Clmass;
    if(atomic_Z==35)mass=Brmass;
    if(atomic_Z==53)mass=Imass;
    if(atomic_Z==999)mass=999.0;
    return mass;
}

void psf_NATOM(int atom_id,char *label,int subst_id,char *subst_name,char *atom_name,char *atom_type,double q,double mass,int flag);
void psf_NATOM(int atom_id,char *label,int subst_id,char *subst_name,char *atom_name,char *atom_type,double q,double mass,int flag)
{
    if(q<0.0)
        printf("%10d %-6s %10d %-8s %-7s %-7s%-9.6lf%14.3lf%12d\n",atom_id,label,subst_id,subst_name,atom_name,atom_type,q,mass,flag);
    else
        printf("%10d %-6s %10d %-8s %-7s %-7s %-9.6lf%13.3lf%12d\n",atom_id,label,subst_id,subst_name,atom_name,atom_type,q,mass,flag);
}
int main(int argc,char **argv)
{
    
    FILE *fp,*fpw;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],mol_name[cmax_length],species[sub_length],sys[4];
    char **atom_name,**atom_type,**subst_name;
    int i,success,num_atoms;
    int *atom_id,*subst_id,*atomic_Z;
    double *x,*y,*z;
    
    char command[cmax_length];
    int atoms,bonds,angles,dihedrals,impropers,int_buffer,id1,id2,id3,id4,read_impr,read_tri;
    
    int *nx,*ny,*nz,entries,mode,*mol;
    double cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz;
    
    char word[cmax_length],*token,sp1[cmax_length],sp2[cmax_length],sp3[cmax_length],sp4[cmax_length];
    int n;
    double d1,d2,d3;
    
    double xy,xz,yz,ax,ay,az,bx,by,bz,cx,cy,cz;

	int fourier;

    //-- Start reading the mol2 input file - atoms only! -----------------------
    //--------------------------------------------------------------------------
    
    // open mol2
    getcwd(current_folder,cmax_length);
    sprintf(file_path,"%s/%s",current_folder,argv[1]);
    fp=fopen(file_path,"r");if(fp==NULL){printf("! ... Could not locate %s ... Exiting ... !\n",file_path);exit(-1);}
    
    // locate and read @<TRIPOS>MOLECULE info
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0){success=1;break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>MOLECULE ... Exiting ... !\n");exit(-1);}
    fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%s",mol_name);if(success!=1){printf("! ... Could not resolve system name ... Exiting ... !\n");exit(-1);}
    fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%d",&num_atoms);if(success!=1){printf("! ... Could not resolve atom information ... Exiting ... !\n");exit(-1);}
    // preallocate arrays
    // TRIPOS vectors
    atom_name=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)atom_name[i]=(char*)malloc(sub_length*sizeof(char));
    atom_type=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)atom_type[i]=(char*)malloc(sub_length*sizeof(char));
    subst_name=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)subst_name[i]=(char*)malloc(sub_length*sizeof(char));
    atom_id=(int*)malloc(num_atoms*sizeof(int));
    subst_id=(int*)malloc(num_atoms*sizeof(int));
    x=(double*)malloc(num_atoms*sizeof(double));
    y=(double*)malloc(num_atoms*sizeof(double));
    z=(double*)malloc(num_atoms*sizeof(double));
    // atomic number vector for species identification
    atomic_Z=(int*)malloc(num_atoms*sizeof(int));
    rewind(fp);
    
    // locate and read @<TRIPOS>ATOM info
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0){success=1;break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>ATOM ... Exiting ... !\n");exit(-1);}
    // read atomic info
    for(i=0;i<num_atoms;++i)
    {
        fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\%s\t%d\t%s",&atom_id[i],atom_name[i],&x[i],&y[i],&z[i],atom_type[i],&subst_id[i],subst_name[i]);if(success!=8){printf("! ... Error reading atomic entry %d: %s! Exiting ... !\n",i+1,buffer);exit(-1);}
    }
    rewind(fp);
    
    fclose(fp);
    
    //-- End of mol2 file reading ----------------------------------------------
    //--------------------------------------------------------------------------
    
    // resolve atomic number
    for(i=0;i<num_atoms;++i)resolve_atomic_number(&atomic_Z[i],atom_type[i]);
    
    // header
    printf("PSF EXT\n");
    printf(" \n");
    printf("%10d !NTITLE\n",1);
    printf(" %s\n",mol_name);
    printf(" \n");
    // atomic
    printf("%10d !NATOM\n",num_atoms);
    for(i=0;i<3;++i)sys[i]=mol_name[i];sys[3]='\0'; // quick and dirty...
    for(i=0;i<num_atoms;++i)
    {
        resolve_species(atomic_Z[i],species);
        psf_NATOM(atom_id[i],sys,subst_id[i],"RES",species,atom_name[i],0.0,resolve_mass(atomic_Z[i]),0);
    }
    printf("\n");
    // open lammps data file and resolve topology
    // check if the file contains a triclinic cell
    sprintf(command,"grep 'xy xz yz' %s",argv[2]);
    fp=popen(command,"r");
    sprintf(buffer,"%s","\0");
    fgets(buffer,cmax_length,fp);
    pclose(fp);
    buffer[strcspn(buffer,"\n")]='\0';
    if(strcmp(buffer,"\0")==0){read_tri=0;}else{read_tri=1;}
    // check if the file contains impropers
    sprintf(command,"grep impropers %s",argv[2]);
    fp=popen(command,"r");
    sprintf(buffer,"%s","\0");
    fgets(buffer,cmax_length,fp);
    pclose(fp);
    buffer[strcspn(buffer,"\n")]='\0';
    if(strcmp(buffer,"\0")==0){read_impr=0;}else{read_impr=1;}
    // open file
    sprintf(file_path,"%s/%s",current_folder,argv[2]);
    fp=fopen(file_path,"r");if(fp==NULL){printf("! ... Could not locate %s ... Exiting ... !\n",file_path);exit(-1);}
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);if(num_atoms!=atoms){printf("MISMATCH!!!\n");exit(-1);}
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
    if(read_impr==1){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&impropers);}
    // store supercell
    // ortho
    if(read_tri==0)
    {
        if(read_impr==0)for(i=0;i<6;++i)fgets(buffer,cmax_length,fp);
	if(read_impr==1)for(i=0;i<7;++i)fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&xlo,&xhi);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&ylo,&yhi);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&zlo,&zhi);
        lx=xhi-xlo;ly=yhi-ylo;lz=zhi-zlo;

    }
    // tri
    if(read_tri==1)
    {
        if(read_impr==0)for(i=0;i<6;++i)fgets(buffer,cmax_length,fp);
	if(read_impr==1)for(i=0;i<7;++i)fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&xlo,&xhi);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&ylo,&yhi);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&zlo,&zhi);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&xy,&xz,&yz);
    }
    // prealloc n
    nx=(int*)malloc(atoms*sizeof(int));
    ny=(int*)malloc(atoms*sizeof(int));
    nz=(int*)malloc(atoms*sizeof(int));
    mol=(int*)malloc(atoms*sizeof(int));
    // locate atoms
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Atoms\n")==0)break;fgets(buffer,cmax_length,fp);
    // mode: 1 for molecular; 2 for full
    fgets(buffer,cmax_length,fp);entries=0;for(i=0;i<cmax_length;++i){if(buffer[i]=='\0')break;if(buffer[i]=='\t')++entries;}if(entries==8)mode=1;
    if(mode==1)
    {
        sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d",&int_buffer,&mol[0],&int_buffer,&x[0],&y[0],&z[0],&nx[0],&ny[0],&nz[0]);
        for(i=1;i<atoms;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d",&int_buffer,&mol[i],&int_buffer,&x[i],&y[i],&z[i],&nx[i],&ny[i],&nz[i]);
        }
    }
    if(read_tri==0)
    {
        cell_alpha=90.0;cell_beta=90.0;cell_gamma=90.0;
        cell_a=lx;cell_b=ly;cell_c=lz;
        ax=cell_a;
        ay=0.0;
        az=0.0;
        bx=0.0;
        by=cell_b;
        bz=0.0;
        cx=0.0;
        cy=0.0;
        cz=cell_c;
        //for(i=0;i<atoms;++i){x[i]=x[i]+lx*nx[i];y[i]=y[i]+ly*ny[i];z[i]=z[i]+lz*nz[i];}
        for(i=0;i<atoms;++i){x[i]=x[i]-xlo;y[i]=y[i]-ylo;z[i]=z[i]-zlo;}
    }
    if(read_tri==1)
    {
        cell_a=xhi-xlo;cell_b=sqrt(xy*xy+(yhi-ylo)*(yhi-ylo));cell_c=sqrt((zhi-zlo)*(zhi-zlo)+xz*xz+yz*yz);
        cell_alpha=acos((xy*xz+(yhi-ylo)*yz)/(cell_b*cell_c))*180.0/pi;
        cell_beta=acos(xz/cell_c)*180.0/pi;
        cell_gamma=acos(xy/cell_b)*180.0/pi;
        
        ax=cell_a;
        ay=0.0;
        az=0.0;
        bx=cell_b*cos(cell_gamma*pi/180.0);
        by=cell_b*sin(cell_gamma*pi/180.0);
        bz=0.0;
        cx=cell_c*cos(cell_beta*pi/180.0);
        cy=(cell_b*cell_c*cos(cell_alpha*pi/180.0)-bx*cx)/by;
        cz=cell_c*cell_c-cx*cx-cy*cy;cz=sqrt(cz);
    }
    
    sprintf(file_path,"%s/%s.cell",current_folder,mol_name);
    fpw=fopen(file_path,"w+");
    fprintf(fpw,"@SET cell_ax %lf\n",ax);
    fprintf(fpw,"@SET cell_ay %lf\n",ay);
    fprintf(fpw,"@SET cell_az %lf\n",az);
    fprintf(fpw,"@SET cell_bx %lf\n",bx);
    fprintf(fpw,"@SET cell_by %lf\n",by);
    fprintf(fpw,"@SET cell_bz %lf\n",bz);
    fprintf(fpw,"@SET cell_cx %lf\n",cx);
    fprintf(fpw,"@SET cell_cy %lf\n",cy);
    fprintf(fpw,"@SET cell_cz %lf\n",cz);
    fprintf(fpw,"@SET cell_alpha %lf\n",cell_alpha);
    fprintf(fpw,"@SET cell_beta %lf\n",cell_beta);
    fprintf(fpw,"@SET cell_gamma %lf\n",cell_gamma);
    fclose(fpw);
    
    // locate bonds
    printf("%10d !NBOND\n",bonds);
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Bonds\n")==0)break;fgets(buffer,cmax_length,fp);
    for(i=0;i<bonds;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&int_buffer,&id1,&id2);
        printf("%10d%10d",id1,id2);
        if(i%4==3)printf("\n");
    }
    if((i-1)%4==3)
        printf("\n");
    else
        printf("\n\n");
    
    // locate angles
    printf("%10d !NTHETA\n",angles);
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Angles\n")==0)break;fgets(buffer,cmax_length,fp);
    for(i=0;i<angles;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&int_buffer,&int_buffer,&id1,&id2,&id3);
        printf("%10d%10d%10d",id1,id2,id3);
        if(i%3==2)printf("\n");
    }
    if((i-1)%3==2)
        printf("\n");
    else
        printf("\n\n");
    
    // locate dihedrals
    printf("%10d !NPHI\n",dihedrals);
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Dihedrals\n")==0)break;fgets(buffer,cmax_length,fp);
    for(i=0;i<dihedrals;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d",&int_buffer,&int_buffer,&id1,&id2,&id3,&id4);
        printf("%10d%10d%10d%10d",id1,id2,id3,id4);
        if(i%2==1)printf("\n");
    }
    if((i-1)%2==1)
        printf("\n");
    else
        printf("\n\n");
    
    // impropers
    if(read_impr==0)
    {
        printf("%10d !NIMPHI\n\n",0);
    }
    else
    {
        // locate impropers
        // the issue of the sequence remains to be solved!!
        printf("%10d !NIMPHI\n",impropers);
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Impropers\n")==0)break;fgets(buffer,cmax_length,fp);
        for(i=0;i<impropers;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d",&int_buffer,&int_buffer,&id1,&id2,&id3,&id4);
            printf("%10d%10d%10d%10d",id1,id2,id3,id4);
            if(i%2==1)printf("\n");
        }
        if((i-1)%2==1)
            printf("\n");
        else
            printf("\n\n");
    }
    fclose(fp);
    
    // print rest as zero...
    printf("%10d !NDON\n\n",0);
    printf("%10d !NACC\n\n",0);
    printf("%10d !NNB\n\n",0);
    printf("%10d !NGRP\n\n",0);

    // pdb
    sprintf(file_path,"%s/%s.pdb",current_folder,mol_name);
    fp=fopen(file_path,"w+");
    fprintf(fp,"CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf\n",cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma);
    // according to:
    // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    for(i=0;i<atoms;++i){
        resolve_species(atomic_Z[i],species);
        fprintf(fp,
                "%-6s%5d %-4s%c%-3s%c%4d%c    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s%2s\n",
                "ATOM",
                i+1,
                atom_name[i],
                ' ',
                "RES",
                ' ',
                mol[i],
                ' ',
                x[i]+nx[i]*ax+ny[i]*bx+nz[i]*cx,
                y[i]+nx[i]*ay+ny[i]*by+nz[i]*cy,
                z[i]+nx[i]*az+ny[i]*bz+nz[i]*cz,
                0.0,
                0.0,
                species,
                "");
        //"------||||| ----|--- |----|   --------||||||||--------||||||------          ||--"
    }
    fprintf(fp,"END\n");
    fclose(fp);
    
    // open the instructions file and find FF data
    
    //
    int ignore_out,elements_out;
    sprintf(file_path,"%s/_instruct_WITH_TABS_",current_folder,mol_name);
    fpw=fopen(file_path,"w+");
    sprintf(file_path,"%s/%s",current_folder,argv[3]);
    fp=fopen(file_path,"r");if(fp==NULL){printf("! ... Could not locate %s ... Exiting ... !\n",file_path);exit(-1);}
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        tokenize(buffer,&ignore_out,&elements_out);
        for(i=0;i<cmax_length;++i){if(buffer[i]=='\0')break;if(buffer[i]=='$')buffer[i]='\t';}
        fprintf(fpw,"%s\n",buffer);
    }
    fclose(fpw);
    fclose(fp);
    sprintf(file_path,"%s/%s",current_folder,argv[3]);
    fpw=fopen(file_path,"w+");
    sprintf(file_path,"%s/_instruct_WITH_TABS_",current_folder,mol_name);
    fp=fopen(file_path,"r");
    while(fgets(buffer,cmax_length,fp)!=NULL)fprintf(fpw,"%s",buffer);
    fclose(fp);
    fclose(fpw);
    system("rm _instruct_WITH_TABS_");
    
    // pot
    sprintf(file_path,"%s/%s.pot",current_folder,mol_name);
    fpw=fopen(file_path,"w+");
    
    sprintf(file_path,"%s/%s",current_folder,argv[3]);
    fp=fopen(file_path,"r");if(fp==NULL){printf("! ... Could not locate %s ... Exiting ... !\n",file_path);exit(-1);}
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"bond_style")==0){fprintf(fpw,"BONDS\n");}
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"bond")==0){
            sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%lf\t%lf",word,&int_buffer,sp1,sp2,word,&d1,&d2);
            fprintf(fpw,"%s\t%s\t%lf\t%lf\n",sp1,sp2,d2,d1);
        }
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"angle_style")==0){fprintf(fpw,"ANGLES\n");}
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"angle")==0){
            sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%s\t%lf\t%lf",word,&int_buffer,sp1,sp2,sp3,word,&d1,&d2);
            fprintf(fpw,"%s\t%s\t%s\t%lf\t%lf\n",sp1,sp2,sp3,d2,d1);
        }

	sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"dihedral_style")==0){fprintf(fpw,"DIHEDRALS\n");}
	sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"dihedral")==0){
            sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d",word,&int_buffer,sp1,sp2,sp3,sp4,word,&fourier);

	    token=strtok(buffer,"\t");


           
	while(token!=NULL)
            {
                if(strcmp(token,"fourier")==0)break;
                token=strtok(NULL,"\t");
            }
	 token=strtok(NULL,"\t");

		

	    for(i=0;i<fourier;++i){
		
            	fprintf(fpw,"%s\t%s\t%s\t%s\t",sp1,sp2,sp3,sp4);

		token=strtok(NULL,"\t");d1=atof(token);
            token=strtok(NULL,"\t");d2=atof(token);
            token=strtok(NULL,"\t");n=atoi(token);
            fprintf(fpw,"%lf\t%d\t%lf\n",d1,n,d2);

		
		
	    }
        }

	sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"improper_style")==0){fprintf(fpw,"IMPROPER\n");}
	sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"improper")==0){
		sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%lf",word,&int_buffer,sp1,sp2,sp3,sp4,word,&d1);
		fprintf(fpw,"%s\t%s\t%s\t%s\t%lf\t%d\t%lf\n",sp1,sp2,sp3,sp4,2.0*d1,0,0.0);

	}

	/*
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"dihedral_style")==0){fprintf(fpw,"DIHEDRALS\n");}
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"dihedral")==0){
            sscanf(buffer,"%s\t%d\t%s\t%s\t%s\t%s",word,&int_buffer,sp1,sp2,sp3,sp4);
            fprintf(fpw,"%s\t%s\t%s\t%s\t",sp1,sp2,sp3,sp4);
            token=strtok(buffer,"\t");
            while(token!=NULL)
            {
                if(strcmp(token,"charmm")==0)break;
                token=strtok(NULL,"\t");
            }
            token=strtok(NULL,"\t");
            fprintf(fpw,"%s\t",token);
            token=strtok(NULL,"\t");
            fprintf(fpw,"%s\t",token);
            token=strtok(NULL,"\t");
            fprintf(fpw,"%s",token);
            //printf("\n");
        }
	*/
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"pair_style")==0){fprintf(fpw,"NONBONDED\n");}
        sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"pair_coeff_amber")==0){
            sscanf(buffer,"%s\t%s\t%lf\t%lf",word,sp1,&d1,&d2);
            fprintf(fpw,"%s\t%lf\t%lf\t%lf\n",sp1,0.0,d2,d1);
        }
    }
    fprintf(fpw,"END\n");
    fclose(fpw);
    fclose(fp);
    
    /*
     token=strtok(buffer,"\t");
     // 1
     token=strtok(NULL,"\t");
     // 2
     token=strtok(NULL,"\t");
     while(token!=NULL)
     {
     
     
     token=strtok(NULL,"\t");
     }
     getchar();
     */
    
    // free
    for(i=0;i<num_atoms;++i){free(atom_name[i]);free(atom_type[i]);free(subst_name[i]);}
    free(atom_name);free(atom_type);free(subst_name);
    free(atom_id);free(subst_id);free(atomic_Z);
    
    free(nx);free(ny);free(nz);free(mol);
    free(x);free(y);free(z);
    
    return 0;
}

void tokenize(char *buffer,int *ignore_out,int *elements_out)
{
    ////////////////////////////////////////////////////////////////////////////////////
    // This function reads a line, tokenizes and merges consecutive white spaces and  //
    // tabs, rewrites each line in a tab-delimited format and returns the number      //
    // of tabs.                                                                       //
    ////////////////////////////////////////////////////////////////////////////////////
    // final revision: 30/04/2015 OGZ
    //
    char second_buffer[cmax_length];
    int j,k,l;
    int end,ignore=0;
    int elements;
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
    *ignore_out=ignore;
    *elements_out=elements;
    //
}

