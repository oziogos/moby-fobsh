
#define cmax_length 1000
#define sub_length 10
#define ljscale 1.122462048309373 // 2^(1/6)
#define pi 3.14159265359

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<math.h>

// OS check: important for libraries, system operations and paths
#ifdef _WIN32
#define OS 1
#elif _WIN64
#define OS 1
#elif __linux__
#define OS 2
#elif __APPLE__
#define OS 3
#elif __MACH__
#define OS 3
#else
#define OS 0
#endif

#define species_length 10
#define mol2_length 10

void eraseBAD(char mode,int argc, char **argv);

void alterA(int argc,char *argv);
void alterB(int argc,char *argv);
void alterD(int argc,char *argv);
void label(int argc,char *argv);
void mtol(int argc, char *argv1, char *argv2);
void topo(int argc, char *argv);

void tokenize(char *buffer,int *ignore_out,int *elements_out);
void select_from_tok(char *tok,int N_toks,int target,char *res);
void inline_from_select(char **words,int W,char *res);
void find_unique_1D(int N, int *array,int **out);

int one_two_init(int *B1, int *B2, int atoms, int bonds);
void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);

void resolve_species(int atomic_M,char *species);

double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz);
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz)
{
    double dot,norm,res;
    dot=Ax*Bx+Ay*By+Az*Bz;
    norm=Ax*Ax+Ay*Ay+Az*Az;
    res=dot/norm;
    return res;
}
void tri_wrap(double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double x,double y,double z,int *nx,int *ny,int *nz);
void tri_wrap(double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double x,double y,double z,int *nx,int *ny,int *nz)
{
    // Based on:
    // Tuckerman, M. Statistical Mechanics: Theory and Molecular Simulation; Oxford University Press, 2010.
    // Appendix B; p. 655
    
    double V,norm,ax_hat,ay_hat,az_hat,bx_hat,by_hat,bz_hat,cx_hat,cy_hat,cz_hat;
    double a0x_hat=1.0,a0y_hat=0.0,a0z_hat=0.0;
    double b0x_hat=0.0,b0y_hat=1.0,b0z_hat=0.0;
    double c0x_hat=0.0,c0y_hat=0.0,c0z_hat=1.0;
    double projaa0,projab0,projac0,projba0,projbb0,projbc0,projca0,projcb0,projcc0;
    double denom,a_norm_REF,b_norm_REF,c_norm_REF,nom1,nom2,nom3,a,b,c,a_sign,b_sign,c_sign,a_norm,b_norm,c_norm;
    
    V=ax*by*cz;
    
    // normalize
    norm=ax*ax+ay*ay+az*az;norm=sqrt(norm);
    ax_hat=ax/norm;
    ay_hat=ay/norm;
    az_hat=az/norm;
    norm=bx*bx+by*by+bz*bz;norm=sqrt(norm);
    bx_hat=bx/norm;
    by_hat=by/norm;
    bz_hat=bz/norm;
    norm=cx*cx+cy*cy+cz*cz;norm=sqrt(norm);
    cx_hat=cx/norm;
    cy_hat=cy/norm;
    cz_hat=cz/norm;
    
    // calculate projections
    projaa0=proj(ax_hat,ay_hat,az_hat,a0x_hat,a0y_hat,a0z_hat);
    projab0=proj(ax_hat,ay_hat,az_hat,b0x_hat,b0y_hat,b0z_hat);
    projac0=proj(ax_hat,ay_hat,az_hat,c0x_hat,c0y_hat,c0z_hat);
    projba0=proj(bx_hat,by_hat,bz_hat,a0x_hat,a0y_hat,a0z_hat);
    projbb0=proj(bx_hat,by_hat,bz_hat,b0x_hat,b0y_hat,b0z_hat);
    projbc0=proj(bx_hat,by_hat,bz_hat,c0x_hat,c0y_hat,c0z_hat);
    projca0=proj(cx_hat,cy_hat,cz_hat,a0x_hat,a0y_hat,a0z_hat);
    projcb0=proj(cx_hat,cy_hat,cz_hat,b0x_hat,b0y_hat,b0z_hat);
    projcc0=proj(cx_hat,cy_hat,cz_hat,c0x_hat,c0y_hat,c0z_hat);
    
    // projection denominator
    denom=(-projac0*projbb0*projca0+projab0*projbc0*projca0+projac0*projba0*projcb0-projaa0*projbc0*projcb0-projab0*projba0*projcc0+projaa0*projbb0*projcc0);
    
    // oblique supercell vector lengths
    a_norm_REF=ax*ax;a_norm_REF=sqrt(a_norm_REF);
    b_norm_REF=bx*bx+by*by;b_norm_REF=sqrt(b_norm_REF);
    c_norm_REF=cx*cx+cy*cy+cz*cz;c_norm_REF=sqrt(c_norm_REF);
    
    // projections
    nom1=-1.0*(projbc0*projcb0*x-projbb0*projcc0*x-projbc0*projca0*y+projba0*projcc0*y+projbb0*projca0*z-projba0*projcb0*z);
    nom2=(projac0*projcb0*x-projab0*projcc0*x-projac0*projca0*y+projaa0*projcc0*y+projab0*projca0*z-projaa0*projcb0*z);
    nom3=-1.0*(projac0*projbb0*x-projab0*projbc0*x-projac0*projba0*y+projaa0*projbc0*y+projab0*projba0*z-projaa0*projbb0*z);
    
    // coordinates with respect to oblique system
    a=nom1/denom;
    b=nom2/denom;
    c=nom3/denom;
    
    // +/- direction
    if((a*ax_hat-0.0)*(ax-0.0)+(a*ay_hat-0.0)*(0.0-0.0)+(a*az_hat-0.0)*(0.0-0.0) < 0.0 ){a_sign=-1.0;}else{a_sign=+1.0;}
    if((b*bx_hat-0.0)*(bx-0.0)+(b*by_hat-0.0)*(by-0.0)+(b*bz_hat-0.0)*(0.0-0.0) < 0.0 ){b_sign=-1.0;}else{b_sign=+1.0;}
    if((c*cx_hat-0.0)*(cx-0.0)+(c*cy_hat-0.0)*(cy-0.0)+(c*cz_hat-0.0)*(cz-0.0) < 0.0 ){c_sign=-1.0;}else{c_sign=+1.0;}
    
    // image calculation
    a_norm=(a*ax_hat)*(a*ax_hat)+(a*ay_hat)*(a*ay_hat)+(a*az_hat)*(a*az_hat);a_norm=sqrt(a_norm);
    b_norm=(b*bx_hat)*(b*bx_hat)+(b*by_hat)*(b*by_hat)+(b*bz_hat)*(b*bz_hat);b_norm=sqrt(b_norm);
    c_norm=(c*cx_hat)*(c*cx_hat)+(c*cy_hat)*(c*cy_hat)+(c*cz_hat)*(c*cz_hat);c_norm=sqrt(c_norm);
    *nx=(int)floor(a_norm/a_norm_REF);if(a_sign<0.0)*nx=-(*nx)-1;
    *ny=(int)floor(b_norm/b_norm_REF);if(b_sign<0.0)*ny=-(*ny)-1;
    *nz=(int)floor(c_norm/c_norm_REF);if(c_sign<0.0)*nz=-(*nz)-1;
}


int main(int argc,char **argv)
{
    FILE *fp,*fpw,*fpr;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],word[cmax_length],command[cmax_length];
    char **juice,species[sub_length];
    char mol2_filename[cmax_length],typestr[cmax_length];
    int start=0,stop=0,i,j,k,lines,backup=0;
    int ignore_out,elements_out,W,atoms;
    char **words;
    double d_buffer,xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz,x,y,z,q;
    int id,mol,type,ix,iy,iz;
    int B_type,topo_bond_types,topo_B_id;
    char B1_species[sub_length],B2_species[sub_length],topo_B1_species[sub_length],topo_B2_species[sub_length];
    int A_type,topo_angle_types,topo_A_id;
    char A1_species[sub_length],A2_species[sub_length],A3_species[sub_length],topo_A1_species[sub_length],topo_A2_species[sub_length],topo_A3_species[sub_length];
    int D_type,topo_dihedral_types,topo_D_id;
    char D1_species[sub_length],D2_species[sub_length],D3_species[sub_length],D4_species[sub_length],topo_D1_species[sub_length],topo_D2_species[sub_length],topo_D3_species[sub_length],topo_D4_species[sub_length];
    int int_buffer;
    int *map,bonds,angles,dihedrals,c0,c1,c2,c3,c4,**fu_out,*map_unique;
    int atom_types,hybrid,*search_vector;
    double d1,d2,d3,d4;
    char table1[cmax_length],table2[cmax_length];
    int pair_type=0;
    char **global_species;
    int iflag=0;
    int improper_types;
    int *B1,*B2,max_1_2,**one_two;
    char **I1_type,**I2_type,**I3_type,**I4_type;
    int *I_type_init;
    char **mol2_species;
    int impropers=0,improper_types_final;
    int one_five=0;
    int *D1,*D2,*D3,*D4;
    int one,four,five,onefive0,*myonefiveL0,*myonefiveR0,onefive,*myonefiveL,*myonefiveR;
    int one_five_types,one_five_types_final,*one_five_type_init;
    char **one_five1_type,**one_five2_type;
    int one_five_final;
    int bond_types;
    double K0,K1,K2,K3,K4;
    char boundx[cmax_length],boundy[cmax_length],boundz[cmax_length],name[cmax_length],cell_type[cmax_length];
    
    int success;
    double cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma;
    double ax,ay,az,bx,by,bz,cx,cy,cz,xy,xz,yz;
    
    int l,fourier_n,fourier_count;
    double fourier_K,fourier_d;

    //
    getcwd(current_folder,cmax_length);
    //
    /*
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"topo",current_folder,"topo");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"label",current_folder,"label");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"alterB",current_folder,"alterB");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"alterA",current_folder,"alterA");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"alterD",current_folder,"alterD");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"mol2_deleteH",current_folder,"mol2_deleteH");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"mol2_rename",current_folder,"mol2_rename");system(command);
    sprintf(command,"gcc -o %s/src/%s %s/src/%s.c -lm -O3",current_folder,"mtol",current_folder,"mtol");system(command);
    
    sprintf(command,"cp %s/src/%s .",current_folder,"topo");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"label");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"alterB");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"alterA");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"alterD");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"mol2_deleteH");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"mol2_rename");system(command);
    sprintf(command,"cp %s/src/%s .",current_folder,"mtol");system(command);
    */
    //
    
    // read the instruction file and store to memory
    sprintf(file_path,"%s/%s",current_folder,argv[1]);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    i=0;while(fgets(buffer,cmax_length,fp)!=NULL){i=i+1;sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"begin")==0){start=i;break;}}rewind(fp);
    i=0;while(fgets(buffer,cmax_length,fp)!=NULL){i=i+1;sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"end")==0){stop=i;break;}}rewind(fp);
    lines=stop-start+1;
    juice=(char**)malloc(lines*sizeof(char*));for(i=0;i<lines;++i)juice[i]=(char*)malloc(cmax_length*sizeof(char));
    i=0;while(fgets(buffer,cmax_length,fp)!=NULL){i=i+1;if(i==start)break;}
    i=0;sprintf(juice[i],"%s",buffer);
    for(i=1;i<lines;++i){fgets(buffer,cmax_length,fp);sprintf(juice[i],"%s",buffer);}
    fclose(fp);
    // look for 'source': store mol2 filename and check if it is present
    for(i=0;i<lines;++i){sprintf(word,"%s","");sscanf(juice[i],"%s",word);if(strcmp(word,"source")==0){sscanf(juice[i],"%s\t%s",word,mol2_filename);break;}}
    sprintf(file_path,"%s/%s.mol2",current_folder,mol2_filename);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}fclose(fp);
    // backup the mol2 file
    /*
    for(i=0;i<lines;++i){sprintf(word,"%s","");sscanf(juice[i],"%s",word);if(strcmp(word,"backup")==0){backup=1;break;}}
    if(backup==1){sprintf(buffer,"cp %s.mol2 %s.mol2.backup",mol2_filename,mol2_filename);system(buffer);}
    */
    //
    for(i=0;i<lines;++i){sprintf(word,"%s","");sscanf(juice[i],"%s",word);if(strcmp(word,"name")==0){sscanf(juice[i],"%s\t%s",word,name);break;}}
    //sprintf(buffer,"cp %s.mol2 result.mol2",mol2_filename);system(buffer);
    //sprintf(mol2_filename,"%s","result");
    sprintf(buffer,"cp %s.mol2 %s.mol2",mol2_filename,name);system(buffer);
    sprintf(mol2_filename,"%s",name);
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"rename")==0)
        {
            tokenize(juice[i],&ignore_out,&elements_out);
            W=elements_out+1+1;
            words=(char**)malloc(W*sizeof(char*));
            for(j=0;j<W;++j)words[j]=(char*)malloc(cmax_length*sizeof(char));
            sprintf(words[0],"%s","./mol2_rename");
            sprintf(words[1],"%s.mol2",mol2_filename);
            for(j=2;j<W-1;++j)
                select_from_tok(juice[i],elements_out,j-1,words[j]);
            sprintf(words[W-1],"> %s.mol2.temp ; mv %s.mol2.temp %s.mol2",mol2_filename,mol2_filename,mol2_filename);
            inline_from_select(words,W,command);
            system(command);
            for(j=0;j<W;++j)free(words[j]);free(words);
        }
    }
    /*
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"rename_CP")==0)
        {
            system("./rename_CP result.mol2 > test.mol2 ; mv test.mol2 result.mol2");
        }
    }
    */
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"deleteH")==0)
        {
            sscanf(juice[i],"%s\t%s",word,species);
            sprintf(command,"./mol2_deleteH %s.mol2 H %s > %s.mol2.temp ; mv %s.mol2.temp %s.mol2",mol2_filename,species,mol2_filename,mol2_filename,mol2_filename);
            system(command);
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"mol2tolammps")==0)
        {
            sscanf(juice[i],"%s\t%s",word,typestr);
            if(strcmp(typestr,"full")==0){
                //sprintf(command,"./mtol %s 2",mol2_filename);
                //system(command);
                mtol(3,mol2_filename,"2");
            }
            if(strcmp(typestr,"molecular")==0){
                //sprintf(command,"./mtol %s 1",mol2_filename);
                //system(command);
                mtol(3,mol2_filename,"1");
            }
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"charge")==0)
        {
            sscanf(juice[i],"%s\t%s\t%lf",word,buffer,&d_buffer);
            if(OS==3)sprintf(command,"sed -i .bak 's/%s_charge/%lf/g' init.lammps",buffer,d_buffer);
            else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s_charge/%lf/g' init.lammps",buffer,d_buffer);
            printf("%s\n",command);
            system(command);
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"mass")==0)
        {
            sscanf(juice[i],"%s\t%s\t%lf",word,buffer,&d_buffer);
            if(OS==3)sprintf(command,"sed -i .bak 's/%s_mass/%lf/g' init.lammps",buffer,d_buffer);
            else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s_mass/%lf/g' init.lammps",buffer,d_buffer);
            printf("%s\n",command);
            system(command);
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"boundary")==0)
        {
            sscanf(juice[i],"%s\t%s\t%s\t%s",word,boundx,boundy,boundz);
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"supercell")==0)
        {
            sprintf(cell_type,"%s","");
            sscanf(juice[i],"%s\t%s",word,cell_type);
            
            if(strcmp(cell_type,"ortho")==0)
            {
                //sscanf(juice[i],"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%s",word,&xlo,&xhi,&ylo,&yhi,&zlo,&zhi,typestr);
                //sscanf(juice[i],"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",word,&xlo,&xhi,&ylo,&yhi,&zlo,&zhi);
                sscanf(juice[i],"%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",word,word,&xlo,&xhi,&ylo,&yhi,&zlo,&zhi);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","xmin",xlo);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","xmin",xlo);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","xmax",xhi);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","xmax",xhi);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","ymin",ylo);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","ymin",ylo);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","ymax",yhi);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","ymax",yhi);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","zmin",zlo);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","zmin",zlo);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","zmax",zhi);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","zmax",zhi);
                printf("%s\n",command);
                system(command);
                
                sprintf(file_path,"%s/%s",current_folder,"init.lammps");
                fp=fopen(file_path,"r");
                sprintf(file_path,"%s/%s",current_folder,"init.lammps.bak");
                fpw=fopen(file_path,"w+");
                if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&atoms);
                while(fgets(buffer,cmax_length,fp)!=NULL)
                {
                    fprintf(fpw,"%s",buffer);
                    if(strcmp(buffer,"Atoms\n")==0)break;
                }
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
                
                lx=xhi-xlo;
                ly=yhi-ylo;
                lz=zhi-zlo;
                if(strcmp(typestr,"full")==0)
                {
                    for(j=0;j<atoms;++j)
                    {
                        fgets(buffer,cmax_length,fp);
                        sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf",&id,&mol,&type,&q,&x,&y,&z);
                        ix=0;iy=0;iz=0;
                        if(x>xhi){ix=(int)floor((x-xlo)/lx);}
                        if(x<xlo){ix=-(int)floor((xhi-x)/lx);}
                        if(y>yhi){iy=(int)floor((y-ylo)/ly);}
                        if(y<ylo){iy=-(int)floor((yhi-y)/ly);}
                        if(z>zhi){iz=(int)floor((z-zlo)/lz);}
                        if(z<zlo){iz=-(int)floor((zhi-z)/lz);}
                        fprintf(fpw,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",id,mol,type,q,x-(double)ix*lx,y-(double)iy*ly,z-(double)iz*lz,ix,iy,iz);
                    }
                }
                if(strcmp(typestr,"molecular")==0)
                {
                    for(j=0;j<atoms;++j)
                    {
                        fgets(buffer,cmax_length,fp);
                        sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf",&id,&mol,&type,&x,&y,&z);
                        ix=0;iy=0;iz=0;
                        if(x>xhi){ix=(int)floor((x-xlo)/lx);}
                        if(x<xlo){ix=-(int)floor((xhi-x)/lx);}
                        if(y>yhi){iy=(int)floor((y-ylo)/ly);}
                        if(y<ylo){iy=-(int)floor((yhi-y)/ly);}
                        if(z>zhi){iz=(int)floor((z-zlo)/lz);}
                        if(z<zlo){iz=-(int)floor((zhi-z)/lz);}
                        fprintf(fpw,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",id,mol,type,x-(double)ix*lx,y-(double)iy*ly,z-(double)iz*lz,ix,iy,iz);
                    }
                }
                fclose(fpw);
                fclose(fp);
                
                system("mv init.lammps.bak init.lammps");
            }
            else if(strcmp(cell_type,"mol2")==0)
            {
                // open mol2 and read CRYSIN
                sprintf(file_path,"%s/%s.mol2",current_folder,name);fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
                success=0;
                while(fgets(buffer,cmax_length,fp)!=NULL)
                {
                    if(strcmp(buffer,"@<TRIPOS>CRYSIN\n")==0){success=1;printf("$ ... Located @<TRIPOS>CRYSIN ...\n");break;}
                }
                if(success==0){printf("! ... Unable to locate @<TRIPOS>CRYSIN ... Exiting ... !\n");exit(-1);}
                // read CRYSIN
                fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf",&cell_a,&cell_b,&cell_c,&cell_alpha,&cell_beta,&cell_gamma);
                printf("$ ... Cell info ...\n$ ... [ a b c ] = [ %lf %lf %lf ] ...\n$ ... [ alpha beta gamma ] = [ %lf %lf %lf ] ...\n",cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma);
                fclose(fp);
                // lammps tri data
                ax=cell_a;
                ay=0.0;
                az=0.0;
                bx=cell_b*cos(cell_gamma*pi/180.0);
                by=cell_b*sin(cell_gamma*pi/180.0);
                bz=0.0;
                cx=cell_c*cos(cell_beta*pi/180.0);
                cy=(cell_b*cell_c*cos(cell_alpha*pi/180.0)-bx*cx)/by;
                cz=cell_c*cell_c-cx*cx-cy*cy;cz=sqrt(cz);
                //
                xlo=0.0;ylo=0.0;zlo=0.0;
                xhi=ax;xy=bx;yhi=by;
                xz=cx;yz=cy;zhi=cz;
                /*
                printf("%lf\t%lf\n",xlo,xhi);
                printf("%lf\t%lf\n",ylo,yhi);
                printf("%lf\t%lf\n",zlo,zhi);
                printf("%lf\t%lf\t%lf\n",xy,xz,yz);
                */
                
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","xmin",xlo);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","xmin",xlo);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","xmax",xhi);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","xmax",xhi);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","ymin",ylo);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","ymin",ylo);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","ymax",yhi);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","ymax",yhi);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","zmin",zlo);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","zmin",zlo);
                printf("%s\n",command);
                system(command);
                if(OS==3)sprintf(command,"sed -i .bak 's/%s/%lf/g' init.lammps","zmax",zhi);
                else if(OS==2 || OS==1)sprintf(command,"sed -i 's/%s/%lf/g' init.lammps","zmax",zhi);
                printf("%s\n",command);
                system(command);
                
                sprintf(file_path,"%s/%s",current_folder,"init.lammps");
                fp=fopen(file_path,"r");
                sprintf(file_path,"%s/%s",current_folder,"init.lammps.bak");
                fpw=fopen(file_path,"w+");
                if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&atoms);
                while(fgets(buffer,cmax_length,fp)!=NULL)
                {
                    fprintf(fpw,"%s",buffer);
                    if(strcmp(buffer,"Atoms\n")==0)break;
                }
                fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
                
                if(strcmp(typestr,"full")==0)
                {
                    for(j=0;j<atoms;++j)
                    {
                        fgets(buffer,cmax_length,fp);
                        sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf",&id,&mol,&type,&q,&x,&y,&z);
                        ix=0;iy=0;iz=0;
                        
                        tri_wrap(ax,ay,az,bx,by,bz,cx,cy,cz,x,y,z,&ix,&iy,&iz);
                        /*
                         // periodic wrapping
                         x[i]=x[i]-nx[i]*ax-ny[i]*bx-nz[i]*cx;
                         y[i]=y[i]-nx[i]*(0.0)-ny[i]*by-nz[i]*cy;
                         z[i]=z[i]-nx[i]*(0.0)-ny[i]*(0.0)-nz[i]*cz;
                         */
                        fprintf(fpw,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",id,mol,type,q,
                                x-(double)ix*ax-(double)iy*bx-(double)iz*cx,
                                y-(double)iy*by-(double)iz*cy,
                                z-(double)iz*cz,
                                ix,
                                iy,
                                iz);
                    }
                }
                if(strcmp(typestr,"molecular")==0)
                {
                    for(j=0;j<atoms;++j)
                    {
                        fgets(buffer,cmax_length,fp);
                        sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf",&id,&mol,&type,&x,&y,&z);
                        ix=0;iy=0;iz=0;
                        
                        tri_wrap(ax,ay,az,bx,by,bz,cx,cy,cz,x,y,z,&ix,&iy,&iz);
                        /*
                         // periodic wrapping
                         x[i]=x[i]-nx[i]*ax-ny[i]*bx-nz[i]*cx;
                         y[i]=y[i]-nx[i]*(0.0)-ny[i]*by-nz[i]*cy;
                         z[i]=z[i]-nx[i]*(0.0)-ny[i]*(0.0)-nz[i]*cz;
                         */
                        fprintf(fpw,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",id,mol,type,
                                x-(double)ix*ax-(double)iy*bx-(double)iz*cx,
                                y-(double)iy*by-(double)iz*cy,
                                z-(double)iz*cz,
                                ix,
                                iy,
                                iz);
                    }
                }
                
                fclose(fpw);
                fclose(fp);
                
                system("mv init.lammps.bak init.lammps");
            }
            
            
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"topo")==0)
        {
            //sprintf(command,"./topo %s",mol2_filename);
            //system(command);
            topo(2,mol2_filename);
        }
    }
    //



//---------------- erase - - erase ---------------------------------------

    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"eraseB")==0)
        {
	    sscanf(juice[i],"%s\t%d",word,&c0);//printf("%s",juice[i]);
            char **dummy_argv;
	    dummy_argv=(char**)malloc(c0*sizeof(char*));
	    for(j=0;j<c0;++j)dummy_argv[j]=(char*)malloc(sub_length*sizeof(char));
	    for(j=0;j<c0;++j)
	    {
		++i;sscanf(juice[i],"%s",dummy_argv[j]);
		}
	    printf("%d\n",c0);for(j=0;j<c0;++j)printf("%s\n",dummy_argv[j]);
            eraseBAD('B',c0,dummy_argv);
	    system("mv topo.log.erased topo.log");
	    system("mv topo.out.erased topo.out");
	    for(j=0;j<c0;++j)free(dummy_argv[j]);free(dummy_argv);
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"eraseA")==0)
        {
	    sscanf(juice[i],"%s\t%d",word,&c0);//printf("%s",juice[i]);
            char **dummy_argv;
	    dummy_argv=(char**)malloc(c0*sizeof(char*));
	    for(j=0;j<c0;++j)dummy_argv[j]=(char*)malloc(sub_length*sizeof(char));
	    for(j=0;j<c0;++j)
	    {
		++i;sscanf(juice[i],"%s",dummy_argv[j]);
		}
	    printf("%d\n",c0);for(j=0;j<c0;++j)printf("%s\n",dummy_argv[j]);
            eraseBAD('A',c0,dummy_argv);
	    system("mv topo.log.erased topo.log");
	    system("mv topo.out.erased topo.out");
	    for(j=0;j<c0;++j)free(dummy_argv[j]);free(dummy_argv);
        }
    }
    //	
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"eraseD")==0)
        {
	    sscanf(juice[i],"%s\t%d",word,&c0);//printf("%s",juice[i]);
            char **dummy_argv;
	    dummy_argv=(char**)malloc(c0*sizeof(char*));
	    for(j=0;j<c0;++j)dummy_argv[j]=(char*)malloc(sub_length*sizeof(char));
	    for(j=0;j<c0;++j)
	    {
		++i;sscanf(juice[i],"%s",dummy_argv[j]);
		}
	    printf("%d\n",c0);for(j=0;j<c0;++j)printf("%s\n",dummy_argv[j]);
            eraseBAD('D',c0,dummy_argv);
	    system("mv topo.log.erased topo.log");
	    system("mv topo.out.erased topo.out");
	    for(j=0;j<c0;++j)free(dummy_argv[j]);free(dummy_argv);
        }
    }
    //

//---------------- erase - - erase ----------------------------------- END

    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"bond")==0)
        {
            sscanf(juice[i],"%s\t%d\t%s\t%s",word,&B_type,B1_species,B2_species);
            sprintf(file_path,"%s/%s",current_folder,"topo.log.bak");
            fpw=fopen(file_path,"w+");
            sprintf(file_path,"%s/%s",current_folder,"topo.log");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            
            
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
                if(strcmp(buffer,"Bonds\n")==0)
                    break;}
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&topo_bond_types);fprintf(fpw,"%s",buffer);
            for(j=0;j<topo_bond_types;++j)
            {
                fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s",&topo_B_id,topo_B1_species,topo_B2_species);
                if(
                   (strcmp(B1_species,topo_B1_species)==0 && strcmp(B2_species,topo_B2_species)==0)||(strcmp(B1_species,topo_B2_species)==0 && strcmp(B2_species,topo_B1_species)==0))
                {
                    fprintf(fpw,"[%d]\t%s\t%s\t\t\t\t\t#%d\n",topo_B_id,topo_B1_species,topo_B2_species,B_type);
                }
                else
                {
                    fprintf(fpw,"%s",buffer);
                }
            }
            
            while(fgets(buffer,cmax_length,fp)!=NULL)
                fprintf(fpw,"%s",buffer);
            
            fclose(fp);
            fclose(fpw);
            
            system("mv topo.log.bak topo.log");
            
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"angle")==0)
        {
            sscanf(juice[i],"%s\t%d\t%s\t%s\t%s",word,&A_type,A1_species,A2_species,A3_species);
            sprintf(file_path,"%s/%s",current_folder,"topo.log.bak");
            fpw=fopen(file_path,"w+");
            sprintf(file_path,"%s/%s",current_folder,"topo.log");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            
            
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
                if(strcmp(buffer,"Angles\n")==0)
                    break;}
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&topo_angle_types);fprintf(fpw,"%s",buffer);
            for(j=0;j<topo_angle_types;++j)
            {
                fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s",&topo_A_id,topo_A1_species,topo_A2_species,topo_A3_species);
                if(
                   (strcmp(A1_species,topo_A1_species)==0 && strcmp(A2_species,topo_A2_species)==0 && strcmp(A3_species,topo_A3_species)==0)
                   ||
                   (strcmp(A1_species,topo_A3_species)==0 && strcmp(A2_species,topo_A2_species)==0 && strcmp(A3_species,topo_A1_species)==0)
                   )
                {
                    fprintf(fpw,"[%d]\t%s\t%s\t%s\t\t\t\t\t#%d\n",topo_A_id,topo_A1_species,topo_A2_species,topo_A3_species,A_type);
                }
                else
                {
                    fprintf(fpw,"%s",buffer);
                }
            }
            
            while(fgets(buffer,cmax_length,fp)!=NULL)
                fprintf(fpw,"%s",buffer);
            
            fclose(fp);
            fclose(fpw);
            
            system("mv topo.log.bak topo.log");
            
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"dihedral")==0)
        {
            sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s",word,&D_type,D1_species,D2_species,D3_species,D4_species);
            sprintf(file_path,"%s/%s",current_folder,"topo.log.bak");
            fpw=fopen(file_path,"w+");
            sprintf(file_path,"%s/%s",current_folder,"topo.log");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            
            
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
                if(strcmp(buffer,"Dihedrals\n")==0)
                    break;}
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&topo_dihedral_types);fprintf(fpw,"%s",buffer);
            for(j=0;j<topo_dihedral_types;++j)
            {
                fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s",&topo_D_id,topo_D1_species,topo_D2_species,topo_D3_species,topo_D4_species);
                if(
                   (strcmp(D1_species,topo_D1_species)==0 && strcmp(D2_species,topo_D2_species)==0 && strcmp(D3_species,topo_D3_species)==0 && strcmp(D4_species,topo_D4_species)==0)
                   ||
                   (strcmp(D1_species,topo_D4_species)==0 && strcmp(D2_species,topo_D3_species)==0 && strcmp(D3_species,topo_D2_species)==0 && strcmp(D4_species,topo_D1_species)==0)
                   )
                {
                    fprintf(fpw,"[%d]\t%s\t%s\t%s\t%s\t\t\t\t\t#%d\n",topo_D_id,topo_D1_species,topo_D2_species,topo_D3_species,topo_D4_species,D_type);
                }
                else
                {
                    fprintf(fpw,"%s",buffer);
                }
            }
            
            while(fgets(buffer,cmax_length,fp)!=NULL)
                fprintf(fpw,"%s",buffer);
            
            fclose(fp);
            fclose(fpw);
            
            system("mv topo.log.bak topo.log");
            
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"groupFF")==0)
        {
            //system("./label topo.log");
            label(2,"topo.log");
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"alterB")==0)
        {
            //system("./alterB bonds");
            alterB(2,"bonds");
            system("mv topo.log.altered topo.log");
            system("mv topo.out.altered topo.out");
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"alterA")==0)
        {
            //system("./alterA angles");
            alterA(2,"angles");
            system("mv topo.log.altered topo.log");
            system("mv topo.out.altered topo.out");
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"alterD")==0)
        {
            //system("./alterD dihedrals");
            alterD(2,"dihedrals");
            system("mv topo.log.altered topo.log");
            system("mv topo.out.altered topo.out");
        }
    }
    //
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"mapFF")==0)
        {
            //
            sprintf(file_path,"%s/%s",current_folder,"topo.log");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Bonds\n")==0)break;
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&int_buffer,&topo_bond_types);
            map=(int*)malloc(topo_bond_types*sizeof(int));
            for(j=0;j<topo_bond_types;++j)
            {
                fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s\t%s\t#%d",word,word,word,&map[j]);
            }
            fclose(fp);
            
            map_unique=(int*)malloc(topo_bond_types*sizeof(int));
            fu_out=(int**)malloc(topo_bond_types*sizeof(int*));
            for(j=0;j<topo_bond_types;++j)fu_out[j]=(int*)malloc(3*sizeof(int));
            find_unique_1D(topo_bond_types,map,fu_out);
            k=-1;
            for(j=0;j<topo_bond_types;++j)if(fu_out[j][2]!=-1){k=k+1;map_unique[k]=map[fu_out[j][0]-1];}
            for(j=0;j<topo_bond_types;++j)free(fu_out[j]);free(fu_out);
            
            sprintf(file_path,"%s/%s",current_folder,"topo.out");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            sprintf(file_path,"%s/%s",current_folder,"topo.out.bak");
            fpw=fopen(file_path,"w+");
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&bonds);
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&angles);
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&dihedrals);
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
                if(strcmp(buffer,"Bonds\n")==0)break;
            }
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
            for(j=0;j<bonds;++j)
            {
                fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&c0,&c1,&c2);
                fprintf(fpw,"%d\t%d\t%d\t%d\n",int_buffer,map_unique[c0-1],c1,c2);
                //fprintf(fpw,"%s",buffer);
            }
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
            }
            fclose(fpw);
            fclose(fp);
            free(map);free(map_unique);
            system("mv topo.out.bak topo.out");
            
            //
            sprintf(file_path,"%s/%s",current_folder,"topo.log");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Angles\n")==0)break;
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&int_buffer,&topo_angle_types);
            map=(int*)malloc(topo_angle_types*sizeof(int));
            for(j=0;j<topo_angle_types;++j)
            {
                fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s\t%s\t%s\t#%d",word,word,word,word,&map[j]);
            }
            fclose(fp);
            
            map_unique=(int*)malloc(topo_angle_types*sizeof(int));
            fu_out=(int**)malloc(topo_angle_types*sizeof(int*));
            for(j=0;j<topo_angle_types;++j)fu_out[j]=(int*)malloc(3*sizeof(int));
            find_unique_1D(topo_angle_types,map,fu_out);
            k=-1;
            for(j=0;j<topo_angle_types;++j)if(fu_out[j][2]!=-1){k=k+1;map_unique[k]=map[fu_out[j][0]-1];}
            for(j=0;j<topo_angle_types;++j)free(fu_out[j]);free(fu_out);
            
            sprintf(file_path,"%s/%s",current_folder,"topo.out");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            sprintf(file_path,"%s/%s",current_folder,"topo.out.bak");
            fpw=fopen(file_path,"w+");
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&bonds);
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&angles);
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&dihedrals);
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
                if(strcmp(buffer,"Angles\n")==0)break;
            }
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
            for(j=0;j<angles;++j)
            {
                fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&int_buffer,&c0,&c1,&c2,&c3);
                fprintf(fpw,"%d\t%d\t%d\t%d\t%d\n",int_buffer,map_unique[c0-1],c1,c2,c3);
                //fprintf(fpw,"%s",buffer);
            }
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
            }
            fclose(fpw);
            fclose(fp);
            free(map);free(map_unique);
            system("mv topo.out.bak topo.out");
            //
            
            sprintf(file_path,"%s/%s",current_folder,"topo.log");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Dihedrals\n")==0)break;
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&int_buffer,&topo_dihedral_types);
            map=(int*)malloc(topo_dihedral_types*sizeof(int));
            for(j=0;j<topo_dihedral_types;++j)
            {
                fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s\t%s\t%s\t%s\t#%d",word,word,word,word,word,&map[j]);
            }
            fclose(fp);
            
            map_unique=(int*)malloc(topo_dihedral_types*sizeof(int));
            fu_out=(int**)malloc(topo_dihedral_types*sizeof(int*));
            for(j=0;j<topo_dihedral_types;++j)fu_out[j]=(int*)malloc(3*sizeof(int));
            find_unique_1D(topo_dihedral_types,map,fu_out);
            k=-1;
            for(j=0;j<topo_dihedral_types;++j)if(fu_out[j][2]!=-1){k=k+1;map_unique[k]=map[fu_out[j][0]-1];}
            for(j=0;j<topo_dihedral_types;++j)free(fu_out[j]);free(fu_out);
            
            sprintf(file_path,"%s/%s",current_folder,"topo.out");
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            sprintf(file_path,"%s/%s",current_folder,"topo.out.bak");
            fpw=fopen(file_path,"w+");
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&bonds);
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&angles);
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);sscanf(buffer,"%d",&dihedrals);
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
                if(strcmp(buffer,"Dihedrals\n")==0)break;
            }
            fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
            for(j=0;j<dihedrals;++j)
            {
                fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d",&int_buffer,&c0,&c1,&c2,&c3,&c4);
                fprintf(fpw,"%d\t%d\t%d\t%d\t%d\t%d\n",int_buffer,map_unique[c0-1],c1,c2,c3,c4);
                //fprintf(fpw,"%s",buffer);
            }
            while(fgets(buffer,cmax_length,fp)!=NULL){
                fprintf(fpw,"%s",buffer);
            }
            fclose(fpw);
            fclose(fp);
            free(map);free(map_unique);
            system("mv topo.out.bak topo.out");
            //
        }
    }
    //
    sprintf(file_path,"%s/%s.in",current_folder,mol2_filename);
    fpw=fopen(file_path,"w+");
    sprintf(file_path,"%s/topo.species",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    
    fprintf(fpw,"\nvariable sim_name index %s\n",name);
    fprintf(fpw,"variable E_tol equal 1.0e-14\nvariable F_tol equal 1.0e-12\nvariable N_iter equal 1000000\nvariable N_eval equal 2500000\n");
    fprintf(fpw,"variable console_step_min equal 1\n\n");
    fprintf(fpw,"units real\natom_style %s\nboundary %s %s %s\nread_data %s.lammps\n\n",typestr,boundx,boundy,boundz,mol2_filename);
    
    atom_types=0;
    fprintf(fpw,"# atom types\n");
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        atom_types=atom_types+1;
        fprintf(fpw,"# %d --> %s",atom_types,buffer);
    }
    rewind(fp);
    //fprintf(fpw,"group my_atoms type ");
    //for(j=0;j<atom_types;++j)fprintf(fpw,"%d ",j+1);fprintf(fpw,"\n");
    fclose(fp);
    
    fprintf(fpw,"\n");
    
    fprintf(fpw,"# bonds\n");
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"bond_style")==0)
        {
            fprintf(fpw,"%s",juice[i]);
            sscanf(juice[i],"%s\t%s",word,buffer);
            if(strcmp(buffer,"hybrid")!=0)
                hybrid=0;
            else
                hybrid=1;
            break;
        }
    }
    
    sprintf(file_path,"%s/%s",current_folder,"topo.log");
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Bonds\n")==0)break;
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&topo_bond_types);
    fclose(fp);
    
    search_vector=(int*)malloc(topo_bond_types*sizeof(int));
    for(i=0;i<topo_bond_types;++i)search_vector[i]=0;
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"bond")==0)
        {
            sscanf(juice[i],"%s\t%d\t%s\t%s\t%s",word,&c0,B1_species,B2_species,buffer);
            if(search_vector[c0-1]==0)
            {
                search_vector[c0-1]=1;
                if(strcmp(buffer,"harmonic")==0)
                {
                    sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%lf\t%lf",word,&c0,B1_species,B2_species,buffer,&d1,&d2);
                    fprintf(fpw,"# %s-%s\n",B1_species,B2_species);
                    fprintf(fpw,"bond_coeff\t%d\t%lf\t%lf\n",c0,d2,d1);
                }
            }
        }
    }
    
    free(search_vector);
    
    //
    
    fprintf(fpw,"\n");
    
    fprintf(fpw,"# angles\n");
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"angle_style")==0)
        {
            fprintf(fpw,"%s",juice[i]);
            sscanf(juice[i],"%s\t%s",word,buffer);
            if(strcmp(buffer,"hybrid")!=0)
                hybrid=0;
            else
                hybrid=1;
            break;
        }
    }
    
    sprintf(file_path,"%s/%s",current_folder,"topo.log");
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Angles\n")==0)break;
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&topo_angle_types);
    fclose(fp);
    
    search_vector=(int*)malloc(topo_angle_types*sizeof(int));
    for(i=0;i<topo_angle_types;++i)search_vector[i]=0;
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"angle")==0)
        {
            sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s",word,&c0,A1_species,A2_species,A3_species,buffer);
            if(search_vector[c0-1]==0)
            {
                search_vector[c0-1]=1;
                if(strcmp(buffer,"harmonic")==0)
                {
                    sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%lf\t%lf",word,&c0,A1_species,A2_species,A3_species,buffer,&d1,&d2);
                    fprintf(fpw,"# %s-%s-%s\n",A1_species,A2_species,A3_species);
                    fprintf(fpw,"angle_coeff\t%d\t%lf\t%lf\n",c0,d2,d1);
                }
            }
        }
    }
    
    free(search_vector);
    
    //
    
    fprintf(fpw,"\n");
    
    fprintf(fpw,"# dihedrals\n");
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"dihedral_style")==0)
        {
            fprintf(fpw,"%s",juice[i]);
            sscanf(juice[i],"%s\t%s",word,buffer);
            if(strcmp(buffer,"hybrid")!=0)
                hybrid=0;
            else
                hybrid=1;
            break;
        }
    }
    
    sprintf(file_path,"%s/%s",current_folder,"topo.log");
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Dihedrals\n")==0)break;
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&topo_dihedral_types);
    fclose(fp);
    
    search_vector=(int*)malloc(topo_dihedral_types*sizeof(int));
    for(i=0;i<topo_dihedral_types;++i)search_vector[i]=0;
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"dihedral")==0)
        {
            sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer);
            if(search_vector[c0-1]==0)
            {
                search_vector[c0-1]=1;
                if(strcmp(buffer,"table")==0)
                {
                    sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,table1,table2);
                    fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                    if(hybrid==0)fprintf(fpw,"dihedral_coeff\t%d\t%s\t%s\n",c0,table1,table2);
                    if(hybrid==1)fprintf(fpw,"dihedral_coeff\t%d\t%s\t%s\t%s\n",c0,buffer,table1,table2);
                }
                if(strcmp(buffer,"opls")==0)
                {
                    sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,&K0,&K1,&K2,&K3);
                    fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                    if(hybrid==0)fprintf(fpw,"dihedral_coeff\t%d\t%lf\t%lf\t%lf\t%lf\n",c0,K0,K1,K2,K3);
                    if(hybrid==1)fprintf(fpw,"dihedral_coeff\t%d\t%s\t%lf\t%lf\t%lf\t%lf\n",c0,buffer,K0,K1,K2,K3);
                }
                if(strcmp(buffer,"multi/harmonic")==0)
                {
                    sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lf",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,&K0,&K1,&K2,&K3,&K4);
                    fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                    if(hybrid==0)fprintf(fpw,"dihedral_coeff\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",c0,K0,K1,K2,K3,K4);
                    if(hybrid==1)fprintf(fpw,"dihedral_coeff\t%d\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\n",c0,buffer,K0,K1,K2,K3,K4);
                }
		if(strcmp(buffer,"fourier")==0)
                {
                    sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,&fourier_count);
                    fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                    fprintf(fpw,"dihedral_coeff\t%d\t%d\t",c0,fourier_count);

                    tokenize(juice[i],&ignore_out,&elements_out);
                    
                    W=elements_out;
                    words=(char**)malloc(W*sizeof(char*));
                    for(j=0;j<W;++j)words[j]=(char*)malloc(cmax_length*sizeof(char));
                    for(j=0;j<W;++j)
                        select_from_tok(juice[i],elements_out,j,words[j]);
                    
                    l=7;
                    for(j=0;j<fourier_count;++j)
                    {
                            l=l+1;
                            sscanf(words[l],"%lf",&fourier_K);
                            l=l+1;
                            sscanf(words[l],"%lf",&fourier_d);
                            l=l+1;
                            sscanf(words[l],"%d",&fourier_n);
                        fprintf(fpw,"%lf\t%d\t%lf\t",fourier_K,fourier_n,fourier_d);
                    }
                    fprintf(fpw,"\n");
                    
                    for(j=0;j<W;++j)free(words[j]);free(words);                    
                }
            }
        }
    }
    
    free(search_vector);
    
    //
    
    
    
    //
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"improper_style")==0)
        {
            iflag=1;
            fprintf(fpw,"\n");
            fprintf(fpw,"# impropers\n");
            fprintf(fpw,"%s",juice[i]);
            sscanf(juice[i],"%s\t%s",word,buffer);
            if(strcmp(buffer,"hybrid")!=0)
                hybrid=0;
            else
                hybrid=1;
            break;
        }
    }
    
    if(iflag==1)
    {
        
        improper_types=0;
        for(i=0;i<lines;++i){
            sprintf(word,"%s","");sscanf(juice[i],"%s",word);
            if(strcmp(word,"improper")==0)
            {
                improper_types=improper_types+1;
            }
        }
        
        search_vector=(int*)malloc(improper_types*sizeof(int));
        for(i=0;i<improper_types;++i)search_vector[i]=0;
        
        for(i=0;i<lines;++i){
            sprintf(word,"%s","");sscanf(juice[i],"%s",word);
            if(strcmp(word,"improper")==0)
            {
                sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer);
                if(search_vector[c0-1]==0)
                {
                    search_vector[c0-1]=1;
                    if(strcmp(buffer,"ring")==0)
                    {
                        sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,&d1,&d2,&d3);
                        fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                        fprintf(fpw,"improper_coeff\t%d\t%lf\t%lf\t%lf\n",c0,d1,d2,d3);
                    }
                    if(strcmp(buffer,"harmonic")==0)
                    {
                        sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,&d1,&d2);
                        fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                        fprintf(fpw,"improper_coeff\t%d\t%lf\t%lf\n",c0,d1,d2);
                    }
                    if(strcmp(buffer,"cvff")==0)
                    {
                        sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf",word,&c0,D1_species,D2_species,D3_species,D4_species,buffer,&d1,&d2,&d3);
                        fprintf(fpw,"# %s-%s-%s-%s\n",D1_species,D2_species,D3_species,D4_species);
                        fprintf(fpw,"improper_coeff\t%d\t%lf\t%d\t%d\n",c0,d1,(int)d2,(int)d3);
                    }
                }
            }
        }
        
        free(search_vector);
        
        //sprintf(file_path,"%s/%s",current_folder,"result.mol2");
        sprintf(file_path,"%s/%s.mol2",current_folder,name);
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}

        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms,&bonds);
        mol2_species=(char**)malloc(atoms*sizeof(char*));
        for(i=0;i<atoms;++i)mol2_species[i]=(char*)malloc(sub_length*sizeof(char));
        B1=(int*)malloc(bonds*sizeof(int));
        B2=(int*)malloc(bonds*sizeof(int));
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
        for(i=0;i<atoms;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%s",&int_buffer,mol2_species[i]);
        }
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>BOND\n")==0)break;
        for(i=0;i<bonds;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d",&int_buffer,&B1[i],&B2[i]);
        }
        fclose(fp);
        
        max_1_2=one_two_init(B1,B2,atoms,bonds);
        one_two=(int**)malloc(atoms*sizeof(int*));
        for(i=0;i<atoms;++i)one_two[i]=(int*)malloc(max_1_2*sizeof(int));
        one_two_build(one_two,B1,B2,atoms,bonds,max_1_2);
        
        I1_type=(char**)malloc(improper_types*sizeof(char*));for(i=0;i<improper_types;++i)I1_type[i]=(char*)malloc(sub_length*sizeof(char));
        I2_type=(char**)malloc(improper_types*sizeof(char*));for(i=0;i<improper_types;++i)I2_type[i]=(char*)malloc(sub_length*sizeof(char));
        I3_type=(char**)malloc(improper_types*sizeof(char*));for(i=0;i<improper_types;++i)I3_type[i]=(char*)malloc(sub_length*sizeof(char));
        I4_type=(char**)malloc(improper_types*sizeof(char*));for(i=0;i<improper_types;++i)I4_type[i]=(char*)malloc(sub_length*sizeof(char));
        I_type_init=(int*)malloc(improper_types*sizeof(int));
        
        j=-1;
        for(i=0;i<lines;++i){
            sprintf(word,"%s","");sscanf(juice[i],"%s",word);
            if(strcmp(word,"improper")==0)
            {
                j=j+1;
                sscanf(juice[i],"%s\t%d\t%s\t%s\t%s\t%s",word,&I_type_init[j],I1_type[j],I2_type[j],I3_type[j],I4_type[j]);
                
            }
        }
        
        for(i=0;i<improper_types;++i)printf("%d\t%s\t%s\t%s\t%s\n",I_type_init[i],I1_type[i],I2_type[i],I3_type[i],I4_type[i]);
        
        improper_types_final=I_type_init[improper_types-1];
        
        impropers=0;
        for(i=0;i<atoms;++i)
        {
            if(one_two[i][0]!=0 && one_two[i][1]!=0 && one_two[i][2]!=0)
            {
                
                for(j=0;j<improper_types;++j)
                {
                    if(
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I4_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I3_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I4_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I2_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I3_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I2_type[j])==0)
                       )
                    {
                        impropers=impropers+1;
                        break;
                    }
                }
                
            }
        }
        
        sprintf(file_path,"%s/%s",current_folder,"impropers.out");
        fp=fopen(file_path,"w+");
        
        fprintf(fp,"%d impropers\n",impropers);
        fprintf(fp,"\n");
        fprintf(fp,"%d improper types\n",improper_types_final);
        fprintf(fp,"\n");
        fprintf(fp,"Impropers\n");
        fprintf(fp,"\n");
        
        k=0;
        for(i=0;i<atoms;++i)
        {
            if(one_two[i][0]!=0 && one_two[i][1]!=0 && one_two[i][2]!=0)
            {
                
                for(j=0;j<improper_types;++j)
                {
                    if(
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I4_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I3_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I4_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I2_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I2_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I3_type[j])==0)
                       ||
                       (strcmp(mol2_species[i],I1_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][0]-1],I4_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][1]-1],I3_type[j])==0 &&
                        strcmp(mol2_species[one_two[i][2]-1],I2_type[j])==0)
                       )
                    {
                        k=k+1;
                        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",k,I_type_init[j],i+1,one_two[i][0],one_two[i][1],one_two[i][2]);
                        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t\t\t#\t%s\t%s\t%s\t%s\n",k,I_type_init[j],one_two[i][0],i+1,one_two[i][1],one_two[i][2],
                        //        mol2_species[one_two[i][0]-1],mol2_species[i+1-1],mol2_species[one_two[i][1]-1],mol2_species[one_two[i][2]-1]);
                        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",k,I_type_init[j],one_two[i][0],i+1,one_two[i][1],one_two[i][2]);
                        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",k,I_type_init[j],i+1,one_two[i][0],one_two[i][1],one_two[i][2]);
                        break;
                    }
                }
                
            }
        }
        
        fclose(fp);
        
        //for(i=0;i<atoms;++i){for(j=0;j<max_1_2;++j)printf("%d\t",one_two[i][j]);printf("\n");}
        
        free(B1);free(B2);
        for(i=0;i<atoms;++i)free(one_two[i]);free(one_two);
        free(I_type_init);
        for(i=0;i<improper_types;++i)free(I1_type[i]);free(I1_type);
        for(i=0;i<improper_types;++i)free(I2_type[i]);free(I2_type);
        for(i=0;i<improper_types;++i)free(I3_type[i]);free(I3_type);
        for(i=0;i<improper_types;++i)free(I4_type[i]);free(I4_type);
        for(i=0;i<atoms;++i)free(mol2_species[i]);free(mol2_species);
    }
    
    //
    
    fprintf(fpw,"\n");
    
    fprintf(fpw,"# pair\n");
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"pair_style")==0)
        {
            fprintf(fpw,"%s",juice[i]);
            sscanf(juice[i],"%s\t%s",word,buffer);
            if(strcmp(buffer,"lj/cut/coul/cut")==0)
                pair_type=1;
            break;
        }
    }
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"LJ_amber")==0){pair_type=2;break;}
            }
    
    global_species=(char**)malloc(atom_types*sizeof(char*));
    for(i=0;i<atom_types;++i)global_species[i]=(char*)malloc(sub_length*sizeof(char));
    
    sprintf(file_path,"%s/topo.species",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    for(i=0;i<atom_types;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%s",global_species[i]);printf("from topo.species:\t%s\n",global_species[i]);
    }
    fclose(fp);
    
    if(pair_type==1)
    {
        for(i=0;i<atom_types;++i)
        {
            for(j=i;j<atom_types;++j)
            {
                
                for(k=0;k<lines;++k)
                {
                    sprintf(word,"%s","");sscanf(juice[k],"%s",word);
                    if(strcmp(word,"pair_coeff")==0)
                    {
                        sscanf(juice[k],"%s\t%s\t%s\t%lf\t%lf",word,B1_species,B2_species,&d1,&d2);
                        if((strcmp(B1_species,global_species[i])==0 && strcmp(B2_species,global_species[j])==0)||(strcmp(B1_species,global_species[j])==0 && strcmp(B2_species,global_species[i])==0))
                        {
                            fprintf(fpw,"pair_coeff %d %d %lf %lf # %s-%s\n",i+1,j+1,d1,d2,global_species[i],global_species[j]);
                            break;
                            
                        }
                    }
                }
            }
        }
    }
    
    if(pair_type==2)
    {
        
        for(i=0;i<atom_types;++i)
        {
            for(j=i;j<atom_types;++j)
            {
                // locate i
                for(k=0;k<lines;++k)
                {
                    sprintf(word,"%s","");sscanf(juice[k],"%s",word);
                    if(strcmp(word,"pair_coeff_amber")==0)
                    {
                        sscanf(juice[k],"%s\t%s\t%lf\t%lf",word,species,&d1,&d2);   // R* and epsilon for i
                        if(strcmp(species,global_species[i])==0)break;
                    }
                }
                // locate j
                for(k=0;k<lines;++k)
                {
                    sprintf(word,"%s","");sscanf(juice[k],"%s",word);
                    if(strcmp(word,"pair_coeff_amber")==0)
                    {
                        sscanf(juice[k],"%s\t%s\t%lf\t%lf",word,species,&d3,&d4);   // R* and epsilon for j
                        if(strcmp(species,global_species[j])==0)break;
                    }
                }
                // output
                fprintf(fpw,"pair_coeff %d %d %lf %lf # %s-%s\n",i+1,j+1,sqrt(d2*d4),(d1+d3)/ljscale,global_species[i],global_species[j]);
            }
        }
        
        
    }
    
    
    for(i=0;i<atom_types;++i)free(global_species[i]);free(global_species);
    
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"special")==0)
        {
            fprintf(fpw,"\n");
            fprintf(fpw,"%s",juice[i+1]);
            break;
        }
    }
    
    /*
     # write preview
     dump xyz_coords all xyz 1 animation_CG.xyz
     dump_modify xyz_coords element C H
     */
    fprintf(fpw,"\n# write preview\n");
    fprintf(fpw,"dump xyz_coords all xyz 1 %s.xyz\n",name);
    fprintf(fpw,"dump_modify xyz_coords element");
    sprintf(file_path,"%s/init.lammps",current_folder);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Masses\n")==0)break;
    fgets(buffer,cmax_length,fp);
    for(i=0;i<atom_types;++i){
        fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%s\t%lf",word,&d1);
        resolve_species((int)floor(d1),species);
        fprintf(fpw," %s",species);
    }
    fclose(fp);
    fprintf(fpw,"\ndump_modify xyz_coords sort id\n\n");
    
    /*
     # -- minimization ------------------------------------------------------
     thermo ${console_step_min}
     thermo_style custom step pe emol evdwl ecoul elong ebond eangle edihed eimp
     min_style cg
     minimize ${E_tol} ${F_tol} ${N_iter} ${N_eval}
     reset_timestep 0
     
     # write restart
     write_restart resume_${sim_name}
     */
    fprintf(fpw,"# -- minimization ------------------------------------------------------\n");
    fprintf(fpw,"thermo ${console_step_min}\n");
    fprintf(fpw,"thermo_style custom step pe emol evdwl ecoul elong ebond eangle edihed eimp\n");
    fprintf(fpw,"min_style cg\n");
    fprintf(fpw,"minimize ${E_tol} ${F_tol} ${N_iter} ${N_eval}\n");
    fprintf(fpw,"reset_timestep 0\n\n");
    
    fprintf(fpw,"# write restart\n");
    fprintf(fpw,"write_restart resume_${sim_name}\n\n");
    
    
    fclose(fpw);
    
    //
    /*
    system("head -3 init.lammps > result.lammps");
    system("head -3 topo.out >> result.lammps");
    if(impropers>0)system("head -1 impropers.out >> result.lammps");
    system("echo >> result.lammps");
    system("sed -n '5,5p' init.lammps >> result.lammps");
    system("sed -n '5,7p' topo.out >> result.lammps");
    if(impropers>0)system("sed -n '3,3p' impropers.out >> result.lammps");
    system("echo >> result.lammps");
    system("sed -n '7,$p' init.lammps >> result.lammps");
    system("sed -n '8,$p' topo.out >> result.lammps");
    if(impropers>0)system("sed -n '4,$p' impropers.out >> result.lammps");
    */
    
    sprintf(command,"head -3 init.lammps > %s.lammps",name);system(command);
    sprintf(command,"head -3 topo.out >> %s.lammps",name);system(command);
    if(impropers>0){
        sprintf(command,"head -1 impropers.out >> %s.lammps",name);system(command);
        }
    sprintf(command,"echo >> %s.lammps",name);system(command);
    sprintf(command,"sed -n '5,5p' init.lammps >> %s.lammps",name);system(command);
    sprintf(command,"sed -n '5,7p' topo.out >> %s.lammps",name);system(command);
    if(impropers>0){
        sprintf(command,"sed -n '3,3p' impropers.out >> %s.lammps",name);system(command);
        }
    sprintf(command,"echo >> %s.lammps",name);system(command);
    sprintf(command,"sed -n '7,$p' init.lammps >> %s.lammps",name);system(command);
    sprintf(command,"sed -n '8,$p' topo.out >> %s.lammps",name);system(command);
    if(impropers>0){
        sprintf(command,"sed -n '4,$p' impropers.out >> %s.lammps",name);system(command);
        }
    
    // add tri data
    if(strcmp(cell_type,"ortho")!=0){
    sprintf(file_path,"%s/%s.lammps",current_folder,name);fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    sprintf(file_path,"%s/%s.lammps.bak",current_folder,name);fpw=fopen(file_path,"w+");
    if(impropers>0)
    {
        for(i=0;i<17;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
        fprintf(fpw,"%lf %lf %lf xy xz yz\n",xy,xz,yz);
        while(fgets(buffer,cmax_length,fp)!=NULL)fprintf(fpw,"%s",buffer);
    }
    else
    {
        for(i=0;i<15;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
        fprintf(fpw,"%lf %lf %lf xy xz yz\n",xy,xz,yz);
        while(fgets(buffer,cmax_length,fp)!=NULL)fprintf(fpw,"%s",buffer);
    }
    fclose(fpw);
    fclose(fp);
    sprintf(command,"mv %s.lammps.bak %s.lammps",name,name);system(command);
    }
    
    //
    /*
    for(i=0;i<lines;++i){
        sprintf(word,"%s","");sscanf(juice[i],"%s",word);
        if(strcmp(word,"one_five")==0)
        {
            one_five=1;
            break;
        }
    }
    
    if(one_five==1)
    {
        sprintf(file_path,"%s/%s",current_folder,"result.lammps");
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
        B1=(int*)malloc(bonds*sizeof(int));
        B2=(int*)malloc(bonds*sizeof(int));
        D1=(int*)malloc(dihedrals*sizeof(int));
        D2=(int*)malloc(dihedrals*sizeof(int));
        D3=(int*)malloc(dihedrals*sizeof(int));
        D4=(int*)malloc(dihedrals*sizeof(int));
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Bonds\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<bonds;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&int_buffer,&B1[i],&B2[i]);
        }
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"Dihedrals\n")==0)break;
        fgets(buffer,cmax_length,fp);
        for(i=0;i<dihedrals;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d",&int_buffer,&int_buffer,&D1[i],&D2[i],&D3[i],&D4[i]);
        }
        fclose(fp);
        
        ////
        
        k=0;
        for(i=0;i<dihedrals;++i)
        {
            four=D4[i];
            for(j=0;j<bonds;++j)
            {
                if(B1[j]==four && B2[j]!=D3[i]){k=k+1;}
                if(B2[j]==four && B1[j]!=D3[i]){k=k+1;}
            }
            one=D1[i];
            for(j=0;j<bonds;++j)
            {
                if(B1[j]==one && B2[j]!=D2[i]){k=k+1;}
                if(B2[j]==one && B1[j]!=D2[i]){k=k+1;}
            }
        }
        onefive0=k;
        
        myonefiveL0 = (int*)malloc(onefive0*sizeof(int));
        myonefiveR0 = (int*)malloc(onefive0*sizeof(int));
        
        k=0;
        for(i=0;i<dihedrals;++i)
        {
            one=D1[i];four=D4[i];
            for(j=0;j<bonds;++j)
            {
                if(B1[j]==four && B2[j]!=D3[i]){k=k+1;five=B2[j];myonefiveL0[k-1]=one;myonefiveR0[k-1]=five;}
                if(B2[j]==four && B1[j]!=D3[i]){k=k+1;five=B1[j];myonefiveL0[k-1]=one;myonefiveR0[k-1]=five;}
            }
            one=D4[i];four=D1[i];
            for(j=0;j<bonds;++j)
            {
                if(B1[j]==four && B2[j]!=D2[i]){k=k+1;five=B2[j];myonefiveL0[k-1]=one;myonefiveR0[k-1]=five;}
                if(B2[j]==four && B1[j]!=D2[i]){k=k+1;five=B1[j];myonefiveL0[k-1]=one;myonefiveR0[k-1]=five;}
            }
        }
        
        for(i=0;i<onefive0-1;++i)
        {
            for(j=i+1;j<onefive0;++j)
            {
                if((myonefiveL0[i]==myonefiveL0[j] && myonefiveR0[i]==myonefiveR0[j])||(myonefiveL0[i]==myonefiveR0[j] && myonefiveR0[i]==myonefiveL0[j]))
                {myonefiveL0[j]=-1;myonefiveR0[j]=-1;}
            }
        }
        k=0;
        for(i=0;i<onefive0;++i)
            if(myonefiveL0[i]!=-1)
            {k=k+1;}
        onefive=k;
        
        myonefiveL = (int*)malloc(onefive*sizeof(int));
        myonefiveR = (int*)malloc(onefive*sizeof(int));
        
        k=-1;
        for(i=0;i<onefive0;++i)
            if(myonefiveL0[i]!=-1)
            {k=k+1;myonefiveL[k]=myonefiveL0[i];myonefiveR[k]=myonefiveR0[i];}
        
        //for(i=0;i<onefive;++i)printf("~~[%d]\t%d\t%d\n",i+1,myonefiveL[i],myonefiveR[i]);

        one_five_types=0;
        for(i=0;i<lines;++i){
            sprintf(word,"%s","");sscanf(juice[i],"%s",word);
            if(strcmp(word,"one_five_coeff")==0)
            {
                one_five_types=one_five_types+1;
            }
        }
        
        one_five1_type=(char**)malloc(one_five_types*sizeof(char*));for(i=0;i<one_five_types;++i)one_five1_type[i]=(char*)malloc(sub_length*sizeof(char));
        one_five2_type=(char**)malloc(one_five_types*sizeof(char*));for(i=0;i<one_five_types;++i)one_five2_type[i]=(char*)malloc(sub_length*sizeof(char));
        one_five_type_init=(int*)malloc(one_five_types*sizeof(int));
        
        j=-1;
        for(i=0;i<lines;++i){
            sprintf(word,"%s","");sscanf(juice[i],"%s",word);
            if(strcmp(word,"one_five_coeff")==0)
            {
                j=j+1;
                sscanf(juice[i],"%s\t%d\t%s\t%s",word,&one_five_type_init[j],one_five1_type[j],one_five2_type[j]);
                
            }
        }
        
        for(i=0;i<one_five_types;++i)printf("%d\t%s\t%s\n",one_five_type_init[i],one_five1_type[i],one_five2_type[i]);
        
        one_five_types_final=one_five_type_init[one_five_types-1];
        
        sprintf(file_path,"%s/%s",current_folder,"result.mol2");
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms,&bonds);
        mol2_species=(char**)malloc(atoms*sizeof(char*));
        for(i=0;i<atoms;++i)mol2_species[i]=(char*)malloc(sub_length*sizeof(char));
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
        for(i=0;i<atoms;++i)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%s",&int_buffer,mol2_species[i]);
        }
        fclose(fp);
        
        one_five_final=0;
        for(i=0;i<onefive;++i)
        {
            for(j=0;j<one_five_types;++j)
            {
                if((strcmp(mol2_species[myonefiveL[i]-1],one_five1_type[j])==0 && strcmp(mol2_species[myonefiveR[i]-1],one_five2_type[j])==0)
                   ||
                   (strcmp(mol2_species[myonefiveL[i]-1],one_five2_type[j])==0 && strcmp(mol2_species[myonefiveR[i]-1],one_five1_type[j])==0)
                   )
                {
                    one_five_final=one_five_final+1;
                    break;
                }
            }
        }
        
        sprintf(file_path,"%s/%s",current_folder,"one_five.out");
        fp=fopen(file_path,"w+");
        
        fprintf(fp,"%d one five pairs\n\n",one_five_final);
        fprintf(fp,"%d one five pair types\n\n",one_five_types_final);
        fprintf(fp,"Pairs\n\n");
        
        k=0;
        for(i=0;i<onefive;++i)
        {
            for(j=0;j<one_five_types;++j)
            {
                if((strcmp(mol2_species[myonefiveL[i]-1],one_five1_type[j])==0 && strcmp(mol2_species[myonefiveR[i]-1],one_five2_type[j])==0)
                   ||
                   (strcmp(mol2_species[myonefiveL[i]-1],one_five2_type[j])==0 && strcmp(mol2_species[myonefiveR[i]-1],one_five1_type[j])==0)
                   )
                {
                    k=k+1;
                    fprintf(fp,"%d\t%d\t%d\t%d\t\t\t#\t%s-%s\n",k,one_five_type_init[j],myonefiveL[i],myonefiveR[i],mol2_species[myonefiveL[i]-1],mol2_species[myonefiveR[i]-1]);
                    break;
                }
            }
        }
        
        fclose(fp);
        
        ////
        for(i=0;i<one_five_types;++i)free(one_five1_type[i]);free(one_five1_type);
        for(i=0;i<one_five_types;++i)free(one_five2_type[i]);free(one_five2_type);
        free(one_five_type_init);
        free(myonefiveL0);free(myonefiveR0);
        free(myonefiveL);free(myonefiveR);
        free(B1);free(B2);free(D1);free(D2);free(D3);free(D4);
        for(i=0;i<atoms;++i)free(mol2_species[i]);free(mol2_species);
        
        //
        
        sprintf(file_path,"%s/%s",current_folder,"result.lammps");
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        
        sprintf(file_path,"%s/%s",current_folder,"result_one_five.lammps");
        fpw=fopen(file_path,"w+");
        
        
        
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
        fprintf(fpw,"%d bonds\n",bonds+one_five_final);
        
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        if(impropers>0){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bond_types);
        fprintf(fpw,"%d bond types\n",bond_types+one_five_types_final);
        
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            fprintf(fpw,"%s",buffer);
            if(strcmp(buffer,"Bonds\n")==0)break;
        }
        fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);
        for(i=0;i<bonds;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
        
        sprintf(file_path,"%s/%s",current_folder,"one_five.out");
        fpr=fopen(file_path,"r");
        if(fpr==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fpr)!=NULL)
            if(strcmp(buffer,"Pairs\n")==0)
                break;
        fgets(buffer,cmax_length,fpr);
        for(i=0;i<one_five_final;++i)
        {
            fgets(buffer,cmax_length,fpr);sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&c0,&c1,&c2);
            fprintf(fpw,"%d\t%d\t%d\t%d\n",bonds+int_buffer,bond_types+c0,c1,c2);
        }
        fclose(fpr);
        
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            fprintf(fpw,"%s",buffer);
        }
        
        fclose(fpw);
        
        fclose(fp);
        
        //
        
    }
    */
    
    //
    for(i=0;i<lines;++i)free(juice[i]);free(juice);
    
    //
    return 0;
}

void tokenize(char *buffer,int *ignore_out,int *elements_out);
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
void select_from_tok(char *tok,int N_toks,int target,char *res);
void select_from_tok(char *tok,int N_toks,int target,char *res)
{
    char word[cmax_length];
    int *pos,j,k,w0,w1;
    
    pos=(int*)malloc(N_toks*sizeof(int));
    k=-1;
    for(j=0;j<cmax_length;++j)
    {
        if(tok[j]=='\0')break;
        if(tok[j]=='$'){k=k+1;pos[k]=j;}
    }
    if(target==0)
        w0=0;
    else
        w0=pos[target-1]+1;
    w1=pos[target]-1;
    k=-1;
    for(j=w0;j<=w1;++j)
    {
        k=k+1;
        word[k]=tok[j];
    }
    word[k+1]='\0';
    sprintf(res,"%s",word);
    free(pos);

}
void inline_from_select(char **words,int W,char *res);
void inline_from_select(char **words,int W,char *res)
{
    char word[cmax_length];
    int j,k,l;
    
    l=-1;
    for(j=0;j<W;++j)
    {
        for(k=0;k<cmax_length;++k)
        {
            if(words[j][k]=='\0')break;
            l=l+1;
            word[l]=words[j][k];
        }
        l=l+1;
        word[l]=' ';
    }
    l=l+1;
    word[l]='\0';
    sprintf(res,"%s",word);
}

void find_unique_1D(int N, int *array,int **out);
void find_unique_1D(int N, int *array,int **out)
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
            if(array[i]==array[j])final[j]=-1;
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
                    if(array[i]==array[j])
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
int one_two_init(int *B1, int *B2, int atoms, int bonds);
int one_two_init(int *B1, int *B2, int atoms, int bonds)
{
    int i,max_1_2,*one_two_counter_array;
    one_two_counter_array=(int*)malloc(atoms*sizeof(int));
    for(i=0;i<atoms;++i)one_two_counter_array[i]=0;
    for(i=0;i<bonds;++i)
    {
        one_two_counter_array[B1[i]-1]=one_two_counter_array[B1[i]-1]+1;
        one_two_counter_array[B2[i]-1]=one_two_counter_array[B2[i]-1]+1;
    }
    max_1_2=1;
    for(i=0;i<atoms;++i)
        if(one_two_counter_array[i]>max_1_2)max_1_2=one_two_counter_array[i];
    free(one_two_counter_array);
    return max_1_2;
}

void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);
void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2)
{
    int i,j;
    for(i=0;i<atoms;++i)
        for(j=0;j<max_1_2;++j)
            one_two[i][j]=0;
    for(i=0;i<bonds;++i)
    {
        j=0;
        while(one_two[B1[i]-1][j]!=0)
            j=j+1;
        one_two[B1[i]-1][j]=B2[i];
        j=0;
        while(one_two[B2[i]-1][j]!=0)
            j=j+1;
        one_two[B2[i]-1][j]=B1[i];
    }
}
void resolve_species(int atomic_M,char *species);
void resolve_species(int atomic_M,char *species)
{
    if(atomic_M==1)sprintf(species,"%s","H");
    if(atomic_M==12)sprintf(species,"%s","C");
    if(atomic_M==14)sprintf(species,"%s","N");
    if(atomic_M==16)sprintf(species,"%s","O");
}

// alterA
void alterA(int argc,char *argv);
void alterA(int argc,char *argv)
{
    FILE *fp;
    
    // fp related
    char current_folder[cmax_length],file_path[cmax_length];
    
    // dummies
    char *fcheck;
    char buffer[cmax_length],temp_species[species_length];
    int int_buffer;
    
    //
    int i,j,k;
    int bonds,angles,dihedrals,length;
    int bonds_N,angles_N,dihedrals_N;
    int B_lines,A_lines,D_lines,lines;
    int flag,temp;
    int rows,cols,unique_count;
    int *sites,*init,*final;
    int *store1,*store2;
    int *init0;
    int *unique,*unique_ID;
    int *bonds_col_1,*bonds_col_2,*bonds_col_3;
    int *angles_col_1,*angles_col_2,*angles_col_3,*angles_col_4;
    int *dihedrals_col_1,*dihedrals_col_2,*dihedrals_col_3,*dihedrals_col_4,*dihedrals_col_5;
    int **info;
    char **B1,**B2;
    char **A1,**A2,**A3;
    char **D1,**D2,**D3,**D4;
    
    //
    int *B_map;
    int *A_map;
    int *D_map;
    
    //------------------------------------------------------------------
    // print exe info
    if(argc==1)
    {
        printf(
               "\n*Reads topo.out and topo.log and unifies the selected angles\n"
               "*The rules are read from <filename>\n"
               "*Generates topo.out.altered and topo.log.altered\n\n"
               "./alterA <filename>\n\n"
               );
        exit(-1);
    }
    
    // working directory
    fcheck=getcwd(current_folder,cmax_length);
    
    // read the rules from file and store in 2d matrix
    sprintf(file_path,"%s/%s",current_folder,argv);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    // find the number of cols
    rows=0;cols=2;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(buffer[0]=='[')
        {
            rows=rows+1;
            sscanf(buffer,"[%d]",&int_buffer);
            if(int_buffer>cols){cols=int_buffer;}
        }
    }
    
    rewind(fp);
    
    // matrix to store the alteration rules
    info=(int**)malloc(rows*sizeof(int*));
    for(i=0;i<rows;++i){info[i]=(int*)malloc(cols*sizeof(int));}
    // initialization
    for(i=0;i<rows;++i){for(j=0;j<cols;++j){info[i][j]=0;}}
    // populate
    i=-1;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(buffer[0]=='[')
        {
            i=i+1;
            sscanf(buffer,"[%d]",&int_buffer);
            for(j=0;j<int_buffer;++j)
            {
                fcheck=fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d",&info[i][j]);
            }
        }
    }
    fclose(fp);
    
    // open the log file created by topo
    sprintf(file_path,"%s/topo.log",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&bonds,&B_lines);			// read bonds
    
    //init=(int*)malloc(B_lines*sizeof(int));
    //final=(int*)malloc(B_lines*sizeof(int));
    store1=(int*)malloc(B_lines*sizeof(int));
    
    B1=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B1[i]=(char*)malloc(species_length*sizeof(char));}
    B2=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B2[i]=(char*)malloc(species_length*sizeof(char));}
    B_map=(int*)malloc(B_lines*sizeof(int));
    //for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s",&init[i],B1[i],B2[i]);}
    for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t#%d",&store1[i],B1[i],B2[i],&B_map[i]);}
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&angles,&A_lines);			// read angles
    
    init=(int*)malloc(A_lines*sizeof(int));
    final=(int*)malloc(A_lines*sizeof(int));
    
    A1=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A1[i]=(char*)malloc(species_length*sizeof(char));}
    A2=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A2[i]=(char*)malloc(species_length*sizeof(char));}
    A3=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A3[i]=(char*)malloc(species_length*sizeof(char));}
    A_map=(int*)malloc(A_lines*sizeof(int));
    for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t#%d",&init[i],A1[i],A2[i],A3[i],&A_map[i]);}
    //for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s",&int_buffer,A1[i],A2[i],A3[i]);}
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&dihedrals,&D_lines);		// read dihedrals
    
    //init=(int*)malloc(D_lines*sizeof(int));
    //final=(int*)malloc(D_lines*sizeof(int));
    store2=(int*)malloc(D_lines*sizeof(int));
    
    D1=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D1[i]=(char*)malloc(species_length*sizeof(char));}
    D2=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D2[i]=(char*)malloc(species_length*sizeof(char));}
    D3=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D3[i]=(char*)malloc(species_length*sizeof(char));}
    D4=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D4[i]=(char*)malloc(species_length*sizeof(char));}
    D_map=(int*)malloc(D_lines*sizeof(int));
    //for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s",&init[i],D1[i],D2[i],D3[i],D4[i]);}
    for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s\t#%d",&store2[i],D1[i],D2[i],D3[i],D4[i],&D_map[i]);}
    
    fclose(fp);
    
    //lines=B_lines;
    lines=A_lines;
    //lines=D_lines;
    
    // save the initial numbering
    init0=(int*)malloc(lines*sizeof(int));
    for(i=0;i<lines;++i){init0[i]=init[i];}
    
    //----------------------------------------------------------------------
    // the first application
    
    // find the number of elements to unify
    length=0;
    for(i=0;i<cols;++i){if(info[0][i]!=0){length=length+1;}}
    sites=(int*)malloc(length*sizeof(int));
    // parse the info from info to sites
    for(i=0;i<length;++i){sites[i]=info[0][i];}
    
    // bubble sort the info
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<length;++i)
        {
            if(sites[i]<sites[i-1])
            {
                temp=sites[i-1];
                sites[i-1]=sites[i];
                sites[i]=temp;
                
                flag=1;
            }
        }
    }
    
    // initialize final array
    for(i=0;i<lines;++i){final[i]=init[i];}
    // alter entries
    for(i=0;i<lines;++i){for(j=0;j<length;++j){if(init[i]==sites[j]){final[i]=sites[0];}}}
    // bubble sort
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<lines;++i)
        {
            if(final[i]<final[i-1])
            {
                temp=final[i-1];
                final[i-1]=final[i];
                final[i]=temp;
                
                temp=init[i-1];
                init[i-1]=init[i];
                init[i]=temp;
                
                // initials
                temp=init0[i-1];
                init0[i-1]=init0[i];
                init0[i]=temp;
                
                //sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
                //sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
                //sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
                //sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);
                
                sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
                sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
                sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);
                temp=A_map[i-1];
                A_map[i-1]=A_map[i];
                A_map[i]=temp;
                
                //sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
                //sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);
                
                flag=1;
            }
        }
    }
    
    // normalize (the code is taken from the erase scripts)
    
    unique=(int*)malloc(lines*sizeof(int));
    unique_ID=(int*)malloc(lines*sizeof(int));
    
    for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
    
    unique_count=0;	// counter for the unique appearances
    for(i=0;i<lines;++i)
    {
        if(unique[i]!=0)
        {
            for(j=i+1;j<lines;++j)
            {
                if(unique[i]==unique[j]){unique[j]=0;}
            }
            unique_count=unique_count+1;
            unique_ID[i]=unique_count;
        }
    }
    
    for(i=0;i<lines;++i)
    {
        //if(final[i]!=0)
        //{
        for(j=0;j<lines;++j)
        {
            if(final[i]==unique[j])
            {
                final[i]=unique_ID[j];
            }
        }
        //}
    }
    
    free(sites);
    
    //----------------------------------------------------------------------
    
    // the rest applications of the rules
    
    for(k=1;k<rows;++k)
    {
        
        //----------------------------------------------------------------------
        // switch!
        for(i=0;i<lines;++i){init[i]=final[i];}
        
        length=0;
        for(i=0;i<cols;++i){if(info[k][i]!=0){length=length+1;}}
        
        sites=(int*)malloc(length*sizeof(int));
        
        for(i=0;i<length;++i){sites[i]=info[k][i];}
        
        for(i=0;i<lines;++i)
        {
            for(j=0;j<length;++j)
            {
                if(sites[j]==init0[i]){sites[j]=final[i];}
            }
        }
        
        // bubble sort
        flag=1;
        while(flag==1)
        {
            flag=0;
            for(i=1;i<length;++i)
            {
                if(sites[i]<sites[i-1])
                {
                    temp=sites[i-1];
                    sites[i-1]=sites[i];
                    sites[i]=temp;
                    
                    flag=1;
                }
            }
        }
        
        // initialize final array
        for(i=0;i<lines;++i){final[i]=init[i];}
        // alter entries
        for(i=0;i<lines;++i){for(j=0;j<length;++j){if(init[i]==sites[j]){final[i]=sites[0];}}}
        // bubble sort
        flag=1;
        while(flag==1)
        {
            flag=0;
            for(i=1;i<lines;++i)
            {
                if(final[i]<final[i-1])
                {
                    temp=final[i-1];
                    final[i-1]=final[i];
                    final[i]=temp;
                    
                    temp=init[i-1];
                    init[i-1]=init[i];
                    init[i]=temp;
                    
                    // initial
                    temp=init0[i-1];
                    init0[i-1]=init0[i];
                    init0[i]=temp;
                    
                    //sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);
                    
                    sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
                    sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
                    sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);
                    temp=A_map[i-1];
                    A_map[i-1]=A_map[i];
                    A_map[i]=temp;
                    
                    //sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);
                    
                    flag=1;
                }
            }
        }
        
        // normalize
        
        for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
        
        unique_count=0;	// counter for the unique appearances
        for(i=0;i<lines;++i)
        {
            if(unique[i]!=0)
            {
                for(j=i+1;j<lines;++j)
                {
                    if(unique[i]==unique[j]){unique[j]=0;}
                }
                unique_count=unique_count+1;
                unique_ID[i]=unique_count;
            }
        }
        
        for(i=0;i<lines;++i)
        {
            //if(final[i]!=0)
            //{
            for(j=0;j<lines;++j)
            {
                if(final[i]==unique[j])
                {
                    final[i]=unique_ID[j];
                }
            }
            //}
        }
        
        free(sites);
        
        //----------------------------------------------------------------------
        
    }
    
    // write altered log file
    sprintf(file_path,"%s/topo.log.altered",current_folder);
    fp=fopen(file_path,"w+");
    /*
     fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
     for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",i+1,B1[i],B2[i]);}
     fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
     for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",i+1,A1[i],A2[i],A3[i]);}
     fprintf(fp,"Dihedrals\n%d(%d) types\n",final[D_lines-1],D_lines);
     for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",final[i],D1[i],D2[i],D3[i],D4[i]);}
     */
    fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
    for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t#%d\n",store1[i],B1[i],B2[i],B_map[i]);}
    fprintf(fp,"Angles\n%d(%d) types\n",final[A_lines-1],A_lines);
    for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t#%d\n",final[i],A1[i],A2[i],A3[i],A_map[i]);}
    fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
    for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\t#%d\n",store2[i],D1[i],D2[i],D3[i],D4[i],D_map[i]);}
    /*
     fprintf(fp,"Bonds\n%d(%d) types\n",final[B_lines-1],B_lines);
     for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",final[i],B1[i],B2[i]);}
     fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
     for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",i+1,A1[i],A2[i],A3[i]);}
     fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
     for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",i+1,D1[i],D2[i],D3[i],D4[i]);}
     */
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    // print to console the switching info
    for(i=0;i<lines;++i){printf("%d\t%d\n",init0[i],final[i]);}
    
    // read topo.out
    sprintf(file_path,"%s/topo.out",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds_N);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles_N);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals_N);
    
    bonds_col_1=(int*)malloc(bonds_N*sizeof(int));
    bonds_col_2=(int*)malloc(bonds_N*sizeof(int));
    bonds_col_3=(int*)malloc(bonds_N*sizeof(int));
    angles_col_1=(int*)malloc(angles_N*sizeof(int));
    angles_col_2=(int*)malloc(angles_N*sizeof(int));
    angles_col_3=(int*)malloc(angles_N*sizeof(int));
    angles_col_4=(int*)malloc(angles_N*sizeof(int));
    dihedrals_col_1=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_2=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_3=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_4=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_5=(int*)malloc(dihedrals_N*sizeof(int));
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<bonds_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&bonds_col_1[i],&bonds_col_2[i],&bonds_col_3[i]);
    }
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<angles_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&int_buffer,&angles_col_1[i],&angles_col_2[i],&angles_col_3[i],&angles_col_4[i]);
    }
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<dihedrals_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d\n",&int_buffer,&dihedrals_col_1[i],&dihedrals_col_2[i],&dihedrals_col_3[i],&dihedrals_col_4[i],&dihedrals_col_5[i]);
    }
    
    fclose(fp);
    
    // write altered topo.out	
    sprintf(file_path,"%s/topo.out.altered",current_folder);
    fp=fopen(file_path,"w+");
    fprintf(fp,"%d bonds\n",bonds_N);
    fprintf(fp,"%d angles\n",angles_N);
    fprintf(fp,"%d dihedrals\n",dihedrals_N);
    fprintf(fp,"\n");
    
    fprintf(fp,"%d bond types\n",bonds);
    //fprintf(fp,"%d bond types\n",final[B_lines-1]);
    
    //fprintf(fp,"%d angle types\n",angles);
    fprintf(fp,"%d angle types\n",final[A_lines-1]);
    
    fprintf(fp,"%d dihedral types\n",dihedrals);
    //fprintf(fp,"%d dihedral types\n",final[D_lines-1]);
    
    fprintf(fp,"\nBonds\n");
    fprintf(fp,"\n");
    for(i=0;i<bonds_N;++i)
    {
        // switch
        //for(j=0;j<B_lines;++j){if(bonds_col_1[i]==init0[j]){temp=final[j];}}
        //fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,temp,bonds_col_2[i],bonds_col_3[i]);
        fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,bonds_col_1[i],bonds_col_2[i],bonds_col_3[i]);
    }
    fprintf(fp,"\nAngles\n");
    fprintf(fp,"\n");
    for(i=0;i<angles_N;++i)
    {
        // switch
        for(j=0;j<A_lines;++j){if(angles_col_1[i]==init0[j]){temp=final[j];}}
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,temp,angles_col_2[i],angles_col_3[i],angles_col_4[i]);
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,angles_col_1[i],angles_col_2[i],angles_col_3[i],angles_col_4[i]);
    }
    fprintf(fp,"\nDihedrals\n");
    fprintf(fp,"\n");
    for(i=0;i<dihedrals_N;++i)
    {
        // switch
        //for(j=0;j<D_lines;++j){if(dihedrals_col_1[i]==init0[j]){temp=final[j];}}
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,temp,dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,dihedrals_col_1[i],dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
    }
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    //------------------------------------------------------------------
    // free allocated memory	
    
    free(init);free(final);//free(sites);
    for(i=0;i<D_lines;++i){free(D1[i]);free(D2[i]);free(D3[i]);free(D4[i]);}free(D1);free(D2);free(D3);free(D4);
    for(i=0;i<A_lines;++i){free(A1[i]);free(A2[i]);free(A3[i]);}free(A1);free(A2);free(A3);
    for(i=0;i<B_lines;++i){free(B1[i]);free(B2[i]);}free(B1);free(B2);
    free(bonds_col_1);free(bonds_col_2);free(bonds_col_3);
    free(angles_col_1);free(angles_col_2);free(angles_col_3);free(angles_col_4);
    free(dihedrals_col_1);free(dihedrals_col_2);free(dihedrals_col_3);free(dihedrals_col_4);free(dihedrals_col_5);
    free(store1);free(store2);
    free(init0);
    for(i=0;i<rows;++i){free(info[i]);}free(info);
    free(unique);free(unique_ID);
    
    free(B_map);free(A_map);free(D_map);
    
}

void alterB(int argc,char *argv);
void alterB(int argc,char *argv)
{
    FILE *fp;
    
    // fp related
    char current_folder[cmax_length],file_path[cmax_length];
    
    // dummies
    char *fcheck;
    char buffer[cmax_length],temp_species[species_length];
    int int_buffer;
    
    //
    int i,j,k;
    int bonds,angles,dihedrals,length;
    int bonds_N,angles_N,dihedrals_N;
    int B_lines,A_lines,D_lines,lines;
    int flag,temp;
    int rows,cols,unique_count;
    int *sites,*init,*final;
    int *store1,*store2;
    int *init0;
    int *unique,*unique_ID;
    int *bonds_col_1,*bonds_col_2,*bonds_col_3;
    int *angles_col_1,*angles_col_2,*angles_col_3,*angles_col_4;
    int *dihedrals_col_1,*dihedrals_col_2,*dihedrals_col_3,*dihedrals_col_4,*dihedrals_col_5;
    int **info;
    char **B1,**B2;
    char **A1,**A2,**A3;
    char **D1,**D2,**D3,**D4;
    
    //
    int *B_map;
    int *A_map;
    int *D_map;
    
    //------------------------------------------------------------------
    // print exe info
    if(argc==1)
    {
        printf(
               "\n*Reads topo.out and topo.log and unifies the selected bonds\n"
               "*The rules are read from <filename>\n"
               "*Generates topo.out.altered and topo.log.altered\n\n"
               "./alterB <filename>\n\n"
               );
        exit(-1);
    }
    
    // working directory
    fcheck=getcwd(current_folder,cmax_length);
    
    // read the rules from file and store in 2d matrix
    sprintf(file_path,"%s/%s",current_folder,argv);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    // find the number of cols
    rows=0;cols=2;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(buffer[0]=='[')
        {
            rows=rows+1;
            sscanf(buffer,"[%d]",&int_buffer);
            if(int_buffer>cols){cols=int_buffer;}
        }
    }
    
    rewind(fp);
    
    // matrix to store the alteration rules
    info=(int**)malloc(rows*sizeof(int*));
    for(i=0;i<rows;++i){info[i]=(int*)malloc(cols*sizeof(int));}
    // initialization
    for(i=0;i<rows;++i){for(j=0;j<cols;++j){info[i][j]=0;}}
    // populate
    i=-1;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(buffer[0]=='[')
        {
            i=i+1;
            sscanf(buffer,"[%d]",&int_buffer);
            for(j=0;j<int_buffer;++j)
            {
                fcheck=fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d",&info[i][j]);
            }
        }
    }
    fclose(fp);
    
    // open the log file created by topo
    sprintf(file_path,"%s/topo.log",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&bonds,&B_lines);			// read bonds
    
    init=(int*)malloc(B_lines*sizeof(int));
    final=(int*)malloc(B_lines*sizeof(int));
    
    B1=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B1[i]=(char*)malloc(species_length*sizeof(char));}
    B2=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B2[i]=(char*)malloc(species_length*sizeof(char));}
    B_map=(int*)malloc(B_lines*sizeof(int));
    for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t#%d",&init[i],B1[i],B2[i],&B_map[i]);}
    //for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s",&store1[i],B1[i],B2[i]);}
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&angles,&A_lines);			// read angles
    
    //init=(int*)malloc(A_lines*sizeof(int));
    //final=(int*)malloc(A_lines*sizeof(int));
    store1=(int*)malloc(A_lines*sizeof(int));
    
    A1=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A1[i]=(char*)malloc(species_length*sizeof(char));}
    A2=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A2[i]=(char*)malloc(species_length*sizeof(char));}
    A3=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A3[i]=(char*)malloc(species_length*sizeof(char));}
    A_map=(int*)malloc(A_lines*sizeof(int));
    //for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s",&init[i],A1[i],A2[i],A3[i]);}
    for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t#%d",&store1[i],A1[i],A2[i],A3[i],&A_map[i]);}
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&dihedrals,&D_lines);		// read dihedrals
    
    //init=(int*)malloc(D_lines*sizeof(int));
    //final=(int*)malloc(D_lines*sizeof(int));
    store2=(int*)malloc(D_lines*sizeof(int));
    
    D1=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D1[i]=(char*)malloc(species_length*sizeof(char));}
    D2=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D2[i]=(char*)malloc(species_length*sizeof(char));}
    D3=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D3[i]=(char*)malloc(species_length*sizeof(char));}
    D4=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D4[i]=(char*)malloc(species_length*sizeof(char));}
    D_map=(int*)malloc(D_lines*sizeof(int));
    //for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s",&init[i],D1[i],D2[i],D3[i],D4[i]);}
    for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s\t#%d",&store2[i],D1[i],D2[i],D3[i],D4[i],&D_map[i]);}
    
    fclose(fp);
    
    lines=B_lines;
    //lines=A_lines;
    //lines=D_lines;
    
    // save the initial numbering
    init0=(int*)malloc(lines*sizeof(int));
    for(i=0;i<lines;++i){init0[i]=init[i];}
    
    //----------------------------------------------------------------------
    // the first application
    
    // find the number of elements to unify
    length=0;
    for(i=0;i<cols;++i){if(info[0][i]!=0){length=length+1;}}
    sites=(int*)malloc(length*sizeof(int));
    // parse the info from info to sites
    for(i=0;i<length;++i){sites[i]=info[0][i];}
    
    // bubble sort the info
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<length;++i)
        {
            if(sites[i]<sites[i-1])
            {
                temp=sites[i-1];
                sites[i-1]=sites[i];
                sites[i]=temp;
                
                flag=1;
            }
        }
    }
    
    // initialize final array
    for(i=0;i<lines;++i){final[i]=init[i];}
    // alter entries
    for(i=0;i<lines;++i){for(j=0;j<length;++j){if(init[i]==sites[j]){final[i]=sites[0];}}}
    // bubble sort
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<lines;++i)
        {
            if(final[i]<final[i-1])
            {
                temp=final[i-1];
                final[i-1]=final[i];
                final[i]=temp;
                
                temp=init[i-1];
                init[i-1]=init[i];
                init[i]=temp;
                
                // initials
                temp=init0[i-1];
                init0[i-1]=init0[i];
                init0[i]=temp;
                
                //sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
                //sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
                //sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
                //sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);
                
                //sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
                //sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
                //sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);
                
                sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
                sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);
                temp=B_map[i-1];
                B_map[i-1]=B_map[i];
                B_map[i]=temp;
                
                flag=1;
            }
        }
    }
    
    // normalize (the code is taken from the erase scripts)
    
    unique=(int*)malloc(lines*sizeof(int));
    unique_ID=(int*)malloc(lines*sizeof(int));
    
    for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
    
    unique_count=0;	// counter for the unique appearances
    for(i=0;i<lines;++i)
    {
        if(unique[i]!=0)
        {
            for(j=i+1;j<lines;++j)
            {
                if(unique[i]==unique[j]){unique[j]=0;}
            }
            unique_count=unique_count+1;
            unique_ID[i]=unique_count;
        }
    }
    
    for(i=0;i<lines;++i)
    {
        //if(final[i]!=0)
        //{
        for(j=0;j<lines;++j)
        {
            if(final[i]==unique[j])
            {
                final[i]=unique_ID[j];
            }
        }
        //}
    }
    
    free(sites);
    
    //----------------------------------------------------------------------
    
    // the rest applications of the rules
    
    for(k=1;k<rows;++k)
    {
        
        //----------------------------------------------------------------------
        // switch!
        for(i=0;i<lines;++i){init[i]=final[i];}
        
        length=0;
        for(i=0;i<cols;++i){if(info[k][i]!=0){length=length+1;}}
        
        sites=(int*)malloc(length*sizeof(int));
        
        for(i=0;i<length;++i){sites[i]=info[k][i];}
        
        for(i=0;i<lines;++i)
        {
            for(j=0;j<length;++j)
            {
                if(sites[j]==init0[i]){sites[j]=final[i];}
            }
        }
        
        // bubble sort
        flag=1;
        while(flag==1)
        {
            flag=0;
            for(i=1;i<length;++i)
            {
                if(sites[i]<sites[i-1])
                {
                    temp=sites[i-1];
                    sites[i-1]=sites[i];
                    sites[i]=temp;
                    
                    flag=1;
                }
            }
        }
        
        // initialize final array
        for(i=0;i<lines;++i){final[i]=init[i];}
        // alter entries
        for(i=0;i<lines;++i){for(j=0;j<length;++j){if(init[i]==sites[j]){final[i]=sites[0];}}}
        // bubble sort
        flag=1;
        while(flag==1)
        {
            flag=0;
            for(i=1;i<lines;++i)
            {
                if(final[i]<final[i-1])
                {
                    temp=final[i-1];
                    final[i-1]=final[i];
                    final[i]=temp;
                    
                    temp=init[i-1];
                    init[i-1]=init[i];
                    init[i]=temp;
                    
                    // initial
                    temp=init0[i-1];
                    init0[i-1]=init0[i];
                    init0[i]=temp;
                    
                    //sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);
                    
                    //sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);
                    
                    sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
                    sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);
                    temp=B_map[i-1];
                    B_map[i-1]=B_map[i];
                    B_map[i]=temp;
                    
                    flag=1;
                }
            }
        }
        
        // normalize
        
        for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
        
        unique_count=0;	// counter for the unique appearances
        for(i=0;i<lines;++i)
        {
            if(unique[i]!=0)
            {
                for(j=i+1;j<lines;++j)
                {
                    if(unique[i]==unique[j]){unique[j]=0;}
                }
                unique_count=unique_count+1;
                unique_ID[i]=unique_count;
            }
        }
        
        for(i=0;i<lines;++i)
        {
            //if(final[i]!=0)
            //{
            for(j=0;j<lines;++j)
            {
                if(final[i]==unique[j])
                {
                    final[i]=unique_ID[j];
                }
            }
            //}
        }
        
        free(sites);
        
        //----------------------------------------------------------------------
        
    }
    
    // write altered log file
    sprintf(file_path,"%s/topo.log.altered",current_folder);
    fp=fopen(file_path,"w+");
    /*
     fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
     for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",i+1,B1[i],B2[i]);}
     fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
     for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",i+1,A1[i],A2[i],A3[i]);}
     fprintf(fp,"Dihedrals\n%d(%d) types\n",final[D_lines-1],D_lines);
     for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",final[i],D1[i],D2[i],D3[i],D4[i]);}
     */
    /*
     fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
     for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",store1[i],B1[i],B2[i]);}
     fprintf(fp,"Angles\n%d(%d) types\n",final[A_lines-1],A_lines);
     for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",final[i],A1[i],A2[i],A3[i]);}
     fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
     for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",store2[i],D1[i],D2[i],D3[i],D4[i]);}
     */
    
    fprintf(fp,"Bonds\n%d(%d) types\n",final[B_lines-1],B_lines);
    for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t#%d\n",final[i],B1[i],B2[i],B_map[i]);}
    fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
    for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t#%d\n",store1[i],A1[i],A2[i],A3[i],A_map[i]);}
    fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
    for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\t#%d\n",store2[i],D1[i],D2[i],D3[i],D4[i],D_map[i]);}
    
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    // print to console the switching info
    for(i=0;i<lines;++i){printf("%d\t%d\n",init0[i],final[i]);}
    
    // read topo.out
    sprintf(file_path,"%s/topo.out",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds_N);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles_N);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals_N);
    
    bonds_col_1=(int*)malloc(bonds_N*sizeof(int));
    bonds_col_2=(int*)malloc(bonds_N*sizeof(int));
    bonds_col_3=(int*)malloc(bonds_N*sizeof(int));
    angles_col_1=(int*)malloc(angles_N*sizeof(int));
    angles_col_2=(int*)malloc(angles_N*sizeof(int));
    angles_col_3=(int*)malloc(angles_N*sizeof(int));
    angles_col_4=(int*)malloc(angles_N*sizeof(int));
    dihedrals_col_1=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_2=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_3=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_4=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_5=(int*)malloc(dihedrals_N*sizeof(int));
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<bonds_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&bonds_col_1[i],&bonds_col_2[i],&bonds_col_3[i]);
    }
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<angles_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&int_buffer,&angles_col_1[i],&angles_col_2[i],&angles_col_3[i],&angles_col_4[i]);
    }
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<dihedrals_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d\n",&int_buffer,&dihedrals_col_1[i],&dihedrals_col_2[i],&dihedrals_col_3[i],&dihedrals_col_4[i],&dihedrals_col_5[i]);
    }
    
    fclose(fp);
    
    // write altered topo.out	
    sprintf(file_path,"%s/topo.out.altered",current_folder);
    fp=fopen(file_path,"w+");
    fprintf(fp,"%d bonds\n",bonds_N);
    fprintf(fp,"%d angles\n",angles_N);
    fprintf(fp,"%d dihedrals\n",dihedrals_N);
    fprintf(fp,"\n");
    
    //fprintf(fp,"%d bond types\n",bonds);
    fprintf(fp,"%d bond types\n",final[B_lines-1]);
    
    fprintf(fp,"%d angle types\n",angles);
    //fprintf(fp,"%d angle types\n",final[A_lines-1]);
    
    fprintf(fp,"%d dihedral types\n",dihedrals);
    //fprintf(fp,"%d dihedral types\n",final[D_lines-1]);
    
    fprintf(fp,"\nBonds\n");
    fprintf(fp,"\n");
    for(i=0;i<bonds_N;++i)
    {
        // switch
        for(j=0;j<B_lines;++j){if(bonds_col_1[i]==init0[j]){temp=final[j];}}
        fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,temp,bonds_col_2[i],bonds_col_3[i]);
        //fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,bonds_col_1[i],bonds_col_2[i],bonds_col_3[i]);
    }
    fprintf(fp,"\nAngles\n");
    fprintf(fp,"\n");
    for(i=0;i<angles_N;++i)
    {
        // switch
        //for(j=0;j<A_lines;++j){if(angles_col_1[i]==init0[j]){temp=final[j];}}
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,temp,angles_col_2[i],angles_col_3[i],angles_col_4[i]);
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,angles_col_1[i],angles_col_2[i],angles_col_3[i],angles_col_4[i]);
    }
    fprintf(fp,"\nDihedrals\n");
    fprintf(fp,"\n");
    for(i=0;i<dihedrals_N;++i)
    {
        // switch
        //for(j=0;j<D_lines;++j){if(dihedrals_col_1[i]==init0[j]){temp=final[j];}}
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,temp,dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,dihedrals_col_1[i],dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
    }
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    //------------------------------------------------------------------
    // free allocated memory	
    
    free(init);free(final);//free(sites);
    for(i=0;i<D_lines;++i){free(D1[i]);free(D2[i]);free(D3[i]);free(D4[i]);}free(D1);free(D2);free(D3);free(D4);
    for(i=0;i<A_lines;++i){free(A1[i]);free(A2[i]);free(A3[i]);}free(A1);free(A2);free(A3);
    for(i=0;i<B_lines;++i){free(B1[i]);free(B2[i]);}free(B1);free(B2);
    free(bonds_col_1);free(bonds_col_2);free(bonds_col_3);
    free(angles_col_1);free(angles_col_2);free(angles_col_3);free(angles_col_4);
    free(dihedrals_col_1);free(dihedrals_col_2);free(dihedrals_col_3);free(dihedrals_col_4);free(dihedrals_col_5);
    free(store1);free(store2);
    free(init0);
    for(i=0;i<rows;++i){free(info[i]);}free(info);
    free(unique);free(unique_ID);
    
    free(B_map);free(A_map);free(D_map);
    
}

void alterD(int argc,char *argv);
void alterD(int argc,char *argv)
{
    FILE *fp;
    
    // fp related
    char current_folder[cmax_length],file_path[cmax_length];
    
    // dummies
    char *fcheck;
    char buffer[cmax_length],temp_species[species_length];
    int int_buffer;
    
    //
    int i,j,k;
    int bonds,angles,dihedrals,length;
    int bonds_N,angles_N,dihedrals_N;
    int B_lines,A_lines,D_lines,lines;
    int flag,temp;
    int rows,cols,unique_count;
    int *sites,*init,*final;
    int *store1,*store2;
    int *init0;
    int *unique,*unique_ID;
    int *bonds_col_1,*bonds_col_2,*bonds_col_3;
    int *angles_col_1,*angles_col_2,*angles_col_3,*angles_col_4;
    int *dihedrals_col_1,*dihedrals_col_2,*dihedrals_col_3,*dihedrals_col_4,*dihedrals_col_5;
    int **info;
    char **B1,**B2;
    char **A1,**A2,**A3;
    char **D1,**D2,**D3,**D4;
    
    //
    int *B_map;
    int *A_map;
    int *D_map;
    
    //------------------------------------------------------------------
    // print exe info
    if(argc==1)
    {
        printf(
               "\n*Reads topo.out and topo.log and unifies the selected dihedrals\n"
               "*The rules are read from <filename>\n"
               "*Generates topo.out.altered and topo.log.altered\n\n"
               "./alterD <filename>\n\n"
               );
        exit(-1);
    }
    
    // working directory
    fcheck=getcwd(current_folder,cmax_length);
    
    // read the rules from file and store in 2d matrix
    sprintf(file_path,"%s/%s",current_folder,argv);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    // find the number of cols
    rows=0;cols=2;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(buffer[0]=='[')
        {
            rows=rows+1;
            sscanf(buffer,"[%d]",&int_buffer);
            if(int_buffer>cols){cols=int_buffer;}
        }
    }
    
    rewind(fp);
    
    // matrix to store the alteration rules
    info=(int**)malloc(rows*sizeof(int*));
    for(i=0;i<rows;++i){info[i]=(int*)malloc(cols*sizeof(int));}
    // initialization
    for(i=0;i<rows;++i){for(j=0;j<cols;++j){info[i][j]=0;}}
    // populate
    i=-1;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(buffer[0]=='[')
        {
            i=i+1;
            sscanf(buffer,"[%d]",&int_buffer);
            for(j=0;j<int_buffer;++j)
            {
                fcheck=fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d",&info[i][j]);
            }
        }
    }
    fclose(fp);
    
    // open the log file created by topo
    sprintf(file_path,"%s/topo.log",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&bonds,&B_lines);			// read bonds
    
    //init=(int*)malloc(B_lines*sizeof(int));
    //final=(int*)malloc(B_lines*sizeof(int));
    store2=(int*)malloc(B_lines*sizeof(int));
    
    B1=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B1[i]=(char*)malloc(species_length*sizeof(char));}
    B2=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B2[i]=(char*)malloc(species_length*sizeof(char));}
    B_map=(int*)malloc(B_lines*sizeof(int));
    //for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s",&init[i],B1[i],B2[i]);}
    for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t#%d",&store2[i],B1[i],B2[i],&B_map[i]);}
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&angles,&A_lines);			// read angles
    
    //init=(int*)malloc(A_lines*sizeof(int));
    //final=(int*)malloc(A_lines*sizeof(int));
    store1=(int*)malloc(A_lines*sizeof(int));
    
    A1=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A1[i]=(char*)malloc(species_length*sizeof(char));}
    A2=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A2[i]=(char*)malloc(species_length*sizeof(char));}
    A3=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A3[i]=(char*)malloc(species_length*sizeof(char));}
    A_map=(int*)malloc(A_lines*sizeof(int));
    //for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s",&init[i],A1[i],A2[i],A3[i]);}
    for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t#%d",&store1[i],A1[i],A2[i],A3[i],&A_map[i]);}
    
    fcheck=fgets(buffer,cmax_length,fp);													// skip line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&dihedrals,&D_lines);		// read dihedrals
    
    init=(int*)malloc(D_lines*sizeof(int));
    final=(int*)malloc(D_lines*sizeof(int));
    
    D1=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D1[i]=(char*)malloc(species_length*sizeof(char));}
    D2=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D2[i]=(char*)malloc(species_length*sizeof(char));}
    D3=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D3[i]=(char*)malloc(species_length*sizeof(char));}
    D4=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D4[i]=(char*)malloc(species_length*sizeof(char));}
    D_map=(int*)malloc(D_lines*sizeof(int));
    for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s\t#%d",&init[i],D1[i],D2[i],D3[i],D4[i],&D_map[i]);}
    //for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s",&store2[i],D1[i],D2[i],D3[i],D4[i]);}
    
    fclose(fp);
    
    //lines=B_lines;
    //lines=A_lines;
    lines=D_lines;
    
    // save the initial numbering
    init0=(int*)malloc(lines*sizeof(int));
    for(i=0;i<lines;++i){init0[i]=init[i];}
    
    //----------------------------------------------------------------------
    // the first application
    
    // find the number of elements to unify
    length=0;
    for(i=0;i<cols;++i){if(info[0][i]!=0){length=length+1;}}
    sites=(int*)malloc(length*sizeof(int));
    // parse the info from info to sites
    for(i=0;i<length;++i){sites[i]=info[0][i];}
    
    // bubble sort the info
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<length;++i)
        {
            if(sites[i]<sites[i-1])
            {
                temp=sites[i-1];
                sites[i-1]=sites[i];
                sites[i]=temp;
                
                flag=1;
            }
        }
    }
    
    // initialize final array
    for(i=0;i<lines;++i){final[i]=init[i];}
    // alter entries
    for(i=0;i<lines;++i){for(j=0;j<length;++j){if(init[i]==sites[j]){final[i]=sites[0];}}}
    // bubble sort
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<lines;++i)
        {
            if(final[i]<final[i-1])
            {
                temp=final[i-1];
                final[i-1]=final[i];
                final[i]=temp;
                
                temp=init[i-1];
                init[i-1]=init[i];
                init[i]=temp;
                
                // initials
                temp=init0[i-1];
                init0[i-1]=init0[i];
                init0[i]=temp;
                
                sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
                sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
                sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
                sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);
                temp=D_map[i-1];
                D_map[i-1]=D_map[i];
                D_map[i]=temp;
                
                //sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
                //sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
                //sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);
                
                //sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
                //sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);
                
                flag=1;
            }
        }
    }
    
    // normalize (the code is taken from the erase scripts)
    
    unique=(int*)malloc(lines*sizeof(int));
    unique_ID=(int*)malloc(lines*sizeof(int));
    
    for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
    
    unique_count=0;	// counter for the unique appearances
    for(i=0;i<lines;++i)
    {
        if(unique[i]!=0)
        {
            for(j=i+1;j<lines;++j)
            {
                if(unique[i]==unique[j]){unique[j]=0;}
            }
            unique_count=unique_count+1;
            unique_ID[i]=unique_count;
        }
    }
    
    for(i=0;i<lines;++i)
    {
        //if(final[i]!=0)
        //{
        for(j=0;j<lines;++j)
        {
            if(final[i]==unique[j])
            {
                final[i]=unique_ID[j];
            }
        }
        //}
    }
    
    free(sites);
    
    //----------------------------------------------------------------------
    
    // the rest applications of the rules
    
    for(k=1;k<rows;++k)
    {
        
        //----------------------------------------------------------------------
        // switch!
        for(i=0;i<lines;++i){init[i]=final[i];}
        
        length=0;
        for(i=0;i<cols;++i){if(info[k][i]!=0){length=length+1;}}
        
        sites=(int*)malloc(length*sizeof(int));
        
        for(i=0;i<length;++i){sites[i]=info[k][i];}
        
        for(i=0;i<lines;++i)
        {
            for(j=0;j<length;++j)
            {
                if(sites[j]==init0[i]){sites[j]=final[i];}
            }
        }
        
        // bubble sort
        flag=1;
        while(flag==1)
        {
            flag=0;
            for(i=1;i<length;++i)
            {
                if(sites[i]<sites[i-1])
                {
                    temp=sites[i-1];
                    sites[i-1]=sites[i];
                    sites[i]=temp;
                    
                    flag=1;
                }
            }
        }
        
        // initialize final array
        for(i=0;i<lines;++i){final[i]=init[i];}
        // alter entries
        for(i=0;i<lines;++i){for(j=0;j<length;++j){if(init[i]==sites[j]){final[i]=sites[0];}}}
        // bubble sort
        flag=1;
        while(flag==1)
        {
            flag=0;
            for(i=1;i<lines;++i)
            {
                if(final[i]<final[i-1])
                {
                    temp=final[i-1];
                    final[i-1]=final[i];
                    final[i]=temp;
                    
                    temp=init[i-1];
                    init[i-1]=init[i];
                    init[i]=temp;
                    
                    // initial
                    temp=init0[i-1];
                    init0[i-1]=init0[i];
                    init0[i]=temp;
                    
                    sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
                    sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
                    sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
                    sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);
                    temp=D_map[i-1];
                    D_map[i-1]=D_map[i];
                    D_map[i]=temp;
                    
                    //sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);
                    
                    //sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
                    //sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);
                    
                    flag=1;
                }
            }
        }
        
        // normalize
        
        for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
        
        unique_count=0;	// counter for the unique appearances
        for(i=0;i<lines;++i)
        {
            if(unique[i]!=0)
            {
                for(j=i+1;j<lines;++j)
                {
                    if(unique[i]==unique[j]){unique[j]=0;}
                }
                unique_count=unique_count+1;
                unique_ID[i]=unique_count;
            }
        }
        
        for(i=0;i<lines;++i)
        {
            //if(final[i]!=0)
            //{
            for(j=0;j<lines;++j)
            {
                if(final[i]==unique[j])
                {
                    final[i]=unique_ID[j];
                }
            }
            //}
        }
        
        free(sites);
        
        //----------------------------------------------------------------------
        
    }
    
    // write altered log file
    sprintf(file_path,"%s/topo.log.altered",current_folder);
    fp=fopen(file_path,"w+");
    
    fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
    for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t#%d\n",store2[i],B1[i],B2[i],B_map[i]);}
    fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
    for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t#%d\n",store1[i],A1[i],A2[i],A3[i],A_map[i]);}
    fprintf(fp,"Dihedrals\n%d(%d) types\n",final[D_lines-1],D_lines);
    for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\t#%d\n",final[i],D1[i],D2[i],D3[i],D4[i],D_map[i]);}
    
    /*
     fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
     for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",store1[i],B1[i],B2[i]);}
     fprintf(fp,"Angles\n%d(%d) types\n",final[A_lines-1],A_lines);
     for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",final[i],A1[i],A2[i],A3[i]);}
     fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
     for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",store2[i],D1[i],D2[i],D3[i],D4[i]);}
     */
    /*
     fprintf(fp,"Bonds\n%d(%d) types\n",final[B_lines-1],B_lines);
     for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",final[i],B1[i],B2[i]);}
     fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
     for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",store1[i],A1[i],A2[i],A3[i]);}
     fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
     for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",store2[i],D1[i],D2[i],D3[i],D4[i]);}
     */
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    // print to console the switching info
    for(i=0;i<lines;++i){printf("%d\t%d\n",init0[i],final[i]);}
    
    // read topo.out
    sprintf(file_path,"%s/topo.out",current_folder);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds_N);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles_N);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals_N);
    
    bonds_col_1=(int*)malloc(bonds_N*sizeof(int));
    bonds_col_2=(int*)malloc(bonds_N*sizeof(int));
    bonds_col_3=(int*)malloc(bonds_N*sizeof(int));
    angles_col_1=(int*)malloc(angles_N*sizeof(int));
    angles_col_2=(int*)malloc(angles_N*sizeof(int));
    angles_col_3=(int*)malloc(angles_N*sizeof(int));
    angles_col_4=(int*)malloc(angles_N*sizeof(int));
    dihedrals_col_1=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_2=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_3=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_4=(int*)malloc(dihedrals_N*sizeof(int));
    dihedrals_col_5=(int*)malloc(dihedrals_N*sizeof(int));
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<bonds_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&bonds_col_1[i],&bonds_col_2[i],&bonds_col_3[i]);
    }
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<angles_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&int_buffer,&angles_col_1[i],&angles_col_2[i],&angles_col_3[i],&angles_col_4[i]);
    }
    
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    fcheck=fgets(buffer,cmax_length,fp);	// ignore line
    
    for(i=0;i<dihedrals_N;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d\n",&int_buffer,&dihedrals_col_1[i],&dihedrals_col_2[i],&dihedrals_col_3[i],&dihedrals_col_4[i],&dihedrals_col_5[i]);
    }
    
    fclose(fp);
    
    // write altered topo.out	
    sprintf(file_path,"%s/topo.out.altered",current_folder);
    fp=fopen(file_path,"w+");
    fprintf(fp,"%d bonds\n",bonds_N);
    fprintf(fp,"%d angles\n",angles_N);
    fprintf(fp,"%d dihedrals\n",dihedrals_N);
    fprintf(fp,"\n");
    
    fprintf(fp,"%d bond types\n",bonds);
    //fprintf(fp,"%d bond types\n",final[B_lines-1]);
    
    fprintf(fp,"%d angle types\n",angles);
    //fprintf(fp,"%d angle types\n",final[A_lines-1]);
    
    //fprintf(fp,"%d dihedral types\n",dihedrals);
    fprintf(fp,"%d dihedral types\n",final[D_lines-1]);
    
    fprintf(fp,"\nBonds\n");
    fprintf(fp,"\n");
    for(i=0;i<bonds_N;++i)
    {
        // switch
        //for(j=0;j<B_lines;++j){if(bonds_col_1[i]==init0[j]){temp=final[j];}}
        //fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,temp,bonds_col_2[i],bonds_col_3[i]);
        fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,bonds_col_1[i],bonds_col_2[i],bonds_col_3[i]);
    }
    fprintf(fp,"\nAngles\n");
    fprintf(fp,"\n");
    for(i=0;i<angles_N;++i)
    {
        // switch
        //for(j=0;j<A_lines;++j){if(angles_col_1[i]==init0[j]){temp=final[j];}}
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,temp,angles_col_2[i],angles_col_3[i],angles_col_4[i]);
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,angles_col_1[i],angles_col_2[i],angles_col_3[i],angles_col_4[i]);
    }
    fprintf(fp,"\nDihedrals\n");
    fprintf(fp,"\n");
    for(i=0;i<dihedrals_N;++i)
    {
        // switch
        for(j=0;j<D_lines;++j){if(dihedrals_col_1[i]==init0[j]){temp=final[j];}}
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,temp,dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,dihedrals_col_1[i],dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
    }
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    //------------------------------------------------------------------
    // free allocated memory	
    
    free(init);free(final);//free(sites);
    for(i=0;i<D_lines;++i){free(D1[i]);free(D2[i]);free(D3[i]);free(D4[i]);}free(D1);free(D2);free(D3);free(D4);
    for(i=0;i<A_lines;++i){free(A1[i]);free(A2[i]);free(A3[i]);}free(A1);free(A2);free(A3);
    for(i=0;i<B_lines;++i){free(B1[i]);free(B2[i]);}free(B1);free(B2);
    free(bonds_col_1);free(bonds_col_2);free(bonds_col_3);
    free(angles_col_1);free(angles_col_2);free(angles_col_3);free(angles_col_4);
    free(dihedrals_col_1);free(dihedrals_col_2);free(dihedrals_col_3);free(dihedrals_col_4);free(dihedrals_col_5);
    free(store1);free(store2);
    free(init0);
    for(i=0;i<rows;++i){free(info[i]);}free(info);
    free(unique);free(unique_ID);
    
    free(B_map);free(A_map);free(D_map);
    
}

void label(int argc, char *argv);
void label(int argc, char *argv)
{
    FILE *fp,*fp_write;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length];
    char temp[cmax_length];
    int i,j,bonds;
    int *b1,*b2;
    int max;
    int *count_b;
    
    int angles,*a1,*a2,*count_a;
    int dihedrals,*d1,*d2,*count_d;
    
    getcwd(current_folder,cmax_length);
    
    sprintf(file_path,"%s/%s",current_folder,argv);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
    
    b1=(int*)malloc(bonds*sizeof(int));
    b2=(int*)malloc(bonds*sizeof(int));
    
    for(i=0;i<bonds;++i){
        fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t#%d",&b1[i],temp,temp,&b2[i]);
    }
    
    max=b2[0];
    for(i=1;i<bonds;++i){if(b2[i]>max){max=b2[i];}}
    
    count_b=(int*)malloc(max*sizeof(int));
    for(i=0;i<max;++i){count_b[i]=0;}
    
    for(i=0;i<bonds;++i)
    {
        count_b[b2[i]-1]=count_b[b2[i]-1]+1;
    }
    
    sprintf(file_path,"%s/bonds",current_folder);
    fp_write=fopen(file_path,"w+");
    
    printf("Bonds\n");
    
    for(i=0;i<max;++i)
    {
        if(count_b[i]!=1)
        {
            fprintf(fp_write,"[%d]\n",count_b[i]);
            printf("[%d]\n",count_b[i]);
            for(j=0;j<bonds;++j)
            {
                if(b2[j]==i+1)
                {
                    fprintf(fp_write,"%d\n",b1[j]);
                    printf("%d\n",b1[j]);
                }
            }
        }
    }
    
    fclose(fp_write);
    
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
    
    a1=(int*)malloc(angles*sizeof(int));
    a2=(int*)malloc(angles*sizeof(int));
    
    for(i=0;i<angles;++i){
        fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t#%d",&a1[i],temp,temp,temp,&a2[i]);
    }
    
    max=a2[0];
    for(i=1;i<angles;++i){if(a2[i]>max){max=a2[i];}}
    
    count_a=(int*)malloc(max*sizeof(int));
    for(i=0;i<max;++i){count_a[i]=0;}
    
    for(i=0;i<angles;++i)
    {
        count_a[a2[i]-1]=count_a[a2[i]-1]+1;
    }
    
    sprintf(file_path,"%s/angles",current_folder);
    fp_write=fopen(file_path,"w+");
    
    printf("Angles\n");
    
    for(i=0;i<max;++i)
    {
        if(count_a[i]!=1)
        {
            fprintf(fp_write,"[%d]\n",count_a[i]);
            printf("[%d]\n",count_a[i]);
            for(j=0;j<angles;++j)
            {
                if(a2[j]==i+1)
                {
                    fprintf(fp_write,"%d\n",a1[j]);
                    printf("%d\n",a1[j]);
                }
            }
        }
    }
    
    fclose(fp_write);
    
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
    
    d1=(int*)malloc(dihedrals*sizeof(int));
    d2=(int*)malloc(dihedrals*sizeof(int));
    
    for(i=0;i<dihedrals;++i){
        fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s\t#%d",&d1[i],temp,temp,temp,temp,&d2[i]);
    }
    
    max=d2[0];
    for(i=1;i<dihedrals;++i){if(d2[i]>max){max=d2[i];}}
    
    count_d=(int*)malloc(max*sizeof(int));
    for(i=0;i<max;++i){count_d[i]=0;}
    
    for(i=0;i<dihedrals;++i)
    {
        count_d[d2[i]-1]=count_d[d2[i]-1]+1;
    }
    
    sprintf(file_path,"%s/dihedrals",current_folder);
    fp_write=fopen(file_path,"w+");
    
    printf("Dihedrals\n");
    
    for(i=0;i<max;++i)
    {
        if(count_d[i]!=1)
        {
            fprintf(fp_write,"[%d]\n",count_d[i]);
            printf("[%d]\n",count_d[i]);
            for(j=0;j<dihedrals;++j)
            {
                if(d2[j]==i+1)
                {
                    fprintf(fp_write,"%d\n",d1[j]);
                    printf("%d\n",d1[j]);
                }
            }
        }
    }
    
    fclose(fp_write);
    
    fclose(fp);
    
    free(b1);
    free(b2);
    free(count_b);
    free(a1);
    free(a2);
    free(count_a);
    free(d1);
    free(d2);
    free(count_d);
}

void mtol(int argc, char *argv1, char *argv2);
void mtol(int argc, char *argv1, char *argv2)
{
    FILE *fp;
    
    // fp related
    char current_folder[cmax_length],file_path[cmax_length];
    
    // dummies
    char buffer[cmax_length],c_buffer[mol2_length];
    char *fcheck;
    int int_buffer;
    
    //
    double *x,*y,*z;
    int i,j,mode;
    int atoms;
    int *atomic_ID,*molecule_ID;
    char **species;
    
    // species ID variables
    int flag,unique_species_number;
    int *unique_species_ID;
    char **unique_species;
    
    //------------------------------------------------------------------
    
    // print exe info
    if(argc==1 || argc==2)
    {
        printf(
               "\n*mol2 to LAMMPS data file converter\n"
               "*Supports atomic, molecular and full style\n"
               "*Supercell parameters, masses and charges must be entered manually\n"
               "*The mol2 filename is used as the simulation title\n"
               "*Generates init.lammps\n\n"
               "./mtol <filename>.mol2 <mode>\n"
               "\tmode=0: 'atomic_style'\n"
               "\tmode=1: 'molecular_style'\n"
               "\tmode=2: 'full_style'\n\n"
               );
        exit(-1);
    }
    
    // the mode flag
    // even in mode=2, the program ignores the mol2 charge info
    mode=atoi(argv2);
    if(mode!=0 && mode!=1 && mode!=2){printf("Invalid mode value\n");exit(-1);}
    
    // read working directory
    fcheck=getcwd(current_folder,cmax_length);
    
    // import data from .mol2 file
    sprintf(file_path,"%s/%s.mol2",current_folder,argv1);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
    
    for(i=0;i<2;++i){fcheck=fgets(buffer,cmax_length,fp);}	// ignore first and second line
    fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms,&int_buffer);	// the bonds info is ignored; topo.c deals with topology
    while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
    //for(i=0;i<4;++i){fcheck=fgets(buffer,cmax_length,fp);}	// ignore lines...
    
    // preallocations
    x=(double*)malloc(atoms*sizeof(double));	// coords
    y=(double*)malloc(atoms*sizeof(double));
    z=(double*)malloc(atoms*sizeof(double));
    molecule_ID=(int*)malloc(atoms*sizeof(int));
    species=(char**)malloc(atoms*sizeof(char*));// species - reads the left species label
    for (i=0;i<atoms;++i){species[i]=(char*)malloc(species_length*sizeof(char));}
    
    // read species and coords data
    for(i=0;i<atoms;++i)
    {
        fcheck=fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d",&int_buffer,species[i],&x[i],&y[i],&z[i],c_buffer,&molecule_ID[i]);
    }
    
    fclose(fp);
    
    // assign an integer counter starting from 1 for every different species present in the dat file
    atomic_ID=(int*)malloc(atoms*sizeof(int));			// the array for the new species representation based on integers
    unique_species=(char**)malloc(atoms*sizeof(char*));	// the array to store the unique species
    for (i=0;i<atoms;++i){unique_species[i]=(char*)malloc(species_length*sizeof(char));}
    
    // identify unique secies
    sprintf(unique_species[0],"%s",species[0]);	// the simplest case: only one species present - assign identifier and store
    unique_species_number=1;					// initialize the counter
    // check the rest species
    for(i=1;i<atoms;++i)
    {
        flag=0;	// flag==0 means that the examined species is unique
        for(j=0;j<unique_species_number;++j)
        {
            if(strcmp(species[i],unique_species[j])==0){flag=1;}
        }
        if(flag==0)
        {
            sprintf(unique_species[unique_species_number],"%s",species[i]);	// store species
            unique_species_number=unique_species_number+1;					// +1 for the counter
        }
    }
    // store the unique species ids in an array
    unique_species_ID=(int*)malloc(unique_species_number*sizeof(int));
    for(i=0;i<unique_species_number;++i){unique_species_ID[i]=i+1;}
    
    // populate the atomic_ID array
    for(i=0;i<atoms;++i)
    {
        for(j=0;j<unique_species_number;++j)
        {
            if(strcmp(species[i],unique_species[j])==0)
            {
                //atomic_ID[i]=j+1;
                atomic_ID[i]=unique_species_ID[j];
            }
        }
    }
    
    // output
    sprintf(file_path,"%s/init.lammps",current_folder);
    fp=fopen(file_path,"w+");
    
    fprintf(fp,"%s\n\n",argv1);
    fprintf(fp,"%d atoms\n\n",atoms);
    fprintf(fp,"%d atom types\n\n",unique_species_number);
    fprintf(fp,"%s %s xlo xhi\n","xmin","xmax");
    fprintf(fp,"%s %s ylo yhi\n","ymin","ymax");
    fprintf(fp,"%s %s zlo zhi\n\n","zmin","zmax");
    fprintf(fp,"Masses\n\n");
    for(i=0;i<unique_species_number;++i){fprintf(fp,"%d\t%s_mass\n",unique_species_ID[i],unique_species[i]);}
    fprintf(fp,"\n");
    fprintf(fp,"Atoms\n\n");
    if(mode==0)
    {
        for(i=0;i<atoms;++i){fprintf(fp,"%d\t%d\t%lf\t%lf\t%lf\n",i+1,atomic_ID[i],x[i],y[i],z[i]);}
    }
    if(mode==1)
    {
        for(i=0;i<atoms;++i){fprintf(fp,"%d\t%d\t%d\t%lf\t%lf\t%lf\n",i+1,molecule_ID[i],atomic_ID[i],x[i],y[i],z[i]);}
    }
    if(mode==2)
    {
        for(i=0;i<atoms;++i){fprintf(fp,"%d\t%d\t%d\t%s_charge\t%lf\t%lf\t%lf\n",i+1,molecule_ID[i],atomic_ID[i],species[i],x[i],y[i],z[i]);}
    }
    
    fclose(fp);
    printf("\nGenerated %s\n",file_path);
    
    //------------------------------------------------------------------
    // free memory
    
    for (i=0;i<atoms;++i){free(species[i]);}free(species);
    for (i=0;i<atoms;++i){free(unique_species[i]);}free(unique_species);
    free(atomic_ID);
    free(unique_species_ID);
    free(molecule_ID);
    free(x);free(y);free(z);
    
}

void topo(int argc, char *argv);
void topo(int argc, char *argv)
{
    FILE *fp;
    
    // fp related
    char current_folder[cmax_length],file_path[cmax_length];
    
    // dummies
    int int_buffer;
    char buffer[cmax_length],word[cmax_length];
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
    sprintf(file_path,"%s/%s.mol2",current_folder,argv);
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
        //sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s",&int_buffer,word,&x[i],&y[i],&z[i],species[i]);
        sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf",&int_buffer,species[i],&x[i],&y[i],&z[i]);
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
}

void eraseBAD(char mode,int argc, char **argv);
void eraseBAD(char mode,int argc, char **argv)
{
	FILE *fp;
	
	// fp related
	char current_folder[cmax_length],file_path[cmax_length];
	
	// dummies
	char *fcheck;
	char buffer[cmax_length],temp_species[species_length];
	int int_buffer;
	
	//
	int i,j,k;
	int bonds,angles,dihedrals;
	int B_lines,A_lines,D_lines,lines;
	int bonds_N,angles_N,dihedrals_N;
	int flag,temp,length,value;
	int unique_count,final_N=0,final_lines;
	
	int *unique,*unique_ID;
	int *sites,*init,*final;
	
	char **B1,**B2;
	char **A1,**A2,**A3;
	char **D1,**D2,**D3,**D4;
	
	int *B_ID,*A_ID,*D_ID;
	int *bonds_col_1,*bonds_col_2,*bonds_col_3;
	int *angles_col_1,*angles_col_2,*angles_col_3,*angles_col_4;
	int *dihedrals_col_1,*dihedrals_col_2,*dihedrals_col_3,*dihedrals_col_4,*dihedrals_col_5;
	
	//------------------------------------------------------------------

	// save the arguments in the sites array
	length=argc;
	sites=(int*)malloc(length*sizeof(int));
	for(i=0;i<argc;++i){sites[i]=atoi(argv[i]);printf("%d ****\n",sites[i]);}

	// working directory
	fcheck=getcwd(current_folder,cmax_length);
	
	// open the log file created by topo
	sprintf(file_path,"%s/topo.log",current_folder);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
	
	fcheck=fgets(buffer,cmax_length,fp);													// skip line
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&bonds,&B_lines);			// read bonds
	
	B_ID=(int*)malloc(B_lines*sizeof(int));
	B1=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B1[i]=(char*)malloc(species_length*sizeof(char));}
	B2=(char**)malloc(B_lines*sizeof(char*));for(i=0;i<B_lines;++i){B2[i]=(char*)malloc(species_length*sizeof(char));}	
	for(i=0;i<B_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s",&B_ID[i],B1[i],B2[i]);}					
	
	fcheck=fgets(buffer,cmax_length,fp);													// skip line
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&angles,&A_lines);			// read angles				
	
	A_ID=(int*)malloc(A_lines*sizeof(int));
	A1=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A1[i]=(char*)malloc(species_length*sizeof(char));}
	A2=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A2[i]=(char*)malloc(species_length*sizeof(char));}
	A3=(char**)malloc(A_lines*sizeof(char*));for(i=0;i<A_lines;++i){A3[i]=(char*)malloc(species_length*sizeof(char));}
	for(i=0;i<A_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s",&A_ID[i],A1[i],A2[i],A3[i]);}
	
	fcheck=fgets(buffer,cmax_length,fp);													// skip line
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d(%d)",&dihedrals,&D_lines);		// read dihedrals
	
	D_ID=(int*)malloc(D_lines*sizeof(int));	
	D1=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D1[i]=(char*)malloc(species_length*sizeof(char));}
	D2=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D2[i]=(char*)malloc(species_length*sizeof(char));}
	D3=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D3[i]=(char*)malloc(species_length*sizeof(char));}
	D4=(char**)malloc(D_lines*sizeof(char*));for(i=0;i<D_lines;++i){D4[i]=(char*)malloc(species_length*sizeof(char));}
	for(i=0;i<D_lines;++i){fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"[%d]\t%s\t%s\t%s\t%s",&D_ID[i],D1[i],D2[i],D3[i],D4[i]);}
	
	fclose(fp);

	if(mode=='B')lines=B_lines;
	if(mode=='A')lines=A_lines;
	if(mode=='D')lines=D_lines;

	// save the numbering info
	init=(int*)malloc(lines*sizeof(int));
	final=(int*)malloc(lines*sizeof(int));
	if(mode=='B'){for(i=0;i<lines;++i){init[i]=B_ID[i];final[i]=B_ID[i];}}
	if(mode=='A'){for(i=0;i<lines;++i){init[i]=A_ID[i];final[i]=A_ID[i];}}
	if(mode=='D'){for(i=0;i<lines;++i){init[i]=D_ID[i];final[i]=D_ID[i];}}

	// set equal to zero the elements for deletion
	for(i=0;i<lines;++i)
	{
		for(j=0;j<length;++j)
		{
			if(mode=='B'){if(B_ID[i]==sites[j]){final[i]=0;}}
			if(mode=='A'){if(A_ID[i]==sites[j]){final[i]=0;}}
			if(mode=='D'){if(D_ID[i]==sites[j]){final[i]=0;}}
		}
	}

	// bubble sort
	flag=1;
	while(flag==1)
	{
		flag=0;
		for(i=1;i<lines;++i)
		{
			if(final[i]<final[i-1])
			{
				temp=final[i-1];
				final[i-1]=final[i];
				final[i]=temp;
				
				temp=init[i-1];
				init[i-1]=init[i];
				init[i]=temp;
				
				if(mode=='B'){temp=B_ID[i-1];B_ID[i-1]=B_ID[i];B_ID[i]=temp;}
				if(mode=='A'){temp=A_ID[i-1];A_ID[i-1]=A_ID[i];A_ID[i]=temp;}
				if(mode=='D'){temp=D_ID[i-1];D_ID[i-1]=D_ID[i];D_ID[i]=temp;}
				
				if(mode=='D'){
				sprintf(temp_species,"%s",D1[i-1]);sprintf(D1[i-1],"%s",D1[i]);sprintf(D1[i],"%s",temp_species);
				sprintf(temp_species,"%s",D2[i-1]);sprintf(D2[i-1],"%s",D2[i]);sprintf(D2[i],"%s",temp_species);
				sprintf(temp_species,"%s",D3[i-1]);sprintf(D3[i-1],"%s",D3[i]);sprintf(D3[i],"%s",temp_species);
				sprintf(temp_species,"%s",D4[i-1]);sprintf(D4[i-1],"%s",D4[i]);sprintf(D4[i],"%s",temp_species);}
				
				if(mode=='A'){
				sprintf(temp_species,"%s",A1[i-1]);sprintf(A1[i-1],"%s",A1[i]);sprintf(A1[i],"%s",temp_species);
				sprintf(temp_species,"%s",A2[i-1]);sprintf(A2[i-1],"%s",A2[i]);sprintf(A2[i],"%s",temp_species);
				sprintf(temp_species,"%s",A3[i-1]);sprintf(A3[i-1],"%s",A3[i]);sprintf(A3[i],"%s",temp_species);}
				
				if(mode=='B'){				
				sprintf(temp_species,"%s",B1[i-1]);sprintf(B1[i-1],"%s",B1[i]);sprintf(B1[i],"%s",temp_species);
				sprintf(temp_species,"%s",B2[i-1]);sprintf(B2[i-1],"%s",B2[i]);sprintf(B2[i],"%s",temp_species);}
				
				flag=1;
			}
		}	
	}

	// assign a new numbering for the remaining elements
	unique=(int*)malloc(lines*sizeof(int));
	unique_ID=(int*)malloc(lines*sizeof(int));
	for(i=0;i<lines;++i){unique[i]=final[i];unique_ID[i]=0;}
	
	unique_count=0;	// counter for the unique appearances
	for(i=0;i<lines;++i)
	{
		if(unique[i]!=0)
		{
			for(j=i+1;j<lines;++j)
			{
				if(unique[i]==unique[j]){unique[j]=0;}
			}
			unique_count=unique_count+1;
			unique_ID[i]=unique_count;
		}
	}

	final_lines=0;	// counter for the final remaining elements
	for(i=0;i<lines;++i)
	{
		if(final[i]!=0)
		{
			for(j=0;j<lines;++j)
			{
				if(final[i]==unique[j])
				{
					final[i]=unique_ID[j];
				}
			}
			final_lines=final_lines+1;
		}
	}

	// output
	for(i=0;i<lines;++i){printf("%d\t%d\n",init[i],final[i]);}printf("_*_*_*_*_\n");

	sprintf(file_path,"%s/topo.log.erased",current_folder);
	fp=fopen(file_path,"w+");
	
	if(mode=='B'){
	fprintf(fp,"Bonds\n%d(%d) types\n",unique_count,final_lines);
	for(i=0;i<B_lines;++i){if(final[i]!=0){fprintf(fp,"[%d]\t%s\t%s\n",final[i],B1[i],B2[i]);}}
	fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
	for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",A_ID[i],A1[i],A2[i],A3[i]);}
	fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
	for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",D_ID[i],D1[i],D2[i],D3[i],D4[i]);}}
	
	
	if(mode=='A'){
	fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
	for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",B_ID[i],B1[i],B2[i]);}
	fprintf(fp,"Angles\n%d(%d) types\n",unique_count,final_lines);
	for(i=0;i<A_lines;++i){if(final[i]!=0){fprintf(fp,"[%d]\t%s\t%s\t%s\n",final[i],A1[i],A2[i],A3[i]);}}
	fprintf(fp,"Dihedrals\n%d(%d) types\n",dihedrals,D_lines);
	for(i=0;i<D_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",D_ID[i],D1[i],D2[i],D3[i],D4[i]);}}
	
	
	if(mode=='D'){
	fprintf(fp,"Bonds\n%d(%d) types\n",bonds,B_lines);
	for(i=0;i<B_lines;++i){fprintf(fp,"[%d]\t%s\t%s\n",B_ID[i],B1[i],B2[i]);}
	fprintf(fp,"Angles\n%d(%d) types\n",angles,A_lines);
	for(i=0;i<A_lines;++i){fprintf(fp,"[%d]\t%s\t%s\t%s\n",A_ID[i],A1[i],A2[i],A3[i]);}
	fprintf(fp,"Dihedrals\n%d(%d) types\n",unique_count,final_lines);
	for(i=0;i<D_lines;++i){if(final[i]!=0){fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",final[i],D1[i],D2[i],D3[i],D4[i]);}}}
	
	fclose(fp);
	printf("\nGenerated %s\n",file_path);

	
	// read topo.out
	sprintf(file_path,"%s/topo.out",current_folder);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Cannot locate %s\n",file_path);exit(-1);}		// file check
	
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds_N);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles_N);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals_N);
	
	bonds_col_1=(int*)malloc(bonds_N*sizeof(int));
	bonds_col_2=(int*)malloc(bonds_N*sizeof(int));
	bonds_col_3=(int*)malloc(bonds_N*sizeof(int));
	angles_col_1=(int*)malloc(angles_N*sizeof(int));
	angles_col_2=(int*)malloc(angles_N*sizeof(int));
	angles_col_3=(int*)malloc(angles_N*sizeof(int));
	angles_col_4=(int*)malloc(angles_N*sizeof(int));
	dihedrals_col_1=(int*)malloc(dihedrals_N*sizeof(int));
	dihedrals_col_2=(int*)malloc(dihedrals_N*sizeof(int));
	dihedrals_col_3=(int*)malloc(dihedrals_N*sizeof(int));
	dihedrals_col_4=(int*)malloc(dihedrals_N*sizeof(int));
	dihedrals_col_5=(int*)malloc(dihedrals_N*sizeof(int));
	
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
	fcheck=fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
	
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	
	for(i=0;i<bonds_N;++i)
	{
		fcheck=fgets(buffer,cmax_length,fp);
		sscanf(buffer,"%d\t%d\t%d\t%d",&int_buffer,&bonds_col_1[i],&bonds_col_2[i],&bonds_col_3[i]);
	}
	
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	
	for(i=0;i<angles_N;++i)
	{
		fcheck=fgets(buffer,cmax_length,fp);
		sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&int_buffer,&angles_col_1[i],&angles_col_2[i],&angles_col_3[i],&angles_col_4[i]);
	}
	
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	fcheck=fgets(buffer,cmax_length,fp);	// ignore line
	
	for(i=0;i<dihedrals_N;++i)
	{
		fcheck=fgets(buffer,cmax_length,fp);
		sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d\n",&int_buffer,&dihedrals_col_1[i],&dihedrals_col_2[i],&dihedrals_col_3[i],&dihedrals_col_4[i],&dihedrals_col_5[i]);
	}
	
	fclose(fp);

	if(mode=='B'){
	for(i=0;i<bonds_N;++i)
	{
		for(j=0;j<lines;++j)
		{
			if(bonds_col_1[i]==init[j]){value=final[j];break;}
		}
		
		if(value!=0)
		{
			final_N=final_N+1;
		}
	}}

	if(mode=='A'){
	for(i=0;i<angles_N;++i)
	{
		for(j=0;j<lines;++j)
		{
			if(angles_col_1[i]==init[j]){value=final[j];break;}
		}
		
		if(value!=0)
		{
			final_N=final_N+1;
		}
	}}

	if(mode=='D'){
	for(i=0;i<dihedrals_N;++i)
	{
		for(j=0;j<lines;++j)
		{
			if(dihedrals_col_1[i]==init[j]){value=final[j];break;}
		}
		
		if(value!=0)
		{
			final_N=final_N+1;
		}
	}}

	// write altered topo.out	
	sprintf(file_path,"%s/topo.out.erased",current_folder);
	fp=fopen(file_path,"w+");

	if(mode=='B'){
	fprintf(fp,"%d bonds\n",final_N);
	fprintf(fp,"%d angles\n",angles_N);
	fprintf(fp,"%d dihedrals\n",dihedrals_N);}
	
	if(mode=='A'){
	fprintf(fp,"%d bonds\n",bonds_N);
	fprintf(fp,"%d angles\n",final_N);
	fprintf(fp,"%d dihedrals\n",dihedrals_N);}

	if(mode=='D'){
	fprintf(fp,"%d bonds\n",bonds_N);
	fprintf(fp,"%d angles\n",angles_N);
	fprintf(fp,"%d dihedrals\n",final_N);}
	
	fprintf(fp,"\n");

	if(mode=='B'){
	fprintf(fp,"%d bond types\n",unique_count);
	fprintf(fp,"%d angle types\n",angles);
	fprintf(fp,"%d dihedral types\n",dihedrals);
	}

	if(mode=='A'){
	fprintf(fp,"%d bond types\n",bonds);
	fprintf(fp,"%d angle types\n",unique_count);
	fprintf(fp,"%d dihedral types\n",dihedrals);
	}

	if(mode=='D'){
	fprintf(fp,"%d bond types\n",bonds);
	fprintf(fp,"%d angle types\n",angles);
	fprintf(fp,"%d dihedral types\n",unique_count);
	}

	if(mode=='B')
	{
		fprintf(fp,"\nBonds\n");
		fprintf(fp,"\n");
		k=0;
		for(i=0;i<bonds_N;++i)
		{
			for(j=0;j<lines;++j)
			{
				if(bonds_col_1[i]==init[j]){value=final[j];break;}
			}
		
			if(value!=0)
			{
				k=k+1;
				fprintf(fp,"%d\t%d\t%d\t%d\n",k,value,bonds_col_2[i],bonds_col_3[i]);
			}
		}
		fprintf(fp,"\nAngles\n");
		fprintf(fp,"\n");
		for(i=0;i<angles_N;++i)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,angles_col_1[i],angles_col_2[i],angles_col_3[i],angles_col_4[i]);
		}
		fprintf(fp,"\nDihedrals\n");
		fprintf(fp,"\n");
		for(i=0;i<dihedrals_N;++i)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,dihedrals_col_1[i],dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
		}
	}

	if(mode=='A')
	{
		fprintf(fp,"\nBonds\n");
		fprintf(fp,"\n");
		for(i=0;i<bonds_N;++i)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,bonds_col_1[i],bonds_col_2[i],bonds_col_3[i]);
		}
		fprintf(fp,"\nAngles\n");
		fprintf(fp,"\n");
		k=0;
		for(i=0;i<angles_N;++i)
		{
			for(j=0;j<A_lines;++j)
			{
				if(angles_col_1[i]==init[j]){value=final[j];break;}
			}
		
			if(value!=0)
			{
				k=k+1;
				fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",k,value,angles_col_2[i],angles_col_3[i],angles_col_4[i]);
			}
		}
		fprintf(fp,"\nDihedrals\n");
		fprintf(fp,"\n");
		for(i=0;i<dihedrals_N;++i)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",i+1,dihedrals_col_1[i],dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
		}
	}

	if(mode=='D')
	{
		fprintf(fp,"\nBonds\n");
		fprintf(fp,"\n");
		for(i=0;i<bonds_N;++i)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\n",i+1,bonds_col_1[i],bonds_col_2[i],bonds_col_3[i]);
		}
		fprintf(fp,"\nAngles\n");
		fprintf(fp,"\n");
		for(i=0;i<angles_N;++i)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",i+1,angles_col_1[i],angles_col_2[i],angles_col_3[i],angles_col_4[i]);
		}
		fprintf(fp,"\nDihedrals\n");
		fprintf(fp,"\n");
		k=0;
		for(i=0;i<dihedrals_N;++i)
		{
			for(j=0;j<lines;++j)
			{
				if(dihedrals_col_1[i]==init[j]){value=final[j];break;}
			}
		
			if(value!=0)
			{
				k=k+1;
				fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",k,value,dihedrals_col_2[i],dihedrals_col_3[i],dihedrals_col_4[i],dihedrals_col_5[i]);
			}
		}
	}

	fclose(fp);
	printf("\nGenerated %s\n",file_path);
	
	//------------------------------------------------------------------
	// free memory
	
	free(unique);free(unique_ID);
	free(sites);
	for(i=0;i<D_lines;++i){free(D1[i]);free(D2[i]);free(D3[i]);free(D4[i]);}free(D1);free(D2);free(D3);free(D4);
	for(i=0;i<A_lines;++i){free(A1[i]);free(A2[i]);free(A3[i]);}free(A1);free(A2);free(A3);
	for(i=0;i<B_lines;++i){free(B1[i]);free(B2[i]);}free(B1);free(B2);
	free(B_ID);free(A_ID);free(D_ID);
	free(init);free(final);
	free(bonds_col_1);free(bonds_col_2);free(bonds_col_3);
	free(angles_col_1);free(angles_col_2);free(angles_col_3);free(angles_col_4);
	free(dihedrals_col_1);free(dihedrals_col_2);free(dihedrals_col_3);free(dihedrals_col_4);free(dihedrals_col_5);
}
