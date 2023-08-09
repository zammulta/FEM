#include <stdio.h>
int main () {
    int nd;
    float ds, dss, u[100], sol[100], ua[100];
    float kmatrix[100][100]={0}, F[100]={0};

    //example input coor & cny
    //float coor[nd]={0,0.25,0.6,1};
    //float coor[9]={0,0.3,0.125,0.55,0.45,0.75,0.6,1,0.9};
    //int cny[2][nd-1]={{1,2,3},{2,3,4}};
    //int cny[2][8]={{1,3,2,5,4,7,6,9},{3,2,5,4,7,6,9,8}};

    //example input boundary condition
    //int lbc=0;
    //float vlbc=0;
    //int rbc=7;
    //float vrbc=0;

    printf("Enter the number of nodes: "); scanf("%d", &nd);
    float coor[nd];
    float c;
    for(int i=0; i<nd; i++){
        printf("Enter coor[%d]: ", i+1); scanf("%f", &c);
        coor[i]=c;
    }
    int cny[2][nd-1];
    int cn;
    for(int i=0; i<2; i++){
        for(int j=0; j<nd-1; j++){
            printf("Enter connectivity (0: left, 1: right)[%d][%d]: ", i,j); scanf("%d", &cn);
            cny[i][j]=cn;
        }
    }
    int lbc;
    float vlbc;
    printf("Enter node of left boundary condition (start from 0: example enter 0 for u(1)): "); scanf("%d", &lbc);
    printf("Enter value of left boundary condition: "); scanf("%f", &vlbc);
    int rbc;
    float vrbc;
    printf("Enter node of right boundary condition (start from 0: example enter 3 for u(4))): "); scanf("%d", &rbc);
    printf("Enter value of right boundary condition: "); scanf("%f", &vrbc);



    //Matrix stiffness
    for(int i = 0; i < nd-1;i++){
        ds=coor[cny[1][i]-1]-coor[cny[0][i]-1];
        dss=1/ds;
        kmatrix[cny[0][i]-1][cny[0][i]-1] += dss;
        kmatrix[cny[0][i]-1][cny[1][i]-1] += (-dss);
        kmatrix[cny[1][i]-1][cny[0][i]-1] += (-dss);
        kmatrix[cny[1][i]-1][cny[1][i]-1] += dss;
        //F
        F[cny[0][i]-1] += -ds/2;
        F[cny[1][i]-1] += -ds/2;

    }

    //combine K and F matrix
    int m = 0;
    float GE[nd-2][nd-1];
    for(int i=0;i<nd;i++){
        int n = 0;
        if(i == lbc){
        }
        else if(i==rbc){
        }
        else{
            for(int j=0;j<nd;j++){
                if(j==lbc){
                    }
                else if(j==rbc){
                }
                else{
                    GE[m][n]=kmatrix[i][j];
                    n+=1;
                }
            }
            GE[m][nd-2]=F[i];
            m+=1;
        }
    }

    printf("Matrix stiffness:\n");
    for(int m=0; m<nd; m++) {
        for(int n=0;n<nd;n++) {
            printf("%f ", kmatrix[m][n]);
        }
        printf("\n");
    }

    printf("\n");
    printf("Matrix F :\n");
    for(int z=0; z<nd; z++){
        printf("%f \n", F[z]);
    }
    printf("\n");

    printf("Reduced stiffness matix according to boundary condition: \n ");
    for(int i=0; i<nd-2; i++) {
        for(int j=0;j<nd-1;j++) {
            printf("%f ", GE[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    //Gaussian Elimination
    float cons;
    for (int i=0; i<nd-2; i++) {
        for (int j = i + 1; j < nd-2; j++) {
            cons = GE[j][i] / GE[i][i];
            for (int k = i; k < nd-1; k++) {
                GE[j][k] = GE[j][k] - cons * GE[i][k];
            }
        }
    }

    printf("Gaussian Elimination: \n ");
    for(int i=0; i<nd-2; i++) {
        for(int j=0;j<nd-1;j++) {
            printf("%f ", GE[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    //Find the u
    for(int i=nd-3; i>=0; i--){
        float a=0;
        for(int j=i+1; j<nd-2; j++){
            a += GE[i][j]*u[j];
        }
        u[i]=(GE[i][nd-2]-a)/GE[i][i];
    }


    printf("One-dimensional linear elements with Gaussian elimination solver: \n");
    int p=0;
    for(int i=0; i<nd; i++){
        if(i == lbc){
            sol[i]=vlbc; //left boundary condition
            printf("u[%d]= %f\n", i+1, sol[i]);
        }
        else if(i == rbc){
            sol[i]=vrbc; //right boundary condition
            printf("u[%d]= %f\n", i+1, sol[i]);
        }
        else{
            sol[i]=u[p];
            printf("u[%d]= %f\n", i+1, sol[i]);
            p += 1;
        }
    }
    printf("\n");

    //Analytical solution (ua)
    printf("Analytical solution: \n");
    for(int i=0; i<nd; i++){
        ua[i]=(coor[i]*coor[i]-coor[i])/2;
        printf("u[%d]: %f\n", i+1, ua[i]);
    }
    return 0;
}
