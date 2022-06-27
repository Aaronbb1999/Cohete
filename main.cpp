#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define G 6.67E-11
#define Mt 5.9736E24
#define Ml 0.07349E24
#define dtl 3.844E8
#define w 2.6617E-6
#define Rt 6.378160E6
#define Rl 1.7374E6
#define m 4.5932E4
#define pi 3.14159265359

//Reescalado
#define delta (1.0*G*Mt/(powf(dtl,3)))
#define mu (1.0*Ml/Mt)

void Kutta (double kutta[], double vector[], double tiempo, double h, double aux2);

int main (void)
{
    //Definicion de variables, parametros y ficheros
    int i,j,l;
    double t=0.0, h=3.0;
    int N=4, n=3, iter=500;
    FILE *f, *fcons;

    f=fopen("posiciones.txt","w");
    fcons=fopen("conservacion.txt","w");

    double x[n], y[n]; //Posicion de los 3 cuerpos
    double mov[N]; //Variables de movimiento del cohete
    double k1[N], k2[N], k3[N], k4[N]; //Coeficientes de Kutta para cada variable de las ecuaciones
    

    //Posiciones iniciales Tierra y Luna.
    x[0]=0;
    y[0]=0;

    x[1]=1;
    y[1]=0;

    //Condicion del despegue del cohete
    double angulo, vi;
    vi=6800;
    angulo=pi/12; 

    mov[0]=Rt/dtl; //r
    mov[1]=0; //phi
    mov[2]=vi*cos(angulo)/dtl; //Pr
    mov[3]=mov[0]*vi*sin(angulo)/dtl; //Pphi

    x[2]=mov[0]*cos(mov[1]);
    y[2]=mov[0]*sin(mov[1]); //En cartesianas

    //Añadimos condiciones iniciales al fichero
    for(j=0; j<n; j++)
    {
        fprintf(f, "%lf, %lf\n", x[j], y[j]);
    }
    fprintf(f,"\n");

    //Conservacion de Hprima
    double H, factor;
    double r, phi, pr, pphi;

    r=(mov[0])*dtl;
    phi=mov[1];
    pr=(mov[2])*m*dtl;
    pphi=(mov[3])*m*dtl*dtl;

    factor=powf(r*r+dtl*dtl-2*r*dtl*cos(-w*t+phi),1/2);
    H= 0.5*((pr*pr)/m+(pphi*pphi)/(m*r*r))-(G*Mt*m/r)-G*Ml*m/(factor)-w*pphi;
    fprintf(fcons, "%lf\n",H);


    //Bucle principal
    double rprima, aux1;
    double aux[N];

    for(j=0;j<iter;j++)
    {
        for(i=0;i<250;i++)
        {
            //Calculo de k1
            rprima=sqrt(1+mov[0]*mov[0]-2*mov[0]*cos(mov[1]-w*t));
            aux1=powf(rprima,3);
            Kutta(k1,mov,t,h,aux1); 
        
            
            //Calculo de k2
            t=t+h/2;
            for(l=0;l<N;l++)
            {
               aux[l]=mov[l]+k1[l]/2;
            }   
            rprima=sqrt(1+mov[0]*mov[0]-2*mov[0]*cos(mov[1]-w*t));
            aux1=powf(rprima,3);
            Kutta(k2,aux,t,h,aux1);


            //Calculo de k3
            for(l=0;l<N;l++)
            {
               aux[l]=mov[l]+k2[l]/2;
            }   
            Kutta(k3,aux,t,h,aux1);

            //Calculo de k4
            t=t+h/2;
            for(l=0;l<N;l++)
            {
               aux[l]=mov[l]+k3[l]/2;
            }   
            rprima=sqrt(1+mov[0]*mov[0]-2*mov[0]*cos(mov[1]-w*t));
            aux1=powf(rprima,3);
            Kutta(k4,aux,t,h,aux1);

            
            //Nuevas coordenadas
            for(l=0;l<N;l++)
            {
                mov[l]=mov[l]+(k1[l]+2*(k2[l]+k3[l])+k4[l])/6.0;
            }
        

        }

        //Nuevas coordenadas de la luna
        x[1]=cos(w*t);
        y[1]=sin(w*t);

        //Nuevas coordenadas del cohete
        x[2]=mov[0]*cos(mov[1]); //coordenada x del cohete
        y[2]=mov[0]*sin(mov[1]); //coordenada y del cohete

        for (l = 0; l < n; l++)
        {
            fprintf(f, "%lf, %lf\n", x[l], y[l]);
        }
        fprintf(f,"\n");
        

        //Conservacion de Hprima
        r=(mov[0])*dtl;
        phi=mov[1];
        pr=(mov[2])*m*dtl;
        pphi=(mov[3])*m*dtl*dtl;

        factor=powf(r*r+dtl*dtl-2*r*dtl*cos(-w*t+phi),1/2);
        H= 0.5*((pr*pr)/m+(pphi*pphi)/(m*r*r))-(G*Mt*m/r)-G*Ml*m/(factor)-w*pphi;
        fprintf(fcons, "%lf\n",H);
    }

    fclose(f);
    fclose(fcons);

    return 0;
}

    void Kutta (double kutta[], double vector[], double tiempo, double h, double aux2)
    {
        kutta[0]=h*vector[2];
        kutta[1]=h*vector[3]/(vector[0]*vector[0]);
        kutta[2]=h*vector[3]*vector[3]*(1/powf(vector[0],3))-delta*(1/(vector[0]*vector[0]))-delta*(mu/aux2)*(vector[0]-cos(vector[1]-w*tiempo));
        kutta[3]=-h*delta*mu*vector[0]*sin(vector[1]-w*tiempo)/aux2;

        return;

    }
