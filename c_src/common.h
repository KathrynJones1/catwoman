/* catwoman: a batman extension to generate morning/evening terminator transit lightcurves
 * Copyright (C) 2019 Kathryn Jones & Néstor Espinoza
 *
 * This program incorporates a modified version of the batman package: fast computation of exoplanet transit light curves
 * Copyright (C) 2015 Laura Kreidberg
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

//#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
//#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Must be defined in the C file that includes this header. */
inline double intensity(double x, double* args);


double area(double d, double x, double R, double theta);

double find_theta(double phi, double d, double b, double x);

double area(double d, double x, double R, double theta)
{
        /*
        Returns area of an overlapping semi-circle with radii R and circle x; separated by a distance d
        */ 
	//define and initialise variables
	int sit = 12;
        double A = 1, x_minpos=1,x_minneg=1, a1=1,b1=1,a2=1,b2=1,a3=1,b3=1;

	A = (x*x-R*R-d*d)/(2*d);
        x_minpos = (-1)*d*cos(theta)+sqrt(fabs(x*x-d*d*sin(theta)*sin(theta)));
	x_minneg = (-1)*d*cos(theta)-sqrt(fabs(x*x-d*d*sin(theta)*sin(theta)));
	a1 = (-1)*d*cos(theta) + sqrt(fabs(x*x - d*d*sin(theta)*sin(theta)));
        b1 = 0.0;
        b2 = -A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A));  //theta should be left raw here (no fabs)
	a2 = b2*tan(theta)+A/cos(theta);      //theta should be left raw here (no fabs)
        a3 = (-1)*d*cos(theta) + x;
        b3 = d*sin(theta);   //theta should be left raw here (no fabs)
	double lim = pow(10,-9);

	if ((fabs(d * sin(theta)) < x)&&(fabs(fabs(d*sin(theta))-x)>=lim)) {   //circle intersects with y=0 line    
               	if ((fabs(x_minpos)<=R||(fabs(fabs(x_minpos)-R)<lim)) && fabs(x_minneg)>R && fabs(fabs(x_minneg)-R)>=lim) {  //if there is only one intersection on base semi-circle
                        if (theta>=0||(fabs(theta)<lim)){                                                                   //if theta is >=0
				if (b3>=b2||(fabs(b3-b2)<lim)){
                               		sit = 1;    //(a-1)
                       		} else if ((b3<b2)&&(fabs(b3-b2)>=lim) && (b3>=b1||(fabs(b3-b1)<lim)) && (a2>=a1||(fabs(a2-a1)<lim))){
                               		sit = 2;   //(a-2)
                       		} else if ((b3<b2)&&(fabs(b3-b2)>=lim) && (b3>=b1||(fabs(b3-b1)<lim)) && (a2<a1)&&(fabs(a1-a2)>=lim)){
                               		sit = 3;    //(a-3)
				}
			} else if (theta<0&&(fabs(theta)>=lim)){ 
                       		if ((b3<=b1||(fabs(b3-b1)<lim)) && (a2<=a1||(fabs(a2-a1)<lim))){
                               		sit = 4;  //(b)
                	       	} 
			}
               	}else if ((fabs(x_minpos)<=R||(fabs(fabs(x_minpos)-R)<lim)) && (fabs(x_minneg)<=R||(fabs(fabs(x_minneg)-R)<lim))){ //if there are two intersections on base semi-circle

                       	if ((fabs(R)>=fabs(A)||fabs(fabs(R)-fabs(A))<lim)&&((((-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))<lim)) && (((-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A))))<lim))) {  //if two intersections on upper part of semi-circle 

                               	sit = 7;  //(c-3) 

                       	} else if ((fabs(A)>fabs(R) && fabs(fabs(A)-fabs(R))>=lim )||(((-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))<-lim) && ((-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A)))<-lim))) { //if no intersections on upper part of semi-circle
				sit = 5; //(c-1)
			}
               	}else{ //if (fabs(x_minpos)>R && (fabs(x_minpos)-R)>=lim && fabs(x_minneg)>R && fabs(x_minneg)-R>=lim) {  //if no intersections on base of semi-circle

			if ((fabs(R)>=fabs(A)||fabs(fabs(R)-fabs(A))<lim)&&((((-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))<lim)) && (((-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A))))<lim))) {  //if two intersections on upper part of semi-circle
				if (theta>=0||fabs(theta)<lim){
					sit = 9;   //intersection of two circles
				} else if (theta<0  &&fabs(theta)>=lim){                               		
					sit = 6;   //(c-2)                
				}
			} else if ((fabs(A)>fabs(R) && fabs(fabs(A)-fabs(R))>=lim )||(((-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))<-lim) && ((-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A)))<-lim))) { //if no intersections of upper part of semi-circle                         
				if ((d<x)&& fabs(d-x)>=lim) {
					sit = 11; //semi-circle completely inside circle
				} else if ((d>x)||(fabs(d-x)>=lim)) {
					sit = 8;  //no intersection
				}
			}
		}
       } else if (theta<0){ 
	      	sit = 8;  //no intersection
       } else if ((fabs(R)>=fabs(A)||fabs(fabs(R)-fabs(A))<lim)&&((((-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin(theta)+cos(theta)*sqrt(fabs(R*R-A*A)))<lim)) && (((-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A)))>=0)||(fabs(-A*sin(theta)-cos(theta)*sqrt(fabs(R*R-A*A))))<lim))){ // if two intersections in upper part of semi-circle	
		sit = 9;  //intersection of two circles
        } else if ((d+x<=R)||(fabs(d+x-R)<lim)) {
                sit = 10;  //circle completely inside semi-circle
	} else if ((d>=x+R)||(fabs(d-(x+R))<lim)) {
		sit = 8; //no intersection
	}
	
	if ((d==0||fabs(d)<lim) && (x==R||fabs(x-R)<lim)){
                sit = 11 ; //semi-circle completely overlaps with circle
        }
	
 	
	switch(sit) {
                case 1:{   //(a-1)
			double h_a2 = asin(a2/R) + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
			double f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a1+d*cos(theta))*(a1+d*cos(theta))))+asin((a1+d*cos(theta))/x);
                        double f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a2+d*cos(theta))*(a2+d*cos(theta))))+asin((a2+d*cos(theta))/x);
                        double delta_f = f_a1 - f_a2;
			double ret = (M_PI*R*R/2 - (x*x/2)*delta_f - d*sin(fabs(theta))*(a2-a1) - (R*R/2)*(M_PI/2-h_a2));
			if (isnormal(ret)||ret == 0){
				return ret;
			} else{
				double c1 = 1.0;
                        	double c4 = 1.0;
				double c6 = 1.0;

				if ((a2/R)<=1&&(a2/R)>=-1){
                                        c1 = asin(a2/R);
                                } else if  ((a2/R)>1) {
                                        c1 = M_PI/2;
                                } else if ((a2/R)<-1) {
                                        c1 = -M_PI/2;
                                }


				if (((a1+d*cos(theta))/x)<=1&&((a1+d*cos(theta))/x)>=-1){
                                        c4 = asin((a1+d*cos(theta))/x);
                    	        } else if  (((a1+d*cos(theta))/x)>1) {
                                        c4 = M_PI/2;
                                } else if (((a1+d*cos(theta))/x)<-1) {
                                        c4 = -M_PI/2;
                                }

				
				if (((a2+d*cos(theta))/x)<=1&&((a2+d*cos(theta))/x)>=-1){
                                        c6 = asin((a2+d*cos(theta))/x);
                                } else if  (((a2+d*cos(theta))/x)>1) {
                                        c6 = M_PI/2;
                                } else if (((a2+d*cos(theta))/x)<-1) {
                                        c6 = -M_PI/2;
                                }


				h_a2 = c1 + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
				f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a1+d*cos(theta))*(a1+d*cos(theta))))+c4;
	                        f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a2+d*cos(theta))*(a2+d*cos(theta))))+c6; //asin((a2+d*cos(theta))/x);
        	                delta_f = f_a1 - f_a2;
				ret = (M_PI*R*R/2 - (x*x/2)*delta_f - d*sin(fabs(theta))*(a2-a1) - (R*R/2)*(M_PI/2-h_a2));
				
				return ret;
			}
                        break; }

                case 2 :{    //(a-2)
			double h_a2 = asin(a2/R) + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
                        double alpha = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));
			
			double ret = ((R*R/2)*(h_a2+M_PI/2)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha)));
                        
			if (isnormal(ret)||ret == 0){
                                return ret;
                        } else{
                        	double c1;
				double c2;

				if ((a2/R)<=1&&(a2/R)>=-1){
                                        c1 = asin(a2/R);
                                } else if  ((a2/R)>1) {
                                        c1 = M_PI/2;
                                } else if ((a2/R)<-1) {
                                        c1 = -M_PI/2;
                                }

                       	
				if ((1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))>=-1 && (1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))<=1){
                                        c2 = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));
                                } else if ((1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))>1) {
                                        c2 = 0.0;
                                } else if ((1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))<-1){
                                        c2 = M_PI;
                                }
	
				
				h_a2 = c1 + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
                        	alpha = c2;
				ret = ((R*R/2)*(h_a2+M_PI/2)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha)));
				
				return ret;
		
			}break;}
                
		case 3 :{   //(a-3)
                        double h_a2 = asin(a2/R) + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
                        double alpha = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));
			
			double ret = ((R*R/2)*(h_a2+M_PI/2)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha)));

                        if (isnormal(ret)||ret == 0){
                                return ret;
                        } else{
                                double c1;
                                double c2;

				if ((a2/R)<=1&&(a2/R)>=-1){
                                        c1 = asin(a2/R);
                                } else if  ((a2/R)>1) {
                                        c1 = M_PI/2;
                                } else if ((a2/R)<-1) {
                                        c1 = -M_PI/2;
                                }


                                if ((1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))>=-1 && (1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))<=1){
                                        c2 = acos(1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x));
                                } else if ((1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))>1) {
                                        c2 = 0.0;
                                } else if ((1-((a2-a1)*(a2-a1)+b2*b2)/(2*x*x))<-1){
                                        c2 = M_PI;
                                }

                                h_a2 = c1 + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
                                alpha = c2;
                                ret = ((R*R/2)*(h_a2+M_PI/2)-(b2*(a2-a1))/2 + (x*x/2)*(alpha - sin(alpha)));

                                return ret;
                        
			}break;}
                
		case 4 :{     //(b) 
                        double f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a1+d*cos(theta))*(a1+d*cos(theta))))+asin((a1+d*cos(theta))/x);
                        double f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a2+d*cos(theta))*(a2+d*cos(theta))))+asin((a2+d*cos(theta))/x);
                        double delta_f = f_a1 - f_a2;
                        double h_a2 = asin(a2/R) + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
			double ret = ((R*R/2)*(h_a2+M_PI/2) + (x*x/2)*delta_f + d*sin(fabs(theta))*(a2-a1));

			if (isnormal(ret)||ret == 0){
                                return ret;
                        } else{
				double c1;
				double c2;
				double c3;
				
				if (((a1+d*cos(theta))/x)<=1&&((a1+d*cos(theta))/x)>=-1){
					c1 = asin((a1+d*cos(theta))/x);
				} else if  (((a1+d*cos(theta))/x)>1) {
					c1 = M_PI/2;
				} else if (((a1+d*cos(theta))/x)<-1) {
					c1 = -M_PI/2;
				}

				if (((a2+d*cos(theta))/x)<=1&&((a2+d*cos(theta))/x)>=-1){
                                        c2 = asin((a2+d*cos(theta))/x);
                                } else if  (((a2+d*cos(theta))/x)>1) {
                                        c2 = M_PI/2;
                                } else if (((a2+d*cos(theta))/x)<-1) {
                                        c2 = -M_PI/2;
                                }

				if ((a2/R)<=1&&(a2/R)>=-1){
                                        c3 = asin(a2/R);
                                } else if  ((a2/R)>1) {
                                        c3 = M_PI/2;
                                } else if ((a2/R)<-1) {
                                        c3 = -M_PI/2;
                                }

				f_a1 = ((a1+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a1+d*cos(theta))*(a1+d*cos(theta))))+c1;
				f_a2 = ((a2+d*cos(theta))/(x*x))*sqrt(fabs(x*x-(a2+d*cos(theta))*(a2+d*cos(theta))))+c2;
				h_a2 = c3 + (a2/R)*sqrt(fabs(1-(a2/R)*(a2/R)));
                        	ret = ((R*R/2)*(h_a2+M_PI/2) + (x*x/2)*delta_f + d*sin(fabs(theta))*(a2-a1));
                        
				return ret;
			
			}break;}
	
                case 5 :{    //(c-1)
                        double s = sqrt(fabs(x*x-d*d*sin(theta)*sin(theta)));
                        double z = d*sin(-(theta));
			double ret = (x*x*acos(z/x) - s*z);
		
			if (isnormal(ret)||ret == 0){
                                return ret;
                        } else{
				double c1;
				
				if ((z/x<=1)&&(z/x>=-1)){
                                        c1 = acos(z/x);
                                }else if (z/x>1){
                                        c1 = 0.0;
                                }else if (z/x<-1){
                                        c1 = M_PI;
                                }

				
				ret = (x*x*c1 - s*z);
				return ret;
			
			}break;}

                case 6 :{  //(c-2)
                        double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(fabs(w));
                        double ret= (A_int - M_PI*R*R/2);
			
			if (isnormal(ret)||ret == 0){
					return ret;
                        } else{
				double c1;
				double c2;
				
				if ((u<=1)&&(u>=-1)){
                                        c1 = acos(u);
                                }else if (u>1){
                                        c1 = 0.0;
                                }else if (u<-1){
                                        c1 = M_PI;
                                }

                                if ((v<=1)&&(v>=-1)){
                                        c2 = acos(v);
                                }else if (v>1){
                                        c2 = 0.0;
                                }else if (v<-1){
                                        c2 = M_PI;
                                }
	
				
				ret = x*x*c1+R*R*c2-0.5*sqrt(fabs(w))- M_PI*R*R/2;
				
				return ret;
				
                        }break;}

                case 7 :{ //(c-3)
                        double s = sqrt(fabs(x*x-d*d*sin(theta)*sin(theta)));
                        double z = d*sin(fabs(theta));
                        double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(fabs(w));
                        double A1 = x*x*acos(z/x) - s*z;
                        double ret=(A_int - A1);
			
			if (isnormal(ret)||ret == 0){
                                return ret;
                        } else{
				double c1;
				double c2;
				double c3;		
				
				if ((u<=1)&&(u>=-1)){
                                        c1 = acos(u);
                                }else if (u>1){
                                        c1 = 0.0;
                                }else if (u<-1){
                                        c1 = M_PI;
                                }
                                
                                if ((v<=1)&&(v>=-1)){
                                        c2 = acos(v);
                                }else if (v>1){
                                        c2 = 0.0;
                                }else if (v<-1){
                                        c2 = M_PI;
                                }
		
				if ((z/x<=1)&&(z/x>=-1)){
                                        c3 = acos(z/x);
                                }else if (z/x>1){
                                        c3 = 0.0;
                                }else if (z/x<-1){
                                        c3 = M_PI;
                                }

				A_int = x*x*c1+R*R*c2-0.5*sqrt(fabs(w));
				A1 = x*x*c3 - s*z;
				ret = (A_int - A1);
				
				return ret;
                        }break;}
                case 8 :{ //no intersection
                        return 0.0;
                        }break;
                case 9 :{ //intersection of two circles
                        double u = (d*d+x*x-R*R)/(2*d*x);
                        double v = (d*d+R*R-x*x)/(2*d*R);
                        double w = (-d+x+R)*(d+x-R)*(d-x+R)*(d+x+R);
                        double A_int = x*x*acos(u)+R*R*acos(v)-0.5*sqrt(fabs(w));
                        
			double ret =  A_int;
			
			if (isnormal(ret)||ret == 0){
                                return ret;
                        } else{
				double c1;	
				double c2;
				
				if ((u<=1)&&(u>=-1)){
					c1 = acos(u);	
				}else if (u>1){
					c1 = 0.0;
				}else if (u<-1){
					c1 = M_PI;
				}
			
				if ((v<=1)&&(v>=-1)){
                                        c2 = acos(v);
                                }else if (v>1){
                                        c2 = 0.0;
                                }else if (v<-1){
                                        c2 = M_PI;
                                }	
				
				A_int = x*x*c1+R*R*c2-0.5*sqrt(fabs(w));
				ret =  A_int;
				return ret;
	
                        }break;}
                case 10 :{ //circle completely inside semi-circle
                        return M_PI*x*x;
                        }break;

        	case 11: { //semi-circle completely inside circle
			return M_PI*R*R/2;
			}break;
		
	        default :
                        return 0.0;
                        break;
        }

} 

double find_theta(double phi, double d, double b, double x)
{
	/* This function finds theta for a given phi, d and b */

	bool before_mid = (x < 0.0);
	
	double lim = pow(10,-9);
	double theta = 0.0;
	if ((b>=0)||(fabs(b)<lim)){	
		if (before_mid){            
			if (phi >= acos(b/d) - lim) {
                        	theta = -(M_PI - asin(b/d) - phi);
                	} else if (phi < acos(b/d) - lim) {
                        	theta = -phi -asin(b/d);
                	}
                
		} else {
			if (phi <= -acos(b/d) + lim) {
                		theta = -(M_PI+phi - asin(b/d));
                	} else if (phi > -acos(b/d) + lim) {
                        	theta = -(-phi +asin(b/d));
                	}
		}
	} else {
		if (before_mid){ 
			if (phi > -acos(fabs(b)/d) + lim) { 
				theta = -(phi - asin(fabs(b)/d));
			} else if (phi <= -acos(fabs(b)/d) + lim) { 
				theta = M_PI - asin(fabs(b)/d) + phi;
			}
		} else {
			if (phi < acos(fabs(b)/d) - lim) { 
				theta = -(-phi - asin(fabs(b)/d));
			} else if (phi >= acos(fabs(b)/d) - lim) { 
				theta = M_PI - asin(fabs(b)/d) - phi;
                        }

		}
	}
	return theta;
}

void calc_limb_darkening(double* f_array, double* d_array, int N, double rprs, double fac, int nthreads, double* intensity_args, double* phi_array, double* b_array, double* x_array, double rp2, bool twoc)
{
	/*
		This function takes an array of sky distances (d_array) of length N, computes stellar intensity by calling intensity with
		intensity_args, and puts the results in f_array.  To use this function, include this file in a .c file and implement
		the intensity function within that .c file.

		The proper way of implementing this function is to accept a pointer to the intensity function.  Unfortunately, few
		compilers that implement OpenACC support function pointers, so this approach is not yet possible.
	*/

	#if defined (_OPENMP) && !defined(_OPENACC)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif
	#if defined (_OPENACC)
	#pragma acc parallel loop copyout(f_array[:N]) present(intensity_args)
	#elif defined (_OPENMP)
	#pragma omp parallel for
	#endif
	for(int i = 0; i < N; i++)
	{
		double b = b_array[i];
		double phi = phi_array[i];

		double d = d_array[i];
		double x = 0.1;
		double x_in = MAX(MIN(d - rp2, d - rprs), 0.);						//lower bound for integration
		double x_out = MIN(MAX(d + rp2, d + rprs), 1.0);					//upper bound for integration
		if(x_in >= 1.) f_array[i] = 1.0;					//flux = 1. if the planet is not transiting
		else if(x_out - x_in < 1.e-9) f_array[i] = 1.0;				//pathological case	
		else
		{
			double delta = 0.;						//variable to store the integrated intensity, \int I dA
			x = x_in;						//starting radius for integration
			double dx = fac*acos(x); 					//initial step size
			x += dx;						//first step
			double A_i = 0.;						//initial area
		
			double theta=0.0, theta2=0.0, lim = pow(10,-9);
			
			theta = find_theta(phi, d, b, x_array[i]);
			
										//adjusting theta to within the definition of the system so the area() eqns work
			if (theta>=(M_PI/2+lim)){
				theta = M_PI - theta;
			} else if (theta<=(-M_PI/2-lim)){
				theta = -(M_PI + theta);
			}
			
										//Finds theta for the second semi-circle
			theta2 = -theta;

			while(x < x_out)
			{
				double A_f = area(d, x, rprs, theta);				//calculates area of overlapping circles
				if (twoc){
					double A_f2 = area(d, x, rp2, theta2);
					A_f = A_f + A_f2;
					
				}
				double I_temp = intensity(x - dx/2., intensity_args); 	//intensity at the midpoint
				delta += (A_f - A_i)*I_temp;				//increase in transit depth for this integration step
				dx = fac*acos(x);  				//updating step size
				x = x + dx;					//stepping to next element
				A_i = A_f;					//storing area
			}

			dx = x_out - x + dx;  					//calculating change in radius for last step  FIXME
			x = x_out;						//final radius for integration
			double A_f = area(d, x, rprs, theta);					//area for last integration step
			if (twoc){
				double A_f2 = area(d, x, rp2, theta2);
				A_f = A_f + A_f2;
			}
			
			double I_temp = intensity(x - dx/2., intensity_args); 		//intensity at the midpoint
			delta += (A_f - A_i)*I_temp;					//increase in transit depth for this integration step
			f_array[i] = 1.0 - delta;	//flux equals 1 - \int I dA
		}
	}
}
