#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);
//zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__)

    {
        Initialize();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double Particles);    //::DONE
    void Get_Wannier_function(int band);

    void Pushing_points_in_FBZ(int n1, int n2, int & n1_eff, int & n2_eff, int L1_,int  L2_);
    void Get_MaximallyLocalizedWannier_functions(int band);
    void Print_Moire_Potential();
    double TotalDensity();   //::DONE
    double E_QM();   //::DONE
    double NIEnergy(double kx_val, double ky_val);

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE

    void Update_Bloch_States_Using_Projection(Mat_1_Complex_doub & Psi_state_,int space_slices,double rx_min,double d_rx,double ry_min,double d_ry);
    double Gaussian(double rx_, double ry_, double Rx_center, double Ry_center, double std_dev);

    int convert_jm_to_int(string jm_val);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    int ns_, l1_, l2_;
    double kx_, ky_;
    double k_plusx, k_minusx, k_plusy, k_minusy;

    double Wnr_center_x, Wnr_center_y;
    double Wnr_center_1, Wnr_center_2;

    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;

    //real space  effective H params
    int L1_eff, L2_eff;
    Mat_2_Complex_doub Tij;
    Mat_2_Complex_doub Uij;
    Mat_2_Complex_doub Aij;
    Mat_2_Complex_doub Xij;

};


void Hamiltonian::Initialize(){

    ns_=Parameters_.ns;
    l1_=Parameters_.Grid_L1;
    l2_=Parameters_.Grid_L2;

    int space=ns_;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    eigs_saved_.resize(space);

    k_plusx = (-1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    k_plusy = (-1.0/3.0)*(2.0*PI/Parameters_.a_moire);
    k_minusx = (-1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    k_minusy = (1.0/3.0)*(2.0*PI/Parameters_.a_moire);


    //real space  effective H params
    L1_eff=4;L2_eff=4;
    Tij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Tij[i].resize(L1_eff*L2_eff);
    }
    Uij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Uij[i].resize(L1_eff*L2_eff);
    }

    Aij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Aij[i].resize(L1_eff*L2_eff);
    }

	
    Xij.resize(L1_eff*L2_eff);
    for(int i=0;i<L1_eff*L2_eff;i++){
        Xij[i].resize(L1_eff*L2_eff);
    }


    Wnr_center_x = (2.0/sqrt(3.0))*Parameters_.a_moire;
    Wnr_center_y = (0.0)*Parameters_.a_moire;


    Wnr_center_1 = (2.0/3.0)*Parameters_.a_moire;
    Wnr_center_2 = (2.0/3.0)*Parameters_.a_moire;


} // ----------

double Hamiltonian::TotalDensity(){

    double n1=0.0;
    /*
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    */
    return n1;

} // ----------



double Hamiltonian::E_QM(){

    return 0.0;

} // ----------

double Hamiltonian::NIEnergy(double kx_val, double ky_val){

    double energy_;
    //energy_ = -1.0*(Parameters_.RedPlanckConst*Parameters_.RedPlanckConst*(kx_val*kx_val  + ky_val*ky_val))*(0.5/Parameters_.MStar);
    energy_ = -1.0*(((3.809842*1000)/Parameters_.MStar)*(kx_val*kx_val  + ky_val*ky_val));

    return energy_;
}

double Hamiltonian::GetCLEnergy(){

    return 0.0;

} // ----------


int Hamiltonian::convert_jm_to_int(string jm_val){

    int val;
    if(jm_val=="3by2_m3by2"){val=0;}
    if(jm_val=="3by2_3by2"){val=1;}
    if(jm_val=="3by2_m1by2"){val=2;}
    if(jm_val=="3by2_1by2"){val=3;}
    if(jm_val=="1by2_m1by2"){val=4;}
    if(jm_val=="1by2_1by2"){val=5;}
    return val;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<Ham_.n_row();i++) {
        for(int j=0;j<Ham_.n_row();j++) {
            if(
                    abs(Ham_(i,j) - conj(Ham_(j,i)))>0.00001
                    ) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(
                        abs(Ham_(i,j) - conj(Ham_(j,i)))<0.00001
                        ); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}

void Hamiltonian::Print_Moire_Potential(){

string file_out =  "Moire_Potential.txt";
ofstream fl_out(file_out.c_str());


Mat_1_doub VParams;
    VParams.resize(3);
    VParams[0]=Parameters_.V1_param;
    VParams[1]=Parameters_.V2_param;
    VParams[2]=Parameters_.V3_param;



   double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);


    double rx_min, ry_min, rx_max, ry_max, d_rx, d_ry;
    rx_min=-3.0*Parameters_.a_moire;
    ry_min=-3.0*Parameters_.a_moire;
    rx_max=3.0*Parameters_.a_moire;
    ry_max=3.0*Parameters_.a_moire;
    int r_ind;
    int space_slices=200;
    d_rx=(rx_max-rx_min)/(space_slices);
    d_ry=(ry_max-ry_min)/(space_slices);



    Mat_2_int neigh_G_shell_1, neigh_G_shell_2;
    neigh_G_shell_1.resize(3);neigh_G_shell_2.resize(3);
    for(int s=0;s<3;s++){
        neigh_G_shell_1[s].resize(6);
        neigh_G_shell_2[s].resize(6);
        }

        //shell=1
        neigh_G_shell_1[0][0]=1;neigh_G_shell_2[0][0]=0;
        neigh_G_shell_1[0][1]=0;neigh_G_shell_2[0][1]=1;
        neigh_G_shell_1[0][2]=-1;neigh_G_shell_2[0][2]=1;
        neigh_G_shell_1[0][3]=-1;neigh_G_shell_2[0][3]=0;
        neigh_G_shell_1[0][4]=0;neigh_G_shell_2[0][4]=-1;
        neigh_G_shell_1[0][5]=1;neigh_G_shell_2[0][5]=-1;

        //shell=2
        neigh_G_shell_1[1][0]=1;neigh_G_shell_2[1][0]=1;
        neigh_G_shell_1[1][1]=-1;neigh_G_shell_2[1][1]=2;
        neigh_G_shell_1[1][2]=-2;neigh_G_shell_2[1][2]=1;
        neigh_G_shell_1[1][3]=-1;neigh_G_shell_2[1][3]=-1;
        neigh_G_shell_1[1][4]=1;neigh_G_shell_2[1][4]=-2;
        neigh_G_shell_1[1][5]=2;neigh_G_shell_2[1][5]=-1;

        //shell=3
        neigh_G_shell_1[2][0]=2;neigh_G_shell_2[2][0]=0;
        neigh_G_shell_1[2][1]=0;neigh_G_shell_2[2][1]=2;
        neigh_G_shell_1[2][2]=-2;neigh_G_shell_2[2][2]=2;
        neigh_G_shell_1[2][3]=-2;neigh_G_shell_2[2][3]=0;
        neigh_G_shell_1[2][4]=0;neigh_G_shell_2[2][4]=-2;
        neigh_G_shell_1[2][5]=2;neigh_G_shell_2[2][5]=-2;




double k_x_, k_y_;
double dx, dy;
complex<double> Del_r;

for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
  dx= rx_min + rx_ind*d_rx;
  dy= ry_min + ry_ind*d_ry;

  Del_r=0.0;
for (int s=0;s<1;s++){

                        for(int neigh_ind=0;neigh_ind<6;neigh_ind++){

			if(neigh_ind==1 || neigh_ind==3  || neigh_ind==5){
                        k_x_ = neigh_G_shell_1[s][neigh_ind]*b1x_ + neigh_G_shell_2[s][neigh_ind]*b2x_ ;
			k_y_ = neigh_G_shell_1[s][neigh_ind]*b1y_ + neigh_G_shell_2[s][neigh_ind]*b2y_ ;
			Del_r += VParams[s]*(  exp(iota_complex*( (k_x_*dx) + (k_y_*dy) + (Parameters_.Psi_param)   ))   
   					     + exp(-1.0*iota_complex*( (k_x_*dx) + (k_y_*dy) + (Parameters_.Psi_param)   ))	);		
			}
                    
                        }
                     }


fl_out<<dx<<"   "<<dy<<"   "<<Del_r.real()<<"   "<<Del_r.imag()<<endl;

                     }

fl_out<<endl;
         }





}

void Hamiltonian::HTBCreate(){


    Ham_.resize(ns_,ns_);
    double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);

    int Bottom_, Top_;
    Bottom_=0;Top_=1;

    //l1_/2,l2_/2 is the k-point

    Mat_2_int neigh_G_shell_1, neigh_G_shell_2;
    neigh_G_shell_1.resize(3);neigh_G_shell_2.resize(3);
    for(int s=0;s<3;s++){
	neigh_G_shell_1[s].resize(6);
	neigh_G_shell_2[s].resize(6);
	}

	//shell=1
	neigh_G_shell_1[0][0]=1;neigh_G_shell_2[0][0]=0;
	neigh_G_shell_1[0][1]=0;neigh_G_shell_2[0][1]=1;
	neigh_G_shell_1[0][2]=-1;neigh_G_shell_2[0][2]=1;
	neigh_G_shell_1[0][3]=-1;neigh_G_shell_2[0][3]=0;
	neigh_G_shell_1[0][4]=0;neigh_G_shell_2[0][4]=-1;
	neigh_G_shell_1[0][5]=1;neigh_G_shell_2[0][5]=-1;
	
	//shell=2
	neigh_G_shell_1[1][0]=1;neigh_G_shell_2[1][0]=1;
        neigh_G_shell_1[1][1]=-1;neigh_G_shell_2[1][1]=2;
        neigh_G_shell_1[1][2]=-2;neigh_G_shell_2[1][2]=1;
        neigh_G_shell_1[1][3]=-1;neigh_G_shell_2[1][3]=-1;
        neigh_G_shell_1[1][4]=1;neigh_G_shell_2[1][4]=-2;
        neigh_G_shell_1[1][5]=2;neigh_G_shell_2[1][5]=-1;	

	//shell=3
        neigh_G_shell_1[2][0]=2;neigh_G_shell_2[2][0]=0;
        neigh_G_shell_1[2][1]=0;neigh_G_shell_2[2][1]=2;
        neigh_G_shell_1[2][2]=-2;neigh_G_shell_2[2][2]=2;
        neigh_G_shell_1[2][3]=-2;neigh_G_shell_2[2][3]=0;
        neigh_G_shell_1[2][4]=0;neigh_G_shell_2[2][4]=-2;
        neigh_G_shell_1[2][5]=2;neigh_G_shell_2[2][5]=-2;




    Mat_1_doub VParams;
    VParams.resize(3);
    VParams[0]=Parameters_.V1_param;
    VParams[1]=Parameters_.V2_param;
    VParams[2]=Parameters_.V3_param;
   

    double kx_local, ky_local;

    Mat_2_doub val_mat;
    val_mat.resize(l1_);
    for(int n1=0;n1<val_mat.size();n1++){
        val_mat[n1].resize(l2_);
        for(int n2=0;n2<val_mat.size();n2++){
        if(  (n1+n2)<((l1_-1)/2)
          || (n1+n2)>((3*(l1_-1))/2)
         ){
        val_mat[n1][n2]=0.0;
        }
        else{
        val_mat[n1][n2]=1.0;
        }
        }
        }




    int row, col;
    int i1_neigh, i2_neigh;
    for(int i1=0;i1<l1_;i1++){
        for(int i2=0;i2<l2_;i2++){
            kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
            ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
            
                row=Coordinates_.Nbasis(i1, i2, 0);
                

                    //1
                    col = row;
                    Ham_(row,col) += (NIEnergy(kx_local, ky_local) + (0.0*Parameters_.Vz_));


		     for (int s=0;s<1;s++){
		   
			for(int neigh_ind=0;neigh_ind<6;neigh_ind++){
			if(neigh_ind==1 || neigh_ind==3 || neigh_ind==5 ){


			//i1_neigh = (i1 + neigh_G_shell_1[s][neigh_ind] + l1_)%l1_;
                        //i2_neigh = (i2 + neigh_G_shell_2[s][neigh_ind] + l2_)%l2_;

			i1_neigh = i1 + neigh_G_shell_1[s][neigh_ind];
                        i2_neigh = i2 + neigh_G_shell_2[s][neigh_ind];			
			if( (i1_neigh<l1_)  && (i2_neigh<l2_)  && (i1_neigh>=0)  && (i2_neigh>=0)  ){
			col = Coordinates_.Nbasis(i1_neigh, i2_neigh, 0);
                        Ham_(row,col) += VParams[s]*exp(iota_complex*Parameters_.Psi_param);
			}
			}

			if(neigh_ind==0 || neigh_ind==2 || neigh_ind==4 ){

			//i1_neigh = (i1 + neigh_G_shell_1[s][neigh_ind] + l1_)%l1_;
                        //i2_neigh = (i2 + neigh_G_shell_2[s][neigh_ind] + l2_)%l2_;
                        i1_neigh = i1 + neigh_G_shell_1[s][neigh_ind];
                        i2_neigh = i2 + neigh_G_shell_2[s][neigh_ind];
                        if( (i1_neigh<l1_)  && (i2_neigh<l2_)  && (i1_neigh>=0)  && (i2_neigh>=0)  ){
                        col = Coordinates_.Nbasis(i1_neigh, i2_neigh, 0);
                        Ham_(row,col) += VParams[s]*exp(-1.0*iota_complex*Parameters_.Psi_param);
                        }
                        }
		


			}


		     }
            
            
        }
    }



} // ----------

double Hamiltonian::Gaussian(double rx_, double ry_, double Rx_center, double Ry_center, double std_dev){

double val;

val = sqrt(1.0/(PI*std_dev*std_dev))*exp(-1.0*(  ((rx_-Rx_center)*(rx_-Rx_center) + (ry_-Ry_center)*(ry_-Ry_center))/(2.0*std_dev*std_dev)) );

return val;
}


void Hamiltonian::Update_Bloch_States_Using_Projection(Mat_1_Complex_doub & Psi_state_,int space_slices,double r1_min,double d_r1,double r2_min,double d_r2){

Mat_1_Complex_doub phi_state;
phi_state.resize(space_slices*space_slices);


double rx_, ry_;
double std_dev;
complex<double> A, sqrt_Ainv, sqrt_A;

int r_ind;

for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
 rx_ = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
 ry_ = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
		

	std_dev = 0.001*Parameters_.a_moire;
	A += conj(Psi_state_[r_ind])*Gaussian(rx_, ry_, Wnr_center_x, Wnr_center_y, std_dev);


	/*if( ((abs( r1_min+r1_ind*d_r1 - Wnr_center_1 ))<d_r1*0.001) &&  
	    ((abs( r2_min+r2_ind*d_r2 - Wnr_center_2 ))<d_r2*0.001)  
	    ){
	A=conj(Psi_state_[r_ind]);
	}*/
}
}


 for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;

                       // Psi_state_[r_ind] += Ham_(comp,eigen_no)*exp(iota_complex*(kx_local*(rx_min + rx_ind*d_rx) +  ky_local*(ry_min + ry_ind*d_ry) ));

		Psi_state_[r_ind] = Psi_state_[r_ind]*(A)*(1.0/abs(A));
		
}
}


/*for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;
}
}
*/


}


void Hamiltonian::Get_MaximallyLocalizedWannier_functions(int band){



}



void Hamiltonian::Pushing_points_in_FBZ(int n1, int n2, int & n1_eff, int & n2_eff, int L1_,int  L2_){

n1_eff=n1;
n2_eff=n2;


}



void Hamiltonian::Get_Wannier_function(int band){


    double b1x_, b1y_, b2x_, b2y_;
    b1x_=(2.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b1y_=(0.0)*(2.0*PI/Parameters_.a_moire);
    b2x_=(1.0/sqrt(3.0))*(2.0*PI/Parameters_.a_moire);
    b2y_=(1.0)*(2.0*PI/Parameters_.a_moire);

    double kx_local, ky_local;
    complex<double> temp_factor;
    complex<double> checknorm;
    int eigen_no, comp;


    double r1_min, r2_min, r1_max, r2_max, d_r1, d_r2;
     double r1_offset, r2_offset;
//	rx_offset=(1.0/(2.0*sqrt(3.0)))*Parameters_.a_moire;
//	ry_offset=(0.5)*Parameters_.a_moire;
    int space_slices=121;



    r1_min=-3.0*Parameters_.a_moire + Wnr_center_1;
    r2_min=-3.0*Parameters_.a_moire + Wnr_center_2;
    r1_max=3.0*Parameters_.a_moire + Wnr_center_1 ;
    r2_max=3.0*Parameters_.a_moire + Wnr_center_2;
    int r_ind;
    d_r1=(r1_max-r1_min)/(space_slices-1);
    d_r2=(r2_max-r2_min)/(space_slices-1);
    double dis_x, dis_y;
   




//remove later
//double rx_min, ry_min, rx_max, ry_max, d_rx, d_ry;
//double rx_offset, ry_offset;
//---------------


	//For interactions like A (interactioins assisted hoppings)
   assert((L1_eff*Parameters_.a_moire)<(r1_max-r1_min));
   assert((L2_eff*Parameters_.a_moire)<(r2_max-r2_min));



    double qx_min, qy_min, qx_max, qy_max, d_qx, d_qy;
    qx_min=-0.5*PI;//Parameters_.a_moire;
    qy_min=-0.5*PI;//Parameters_.a_moire;
    qx_max=0.5*PI;//Parameters_.a_moire;
    qy_max=0.5*PI;//Parameters_.a_moire;
    int q_slices=50;
    d_qx=(qx_max-qx_min)/(q_slices);
    d_qy=(qy_max-qy_min)/(q_slices);
    double eta_q=0.001;
    double d_sqr=360000;
    double screening=0.0;



    double q_max, d_q,  d_theta;
    q_max=0.1*PI;
    d_q=(q_max)/(q_slices);
    int theta_slices=50;
    d_theta = (2.0*PI)/(theta_slices);



    int L1_,L2_;
    L1_=Parameters_.BZ_L1;
    L2_=Parameters_.BZ_L2;

    int Bottom_, Top_;
    Bottom_=0;Top_=1;

	
    Mat_2_doub val_mat;
    val_mat.resize(space_slices);
    for(int n1=0;n1<val_mat.size();n1++){
	val_mat[n1].resize(space_slices);
	for(int n2=0;n2<val_mat.size();n2++){
	if(  (n1+n2)<((space_slices-1)/2)
	  || (n1+n2)>((3*(space_slices-1))/2)
	 ){
	val_mat[n1][n2]=0.0;
	}
	else{
	val_mat[n1][n2]=1.0;
	}
	}
	}


    



    //For maximally localized Wannier wavefunction
    int NB_=2;
    int L1p2, L2p2;
    L1p2=L1_+2;L2p2=L2_+2;
    Mat_2_Complex_doub U_Bloch_mat; //U_Bloch_mat[k][r]
    U_Bloch_mat.resize((L1p2)*(L2p2));
    for(int i=0;i<(L1p2)*(L2p2);i++){
	U_Bloch_mat[i].resize(space_slices*space_slices);
    }

    Mat_2_Complex_doub M_Mat_old, M_Mat_new, q_mat, R_mat, T_mat, A_mat, S_mat;
    M_Mat_old.resize(L1_*L2_);M_Mat_new.resize(L1_*L2_);
    q_mat.resize(L1_*L2_);R_mat.resize(L1_*L2_);
    T_mat.resize(L1_*L2_);A_mat.resize(L1_*L2_);
    S_mat.resize(L1_*L2_);
    for(int i=0;i<L1_*L2_;i++){
	q_mat[i].resize(6);R_mat[i].resize(6);
	T_mat[i].resize(6);A_mat[i].resize(6);
        S_mat[i].resize(6); 
	M_Mat_old[i].resize(6);  //6 neighbours for triangular reciprocal lattice
        M_Mat_new[i].resize(6);
	}
	

    Mat_1_Complex_doub Unit_mat_old, Unit_mat_new, Delta_W;
    Unit_mat_old.resize(L1_*L2_);
    Unit_mat_new.resize(L1_*L2_);
    Delta_W.resize(L1_*L2_);
  
   double Sigma_old, Sigma_new, d_Sigma, wb_, w_, rx_mean_new, ry_mean_new, rx_mean_old, ry_mean_old, r2_mean_new, r2_mean_old;
   double alpha_;
        
   //------------------------------------------

    Mat_1_Complex_doub Wnr_state_, Psi_state_;
        Wnr_state_.resize(space_slices*space_slices);
	for(int i=0;i<Wnr_state_.size();i++){
	Wnr_state_[i]=0.0;
	}
        Psi_state_.resize(space_slices*space_slices);


    Mat_2_Complex_doub Mq_;
    Mq_.resize(q_slices);
    for(int q1=0;q1<q_slices;q1++){
        Mq_[q1].resize(q_slices);
    }

    Mat_2_Complex_doub Mq_SphC, Mmq_SphC;
    Mq_SphC.resize(q_slices);
    Mmq_SphC.resize(q_slices);
    for(int q_i=0;q_i<q_slices;q_i++){
        Mq_SphC[q_i].resize(theta_slices);
	Mmq_SphC[q_i].resize(theta_slices);
    }

    Mat_3_Complex_doub MqR_SphC;
    MqR_SphC.resize(q_slices);
    for(int q_i=0;q_i<q_slices;q_i++){
        MqR_SphC[q_i].resize(theta_slices);
	for(int th_i=0;th_i<theta_slices;th_i++){
	MqR_SphC[q_i][th_i].resize(L1_eff*L2_eff);
	}
    }



    Mat_2_Complex_doub Vq_;
    Vq_.resize(q_slices);
    for(int q1=0;q1<q_slices;q1++){
        Vq_[q1].resize(q_slices);
    }


    int center_=int(L1_eff/2) + L1_eff*int(L2_eff/2);
    int center_neigh;
    //Tij[center_][center_p1]=0.0;
    eigen_no=(l1_*l2_)-1-band;

    //eigen_no=13;

    for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){


            cout<<"doing "<<n1<<"  "<<n2<<endl;

	    int n1_eff, n2_eff;
	
	    //Pushing_points_in_FBZ(n1, n2, n1_eff, n2_eff, L1_, L2_);

	
            kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3.0)*L1_))  +  n2*(1.0/(sqrt(3.0)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
           



	     HTBCreate();
            char Dflag='V';
            Diagonalize(Dflag);

            //------------------


            //Hopping-real space-------
      
      for(int r1_=0;r1_<L1_eff;r1_++){
                for(int r2_=0;r2_<L2_eff;r2_++){
                    center_neigh = (r1_) + L1_eff*(r2_);
                    dis_x = ((sqrt(3.0)/2.0)*(r1_-(L1_eff/2)) +  (sqrt(3.0)/2.0)*(r2_-(L2_eff/2)))*Parameters_.a_moire;
                    dis_y = (-0.5*(r1_-(L1_eff/2)) + 0.5*(r2_-(L2_eff/2)))*Parameters_.a_moire;
                    Tij[center_][center_neigh]+=(1.0/(L1_*L2_))*exp(iota_complex*( kx_*(dis_x) +  ky_*(dis_y) ))*eigs_[eigen_no];
                }
            }

            //--------------------

	complex<double> phase_1=conj(Ham_(0,eigen_no))*(1.0/abs(Ham_(0,eigen_no)));
	complex<double> phase_2=conj(Ham_(0, eigen_no-1))*(1.0/abs(Ham_(0,eigen_no-1)));
 	

            for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
                        Psi_state_[r_ind]=0.0;
                        for(int i1=0;i1<l1_;i1++){
                            for(int i2=0;i2<l2_;i2++){
                                kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
                                ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
                                comp = Coordinates_.Nbasis(i1, i2, 0);
				dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
				dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

                   //             Psi_state_[r_ind] += (phase_1*(Ham_(comp,eigen_no)))*exp(iota_complex*( kx_local*(rx_min + rx_ind*d_rx) +  ky_local*(ry_min + ry_ind*d_ry) ));
			Psi_state_[r_ind] += Ham_(comp,eigen_no)*exp(iota_complex*(kx_local*(dis_x) +  ky_local*(dis_y) ));
//			Psi_state_[r_ind] += 1.0; //exp(iota_complex*(kx_local*(dis_x) +  ky_local*(dis_y) ));	
                            }
                        }                
                }
            }



	Update_Bloch_States_Using_Projection(Psi_state_, space_slices, r1_min, d_r1, r2_min, d_r2);






	complex<double> checknorm2=0.0;
                        for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                            for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                                r_ind=r1_ind + (space_slices)*r2_ind;
                                checknorm2 +=   (sqrt(3.0)/2.0)*d_r1*d_r2*Psi_state_[r_ind]*conj(Psi_state_[r_ind]); //may be normailzation need da_vec = dr1 \crossproduct dr2
                            }
                        }
		for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                            for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                                r_ind=r1_ind + (space_slices)*r2_ind;
				Psi_state_[r_ind]=Psi_state_[r_ind]/sqrt(abs(checknorm2));
	}
	}



/*
	 if(n1==-1 && n2==1){
	complex<double> val_temp3_;
        string file_Bloch_out="Bloch_function.txt";
    ofstream FileBlochOut(file_Bloch_out.c_str());
    FileBlochOut<<"#index r1_ind r2_ind r1_dis r2_dis rx  ry u(r)  ....."<<endl;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
            dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
           val_temp3_ = Psi_state_[r_ind]*exp(-1.0*iota_complex*(kx_*(dis_x) +  ky_*(dis_y) ));
	    FileBlochOut<<r_ind<<"  "<<r1_ind<<"  "<<r2_ind<<"   "<< (r1_min + r1_ind*d_r1)/(Parameters_.a_moire)<<"   "<<(r2_min + r2_ind*d_r2)/(Parameters_.a_moire)<<"   "<<dis_x/(Parameters_.a_moire)<<"   "<<dis_y/(Parameters_.a_moire)<<"   "<<val_temp3_.real()<<"   "<<val_temp3_.imag() <<endl;
        }
        FileBlochOut<<endl;
    }
        }

*/

	



            //-------------------



		 int r1_ind_rotPiby3, r1_ind_rot2Piby3, r2_ind_rotPiby3, r2_ind_rot2Piby3, r_ind_rotPiby3, r_ind_rot2Piby3;
            for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){

		if(val_mat[r1_ind][r2_ind]==0){
                    Wnr_state_[r_ind] += 0.0;
			}        
		else{
		r_ind=r1_ind + (space_slices)*r2_ind;
				
		r2_ind_rotPiby3 = 3*((space_slices-1)/2) - r1_ind -r2_ind;
		r1_ind_rotPiby3 = r2_ind;

		r2_ind_rot2Piby3 = 3*((space_slices-1)/2) - r1_ind_rotPiby3 -r2_ind_rotPiby3;
                r1_ind_rot2Piby3 = r2_ind_rotPiby3;

		r_ind_rotPiby3 =  r1_ind_rotPiby3 + (space_slices)*r2_ind_rotPiby3;
		r_ind_rot2Piby3 =  r1_ind_rot2Piby3 + (space_slices)*r2_ind_rot2Piby3;

    //           Wnr_state_[r_ind] += (1.0/sqrt(L1_*L2_))*((1.0/3.0)*(Psi_state_[r_ind] + Psi_state_[r_ind_rotPiby3]+ Psi_state_[r_ind_rot2Piby3]))*val_mat[r1_ind][r2_ind];
		Wnr_state_[r_ind] += (1.0/sqrt(L1_*L2_))*((1.0/1.0)*(Psi_state_[r_ind_rot2Piby3]))*val_mat[r1_ind][r2_ind];
		} 
                }
            }



        }
    }






    checknorm=0.0;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            checknorm += (sqrt(3.0)/2.0)*d_r1*d_r2*Wnr_state_[r_ind]*conj(Wnr_state_[r_ind]);
        }
    }

    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            Wnr_state_[r_ind] = Wnr_state_[r_ind]*(1.0/sqrt(abs(checknorm)));
        }
    }


 
   


     double rx_mean_temp, ry_mean_temp, rsqr_mean_temp;
     rx_mean_temp=0.0; ry_mean_temp=0.0;rsqr_mean_temp=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
    	     dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
             dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
            rx_mean_temp += (sqrt(3.0)/2.0)*d_r1*d_r2*abs(Wnr_state_[r_ind])*abs(Wnr_state_[r_ind])*(dis_x);
            ry_mean_temp += (sqrt(3.0)/2.0)*d_r1*d_r2*abs(Wnr_state_[r_ind])*abs(Wnr_state_[r_ind])*(dis_y);
	    rsqr_mean_temp += (sqrt(3.0)/2.0)*d_r1*d_r2*abs(Wnr_state_[r_ind])*abs(Wnr_state_[r_ind])*
			((dis_x)*(dis_y) +
			(dis_x)*(dis_y)	);
        }
    }

    cout<<"rx_mean (using Wr) = "<<rx_mean_temp/Parameters_.a_moire<<endl;
    cout<<"ry_mean (using Wr) = "<<ry_mean_temp/Parameters_.a_moire<<endl;
    cout<<"rsqr_mean (using Wr) = "<<rsqr_mean_temp/(Parameters_.a_moire*Parameters_.a_moire)<<endl;	

    checknorm=0.0;
    string file_Wnr_out="Wannier_functions_band" + to_string(band) +".txt";
    ofstream FileWNROut(file_Wnr_out.c_str());
    FileWNROut<<"#index r1_ind r2_ind r1_dis r2_dis rx  ry W_r  ....."<<endl;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
	    dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
            dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
            FileWNROut<<r_ind<<"  "<<r1_ind<<"  "<<r2_ind<<"   "<< (r1_min + r1_ind*d_r1)/(Parameters_.a_moire)<<"   "<<(r2_min + r2_ind*d_r2)/(Parameters_.a_moire)<<"   "<<dis_x/(Parameters_.a_moire)<<"   "<<dis_y/(Parameters_.a_moire)<<"   "<<Wnr_state_[r_ind].real()<<"   "<<Wnr_state_[r_ind].imag() <<endl;
            checknorm += (sqrt(3.0)/2.0)*d_r1*d_r2*Wnr_state_[r_ind]*conj(Wnr_state_[r_ind]);
        }
        FileWNROut<<endl;
    }

    FileWNROut<<"#a_moire = "<<Parameters_.a_moire<<endl;
    FileWNROut<<"#norm = "<<checknorm.real()<<"  "<<checknorm.imag()<<endl;
    cout<<"r2_max = "<<r2_max<<endl;
    cout<<"r1_max = "<<r1_max<<endl;
	

    //cout<<"Tij[0][0+a1] = "<<Tij[center_][center_p1]<<endl;

    cout<<"--------------------------------------Tij[center][neigh]--------------------------------------"<<endl;
    for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Tij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;




//assert(false);

	for(int n1=0;n1<L1p2;n1++){
	for(int n2=0;n2<L2p2;n2++){

	for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                   r_ind=r1_ind + (space_slices)*r2_ind;
                   U_Bloch_mat[n1+(n2*L1p2)][r_ind]=0.0;
        }
        }

            cout<<"doing "<<n1<<"  "<<n2<<endl;

            kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3.0)*L1_))  +  n2*(1.0/(sqrt(3.0)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
            HTBCreate();
            char Dflag='V';
            Diagonalize(Dflag);

		   for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
                        Psi_state_[r_ind]=0.0;
                        for(int i1=0;i1<l1_;i1++){
                            for(int i2=0;i2<l2_;i2++){
                                kx_local = kx_ + (-(l1_/2)+i1)*(b1x_) + (-(l2_/2)+i2)*(b2x_);
                                ky_local = ky_ + (-(l1_/2)+i1)*(b1y_) + (-(l2_/2)+i2)*(b2y_);
                                comp = Coordinates_.Nbasis(i1, i2, 0);

				dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
                                dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

                   //             Psi_state_[r_ind] += (phase_1*(Ham_(comp,eigen_no)))*exp(iota_complex*( kx_local*(rx_min + rx_ind*d_rx) +  ky_local*(ry_min + ry_ind*d_ry) ));
                        Psi_state_[r_ind] += Ham_(comp,eigen_no)*exp(iota_complex*(kx_local*(dis_x) +  ky_local*(dis_y) ));

                            }
                        }



                }
            }

        Update_Bloch_States_Using_Projection(Psi_state_, space_slices, r1_min, d_r1, r2_min, d_r2);

        complex<double> checknorm2=0.0;
                        for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                            for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                                r_ind=r1_ind + (space_slices)*r2_ind;
                                checknorm2 +=   (sqrt(3.0)/2.0)*d_r1*d_r2*Psi_state_[r_ind]*conj(Psi_state_[r_ind]);
                            }
                        }
                for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                            for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                                r_ind=r1_ind + (space_slices)*r2_ind;
                                Psi_state_[r_ind]=Psi_state_[r_ind]/sqrt(abs(checknorm2));
        }
        }



	
	for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                   r_ind=r1_ind + (space_slices)*r2_ind;
	dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
        dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
	
                   U_Bloch_mat[n1+(n2*L1p2)][r_ind]=exp(-1.0*iota_complex*( kx_*(dis_x)  + ky_*(dis_y)  ) )*Psi_state_[r_ind];
        }
        }
	}
	}







    //For Uij--------------------------------------------------------------------//

/*
  double q_max, d_q,  d_theta;
    q_max=2.0*PI;
    d_q=(q_max)/(q_slices);
    int theta_slices=200;
    d_theta = (2.0*PI)/(theta_slices);
*/


//----- Maximal localization of Wannier function-----------



cout<<"---------Maximal Localization of Wannier functions started----------"<<endl;


if(false){

//initialization:

string file_BF_out="Bloch_U_functions_band" + to_string(band) +".txt";
    ofstream FileBFOut(file_BF_out.c_str());
    FileBFOut<<"#index rx ry unk  ....."<<endl;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
                                dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

            FileBFOut<<r_ind<<"  "<<(r1_min + r1_ind*d_r1)/(Parameters_.a_moire)<<"   "<<(r2_min + r2_ind*d_r2)/(Parameters_.a_moire)<<"   "<< dis_x/(Parameters_.a_moire)<<"  "<<dis_y/(Parameters_.a_moire)<<"   " <<U_Bloch_mat[3+((3)*L1p2)][r_ind].real()<<"   "<<U_Bloch_mat[3+((3)*L1p2)][r_ind].imag() <<endl;
            
        }
        FileBFOut<<endl;
    }


Mat_1_int b_Neighs_1, b_Neighs_2;
b_Neighs_1.clear();b_Neighs_2.clear();
b_Neighs_1.push_back(1);b_Neighs_2.push_back(1);
b_Neighs_1.push_back(0);b_Neighs_2.push_back(1);
b_Neighs_1.push_back(-1);b_Neighs_2.push_back(0);
b_Neighs_1.push_back(-1);b_Neighs_2.push_back(-1);
b_Neighs_1.push_back(0);b_Neighs_2.push_back(-1);
b_Neighs_1.push_back(1);b_Neighs_2.push_back(0);

int n_ind, n_ind_, n_neigh_ind, n1_neigh, n2_neigh;

for(int n1_=0;n1_<L1_;n1_++){
        for(int n2_=0;n2_<L2_;n2_++){
	int n1, n2;
	n1=n1_;n2=n2_;
	if(n1_==0){
	n1=L1_;
	}
	if(n2_==0){
	n2=L2_;
	}
	n_ind_=n1_+(n2_*L1_);
	n_ind =n1 +(n2 *L1p2);
for(int b_neigh=0;b_neigh<6;b_neigh++){

	n1_neigh = (n1 + b_Neighs_1[b_neigh]+(L1p2))%(L1p2);
        n2_neigh = (n2 + b_Neighs_2[b_neigh]+(L2p2))%(L2p2);
	
	n_neigh_ind = n1_neigh + (n2_neigh*L1p2);

	assert(n_neigh_ind<L1p2*L2p2);


//	cout<<"n, b_neigh, n_neigh_ind = "<<n_ind<<"  "<<b_neigh<<"  "<<n_neigh_ind<<endl;

	M_Mat_old[n_ind_][b_neigh]=0.0;

for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                   r_ind=r1_ind + (space_slices)*r2_ind;

	M_Mat_old[n_ind_][b_neigh]+= (sqrt(3.0)/2.0)*d_r1*d_r2*conj(U_Bloch_mat[n_ind][r_ind])*U_Bloch_mat[n_neigh_ind][r_ind];
	}
	}	
	}


//Unit_mat_old[n_ind]=1.0;
	
}
}



 string file_M_out="M_Mat.txt";
    ofstream FileMOut(file_M_out.c_str());
    FileMOut<<"#index n_ind n1 n2 kx ky M[b]"<<endl;
    for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
           double kx_temp=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3.0)*L1_))  +  n2*(1.0/(sqrt(3.0)*L2_)));
           double ky_temp=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));
 
            FileMOut<<n1+(n2*L1_)<<"  "<<n1<<"  "<<n2<<"  "<<kx_temp<<"  "<<ky_temp<<"  "<< M_Mat_old[n1+(n2*L1_)][2].real()<<"  "<<M_Mat_old[n1+(n2*L1_)][2].imag()<<"  "<<(log(M_Mat_old[n1+(n2*L1_)][2])).imag()<<endl;
            
        }
	FileMOut<<endl;
    }





M_Mat_new=M_Mat_old;
double kx_1=(2.0*PI/Parameters_.a_moire)*(1*(1.0/(sqrt(3.0)*L1_))  +  1*(1.0/(sqrt(3.0)*L2_)));
double ky_1=(2.0*PI/Parameters_.a_moire)*(1*(-1.0/(L1_))  +  1*(1.0/(L2_)));

double b_sqr;
b_sqr= kx_1*kx_1 + ky_1*ky_1 ;
wb_=2.0/(6*b_sqr);
w_=2.0/b_sqr;
alpha_=0.1;

rx_mean_old=0.0;
ry_mean_old=0.0;
r2_mean_old=0.0;
for(int n1=1;n1<L1_;n1++){
        for(int n2=1;n2<L2_;n2++){
        n_ind=n1+(n2*L1_);
for(int b_neigh=0;b_neigh<6;b_neigh++){
	
        double bx_=(2.0*PI/Parameters_.a_moire)*(b_Neighs_1[b_neigh]*(1.0/(sqrt(3.0)*L1_))  +  b_Neighs_2[b_neigh]*(1.0/(sqrt(3.0)*L2_)));
	double by_=(2.0*PI/Parameters_.a_moire)*(b_Neighs_1[b_neigh]*(-1.0/(L1_))  +  b_Neighs_2[b_neigh]*(1.0/(L2_)));

           rx_mean_old +=(-1.0/(L1_*L2_))*wb_*bx_*(log(M_Mat_old[n_ind][b_neigh])).imag();
	   ry_mean_old +=(-1.0/(L1_*L2_))*wb_*by_*(log(M_Mat_old[n_ind][b_neigh])).imag();


	   r2_mean_old +=(1.0/(L1_*L2_))*wb_*(
				1.0*(  1.0  - (abs(M_Mat_old[n_ind][b_neigh])*abs(M_Mat_old[n_ind][b_neigh]))  )   +
				(  ((log(M_Mat_old[n_ind][b_neigh])).imag())*((log(M_Mat_old[n_ind][b_neigh])).imag())    )
					);
}
}
}


cout<<"rx_mean_initial(in a_m) = "<<rx_mean_old/(Parameters_.a_moire)<<endl;
cout<<"ry_mean_initial(in a_m) = "<<ry_mean_old/(Parameters_.a_moire)<<endl;
cout<<"r2_mean_initial(in (a_m)^2) = "<<r2_mean_old/(Parameters_.a_moire*Parameters_.a_moire)<<endl;


assert(true);
if(true){

Sigma_old = r2_mean_old - ((rx_mean_old*rx_mean_old) +  (ry_mean_old*ry_mean_old));
d_Sigma=100.0;

for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
        n_ind=n1+n2*L1_;
Unit_mat_old[n_ind]=1.0;
}}


while(abs(d_Sigma)>=0.0001){

for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
        n_ind=n1+n2*L1_;
	Delta_W[n_ind]=0.0;
for(int b_neigh=0;b_neigh<6;b_neigh++){

        double bx_=(2.0*PI/Parameters_.a_moire)*(b_Neighs_1[b_neigh]*(1.0/(sqrt(3)*L1_))  +  b_Neighs_2[b_neigh]*(1.0/(sqrt(3)*L2_)));
        double by_=(2.0*PI/Parameters_.a_moire)*(b_Neighs_1[b_neigh]*(-1.0/(L1_))  +  b_Neighs_2[b_neigh]*(1.0/(L2_)));

q_mat[n_ind][b_neigh] = (log(M_Mat_new[n_ind][b_neigh])).imag() + bx_*rx_mean_old + by_*ry_mean_old; 

R_mat[n_ind][b_neigh] = 1.0;

T_mat[n_ind][b_neigh] = R_mat[n_ind][b_neigh]*q_mat[n_ind][b_neigh];

A_mat[n_ind][b_neigh] = 0.0;

S_mat[n_ind][b_neigh] = (conj(T_mat[n_ind][b_neigh]) + T_mat[n_ind][b_neigh])*(-0.5*iota_complex);

Delta_W[n_ind]+=(alpha_/w_)*( wb_*(A_mat[n_ind][b_neigh] - S_mat[n_ind][b_neigh]) );


}
//cout<<Delta_W[n_ind]<<endl;
}
}

for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
        n_ind=n1+n2*L1_;
Unit_mat_new[n_ind] =  Unit_mat_old[n_ind]*exp(Delta_W[n_ind]);}}

for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
        n_ind=n1+n2*L1_;
for(int b_neigh=0;b_neigh<6;b_neigh++){

	n1_neigh = (n1 + b_Neighs_1[b_neigh]+L1_)%(L1_);
        n2_neigh = (n2 + b_Neighs_2[b_neigh]+L2_)%(L2_);
        n_neigh_ind = n1_neigh + n2_neigh*L1_;
M_Mat_new[n_ind][b_neigh] = conj(Unit_mat_new[n_ind])*M_Mat_old[n_ind][b_neigh]*Unit_mat_new[n_neigh_ind];
}
}
}


rx_mean_new=0.0;
ry_mean_new=0.0;
r2_mean_new=0.0;
for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
        n_ind=n1+n2*L1_;
for(int b_neigh=0;b_neigh<6;b_neigh++){

        double bx_=(2.0*PI/Parameters_.a_moire)*(b_Neighs_1[b_neigh]*(1.0/(sqrt(3)*L1_))  +  b_Neighs_2[b_neigh]*(1.0/(sqrt(3)*L2_)));
        double by_=(2.0*PI/Parameters_.a_moire)*(b_Neighs_1[b_neigh]*(-1.0/(L1_))  +  b_Neighs_2[b_neigh]*(1.0/(L2_)));

           rx_mean_new +=(-1.0/(L1_*L2_))*wb_*bx_*(log(M_Mat_new[n_ind][b_neigh])).imag();
           ry_mean_new +=(-1.0/(L1_*L2_))*wb_*by_*(log(M_Mat_new[n_ind][b_neigh])).imag();


           r2_mean_new +=(1.0/(L1_*L2_))*wb_*(
                                1.0*(  1.0  - (abs(M_Mat_new[n_ind][b_neigh])*abs(M_Mat_new[n_ind][b_neigh]))  )   +
                                (  ((log(M_Mat_new[n_ind][b_neigh])).imag())*((log(M_Mat_new[n_ind][b_neigh])).imag())    )
                                        );
}
}
}


//M_Mat_old = M_Mat_new;
Unit_mat_old = Unit_mat_new;
rx_mean_old = rx_mean_new;
ry_mean_old = ry_mean_new;
r2_mean_old = r2_mean_new;

Sigma_new= r2_mean_new - ((rx_mean_new*rx_mean_new) +  (ry_mean_new*ry_mean_new));


d_Sigma = (Sigma_new - Sigma_old)/(Parameters_.a_moire*Parameters_.a_moire);

cout<<rx_mean_old/(Parameters_.a_moire)<<"   "<<ry_mean_old/(Parameters_.a_moire)<<"   "<<r2_mean_old/(Parameters_.a_moire*Parameters_.a_moire)<<"   "<<d_Sigma<<endl;

Sigma_old=Sigma_new;
}



for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
//	cout<<(Unit_mat_new[n1+n2*L1_])<<endl;
	
for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
        r_ind=r1_ind + (space_slices)*r2_ind;
        U_Bloch_mat[n1+n2*L1p2][r_ind]=U_Bloch_mat[n1+n2*L1p2][r_ind]*Unit_mat_new[n1+n2*L1_];

//exp(-1.0*iota_complex*( kx_*(rx_min + rx_ind*d_rx) + ky_*(ry_min + ry_ind*d_ry) ) )*Psi_state_[r_ind];
        }
        }
}
}


for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
                        Wnr_state_[r_ind]=0.0;
}
}

//HERE
for(int n1=0;n1<L1_;n1++){
        for(int n2=0;n2<L2_;n2++){
	    kx_=(2.0*PI/Parameters_.a_moire)*(n1*(1.0/(sqrt(3)*L1_))  +  n2*(1.0/(sqrt(3)*L2_)));
            ky_=(2.0*PI/Parameters_.a_moire)*(n1*(-1.0/(L1_))  +  n2*(1.0/(L2_)));

for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                   r_ind=r1_ind + (space_slices)*r2_ind;
dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

                   Psi_state_[r_ind]=U_Bloch_mat[n1+n2*L1p2][r_ind]*exp(1.0*iota_complex*( kx_*(dis_x) + ky_*(dis_y) ) );
}
}

//Update_Bloch_States_Using_Projection(Psi_state_, space_slices, rx_min, d_rx, ry_min, d_ry);
	
            for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;
                        Wnr_state_[r_ind] += (1.0/sqrt(L1_*L2_))*(Psi_state_[r_ind]);
                }
            }

        }
    }



    checknorm=0.0;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            checknorm += (sqrt(3.0)/2.0)*d_r1*d_r2*Wnr_state_[r_ind]*conj(Wnr_state_[r_ind]);
        }
    }

    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
            Wnr_state_[r_ind] = Wnr_state_[r_ind]*(1.0/sqrt(abs(checknorm)));
        }
    }






    checknorm=0.0;
    string file_LWnr_out="LocalizedWannier_functions_band" + to_string(band) +".txt";
    ofstream FileLWNROut(file_LWnr_out.c_str());
    FileLWNROut<<"#index r1 r2 rx ry W_r  ....."<<endl;
    for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
           dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));
            FileLWNROut<<r_ind<<"  "<<(r1_min + r1_ind*d_r1)/(Parameters_.a_moire)<<"   "<<(r2_min + r2_ind*d_r2)/(Parameters_.a_moire)<<"   "<< dis_x/(Parameters_.a_moire)<<"   "<<dis_y/(Parameters_.a_moire)<< "   "<<Wnr_state_[r_ind].real()<<"   "<<Wnr_state_[r_ind].imag() <<endl;
            checknorm += (sqrt(3.0)/2.0)*d_r1*d_r2*Wnr_state_[r_ind]*conj(Wnr_state_[r_ind]);
        }
        FileLWNROut<<endl;

    }

    FileLWNROut<<"#a_moire = "<<Parameters_.a_moire<<endl;
    FileLWNROut<<"#norm = "<<checknorm.real()<<"  "<<checknorm.imag()<<endl;


}
cout<<"---------Maximal Localization of Wannier function done-----------"<<endl;


rx_mean_temp=0.0; ry_mean_temp=0.0;rsqr_mean_temp=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
        for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
            r_ind=r1_ind + (space_slices)*r2_ind;
	 dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

            rx_mean_temp += (sqrt(3.0)/2.0)*d_r1*d_r2*abs(Wnr_state_[r_ind])*abs(Wnr_state_[r_ind])*(dis_x);
            ry_mean_temp += (sqrt(3.0)/2.0)*d_r1*d_r2*abs(Wnr_state_[r_ind])*abs(Wnr_state_[r_ind])*(dis_y);
            rsqr_mean_temp += (sqrt(3.0)/2.0)*d_r1*d_r2*abs(Wnr_state_[r_ind])*abs(Wnr_state_[r_ind])*
                        ((dis_x)*(dis_x) +
                        (dis_y)*(dis_y)   );
        }
    }

    cout<<"rx_mean (using MLWR) = "<<rx_mean_temp/Parameters_.a_moire<<endl;
    cout<<"ry_mean (using MLWR) = "<<ry_mean_temp/Parameters_.a_moire<<endl;
    cout<<"rsqr_mean (using MLWR) = "<<rsqr_mean_temp/(Parameters_.a_moire*Parameters_.a_moire)<<endl;




}
//----------------------------------------



//assert(false);
double Mq_norm, Mq_Sph_norm;
Mq_norm=0.0; Mq_Sph_norm=0.0;

double q_val, theta_val;
double qx_,qy_;
for(int q_ind=0;q_ind<q_slices;q_ind++){
 q_val = q_ind*d_q;
for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
 theta_val = theta_ind*d_theta;

qx_=q_val*cos(theta_val);
qy_=q_val*sin(theta_val);


//( kx_*((1.0/(2.0*sqrt(3.0)))*Parameters_.a_moire) + ky_*( (0.5)*Parameters_.a_moire )
      Mq_SphC[q_ind][theta_ind]=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
      r_ind=r1_ind + (space_slices)*r2_ind;
	 dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
         dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

      Mq_SphC[q_ind][theta_ind] += (sqrt(3.0)/2.0)*d_r1*d_r2*(
      Wnr_state_[r_ind]*conj(Wnr_state_[r_ind]))
                            *exp(iota_complex*(qx_*((dis_x - Wnr_center_x )) + qy_*((dis_y  - Wnr_center_y ))));
                }
            }

	Mq_Sph_norm += sqrt((qx_*qx_) + (qy_*qy_))*d_q*d_theta*abs(Mq_SphC[q_ind][theta_ind]);

}
}
cout << "Mq_Sph_norm = "<<Mq_Sph_norm<<endl;



for(int q_ind=0;q_ind<q_slices;q_ind++){
 q_val = q_ind*d_q;
for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
 theta_val = theta_ind*d_theta;

qx_=q_val*cos(theta_val);
qy_=q_val*sin(theta_val);


//( kx_*((1.0/(2.0*sqrt(3.0)))*Parameters_.a_moire) + ky_*( (0.5)*Parameters_.a_moire )
      Mmq_SphC[q_ind][theta_ind]=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
      r_ind=r1_ind + (space_slices)*r2_ind;
 dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

      Mmq_SphC[q_ind][theta_ind] += (sqrt(3.0)/2.0)*d_r1*d_r2*(
      Wnr_state_[r_ind]*conj(Wnr_state_[r_ind]))
                            *exp(iota_complex*(-1.0*qx_*((dis_x - Wnr_center_x )) - 1.0*qy_*((dis_y  - Wnr_center_y ))));
                }
            }


}
}



double dis_1, dis_2;
int r1_ind_shift, r2_ind_shift;


complex<double> val_temp1_;



for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
	
	if( ((r1_==int((L1_eff/2)+1)) && (r2_==int(L2_eff/2))) ||
	  ((r1_==int((L1_eff/2))) && (r2_==int(L2_eff/2)))
 	 ||   ((r1_==int((L1_eff/2)-1)) && (r2_==int(L2_eff/2))) 
	 ||  ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)+1)))
	 ||  ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)-1)))	
	  ||  ((r1_==int((L1_eff/2)-1)) && (r2_==int((L2_eff/2)+1)))
	 ||  ((r1_==int((L1_eff/2)+1)) && (r2_==int((L2_eff/2)-1)))

  ){

	string Mr_file="Mr_integrand_" + to_string(r1_) +  "_" +  to_string(r2_)  +".txt";
	ofstream MrFILE(Mr_file.c_str());

        center_neigh = (r1_) + L1_eff*(r2_);
//            dis_x = ((sqrt(3.0)/2.0)*(r1_-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
  //          dis_y = (-0.5*(r1_-((1.0*L1_eff)/2.0)) + 0.5*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;



	dis_1 =  (r1_-((1.0*L1_eff)/2.0))*Parameters_.a_moire;
	dis_2 =  (r2_-((1.0*L2_eff)/2.0))*Parameters_.a_moire;


for(int q_ind=0;q_ind<q_slices;q_ind++){
cout<<"q_ind of MqR ("<<q_slices<<")  =  "<<q_ind<<endl;
 q_val = q_ind*d_q;
for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
 theta_val = theta_ind*d_theta;

qx_=q_val*cos(theta_val);
qy_=q_val*sin(theta_val);




//( kx_*((1.0/(2.0*sqrt(3.0)))*Parameters_.a_moire) + ky_*( (0.5)*Parameters_.a_moire)
      MqR_SphC[q_ind][theta_ind][center_neigh]=0.0;
      for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
      for(int r2_ind=0;r2_ind<space_slices;r2_ind++){

	int r1_ind_r = int((r1_ind*d_r1 + (dis_1) + d_r1*0.5)/(d_r1));
        int r2_ind_r = int((r2_ind*d_r2 + (dis_2) + d_r2*0.5)/(d_r2));

/*	if( (rx_ind_l!=rx_ind_r) || (rx_ind!=rx_ind_l) ){
	cout<<"rx_ind_l = "<<rx_ind_l<<endl;
	cout<<"rx_ind_r = "<<rx_ind_r<<endl;
	cout<<"rx_ind = "<<rx_ind<<endl;
	cout<<"dis_x = "<<dis_x<<endl;
	cout<<"dis_y = "<<dis_y<<endl;
	assert(rx_ind_l==rx_ind_r);
	assert(rx_ind==rx_ind_l);}*/

//	int rx_ind_l = rx_ind;
//	int ry_ind_l = ry_ind;
//	int rx_ind_r = rx_ind;
//	int ry_ind_r = ry_ind;	

	dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
        dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));



	if( r1_ind_r >=0 && r1_ind_r<space_slices &&
	    r2_ind_r >=0 && r2_ind_r<space_slices
    ) {


     int r_ind_r=r1_ind_r + (space_slices)*r2_ind_r;
     int r_ind_l=r1_ind + (space_slices)*r2_ind;

      

	MqR_SphC[q_ind][theta_ind][center_neigh] += (sqrt(3.0)/2.0)*d_r1*d_r2*(
      Wnr_state_[r_ind_r]*conj(Wnr_state_[r_ind_l]))
//	*cos((qx_*((dis_x - Wnr_center_x )) + qy_*((dis_y  - Wnr_center_y ))));   
*exp(iota_complex*(qx_*((dis_x - Wnr_center_x)) + qy_*((dis_y  - Wnr_center_y))));


	val_temp1_ = (Wnr_state_[r_ind_r]*conj(Wnr_state_[r_ind_l]))
      *exp(iota_complex*(qx_*((dis_x - Wnr_center_x)) + qy_*((dis_y  - Wnr_center_y))));


	if( q_ind==0 && theta_ind==0 && 
           ( ((r1_==int((L1_eff/2)+1)) && (r2_==int(L2_eff/2)))
	     || ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)+1) )) ||
	     ((r1_==int((L1_eff/2)-1)) && (r2_==int((L2_eff/2)+1) ))
 	   )
	){
	MrFILE<<r1_ind<<"  "<<r2_ind<<"  "<<dis_x/Parameters_.a_moire<<"  "<<dis_y/Parameters_.a_moire<<"  "<<val_temp1_.real()<<"   "<<val_temp1_.imag()<<endl; 
	} 
		}

	else{
	if( q_ind==0 && theta_ind==0 && 
	( ((r1_==int((L1_eff/2)+1)) && (r2_==int(L2_eff/2))) 
            || ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)+1) )) ||
             ((r1_==int((L1_eff/2)-1)) && (r2_==int((L2_eff/2)+1) ))
           )
	){
        MrFILE<<r1_ind<<"  "<<r2_ind<<"  "<<dis_x/Parameters_.a_moire<<"  "<<dis_y/Parameters_.a_moire<<"  "<<0.0<<"   "<<0.0<<endl;
        }
	 }


		}
	if(q_ind==0 && theta_ind==0 && 
	( ((r1_==int((L1_eff/2)+1)) && (r2_==int(L2_eff/2)))
            || ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)+1) )) ||
             ((r1_==int((L1_eff/2)-1)) && (r2_==int((L2_eff/2)+1) ))
           )
	 ){
	MrFILE<<endl;
	}
            }

//        Mq_Sph_norm += sqrt((qx_*qx_) + (qy_*qy_))*d_q*d_theta*abs(Mq_SphC[q_ind][theta_ind]);

}
}
}
}
}


    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;

            Mq_[qx_ind][qy_ind]=0.0;
            for(int r1_ind=0;r1_ind<space_slices;r1_ind++){
                for(int r2_ind=0;r2_ind<space_slices;r2_ind++){
                    r_ind=r1_ind + (space_slices)*r2_ind;

        dis_x = ((r1_min + r1_ind*d_r1)*(sqrt(3.0)/2.0)) + ((r2_min + r2_ind*d_r2)*(sqrt(3.0)/2.0));
        dis_y = ((r1_min + r1_ind*d_r1)*((-1.0)/2.0)) + ((r2_min + r2_ind*d_r2)*((1.0)/2.0));

                    Mq_[qx_ind][qy_ind] += (sqrt(3.0)/2.0)*d_r1*d_r2*(
                                Wnr_state_[r_ind]*conj(Wnr_state_[r_ind])) 
                            *exp(iota_complex*(qx_*((dis_x - Wnr_center_x  )) + qy_*((dis_y  -  Wnr_center_y ))));
                }
            }
		Mq_norm += abs(Mq_[qx_ind][qy_ind])*d_qx*d_qy;
            //Mq_[qx_ind][qy_ind]=1.0;
        }
    }


	cout << "Mq_norm = "<<Mq_norm<<endl;

/*
    double rx_, ry_;
    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;

            Vq_[qx_ind][qy_ind]=0.0;
            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    rx_ = rx_min + rx_ind*d_rx;
                    ry_ = ry_min + ry_ind*d_ry;

                    Vq_[qx_ind][qy_ind] += d_rx*d_ry*((14.3952*1000)/(Parameters_.eps_DE))*
                            ((1.0/(sqrt(rx_*rx_ + ry_*ry_ )+eta_q)) - screening*(1.0/ (sqrt(rx_*rx_ + ry_*ry_ + d_sqr)) ))
                            *exp(iota_complex*(qx_*((rx_)) + qy_*((ry_))));
                }
            }
        }
    }

*/

    double V_;
    for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
            center_neigh = (r1_) + L1_eff*(r2_);



          dis_x = ((sqrt(3.0)/2.0)*(r1_-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
          dis_y = (-0.5*(r1_-((1.0*L1_eff)/2.0)) + 0.5*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
            Uij[center_][center_neigh]=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    Uij[center_][center_neigh]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*abs(Mq_SphC[q_ind][theta_ind])*abs(Mq_SphC[q_ind][theta_ind]))
                            *exp(iota_complex*( qx_*(dis_x) +  qy_*(dis_y) ));
                }
            }
        }
    }




//interaction assisted hoppng Aij
complex<double> val_temp;
int center_temp = (int(L1_eff/2)) + L1_eff*(int(L2_eff/2));
for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            dis_x = ((sqrt(3.0)/2.0)*(r1_-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
            dis_y = (-0.5*(r1_-((1.0*L1_eff)/2.0)) + 0.5*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;


    string AqSphC_file="AqSphC" + to_string(r1_) + "_" +to_string(r2_)+".txt";
    ofstream AqSphCFILE(AqSphC_file.c_str());

    
	

            Aij[center_temp][center_neigh]=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    Aij[center_temp][center_neigh]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*(conj(MqR_SphC[q_ind][theta_ind][center_temp]))*(MqR_SphC[q_ind][theta_ind][center_neigh]));
			    //*cos(-1.0*(qx_*(dis_x*0.5) +  qy_*(dis_y*0.5)));
                            //*exp(-1.0*iota_complex*( qx_*(dis_x*0.5) +  qy_*(dis_y*0.5) ));



val_temp=(1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*(conj(MqR_SphC[q_ind][theta_ind][center_temp]))*(MqR_SphC[q_ind][theta_ind][center_neigh]));
   			   //*cos(-1.0*(qx_*(dis_x*0.5) +  qy_*(dis_y*0.5)));
                            //*exp(-1.0*iota_complex*( qx_*(dis_x*0.5) +  qy_*(dis_y*0.5) ));

	AqSphCFILE<<q_val<<"  "<<theta_val<<"  "<<val_temp.real()<<"   "<<val_temp.imag()<<endl; 
                }
		AqSphCFILE<<endl;
            }
        }
    }






//Exchange interaction Aij
val_temp=0.0;
center_temp = (int(L1_eff/2)) + L1_eff*(int(L2_eff/2));
for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            dis_x = ((sqrt(3.0)/2.0)*(r1_-((1.0*L1_eff)/2.0)) +  (sqrt(3.0)/2.0)*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;
            dis_y = (-0.5*(r1_-((1.0*L1_eff)/2.0)) + 0.5*(r2_-((1.0*L2_eff)/2.0)))*Parameters_.a_moire;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;

            Xij[center_temp][center_neigh]=0.0;
                for(int q_ind=0;q_ind<q_slices;q_ind++){
                 q_val = q_ind*d_q;
                for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
                 theta_val = theta_ind*d_theta;
                qx_=q_val*cos(theta_val);
                qy_=q_val*sin(theta_val);


      //Mq_SphC[q_ind][theta_ind]
                    V_= (2*PI*14.399*1000)/(Parameters_.eps_DE);
                    Xij[center_temp][center_neigh]+= (1.0/(4.0*PI*PI))*(d_q*d_theta)*(V_*(conj(MqR_SphC[q_ind][theta_ind][center_neigh]))*(MqR_SphC[q_ind][theta_ind][center_neigh]));
                            //*cos(-1.0*(qx_*(dis_x*0.5) +  qy_*(dis_y*0.5)));
                            //*exp(-1.0*iota_complex*( qx_*(dis_x*0.5) +  qy_*(dis_y*0.5) ));



                           //*cos(-1.0*(qx_*(dis_x*0.5) +  qy_*(dis_y*0.5)));
                            //*exp(-1.0*iota_complex*( qx_*(dis_x*0.5) +  qy_*(dis_y*0.5) ));

                }
            }
        }
    }








    cout<<endl;
    cout<<"--------------------------------------Uij[center][neigh]--------------------------------------"<<endl;
     for(int r2_=L2_eff-1;r2_>=0;r2_--){
	for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Uij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;



   cout<<endl<<endl;
    cout<<"--------------------------------------Aij[center][neigh]--------------------------------------"<<endl;
    for(int r2_=L2_eff-1;r2_>=0;r2_--){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Aij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;




    cout<<endl<<endl;
    cout<<"--------------------------------------Aij[center][neigh]--------------------------------------"<<endl;
    for(int r2_=L2_eff-1;r2_>=0;r2_--){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Xij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;


    string Mq_file="Mq.txt";
    ofstream MqFILE(Mq_file.c_str());

    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;
            MqFILE<<qx_<<"   "<<qy_<<"   "<<Mq_[qx_ind][qy_ind].real()<<"   "<<Mq_[qx_ind][qy_ind].imag()<<endl;
        }
        MqFILE<<endl;
    }


    string MqSphC_file="MqSphC.txt";
    ofstream MqSphCFILE(MqSphC_file.c_str());

    for(int q_ind=0;q_ind<q_slices;q_ind++){
        q_val = q_ind*d_q;
        for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
        theta_val = theta_ind*d_theta;
            MqSphCFILE<<q_val<<"   "<<theta_val<<"   "<<Mq_SphC[q_ind][theta_ind].real()<<"   "<<Mq_SphC[q_ind][theta_ind].imag()<<endl;
        }
        MqSphCFILE<<endl;
    }



//MqR_SphC[q_ind][theta_ind][center_neigh] 


for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){

        if( ((r1_==int((L1_eff/2)+1)) && (r2_==int(L2_eff/2)))
        ||  ((r1_==int((L1_eff/2))) && (r2_==int(L2_eff/2))) ||
            ((r1_==int((L1_eff/2)-1)) && (r2_==int(L2_eff/2))) 
        ||  ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)+1)))   
        ||  ((r1_==int((L1_eff/2))) && (r2_==int((L2_eff/2)-1)))   
         ||  ((r1_==int((L1_eff/2)-1)) && (r2_==int((L2_eff/2)+1)))
        ||  ((r1_==int((L1_eff/2)+1)) && (r2_==int((L2_eff/2)-1)))

  ){

            center_neigh = (r1_) + L1_eff*(r2_);


    string MqRSphC_file="MqRSphC_site_" + to_string(r1_) +"_"+ to_string(r2_)+".txt";
    ofstream MqRSphCFILE(MqRSphC_file.c_str());

    for(int q_ind=0;q_ind<q_slices;q_ind++){
        q_val = q_ind*d_q;
        for(int theta_ind=0;theta_ind<theta_slices;theta_ind++){
        theta_val = theta_ind*d_theta;
            MqRSphCFILE<<q_val<<"   "<<theta_val<<"   "<<MqR_SphC[q_ind][theta_ind][(r1_) + L1_eff*(r2_)].real()<<"   "<<MqR_SphC[q_ind][theta_ind][r1_ + L1_eff*(r2_)].imag()<<endl;
        }
        MqRSphCFILE<<endl;
    }


}

}}



    //For Uij--------------------------------------------------------------------//
/*
    double qx_,qy_;
    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;

            Mq_[qx_ind][qy_ind]=0.0;
            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    r_ind=rx_ind + (space_slices)*ry_ind;

                    Mq_[qx_ind][qy_ind] += d_rx*d_ry*(
                                Wnr_state_[0][r_ind]*conj(Wnr_state_[0][r_ind]) +
                            Wnr_state_[1][r_ind]*conj(Wnr_state_[1][r_ind]) )
                            *exp(iota_complex*(qx_*((rx_min + rx_ind*d_rx)) + qy_*((ry_min + ry_ind*d_ry))));
                }
            }
            //Mq_[qx_ind][qy_ind]=1.0;
        }
    }



    double rx_, ry_;
    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;

            Vq_[qx_ind][qy_ind]=0.0;
            for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
                for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
                    rx_ = rx_min + rx_ind*d_rx;
                    ry_ = ry_min + ry_ind*d_ry;

                    Vq_[qx_ind][qy_ind] += d_rx*d_ry*((14.3952*1000)/(Parameters_.eps_DE))*
                            ((1.0/(sqrt(rx_*rx_ + ry_*ry_ )+eta_q)) - screening*(1.0/ (sqrt(rx_*rx_ + ry_*ry_ + d_sqr)) ))
                            *exp(iota_complex*(qx_*((rx_)) + qy_*((ry_))));
                }
            }
        }
    }



    //double Vq_;
    for(int r1_=0;r1_<L1_eff;r1_++){
        for(int r2_=0;r2_<L2_eff;r2_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            dis_x = ((sqrt(3.0)/2.0)*(r1_-(L1_eff/2)) +  (sqrt(3.0)/2.0)*(r2_-(L2_eff/2)))*Parameters_.a_moire;
            dis_y = (-0.5*(r1_-(L1_eff/2)) + 0.5*(r2_-(L2_eff/2)))*Parameters_.a_moire;

            //dis_x = ((r1_-(L1_eff/2)))*Parameters_.a_moire;
            //dis_y = ((r2_-(L2_eff/2)))*Parameters_.a_moire;
            Uij[center_][center_neigh]=0.0;
            for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
                qx_=qx_min + qx_ind*d_qx;
                for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
                    qy_=qy_min + qy_ind*d_qy;

                    //Vq_= (2*PI*14.3952*1000)/(Parameters_.eps_DE*(sqrt(qx_*qx_ + qy_*qy_)+eta_q));
                    Uij[center_][center_neigh]+= (1.0/(4.0*PI*PI))*(d_qx*d_qy)*(Vq_[qx_ind][qy_ind]*abs(Mq_[qx_ind][qy_ind])*abs(Mq_[qx_ind][qy_ind]))
                            *exp(iota_complex*( qx_*(dis_x) +  qy_*(dis_y) ));
                }
            }
        }
    }





    cout<<endl;
    cout<<"--------------------------------------Uij[center][neigh]--------------------------------------"<<endl;
    for(int r2_=0;r2_<L2_eff;r2_++){
        for(int r1_=0;r1_<L1_eff;r1_++){
            center_neigh = (r1_) + L1_eff*(r2_);
            cout<<Uij[center_][center_neigh]<<"  ";
        }
        cout<<endl;
    }

    cout<<"----------------------------------------------------------------------------------"<<endl;



    string Vq_file="Vq.txt";
    ofstream VqFILE(Vq_file.c_str());

    for(int qx_ind=0;qx_ind<q_slices;qx_ind++){
        qx_=qx_min + qx_ind*d_qx;
        for(int qy_ind=0;qy_ind<q_slices;qy_ind++){
            qy_=qy_min + qy_ind*d_qy;
            VqFILE<<qx_<<"   "<<qy_<<"   "<<Vq_[qx_ind][qy_ind].real()<<"   "<<Vq_[qx_ind][qy_ind].imag()<<endl;
        }
        VqFILE<<endl;
    }




    string Vr_file="Vr.txt";
    ofstream VrFILE(Vr_file.c_str());

    for(int rx_ind=0;rx_ind<space_slices;rx_ind++){
        for(int ry_ind=0;ry_ind<space_slices;ry_ind++){
            rx_ = rx_min + rx_ind*d_rx;
            ry_ = ry_min + ry_ind*d_ry;

            VrFILE<<rx_<<"   "<<ry_<<"   "<< ((14.3952*1000)/(Parameters_.eps_DE))*((1.0/(sqrt(rx_*rx_ + ry_*ry_ )+eta_q)) - screening*(1.0/ (sqrt(rx_*rx_ + ry_*ry_ + d_sqr)) )) <<endl;
        }
        VrFILE<<endl;
    }


    */
    //---------------------------------------------------------------------------


}

void Hamiltonian::Hoppings(){

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
