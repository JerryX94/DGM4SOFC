// DGM4SOFC.cpp: 主项目文件。
// 阳极质量传输计算引入离散化求解

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>

#define INPUTFILENAME1 "inputie.dat"
#define INPUTFILENAME2 "input.dat"
#define OUTPUTFILENAME "output.dat"
#define ALPHA 1.//Anodic transfer coefficient
#define BETA .5//Cathodic transfer coefficient
#define D 2.E-6//m
#define EACT 1.E5//J/mol
#define ERRMAX 1.E-6
#define F 96485.//C/mol
#define K1 1.92//Anodic Thickness Coefficient
#define K2 1.77//Overpotential Coefficient 1.75?
#define KP .01;//WGSR Coefficient
#define MAXIT 100000//Max Iterature Num.
#define MAXNIE 30//Operation Current Density Input Num.
#define MAXOVERP 1.//Specific Max Overpotential
#define N 4//Component Num.
#define NDISC 5//Discretization Num.
#define NREACT 2//Eletrochemical Reaction Num.
#define PA 1.013E5//Pa
#define PA0 1.013E5//Pa
#define PI 3.1416
#define PORO .45
#define R 8.314//J/mol.K
#define T 1073.15//K
#define TREF 1073.15//K
#define THICKA 1.1E-3//m
#define THICKE 1.E-5//m
#define PODIVTO .1

int InIe(double *IeQueue);
void Init();
void InWg();
double CBg();
double CDKn(int species);
double CDij(int speciesi, int speciesj);
int Iter(double Bg, double DKn[], double Dij[][N]);
void Exha(int speciesi);
void CWg(double Z);
double CMu(double *Xave);
double Err(double *X0, double *X1, double *DX, double P0, double P1, \
	int it);
double CEta();
void Enqueue(int it, int iIE, double OutputQueue[][N + 3]);
void Outp(int nIE, double OutputQueue[][N + 3]);

class Gas{
public:
	double P0;//Pa
	double DP;//Pa
	double DPlist[NDISC];//Pa
	double X0[N];
	double X1[N];
	double Xlist[N][NDISC];
	double TotalNg;//mol/m2.s
	double Ng[N];//mol/m2.s
	double DNg[N];//mol/m2.s
	double Mg[N];//g/mol
	double Vd[N];//cm3/mol
	double ReactOrder[N];
	double J0ref[N];//A/m3
	double ViscCoef[N][7];
};

class React1{
public:
	short Ne;
	short Reactant;
	short Product;
	double w;//Species Current Factor
};

class React2{
public:
	short Switch;
	short Reactant[2];
	short Product[2];
	double Keq;
};

Gas mixture;
React1 elechem[NREACT];
React2 wgsr;
double IE = 0.;

int main(){
	double IeQueue[MAXNIE];
	double OutputQueue[MAXNIE][N + 3];
	double DKn[N];//m2/s
	double Dij[N][N] = { 0 };//m2/s
	int nIE = InIe(IeQueue);
	for (int iIE = 0; iIE < nIE; iIE++){
		IE = IeQueue[iIE];
		Init();
		double Bg = CBg();
		for (int i = 0; i < N; i++)
			DKn[i] = CDKn(i);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++){
				if (i != j) Dij[i][j] = Dij[j][i] = CDij(i, j);
			}
		Enqueue(Iter(Bg, DKn, Dij), iIE, OutputQueue);
	}
	Outp(nIE, OutputQueue);
	return 0;
}

int InIe(double *IeQueue){
	using namespace std;
	ifstream fin(INPUTFILENAME1);
	int nIE;
	fin >> nIE;
	for (int i = 0; i < nIE; i++){
		fin >> IeQueue[i];
		IeQueue[i] *= 10000.;//单位变换
	}
	fin.close();
	return nIE;
}

void Init(){
	using namespace std;
	ifstream fin(INPUTFILENAME2);
	for (int i = 0; i < NREACT; i++)
		fin >> elechem[i].Ne >> elechem[i].Reactant >> elechem[i].Product;
	fin >> wgsr.Switch;
	if (wgsr.Switch){
		fin >> wgsr.Reactant[0] >> wgsr.Reactant[1];
		fin >> wgsr.Product[0] >> wgsr.Product[1];
		InWg();
	}
	mixture.P0 = PA;
	mixture.DP = 0.;
	for (int i = 0; i < N; i++){
		fin >> mixture.X0[i];
		mixture.X1[i] = mixture.X0[i];
	}
	for (int i = 0; i < N; i++){
		mixture.DNg[i] = 0.;
		fin >> mixture.Ng[i];//Input Stoichiometric Coefficient
		short flag = 0;
		for (int j = 0; j < NREACT; j++){
			if ((elechem[j].Reactant == i) || (elechem[j].Product == i)){
				mixture.Ng[i] *= IE / (elechem[j].Ne * F);
				flag++;
				break;
			}
		}
		if (!flag) mixture.Ng[i] = 0.;
	}
	for (int i = 0; i < N; fin >> mixture.Mg[i++]);
	for (int i = 0; i < N; fin >> mixture.Vd[i++]);
	for (int i = 0; i < N; fin >> mixture.ReactOrder[i++]);
	for (int i = 0; i < N; fin >> mixture.J0ref[i++]);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 7; j++)
			fin >> mixture.ViscCoef[i][j];
	double rt = 0.;
	for (int i = 0; i < NREACT; i++){
		rt += mixture.X0[elechem[i].Reactant] * elechem[i].Ne;
	}
	mixture.TotalNg = 0.;
	for (int i = 0; i < NREACT; i++){//初值设置是否正确？
		elechem[i].w = mixture.X0[elechem[i].Reactant] * elechem[i].Ne / rt;
		mixture.Ng[elechem[i].Reactant] *= elechem[i].w;
		mixture.Ng[elechem[i].Product] *= elechem[i].w;
		mixture.TotalNg += mixture.Ng[elechem[i].Reactant] * elechem[i].Ne;
	}
	fin.close();
	return;
}

void InWg(){
	double Z = 1000. / T - 1.;
	wgsr.Keq = -.2935 * Z * Z * Z + .635 * Z * Z + 4.1788 * Z + .3169;
	wgsr.Keq = exp(wgsr.Keq);
}

double CBg(){
	double Bg = (pow(D, 2) / 180.) / (PORO * PORO / pow(1 - PORO, 2));
	//double Bg = 1.7E-10;
	return Bg;
}

double CDKn(int species){
	double eDKn = D / 3. * sqrt(8. * R * T / (PI * mixture.Mg[species]\
		/ 1000.));
	eDKn *= PODIVTO;
	return eDKn;
}

double CDij(int speciesi, int speciesj){
	double k = 1.0E-3;
	double mij = 1. / mixture.Mg[speciesi] + 1. / mixture.Mg[speciesj];
	double nume = k * pow(T, 1.75) * sqrt(mij);
	double deno = pow(mixture.Vd[speciesi], 1. / 3)\
		+ pow(mixture.Vd[speciesj], 1. / 3);
	deno = deno * deno * PA / PA0;
	double eDij = PODIVTO * nume / deno * 1.E-4;
	return eDij;
}

int Iter(double Bg, double DKn[], double Dij[][N]){
	int it = 0;//迭代次数
	double Z = K1 * THICKA / NDISC;//离散等效电极厚度
	for (int k = 0; k < NDISC; k++){
		double DX[N] = { 0 };//摩尔因数差
		double Xave[N];//平均摩尔因数
		double Ngave[N];//平均流量
		double DPnew = mixture.DP / NDISC;//更新压力差
		double DPold = DPnew;
		double Pave;//平均压力
		do{
			if (it > MAXIT){
				it = -1;
				break;
			}
			//水气转换反应，先迭代两次使浓度分布偏差减小后再行实施
			if (wgsr.Switch){
				if (it > 1) CWg(Z);
			}
			//调整权重系数
			for (int i = 0; i < NREACT; i++){
				elechem[i].w = mixture.Ng[elechem[i].Reactant] * elechem[i].Ne\
					/ mixture.TotalNg;
			}
			//获取新值
			for (int i = 0; i < N; i++){
				mixture.X1[i] = mixture.X0[i] + DX[i];
				if (mixture.X1[i] < 0.){
					mixture.X1[0] += mixture.X1[i];
					mixture.X1[i] = 0.;
					Exha(i);
				}
			}
			//调整权重系数
			for (int i = 0; i < NREACT; i++){
				elechem[i].w = mixture.Ng[elechem[i].Reactant] * elechem[i].Ne\
					/ mixture.TotalNg;
			}
			//更新计算指标
			DPold = DPnew;
			for (int i = 0; i < N; i++){
				Xave[i] = .5 * (mixture.X0[i] + mixture.X1[i]);
				Ngave[i] = mixture.Ng[i] - .5 * mixture.DNg[i];
			}
			Pave = mixture.P0 + .5 * DPold;
			double Mu = CMu(Xave);//重新计算粘度
			//计算压力差
			double nume = 0.;
			double deno = 1.;
			for (int i = 0; i < N; i++)
				nume += Ngave[i] / DKn[i];
			for (int i = 0; i < N; i++)
				deno += (Bg / Mu) * Pave * Xave[i] / DKn[i];//已由v3.7修正
			DPnew = -(nume / deno) * Z * R * T;
			//计算浓度差
			for (int i = 1; i < N; i++){
				DX[i] = 0.;
				for (int j = 0; j < N; j++){
					if (j != i){
						nume = Xave[j] * Ngave[i] - Xave[i] * Ngave[j];
						deno = Dij[i][j];
						DX[i] += nume / deno;
					}
				}
				DX[i] = R * T * (DX[i] + Ngave[i] / DKn[i]) + Xave[i]\
					* DPold * (1. + Bg * Pave / (Mu * DKn[i])) / Z;
				DX[i] *= -Z / Pave;
			}
			DX[0] = 0.;
			for (int i = 1; i < N; DX[0] -= DX[i++]);
		} while (Err(mixture.X0, mixture.X1, DX, DPold, DPnew, ++it)\
			> ERRMAX);
		for (int i = 0; i < N; i++){
			mixture.Xlist[i][k] = mixture.X0[i];
			mixture.X0[i] = mixture.X1[i];
		}
		mixture.DPlist[k] = DPnew;
		mixture.DP += DPnew;
	}
	for (int i = 0; i < N; i++)
		mixture.X0[i] = mixture.Xlist[i][0];
	return it;
}

void Exha(int speciesi){
	double relax = .99;
	double dNg = 0.;
	for (int i = 0; i < NREACT; i++)
		if (elechem[i].Reactant == speciesi){
			dNg += mixture.Ng[elechem[i].Reactant] * elechem[i].Ne\
				/ elechem[0].Ne * (1. - relax);
			mixture.Ng[elechem[i].Reactant] *= relax;
			mixture.Ng[elechem[i].Product] *= relax;
		}
	mixture.Ng[elechem[0].Reactant] += dNg;
	mixture.Ng[elechem[0].Product] -= dNg;
	return;
}

void CWg(double Z){
	const int WGSplit = 5;//空间离散节点数
	double r[WGSplit];
	for (int i = 0; i < WGSplit; i++){
		double x1 = mixture.X0[wgsr.Product[0]]\
			+ i * (mixture.X1[wgsr.Product[0]]\
			- mixture.X0[wgsr.Product[0]]) / (WGSplit - 1);
		double x2 = mixture.X0[wgsr.Reactant[0]]\
			+ i * (mixture.X1[wgsr.Reactant[0]]\
			- mixture.X0[wgsr.Reactant[0]]) / (WGSplit - 1);
		double x3 = mixture.X0[wgsr.Reactant[1]]\
			+ i * (mixture.X1[wgsr.Reactant[1]]\
			- mixture.X0[wgsr.Reactant[1]]) / (WGSplit - 1);
		double x4 = mixture.X0[wgsr.Product[1]]\
			+ i * (mixture.X1[wgsr.Product[1]]\
			- mixture.X0[wgsr.Product[1]]) / (WGSplit - 1);
		r[i] = .0171 * exp(-103191. / (R * T))\
			* (x2 * x3 - x1 * x4 / wgsr.Keq) * PA * PA;
	}
	double rtotal = 0.;
	for (int i = 0; i < WGSplit; rtotal += r[i++]);
	rtotal *= KP;
	for (int i = 0; i < 2; i++){
		mixture.DNg[wgsr.Reactant[i]] = -rtotal / WGSplit * Z;
		mixture.DNg[wgsr.Product[i]] = rtotal / WGSplit * Z;
	}
	return;
}

double CMu(double *Xave){
	double MuSp[N] = { 0 };
	double phi[N][N] = { 0 };
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 7; j++)
			MuSp[i] += mixture.ViscCoef[i][j] * pow(T / 1000., j);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (j != i){
				double nume = 1.;
				nume += sqrt(MuSp[i] / MuSp[j]) * pow(mixture.Mg[j]\
					/ mixture.Mg[i], .25);
				nume *= nume;
				double deno = sqrt(8. * (1 + mixture.Mg[i]\
					/ mixture.Mg[j]));
				phi[i][j] = nume / deno;
			}
	double Mum = 0.;
	for (int i = 0; i < N; i++){
		double deno = 1.;
		for (int j = 0; j < N; j++){
			if (j != i) deno += Xave[j] * phi[i][j] / Xave[i];
		}
		Mum += MuSp[i] / deno * 1.E-7;//单位变换：微泊(muP)->Pa.s
	}
	return Mum;
}

double Err(double *X0, double *X1, double *DX, double DP0, double DP1, \
	int it){
	double Err2 = 0.;
	for (int i = 0; i < N; i++){
		double diff = X0[i] + DX[i] - X1[i];
		Err2 += (diff * diff) / (X0[i] * X0[i] + 1.E-10);
	}
	double diff = DP0 - DP1;
	Err2 += (diff * diff) / (DP0 * DP0 + 1.E-10);
	if (!(it % 100000))
		std::cout << "Iteration = " << it << ", Err = " << Err2 << "\n";
	return Err2;
}

double CEta(){
	double TotalEta = 1.E10;
	double SpeciesEta = 0.;
	double OhmEta = 0.;
	double IonCond = 3.34E4 * exp(-10300. / T);
	for (int i = 0; i < NREACT; i++){
		double I0refb = IonCond * R * T\
			/ ((ALPHA + BETA) * elechem[i].Ne * F);
		I0refb *= mixture.J0ref[elechem[i].Reactant];
		I0refb *= pow(PA * mixture.X0[elechem[i].Reactant] / PA0, \
			mixture.ReactOrder[elechem[i].Reactant]);
		I0refb *= pow(PA * mixture.X0[elechem[i].Product] / PA0, \
			mixture.ReactOrder[elechem[i].Product]);
		I0refb = sqrt(I0refb);
		double I0 = I0refb * exp(-EACT / R * (1. / T - 1 / TREF));
		double Icomp = elechem[i].w * IE;
		double sub0, sub = 0.;
		do{
			sub0 = sub;
			sub = exp(-BETA * sub) * mixture.X1[elechem[i].Product]\
				/ mixture.X0[elechem[i].Product] + Icomp\
				/ (I0 * mixture.X0[elechem[i].Reactant]);
			sub = log(sub * mixture.X0[elechem[i].Reactant]\
				/ mixture.X1[elechem[i].Reactant]) / ALPHA;
		} while (abs(sub - sub0) / (sub0 + 1.E-10) > ERRMAX);
		SpeciesEta = sub * R * T / (elechem[i].Ne * F);
		if (SpeciesEta > MAXOVERP) SpeciesEta = MAXOVERP;
		if (TotalEta > SpeciesEta) TotalEta = SpeciesEta;
	}
	OhmEta = IE * THICKE / IonCond;
	TotalEta += OhmEta;
	return TotalEta;
}

void Enqueue(int it, int iIE, double OutputQueue[][N + 3]){
	double Eta = K2 * CEta();
	OutputQueue[iIE][0] = (double)it;
	for (int i = 0; i < N; i++)
		OutputQueue[iIE][i + 1] = mixture.X1[i];
	OutputQueue[iIE][N + 1] = mixture.DP;
	OutputQueue[iIE][N + 2] = Eta;
}

void Outp(int nIE, double OutputQueue[][N + 3]){
	using namespace std;
	ofstream fout(OUTPUTFILENAME);
	/*
	fout << "Iterating " << it << " Times Totally." << endl << endl;
	for (int i = 0; i < N; i++)
	fout << "Species Molar Fraction " << i << \
	" on Electrolyte Interface : " << mixture.X1[i] << endl;
	fout << endl << "Pressure Diffrence : " << mixture.DP << " Pa\n";
	fout << endl << "Concentration Overpotential : " << Eta << " V";
	*/
	for (int k = 0; k < nIE; k++){
		fout << (int)OutputQueue[k][0] << " ";
		for (int i = 1; i < N + 3; fout << OutputQueue[k][i++] << " ");
		fout << endl;
	}
	fout.close();
	return;
}
